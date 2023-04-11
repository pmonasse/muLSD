/*----------------------------------------------------------------------------  
  This code is part of the following publication and was subject
  to peer review:
  "Multiscale line segment detector for robust and accurate SfM" by
  Yohann Salaun, Renaud Marlet, and Pascal Monasse
  ICPR 2016
  
  Copyright (c) 2016 Yohann Salaun <yohann.salaun@imagine.enpc.fr>
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the Mozilla Public License as
  published by the Mozilla Foundation, either version 2.0 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  Mozilla Public License for more details.
  
  You should have received a copy of the Mozilla Public License
  along with this program. If not, see <https://www.mozilla.org/en-US/MPL/2.0/>.

  ----------------------------------------------------------------------------*/

#include "mlsd.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
using namespace std;

/// Minimum number of pixels for a valid cluster
static const int MIN_SIZE_CLUSTER=10;

// --- Cluster class ---

/// Create a cluster from a list of pixels. This constructor needs to find
/// the corresponding rectangle.
Cluster::Cluster(image_double angles, image_double modgrad, double logNT,
                 const vector<point>& d, double t, double p,
                 int idx, int s)
: pixels(d), theta(t), prec(p),
  index(idx), scale(s), merged(false) {
    region2rect(pixels.data(), (int)d.size(), modgrad, t, p, p/M_PI, &rec);
    nfa = rect_nfa(&rec, angles, logNT);
}

/// The rectangle containing the list of pixels is already known (after LSD).
Cluster::Cluster(image_double angles, double logNT,
                 const point* d, int dsize, rect& r, int idx, int s)
: pixels(d,d+dsize), rec(r), theta(r.theta), prec(r.prec),
  index(idx), scale(s), merged(false) {
    nfa = rect_nfa(&rec, angles, logNT);
}

double Cluster::length() const {
    return hypot(rec.x1-rec.x2, rec.y1-rec.y2);
}

Segment Cluster::toSegment() const {
    return Segment(rec.x1 + 0.5, rec.y1 + 0.5, rec.x2 + 0.5, rec.y2 + 0.5,
                   rec.width, rec.p, nfa, scale);
}

void Cluster::setUsed(image_char& used) const {
    for(size_t j=0; j<pixels.size(); j++)
        used->data[pixels[j].x + pixels[j].y*used->xsize] = USED;
}

/// Build the union of the current cluster and the ones with indices in
/// \a indexMerge.
Cluster Cluster::united(const vector<Cluster>& clusters,
                        const set<int>& indexMerge,
                        image_double angles, image_double modgrad,
                        double logNT) const {
    // merge pixel sets
    vector<point> mergedPixels = pixels;
    for(set<int>::iterator it=indexMerge.begin(); it!=indexMerge.end(); ++it) {
        const Cluster& c = clusters[*it];
        mergedPixels.insert(mergedPixels.end(),c.pixels.begin(),c.pixels.end());
    }
    return Cluster(angles, modgrad, logNT, mergedPixels, theta, prec, 
                   clusters.size(), scale);
}

// --- Point2d class: pixel with double coordinates ---

struct Point2d {
    double x;
    double y;
    Point2d(double X, double Y): x(X), y(Y) {}
    Point2d() {}
};

inline void operator+=(Point2d& p, const Point2d& q) {
    p.x += q.x;
    p.y += q.y;
}

inline void operator-=(Point2d& p, const Point2d& q) {
    p.x += q.x;
    p.y += q.y;
}

inline void operator-=(Point2d& p, const point& q) {
    p.x -= q.x;
    p.y -= q.y;
}

inline Point2d operator+(Point2d p, const Point2d& q) {
    p += q;
    return p;
}

inline Point2d operator-(Point2d p, const Point2d& q) {
    p -= q;
    return p;
}

static Point2d operator*(const double d, const Point2d& p) {
    return Point2d(d*p.x, d*p.y);
}

inline double crossProd(const Point2d &p, const Point2d& q) {
    return p.x*q.y - p.y*q.x;
}

static bool insideRect(Point2d pu, Point2d pd,
                       Point2d qu, Point2d qd, const point& t) {
    pu-=t; pd-=t; qu-=t; qd-=t;
    return crossProd(pu, qu)*crossProd(pd, qd) < 0 &&
           crossProd(pu, pd)*crossProd(qu, qd) < 0;
}

Point2d getCenter(const Cluster& c) {
    const rect& r = c.rectangle();
    return Point2d(r.x,r.y);
}

Point2d getSlope(const Cluster& c) {
    const rect& r = c.rectangle();
    double dx=r.x2-r.x1, dy=r.y2-r.y1;
    return (fabs(dx)>fabs(dy))? Point2d(1, dy/dx): Point2d(dx/dy, 1);
}

// --- ROI class ---

// Functor to sort clusters by decreasing NFA value.
class CompareClusters {
    const std::vector<Cluster>& c;
public:
    CompareClusters(const std::vector<Cluster>& v): c(v) {}
    bool operator()(int i, int j) const { return c[i].getNFA()>c[j].getNFA(); }
};

/// Region of interest: a rectangle inside the image.
class ROI {
    static const int CLUSTER_NULL;

    int scale;
    // angles info pointers
    image_double angles, modgrad;
    double logNT;

    // rectangle parameters
    int width; 
    int xMin, xMax, yMin, yMax; ///< Axis-aligned bounding box
    Point2d p_up, p_down, q_up, q_down; ///< Corners of rectangle

    // aligned pixels
    double theta, prec;
    vector<int> definedPixels;
    vector<int> pixelCluster; /// Maps pixel to cluster index
    vector<Cluster> clusters;

    // conversion from 2D point to 1D index
    int pixelToIndex(int x, int y) const {
        return (x-xMin) + (y-yMin)*width;
    }
    int pixelToIndex(const point& p) const {
        return pixelToIndex(p.x,p.y);
    }
    int sizeBB() const { return width*(yMax-yMin+1); }
    bool inBB(int x, int y) const {
        return (xMin<=x && x<=xMax && yMin<=y && y<=yMax);
    }
    bool inBB(const point& p) const {
        return inBB(p.x,p.y);
    }

    void computeRectangle(const Segment& rawSegment);
    void findAlignedPixels(vector<point>& alignedPixels);
    void agregatePixels(const vector<point>& alignedPixels);

    void findIntersect(const Cluster& c, set<int>& inter, bool postLSD) const;

public:
    ROI(const Segment& seg,
        image_double a, image_double m, double lNT, int s);
    ROI(const vector<Cluster>& c,
        image_double a, image_double m, double lNT, int s);

    bool isVoid() const { return clusters.empty(); }
    void mergeClusters(bool postLSD);
    bool filterClusters(vector<Cluster>& filtered, image_char& used,
                        double log_eps, bool postLSD);
};

const int ROI::CLUSTER_NULL=-1;

/// First constructor: a single segment (from previous scale).
ROI::ROI(const Segment& rawSegment, image_double a, image_double m,
         double lNT, int s)
: scale(s), angles(a), modgrad(m), logNT(lNT) {
    // gradient parameters
    theta = rawSegment.angle;
    prec = M_PI*rawSegment.prec;

    computeRectangle(rawSegment);
    pixelCluster = vector<int>(sizeBB(), CLUSTER_NULL);
    vector<point> alignedPixels;
    findAlignedPixels(alignedPixels);
    agregatePixels(alignedPixels);
}

/// Second constructor: the ROI is the full image, clusters already computed.
ROI::ROI(const vector<Cluster>& c, image_double a, image_double m,
         double lNT, int s)
: scale(s), angles(a), modgrad(m), logNT(lNT) {
    width = (int)a->xsize;
    xMin = yMin = 0;
    xMax = a->xsize-1;
    yMax = a->ysize-1;
    clusters = c;
    pixelCluster = vector<int>(sizeBB(), CLUSTER_NULL);
    for(size_t i=0; i<clusters.size(); i++)
        for(size_t j=0; j<clusters[i].getPixels().size(); j++) {
            int idx = pixelToIndex( clusters[i].getPixels()[j] );
            pixelCluster[idx] = clusters[i].getIndex();
        }
}

/// Compute rectangle and bounding box for segment \a seg.
void ROI::computeRectangle(const Segment& seg) {
    // Compute rectangle corners
    Point2d P(seg.x1,seg.y1), Q(seg.x2,seg.y2);
    Point2d delta = (seg.width/2.0) / seg.length * (Q-P);
    p_up   = P+delta;
    p_down = P-delta;
    q_up   = Q+delta;
    q_down = Q-delta;

    // compute min/max along x/y axis
    xMin = max((int)floor(min({p_up.x, p_down.x, q_up.x, q_down.x})),0);
    xMax = min((int)ceil (max({p_up.x, p_down.x, q_up.x, q_down.x})),
               (int)angles->xsize-1);
    yMin = max((int)floor(min({p_up.y, p_down.y, q_up.y, q_down.y})),0);
    yMax = min((int)ceil (max({p_up.y, p_down.y, q_up.y, q_down.y})),
               (int)angles->ysize-1);
    width = xMax-xMin+1;
}

/// Insert in \a alignedPixels the pixels inside \c rect whose gradient
/// direction is aligned. Pixels inside the rectangle with non-null gradient
/// are registered in \c definedPixels. 
void ROI::findAlignedPixels(vector<point>& alignedPixels) {
    for(point p = {xMin,yMin}; p.y<=yMax; p.y++)
        for(p.x=xMin; p.x<=xMax; p.x++)
            if(insideRect(p_up, p_down, q_up, q_down, p)) {
                int i = p.x + p.y*angles->xsize;
                if(angles->data[i] != NOTDEF) {
                    definedPixels.push_back(i);
                    if(angle_diff(angles->data[i], theta) < prec) {
                        alignedPixels.push_back(p);
                        pixelCluster[pixelToIndex(p)] = CLUSTER_NULL;
                    }
                }
            }
}

/// Find 8-connected components of alignedPixels, i.e., the initial clusters.
/// The ones less than \c MIN_SIZE_CLUSTER are just ignored.
void ROI::agregatePixels(const vector<point>& alignedPixels) {
    for(size_t i=0; i<alignedPixels.size(); i++) {
        int index = pixelToIndex(alignedPixels[i]);
        if(pixelCluster[index] != CLUSTER_NULL) continue;

        // initialize new cluster
        vector<point> data(1, alignedPixels[i]);
        const int idCluster = pixelCluster[index] = clusters.size();

        // growing region from the current pixel
        for(size_t j=0; j<data.size(); ++j) {
            point seed = data[j];
            // look inside 8-neighbourhood
            for (int dx=-1; dx<=1; dx++) {
                for (int dy=-1; dy<=1; dy++) {
                    if(dx==0 && dy==0) continue;
                    point p = {seed.x + dx, seed.y + dy};
                    int idx = pixelToIndex(p);
                    if(inBB(p) && pixelCluster[idx]==CLUSTER_NULL) {
                        pixelCluster[idx] = idCluster;
                        data.push_back(p);
                    }
                }
            }
        }
        // suppress clusters with too few pixels
        if(data.size() < MIN_SIZE_CLUSTER)
            for(size_t j=0; j<data.size(); j++)
                pixelCluster[pixelToIndex(data[j])] = CLUSTER_NULL;
        else
            clusters.push_back(Cluster(angles, modgrad, logNT, data, theta,
                                       prec, idCluster, scale));
    }
}

void ROI::findIntersect(const Cluster& c, set<int>& inter, bool postLSD) const {
    const Point2d center = getCenter(c);
    const Point2d step = getSlope(c);
    const double theta = c.getTheta();
    const double prec = c.getPrec();

    // find clusters intersecting with current cluster direction
    for (int s = -1; s <= 1; s += 2){
        const Point2d signed_step = s*step;
        const int cidx = c.getIndex();
        Point2d p = center;
        while(true) {
            p += signed_step;
            int X = floor(p.x + 0.5);
            int Y = floor(p.y + 0.5);
            if(! inBB(X,Y)) break;

            int idx = pixelCluster[pixelToIndex(X,Y)];
            if(idx!=CLUSTER_NULL && idx!=cidx && !clusters[idx].isMerged()){
                if(postLSD && angle_diff(clusters[idx].getTheta(),theta)>prec)
                    continue;
                inter.insert(idx);
                // Post lsd processing: only look for next neighbor?
                if(postLSD) break;
            }
        }
    }
}

void ROI::mergeClusters(bool postLSD) {
    // sort clusters by decreasing NFA
    vector<int> clusterStack(clusters.size());
    std::iota(clusterStack.begin(), clusterStack.end(), 0);
    std::sort(clusterStack.begin(), clusterStack.end(),
              CompareClusters(clusters));

    for (size_t i=0; i<clusterStack.size(); i++) {
        int currIndex = clusterStack[i];
        Cluster& c = clusters[currIndex];
        if( c.isMerged() ) continue;

        set<int> intersect;
        findIntersect(c, intersect, postLSD);
        if( intersect.empty() ) continue;

        // compute merged cluster
        Cluster megaCluster;

        set<int>::const_iterator it=intersect.begin(), end=intersect.end();
        if(!postLSD) {
            megaCluster = c.united(clusters, intersect, angles, modgrad, logNT);
            bool toMerge = megaCluster.getNFA()>c.getNFA();
            for(it=intersect.begin(); it!=end && toMerge; ++it)
                if(megaCluster.getNFA()<=clusters[*it].getNFA())
                    toMerge = false;
                
            if(!toMerge) continue;
            // labelize merged clusters as merged
            for(it=intersect.begin(); it!=end; ++it)
                clusters[*it].setMerged();
        } else {
            bool toMerge = false;
            for(; it!=end && !toMerge; ++it) {
                set<int> temp;
                temp.insert(*it);
                megaCluster = c.united(clusters, temp, angles, modgrad, logNT);
                toMerge = megaCluster.getNFA() > max(c.getNFA(),
                                                     clusters[*it].getNFA());
                // labelize merged clusters as merged
                if(toMerge)
                    clusters[*it].setMerged();
            }
            if(! toMerge) continue;
        }

        c.setMerged();
        clusterStack.push_back(clusters.size());
        clusters.push_back(megaCluster);

        // labelize clustered points with their new label
        for (size_t j=0; j<megaCluster.getPixels().size(); j++) {
            const point& p = megaCluster.getPixels()[j];
            if(inBB(p)) // due to width reduction, there can be some issue
                pixelCluster[pixelToIndex(p)] = megaCluster.getIndex();
        }
    }
}

/// Append valid clusters of the ROI to \a filtered.
/// A valid cluster has not been merged and has a sufficiently large -log(NFA).
/// If no cluster is valid in the ROI, discard all pixels aligned with the ROI,
/// so that the next LSD pass will not detect anything there.
bool ROI::filterClusters(vector<Cluster>& filtered, image_char& used,
                         double log_eps, bool postLSD) {
    bool adding=false;
    for (size_t i=0; i<clusters.size(); i++) {
        // cluster has been merged and has no more meaning
        if(clusters[i].isMerged()) continue;
        // NFA is too low
        if (clusters[i].getNFA() <= log_eps) continue;

        adding = true;
        clusters[i].setUsed(used);
        clusters[i].setIndex(filtered.size());
        filtered.push_back(clusters[i]);
    }

    // if no valid cluster, set pixels to used to avoid useless computation
    if(!postLSD && !adding)
        for(size_t j=0; j<definedPixels.size(); j++)
            used->data[definedPixels[j]] = USED;
    return adding;
}

// --- Free functions ---

vector<Cluster> refineRawSegments(const vector<Segment>& rawSegments,
                                  vector<Segment>& finalLines, int i_scale,
                                  image_double angles, image_double modgrad,
                                  image_char& used,
                                  double logNT, double log_eps) {
    vector<Cluster> clusters;
    for (size_t i=0; i<rawSegments.size(); i++){
        const Segment& seg = rawSegments[i];
        if(seg.scale != i_scale - 1) {
            finalLines.push_back(seg);
            continue;
        }

        ROI roi(seg, angles, modgrad, logNT, i_scale);

        // if no pixels are detected, keep the former one
        if (roi.isVoid()) {
            finalLines.push_back(seg);
            continue;
        }

        // 1. merge greedily aligned clusters that should be (in NFA meaning)
        roi.mergeClusters(false);
        // 2. compute associated lines
        if(! roi.filterClusters(clusters, used, log_eps, false))
            // no lines have been found ==> keep the former one
            finalLines.push_back(seg);
    }
    return clusters;
}

void mergeSegments(vector<Cluster>& clusters,
                   double segment_length_threshold, int i_scale,
                   image_double angles, image_double modgrad, image_char& used,
                   double logNT, double log_eps) {
    // filter refinedLines wrt segment line
    if(segment_length_threshold>0) {
        vector<Cluster> temp;
        for (size_t i=0; i<clusters.size(); i++) {
            if(clusters[i].length() > segment_length_threshold) {
                clusters[i].setIndex(temp.size());
                temp.push_back(clusters[i]);
            }
        }
        clusters = temp;
    }
    // merge greedily aligned clusters that should be merged (in NFA meaning)
    ROI roi(clusters, angles, modgrad, logNT, i_scale);
    roi.mergeClusters(true);
    clusters.clear();
    roi.filterClusters(clusters, used, log_eps, true);
}

void denseGradientFilter(vector<int>& noisyTexture, int w, int h, 
                         const image_double& angles, image_char& used,
                         int xsize, int ysize, int N){
    if(noisyTexture.size() == 0){
        noisyTexture = vector<int>(N, 0);
        const int size_kernel = max(ceil(0.01*(w+h)/4), 5.0); // 0.5% of image size
        const int step = size_kernel/2;
        const double thresh_noise =0.99;
        // integral image for faster computation
        vector<double> integral_image(N);
        for (int i = 0; i < xsize; i++){
            for (int j = 0; j < ysize; j++){
                integral_image[i+j*xsize] = 0;
                if(angles->data[i + j * xsize] != NOTDEF &&
                        used->data[i + j * xsize] == NOTUSED){
                    integral_image[i+j*xsize] = 1;
                }
                integral_image[i+j*xsize] += (i==0)? 0:integral_image[i-1+j*xsize];
                integral_image[i+j*xsize] += (j==0)? 0:integral_image[i+(j-1)*xsize];
                integral_image[i+j*xsize] -= (i==0 || j==0)? 0:integral_image[i-1+(j-1)*xsize];
            }
        }
        // compute scores for each points
        for(int x = 0; x < xsize; x++){
            for(int y = 0; y < ysize; y++){
                int xMin = max(0.,double(x-step)), xMax = min(double(xsize-1), double(x+step)), yMin = max(0., double(y-step)), yMax = min(double(ysize-1), double(y+step));
                int n_gradient = integral_image[xMin + yMin*xsize] + integral_image[xMax + yMax*xsize]
                        -integral_image[xMin + yMax*xsize] - integral_image[xMax + yMin*xsize];
                if(n_gradient > thresh_noise*(xMax-xMin)*(yMax-yMin)){
                    used->data[x + y*xsize] = USED;
                    noisyTexture[x + y*xsize] = 1;
                }
            }
        }
        // expand a bit the noisy texture area
        vector<bool> expandedArea(N, false);
        const int expansion = 10;//10;
        for(int x = 0; x < xsize; x++){
            for(int y = 0; y < ysize; y++){
                if(noisyTexture[x + y*xsize] != 1){continue;}
                for(int xx = -expansion; xx <= expansion; xx++){
                    for(int yy = -expansion; yy <= expansion; yy++){
                        if(x+xx < 0 || x+xx >= xsize || y+yy < 0 || y+yy >= ysize){continue;}
                        expandedArea[x+xx + (y+yy)*xsize] = true;
                    }
                }
            }
        }
        for(int x = 0; x < xsize; x++){
            for(int y = 0; y < ysize; y++){
                if(noisyTexture[x + y*xsize] == 1){continue;}
                if(expandedArea[x + y*xsize]){
                    noisyTexture[x + y*xsize] = 1;
                }
            }
        }
        noisyTexture[0] = xsize;
    }
    else{
        double tempXsize = noisyTexture[0];
        double ratio = tempXsize/xsize;
        for(int x = 0; x < xsize; x++){
            for(int y = 0; y < ysize; y++){
                int X = (x+0.5)*ratio;
                int Y = (y+0.5)*ratio;
                if(noisyTexture[X + Y*tempXsize] == 1){
                    used->data[x + y*xsize] = USED;
                }
            }
        }
    }
}
