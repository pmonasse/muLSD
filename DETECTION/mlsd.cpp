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
using namespace std;

inline void operator+=(Point2d& p, const Point2d& q) {
    p.x += q.x;
    p.y += q.y;
}

inline void operator-=(Point2d& p, const point& q) {
    p.x -= q.x;
    p.y -= q.y;
}

Point2d operator*(const double d, const Point2d& p) {
    return Point2d(d*p.x, d*p.y);
}

inline double crossProd(const Point2d &p, const Point2d &q) {
    return p.x*q.y - p.y*q.x;
}

bool insideRect(Point2d pu, Point2d pd,
                Point2d qu, Point2d qd, const point &t) {
    pu-=t; pd-=t; qu-=t; qd-=t;
    return crossProd(pu, qu)*crossProd(pd, qd) < 0 &&
           crossProd(pu, pd)*crossProd(qu, qd) < 0;
}

// --- NFA_params class ---

NFA_params::NFA_params() {
    computed = false;
}

NFA_params::NFA_params(double t, double p){
    theta = t;
    prec = p;
    computed = false;
}

void NFA_params::computeNFA(rect &rec, image_double angles, double logNT) {
    value = rect_nfa(&rec, angles, logNT);
    computed = true;
}

double NFA_params::getValue() const {
    if(! computed) std::cerr << "ERROR in NFA computation" << std::endl;
    return value;
}

double NFA_params::getTheta() const {
    return theta;
}

double NFA_params::getPrec() const {
    return prec;
}

// --- Cluster class ---
Cluster::Cluster(){}

Cluster::Cluster(image_double angles, image_double modgrad, double logNT,
                 const vector<point>& d, double t, double p,
                 int idx, int s, double nfaSeparate)
: data(d), nfa(t,p), nfa_separated_clusters(nfaSeparate), index(idx), scale(s),
  merged(false) {
    // compute rectangle
    region2rect(data.data(), (int)data.size(), modgrad, t, p, p/M_PI, &rec);
    nfa.computeNFA(rec, angles, logNT);
}

Cluster::Cluster(image_double angles, double logNT,
                 const point* d, int dsize, rect &r, int idx, int s)
: data(d,d+dsize), rec(r), nfa(r.theta,r.prec), nfa_separated_clusters(-1),
  index(idx), scale(s), merged(false) {
    nfa.computeNFA(rec, angles, logNT);
}

Cluster Cluster::mergedCluster(const vector<Cluster> &clusters,
                               const set<int>&indexToMerge,
                               image_double angles, image_double modgrad,
                               double logNT) const {
    // merge pixel sets
    vector<point> mergedData = data;
    double nfaSeparatedClusters = nfa.getValue() - log10(double(rec.n + 1));
    for (set<int>::iterator it = indexToMerge.begin(); it != indexToMerge.end(); it++){
        for (size_t j=0; j<clusters[(*it)].data.size(); j++)
            mergedData.push_back(clusters[(*it)].data[j]);
        nfaSeparatedClusters += clusters[(*it)].nfa.getValue() - log10(double(clusters[(*it)].rec.n + 1));
    }
    return Cluster(angles, modgrad, logNT, mergedData, nfa.getTheta(),
                   nfa.getPrec(), clusters.size(), scale,
                   nfaSeparatedClusters);
}

bool Cluster::isToMerge(image_double angles, double logNT) {
    // check with log_nfa
    int N = rec.width;
    double dx = rec.x1 - rec.x2;
    double dy = rec.y1 - rec.y2;
    int M = sqrt(dx*dx + dy*dy);
    double diff_binom = -(5 / 2 * log10(double(N*M)) - log10(2.0)) + nfa_separated_clusters + log10(double(rec.n + 1));
    double fusion_score = diff_binom - nfa.getValue();

    if (fusion_score > 0) {
        /* try to reduce width */
        rect reduce_rec;
        float delta = 0.5;
        NFA_params nfa_optimized(nfa);
        rect_copy(&rec, &reduce_rec);
        for (int n = 0; n < 5; n++)
        {
            if ((reduce_rec.width - delta) >= 0.5)
            {
                reduce_rec.width -= delta;
                nfa_optimized.computeNFA(reduce_rec, angles, logNT);
                if (nfa_optimized.getValue() > nfa.getValue()){
                    rect_copy(&reduce_rec, &rec);
                    fusion_score = diff_binom - nfa_optimized.getValue();
                    if (fusion_score < 0){ break; }
                }
            }
        }
        if (fusion_score > 0)
            return false;

        // recompute cluster data
        rect_iter * rec_it;
        data.clear();
        for(rec_it=ri_ini(&rec); !ri_end(rec_it); ri_inc(rec_it)) {
            if(rec_it->x >= 0 && rec_it->y >= 0 &&
               rec_it->x<(int)angles->xsize && rec_it->y<(int)angles->ysize) {
                if (isaligned(rec_it->x,rec_it->y, angles,rec.theta,rec.prec)) {
                    point p = {rec_it->x, rec_it->y};
                    data.push_back(p);
                }
            }
        }
        ri_del(rec_it);
    }
    nfa.computeNFA(rec, angles, logNT);
    return true;
}

double Cluster::length() const {
    float dx = rec.x1 - rec.x2;
    float dy = rec.y1 - rec.y2;
    return sqrt(dx*dx + dy*dy);
}

Segment Cluster::toSegment() const {
    return Segment(rec.x1 + 0.5, rec.y1 + 0.5, rec.x2 + 0.5, rec.y2 + 0.5,
                   rec.width, rec.p, nfa.getValue(), scale);
}

bool Cluster::isMerged() const {
    return merged;
}

void Cluster::setMerged() {
    merged = true;
}

double Cluster::getNFA() const {
    return nfa.getValue();
}

int Cluster::getIndex() const {
    return index;
}

int Cluster::getScale() const {
    return scale;
}

double Cluster::getTheta() const {
    return nfa.getTheta();
}

double Cluster::getPrec() const {
    return nfa.getPrec();
}

Point2d Cluster::getCenter() const {
    return Point2d(rec.x, rec.y);
}

Point2d Cluster::getSlope() const {
    double c_dx, c_dy;
    c_dx = rec.x2-rec.x1;
    c_dy = rec.y2-rec.y1;

    if(fabs(c_dx) > fabs(c_dy))
        return Point2d(1, c_dy/c_dx);
    return Point2d(c_dx/c_dy, 1);
}

const vector<point>& Cluster::getData() const {
    return data;
}

void Cluster::setUsed(image_char used) const {
    for(size_t j=0; j<data.size(); j++)
        used->data[data[j].x + data[j].y * used->xsize] = USED;
}

void Cluster::setIndex(int i) {
    index = i;
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
    static const int CLUSTER_NULL = -1;

    int scale;

    // angles info pointers
    image_double angles, modgrad;
    double logNT;

    // rectangle parameters
    int width; 
    int xMin, xMax, yMin, yMax; // Axis-aligned bounding box
    Point2d p_up, p_down, q_up, q_down;

    // aligned pixels
    double theta, prec;
    vector<point> alignedPixels;
    vector<int> definedPixels;
    vector<int> pixelCluster;

    // clusters of pixels
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
    void computeAlignedPixels();

public:
    ROI(const Segment &rawSegment, image_double a, image_double m,
        double lNT, int s);
    ROI(image_double a, image_double m, double lNT,
        const vector<Cluster>& c, int s);

    bool isVoid() const { return alignedPixels.size() <= 1; }
    void agregatePixels();
    void mergeClusters(bool post_lsd);
    bool keepCorrectLines(vector<Cluster> &refinedLines, image_char used,
                          double log_eps, bool post_lsd);
};

void ROI::computeRectangle(const Segment& rawSegment) {
    Point2d P(rawSegment.x1,rawSegment.y1), Q(rawSegment.x2,rawSegment.y2);

    // compute rectangle extremities
    double length = rawSegment.length;
    double dx = Q.x - P.x;
    double dy = Q.y - P.y;
    double dW = (rawSegment.width / 2 + 1) / length;

    p_up   = Point2d(P.x + dW*dy - dx / length, P.y - dW*dx - dy / length);
    p_down = Point2d(P.x - dW*dy - dx / length, P.y + dW*dx - dy / length);
    q_up   = Point2d(Q.x + dW*dy + dx / length, Q.y - dW*dx + dy / length);
    q_down = Point2d(Q.x - dW*dy + dx / length, Q.y + dW*dx + dy / length);

    // compute min/max along x/y axis
    xMin = max((int)floor(min({p_up.x, p_down.x, q_up.x, q_down.x})),0);
    xMax = min((int)ceil (max({p_up.x, p_down.x, q_up.x, q_down.x})),
               (int)angles->xsize-1);
    yMin = max((int)floor(min({p_up.y, p_down.y, q_up.y, q_down.y})),0);
    yMax = min((int)ceil (max({p_up.y, p_down.y, q_up.y, q_down.y})),
               (int)angles->ysize-1);
    width = xMax-xMin+1;
}

void ROI::computeAlignedPixels() {
    pixelCluster = vector<int>(sizeBB(), NOTDEF);
    for(point p = {xMin,yMin}; p.y<=yMax; p.y++) {
        for(p.x=xMin; p.x<=xMax; p.x++) {
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
    }
}

ROI::ROI(const Segment &rawSegment, image_double a, image_double m,
         double lNT, int s) {
    // images parameters
    angles = a;
    modgrad = m;
    logNT = lNT;

    // gradient parameters
    theta = rawSegment.angle;
    prec = M_PI*rawSegment.prec;

    scale = s;

    computeRectangle(rawSegment);
    computeAlignedPixels();
}

ROI::ROI(image_double a, image_double m, double lNT,
         const vector<Cluster>& c, int s) {
    // images parameters
    angles = a;
    modgrad = m;
    logNT = lNT;

    scale = s;

    width = (int)a->xsize;
    xMin = yMin = 0;
    xMax = a->xsize-1;
    yMax = a->ysize-1;
    pixelCluster = vector<int>(sizeBB(), NOTDEF);
    clusters = c;
    for(size_t i=0; i<clusters.size(); i++) {
        for(size_t j=0; j<clusters[i].getData().size(); j++) {
            int idx = pixelToIndex( clusters[i].getData()[j] );
            pixelCluster[idx] = clusters[i].getIndex();
        }
    }
}

void ROI::agregatePixels() {
    for(size_t i=0; i<alignedPixels.size(); i++) {
        int index = pixelToIndex(alignedPixels[i]);

        // if pixel already clusterized or has no angle defined go on
        if (pixelCluster[index] != CLUSTER_NULL){ continue; }

        // initialize new cluster
        vector<point> data(1, alignedPixels[i]);
        pixelCluster[index] = clusters.size();

        // growing region from the current pixel
        for(size_t currIndex=0; currIndex<data.size(); ++currIndex) {
            point seed = data[currIndex];
            // look inside 8-neighbourhood
            for (int dx=-1; dx<=1; dx++) {
                for (int dy=-1; dy<=1; dy++) {
                    if(dx==0 && dy==0) continue;
                    point p = {seed.x + dx, seed.y + dy};
                    int idx = pixelToIndex(p);

                    // add neighbor that have correct gradient directions and is not already in the cluster
                    if (!inBB(p) || pixelCluster[idx]!=CLUSTER_NULL)
                        continue;

                    pixelCluster[idx] = clusters.size();
                    data.push_back(p);
                }
            }
        }
        // suppress clusters with too few pixels
        if (data.size() <= 10){
            for(size_t j=0; j<data.size(); j++)
                pixelCluster[pixelToIndex(data[j])] = NOTDEF;
            continue;
        }

        Cluster c(angles, modgrad, logNT, data, theta, prec, clusters.size(), scale, -1);
        clusters.push_back(c);
    }
}

void ROI::mergeClusters(bool post_lsd) {
    // sort clusters by decreasing NFA
    vector<int> clusterStack(clusters.size());
    std::iota(clusterStack.begin(), clusterStack.end(), 0);
    std::sort(clusterStack.begin(), clusterStack.end(),
              CompareClusters(clusters));

    // merge clusters if necessary
    for (size_t i=0; i<clusterStack.size(); i++) {
        int currIndex = clusterStack[i];
        Cluster& c = clusters[currIndex];

        if( c.isMerged() ) continue;
        // define line of intersection
        const Point2d center = c.getCenter();
        const Point2d step = c.getSlope();
        const double theta = c.getTheta();
        const double prec = c.getPrec();

        // find clusters intersecting with current cluster direction
        set<int> intersect;
        for (int s = -1; s <= 1; s += 2){
            const Point2d signed_step = s*step;
            const int cidx = c.getIndex();
            Point2d p = center;
            while (true){
                p += signed_step;
                int X = floor(p.x + 0.5);
                int Y = floor(p.y + 0.5);
                if(! inBB(X,Y)) break;

                int idx = pixelCluster[pixelToIndex(X,Y)];
                if (idx!=NOTDEF && idx!=cidx && !clusters[idx].isMerged()){
                    if(post_lsd &&
                       angle_diff(clusters[idx].getTheta(),theta) > prec)
                        continue;

                    intersect.insert(idx);
                    // Post lsd processing: only look for next neighbor?
                    if(post_lsd) break;
                }
            }
        }
        // if no cluster intersections
        if( intersect.empty() ) continue;

        // compute merged cluster
        Cluster megaCluster;

        set<int>::const_iterator it=intersect.begin(), end=intersect.end();
        if(!post_lsd) {
            megaCluster = c.mergedCluster(clusters, intersect, angles, modgrad, logNT);
            if (! megaCluster.isToMerge(angles, logNT)) continue;
            // labelize merged clusters as merged
            for (; it!=end; it++)
                clusters[(*it)].setMerged();
        } else {
            bool toMerge = false;
            for(; it!=end && !toMerge; it++) {
                set<int> temp;
                temp.insert(*it);
                megaCluster = c.mergedCluster(clusters, temp, angles, modgrad, logNT);
                toMerge = megaCluster.isToMerge(angles, logNT);
                // labelize merged clusters as merged
                if(toMerge)
                    clusters[(*it)].setMerged();
            }
            if(! toMerge) continue;
        }

        c.setMerged();
        clusterStack.push_back(clusters.size());
        clusters.push_back(megaCluster);

        // labelize clustered points with their new label
        for (size_t j=0; j<megaCluster.getData().size(); j++) {
            const point& p = megaCluster.getData()[j];
            if(inBB(p)) // due to width reduction, there can be some issue
                pixelCluster[pixelToIndex(p)] = megaCluster.getIndex();
        }
    }
}

bool ROI::keepCorrectLines(vector<Cluster> &refinedLines, image_char used,
                           double log_eps, bool post_lsd) {
    const size_t nRefinedLines = refinedLines.size();
    for (size_t i=0; i<clusters.size(); i++) {
        // cluster has been merged and has no more meaning
        if (clusters[i].isMerged()) continue;
        // NFA is too low
        if (clusters[i].getNFA() <= log_eps) continue;

        clusters[i].setUsed(used);
        clusters[i].setIndex(refinedLines.size());
        refinedLines.push_back(clusters[i]);
    }

    // if no lines kept return false and set pixels to used to avoid useless computation
    if (!post_lsd && nRefinedLines == refinedLines.size()) {
        for (size_t j=0; j<definedPixels.size(); j++)
            used->data[definedPixels[j]] = USED;
        return false;
    }
    return true;
}

// --- Free functions ---

vector<Cluster> refineRawSegments(const vector<Segment> &rawSegments,
                                  vector<Segment> &finalLines, int i_scale,
                                  image_double angles, image_double modgrad,
                                  image_char used,
                                  double logNT, double log_eps) {
    vector<Cluster> refinedLines;
    for (size_t i_line=0; i_line<rawSegments.size(); i_line++){
        const Segment& seg = rawSegments[i_line];
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

        // 1. compute clusters
        roi.agregatePixels();
        // 2. merge greedily aligned clusters that should be (in NFA meaning)
        roi.mergeClusters(false);
        // 3. compute associated lines
        if(! roi.keepCorrectLines(refinedLines, used, log_eps, false))
            // no lines have been found ==> keep the former one
            finalLines.push_back(seg);
    }
    return refinedLines;
}

void mergeSegments(vector<Cluster> &refinedLines,
                   double segment_length_threshold, int i_scale,
                   image_double angles, image_double modgrad, image_char& used,
                   double logNT, double log_eps) {
    // filter refinedLines wrt segment line
    if(segment_length_threshold>0) {
        vector<Cluster> temp;
        for (size_t i=0; i<refinedLines.size(); i++) {
            if(refinedLines[i].length() > segment_length_threshold) {
                refinedLines[i].setIndex(temp.size());
                temp.push_back(refinedLines[i]);
            }
        }
        refinedLines = temp;
    }
    // merge greedily aligned clusters that should be merged (in NFA meaning)
    ROI clusters(angles, modgrad, logNT, refinedLines, i_scale);
    clusters.mergeClusters(true);
    refinedLines.clear();
    clusters.keepCorrectLines(refinedLines, used, log_eps, true);
}

void denseGradientFilter(vector<int>& noisyTexture, int w, int h, 
                         const image_double& angles, const image_char& used,
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
