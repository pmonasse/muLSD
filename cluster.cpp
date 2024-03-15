// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file cluster.cpp
 * @brief muLSD: clusters (segment candidates)
 * @author Yohann Salaun <yohann.salaun@imagine.enpc.fr>
 *         Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2016, 2023-2024
 */

#include "cluster.hpp"
#include <queue>
#include <algorithm>
#include <numeric>
#include <cmath>
using namespace std;

/// Minimum number of pixels for a valid cluster
static const int MIN_SIZE_CLUSTER=10;

// --- Cluster class ---

/// Create a cluster from a list of pixels. This constructor needs to find
/// the enclosing rectangle.
Cluster::Cluster(Dimage angles, Dimage modgrad, double logNT,
                 const std::vector<point>& d, double t, double p, int idx)
: pixels(d), index(idx), merged(false) {
    region2rect(pixels.data(), (int)d.size(), modgrad, t, p, p/M_PI, &rec);
    nfa = rect_nfa(&rec, angles, logNT);
}

/// Create a cluster from a list of pixels.
/// The rectangle containing the list of pixels is already known (after LSD).
Cluster::Cluster(Dimage angles, double logNT,
                 const point* d, int dsize, rect& r, int idx)
: pixels(d,d+dsize), rec(r), index(idx), merged(false) {
    nfa = rect_nfa(&rec, angles, logNT);
}

/// Length of enclosing rectangle.
double Cluster::length() const {
    return hypot(rec.x1-rec.x2, rec.y1-rec.y2);
}

/// Convert to Segment.
Segment Cluster::toSegment() const {
    return Segment(rec.x1 + 0.5, rec.y1 + 0.5, rec.x2 + 0.5, rec.y2 + 0.5,
                   rec.width, rec.p, nfa);
}

/// Mark in image \c used the footprint of the cluster (its pixels).
void Cluster::setUsed(Cimage& used) const {
    for(size_t j=0; j<pixels.size(); j++)
        used->data[pixels[j].x + pixels[j].y*used->xsize] = USED;
}

/// Build the union of the current cluster and the ones with index in
/// \a indexMerge.
Cluster Cluster::united(const vector<Cluster>& clusters,
                        const set<int>& indexMerge,
                        Dimage angles, Dimage modgrad,
                        double logNT) const {
    // merge pixel sets
    vector<point> mergedPixels = pixels;
    for(set<int>::iterator it=indexMerge.begin(); it!=indexMerge.end(); ++it) {
        const Cluster& c = clusters[*it];
        mergedPixels.insert(mergedPixels.end(),c.pixels.begin(),c.pixels.end());
    }
    return Cluster(angles, modgrad, logNT, mergedPixels, getTheta(), rec.prec,
                   clusters.size());
}

// --- Point2d class: pixel with double coordinates ---

/// Represent either a point or a vector
struct Point2d {
    double x;
    double y;
    Point2d(double X, double Y): x(X), y(Y) {}
    Point2d() {}
};

// Various operators for compact geometric constructions.

/// Point+vector -> point.
inline void operator+=(Point2d& p, const Point2d& q) {
    p.x += q.x;
    p.y += q.y;
}
inline Point2d operator+(Point2d p, const Point2d& q) {
    p += q;
    return p;
} ///< Point+point -> point
/// Point-vector -> point
inline void operator-=(Point2d& p, const Point2d& q) {
    p.x -= q.x;
    p.y -= q.y;
}
/// Same with integer vector \c q.
inline void operator-=(Point2d& p, const point& q) {
    p.x -= q.x;
    p.y -= q.y;
}
inline Point2d operator-(Point2d p, const Point2d& q) {
    p -= q;
    return p;
} ///< Point - point -> point

/// Scale \c p, considered as vector, by factor \a d.
inline Point2d operator*(double d, const Point2d& p) {
    return Point2d(d*p.x, d*p.y);
}

/// Vector directly orthogonal to vector \c p.
inline Point2d orthogonal(Point2d p) {
    std::swap(p.x,p.y);
    p.y = -p.y;
    return p;
}

/// Third component of cross-product of vectors, determinant of [p q].
inline double crossProd(const Point2d &p, const Point2d& q) {
    return p.x*q.y - p.y*q.x;
}

/// Test if a pixel position \c t is inside a rectangle defined by 4 vertices.
static bool insideRect(Point2d p1, Point2d p2,
                       Point2d q1, Point2d q2, const point& t) {
    p1-=t; p2-=t; q1-=t; q2-=t;
    return crossProd(p1, q1)*crossProd(p2, q2) < 0 &&
           crossProd(p1, p2)*crossProd(q1, q2) < 0;
}

/// Get center of rectangle: weigthed centroid of pixels in the cluster.
static Point2d getCenter(const Cluster& c) {
    const rect& r = c.rectangle();
    return Point2d(r.x,r.y);
}

/// Vector for sampling points along the rectangle direction. It has a step
/// of 1 along one dimension and less than 1 in the other, depending on whether
/// the direction is more vertical or horizontal.
static Point2d getSlope(const Cluster& c) {
    const rect& r = c.rectangle();
    double dx=r.x2-r.x1, dy=r.y2-r.y1;
    return (fabs(dx)>fabs(dy))? Point2d(1, dy/dx): Point2d(dx/dy, 1);
}

// --- ROI class ---

// Functor to sort clusters by decreasing NFA value.
class CompareClusters {
    const vector<Cluster>& c;
public:
    CompareClusters(const vector<Cluster>& v): c(v) {}
    bool operator()(int i, int j) const { return c[i].getNFA()<c[j].getNFA(); }
};

/** \brief Region of interest: a rectangle inside the image. */
/// It is used to merge clusters within the rectangle.
/// There are two constructors:
/// 1. with a \a Segment: ROI is the rectangle defined by the segment.
/// 2. with a collections of clusters: ROI is the full image.
/// Fields #p1, #p2, #q1, #q2 and #definedPixels are used only for case 1.
class ROI {
    static const int CLUSTER_NULL;

    // angles info pointers
    Dimage angles, modgrad;
    double logNT;

    // rectangle parameters
    int width; ///< Logical width of 2D array #pixelCluster
    int xMin, xMax, yMin, yMax; ///< Axis-aligned bounding box
    Point2d p1, p2, q1, q2; ///< Corners of rectangle

    // aligned pixels
    double theta, prec; ///< Direction angle and precision
    vector<int> definedPixels;///< Pixels inside ROI having a direction
    vector<int> pixelCluster; ///< Map pixel to cluster index (logical 2D array)
    vector<Cluster> clusters; ///< Clusters detected inside the ROI

    // conversion from 2D point to 1D index
    int pixelToIndex(int x, int y) const {
        return (x-xMin) + (y-yMin)*width;
    }
    int pixelToIndex(const point& p) const {
        return pixelToIndex(p.x,p.y);
    }
    /// Size of bounding box
    int sizeBB() const { return width*(yMax-yMin+1); }
    /// Test position inside bounding box
    bool inBB(int x, int y) const {
        return (xMin<=x && x<=xMax && yMin<=y && y<=yMax);
    }
    bool inBB(const point& p) const {
        return inBB(p.x,p.y);
    }
    
    void computeRectangle(const Segment& rawSegment);
    void findAlignedPixels(vector<point>& alignedPixels, Cimage used);
    vector<point> findCC(point seed, int idx);
    void aggregatePixels(const vector<point>& alignedPixels);

    set<int> findIntersect(const Cluster& c, bool postLSD) const;
    bool mergeCluster(Cluster& c, const set<int>& inter, Cluster& merged);
public:
    ROI(const Segment& seg,       Dimage a, Dimage m, double lNT, Cimage used);
    ROI(const vector<Cluster>& c, Dimage a, Dimage m, double lNT);

    bool isVoid() const { return clusters.empty(); } ///< No cluster found.
    void setUsed(Cimage& used);
    void mergeClusters(bool postLSD);
    bool filterClusters(vector<Cluster>& filtered, Cimage& used,
                        double log_eps);
};

const int ROI::CLUSTER_NULL=-1;

/// Mark pixels as used in global image. This prevents multiple detections of
/// the same cluster.
void ROI::setUsed(Cimage& used) {
    for(size_t i=0; i<definedPixels.size(); i++)
        used->data[definedPixels[i]] = USED;
}

/// First constructor: a single segment (from previous scale).
/// Called in #refineRawSegments for merging clusters inside a rectangle.
ROI::ROI(const Segment& rawSegment, Dimage a, Dimage m, double lNT, Cimage used)
: angles(a), modgrad(m), logNT(lNT) {
    // gradient parameters
    theta = rawSegment.angle();
    prec = M_PI*rawSegment.prec;

    computeRectangle(rawSegment);
    pixelCluster = vector<int>(sizeBB(), NOTDEF);
    vector<point> alignedPixels;
    findAlignedPixels(alignedPixels, used);
    aggregatePixels(alignedPixels);
}

/// Second constructor: the ROI is the full image, clusters already computed.
/// Used in #mergeClusters.
ROI::ROI(const vector<Cluster>& c, Dimage a, Dimage m, double lNT)
: angles(a), modgrad(m), logNT(lNT) {
    width = (int)a->xsize;
    xMin = yMin = 0;
    xMax = a->xsize-1;
    yMax = a->ysize-1;
    clusters = c;
    pixelCluster = vector<int>(sizeBB(), NOTDEF);
    for(size_t i=0; i<clusters.size(); i++)
        for(size_t j=0; j<clusters[i].getPixels().size(); j++) {
            int idx = pixelToIndex( clusters[i].getPixels()[j] );
            pixelCluster[idx] = clusters[i].getIndex();
        }
}

/// Compute rectangle and bounding box for segment \a seg.
void ROI::computeRectangle(const Segment& seg) {
    // Compute rectangle corners
    Point2d P(seg.x1-1,seg.y1-1), Q(seg.x2-1,seg.y2-1); // -1: cluster->segment involved +0.5, but at scale 1/2.
    Point2d delta = (seg.width/2.0)/seg.length() * orthogonal(Q-P);
    p1 = P+delta;
    p2 = P-delta;
    q1 = Q+delta;
    q2 = Q-delta;

    // compute min/max along x/y axis
    xMin = max((int)floor(min({p1.x, p2.x, q1.x, q2.x})), 0);
    xMax = min((int)ceil (max({p1.x, p2.x, q1.x, q2.x})), (int)angles->xsize-1);
    yMin = max((int)floor(min({p1.y, p2.y, q1.y, q2.y})), 0);
    yMax = min((int)ceil (max({p1.y, p2.y, q1.y, q2.y})), (int)angles->ysize-1);
    width = xMax-xMin+1;
}

/// Insert in \a alignedPixels the pixels inside the ROI whose gradient
/// direction is aligned. Pixels inside the ROI with non-null gradient
/// are registered in field \a definedPixels.
/// As different ROI may overlap, some pixels may have been already used in a
/// previous ROI (checked with image \a used); they are then skipped here.
void ROI::findAlignedPixels(vector<point>& alignedPixels, Cimage used) {
    for(point p = {xMin,yMin}; p.y<=yMax; p.y++)
        for(p.x=xMin; p.x<=xMax; p.x++)
            if(insideRect(p1, p2, q1, q2, p)) {
                int i = p.x + p.y*angles->xsize;
                if(used->data[i] == NOTUSED && angles->data[i] != NOTDEF) {
                    definedPixels.push_back(i);
                    if(angle_diff(angles->data[i], theta) < prec) {
                        alignedPixels.push_back(p);
                        pixelCluster[pixelToIndex(p)] = CLUSTER_NULL;
                    }
                }
            }
}

/// Find 8-connected component of pixels aligned with the ROI direction.
/// The pixels are marked in array #pixelCluster with index \a idx. They
/// were initialized previously with label \c CLUSTER_NULL.
vector<point> ROI::findCC(point seed, int idx) {
    pixelCluster[pixelToIndex(seed)] = idx;
    vector<point> data(1, seed);
    for(size_t i=0; i<data.size(); ++i) {
        seed = data[i];
        // look inside 8-neighbourhood
        for (int dx=-1; dx<=1; dx++)
            for (int dy=-1; dy<=1; dy++) {
                if(dx==0 && dy==0) continue;
                point p = {seed.x + dx, seed.y + dy};
                int j = pixelToIndex(p);
                if(inBB(p) && pixelCluster[j]==CLUSTER_NULL) {
                    pixelCluster[j] = idx;
                    data.push_back(p);
                }
            }
    }
    return data;
}

/// Find 8-connected components of \a alignedPixels, i.e., the initial clusters.
/// The ones with fewer than \c MIN_SIZE_CLUSTER pixels are just ignored.
void ROI::aggregatePixels(const vector<point>& alignedPixels) {
    for(size_t i=0; i<alignedPixels.size(); i++) {
        int index = pixelToIndex(alignedPixels[i]);
        if(pixelCluster[index] != CLUSTER_NULL)
            continue;
        vector<point> data = findCC(alignedPixels[i], clusters.size());
        // suppress clusters with too few pixels: marked with
        // NOTDEF!=CLUSTER_NULL so that the same CC is not extracted later
        if(data.size() < MIN_SIZE_CLUSTER)
            for(size_t j=0; j<data.size(); j++)
                pixelCluster[pixelToIndex(data[j])] = NOTDEF;
        else
            clusters.push_back(Cluster(angles, modgrad, logNT, data, theta,
                                       prec, clusters.size()));
    }
}

/// Return indexes of clusters inside the ROI that meet the supporting line of
/// cluster \a c. The search includes both sides of the segment. If \a postLSD,
/// the search in each direction stops at the closest cluster. It is kept only
/// if its direction is roughly aligned and its NFA is worse. This last
/// criterion allows more significant clusters to absorb lesser ones, but not
/// the contrary.
set<int> ROI::findIntersect(const Cluster& c, bool postLSD) const {
    set<int> inter;
    const Point2d center = getCenter(c);
    const Point2d step = getSlope(c);
    const double theta = c.getTheta();
    const double len  = c.length();
    const int cidx = c.getIndex();
    const double prec = 1.0/len;

    for(int s=-1; s<=1; s += 2) { // Both directions
        const Point2d signed_step = s*step;
        for(Point2d p = center+0.5*len*signed_step; true; p+=signed_step) {
            int X = floor(p.x + 0.5);
            int Y = floor(p.y + 0.5);
            if(! inBB(X,Y)) break;

            int idx = pixelCluster[pixelToIndex(X,Y)];
            if(idx!=NOTDEF && idx!=cidx && !clusters[idx].isMerged()) {
                if(postLSD &&
                   (angle_diff(clusters[idx].getTheta(),theta)>prec ||
                    c.getNFA() < clusters[idx].getNFA()))
                    break;
                inter.insert(idx);
                if(postLSD) break;
            }
        }
    }
    return inter;
}

/// Try to merge cluster \a c with clusters of index in \a inter.
/// Two conditions must be satisfed:
/// 1. The NFA of the merged cluster must be better than the one of \a c.
/// 2. It must be better than the one of each consistuent in \a inter.
/// If successful (return value \c true), \a c and the clusters referred in
/// \a inter are marked as merged and the pixel positions get the new index.
bool ROI::mergeCluster(Cluster& c, const set<int>& inter, Cluster& merged) {
    set<int>::const_iterator it=inter.begin(), end=inter.end();
    merged = c.united(clusters, inter, angles, modgrad, logNT);
    if(merged.getNFA() <= c.getNFA())
        return false;
    for(it=inter.begin(); it!=end; ++it)
        if(merged.getNFA()<=clusters[*it].getNFA())
            return false;

    // labelize merged clusters as merged
    c.setMerged();
    for(it=inter.begin(); it!=end; ++it)
        clusters[*it].setMerged();
    // labelize clustered points with their new label
    for (size_t j=0; j<merged.getPixels().size(); j++) {
        const point& p = merged.getPixels()[j];
        if(inBB(p)) // due to width reduction, there can be some issue
            pixelCluster[pixelToIndex(p)] = merged.getIndex();
    }
    return true;
}

/// Merged clusters inside the ROI. This is done by NFA priority: the most
/// significant one finds the ones it can merge with (hence of worse NFA).
/// Notice that if the merge is successful, the new cluster has the best NFA
/// and is candidate for further merging.
/// If \a postLSD, merge is attempted with a single other cluster at once.
void ROI::mergeClusters(bool postLSD) {
    // sort clusters by decreasing NFA
    CompareClusters cmp(clusters);
    vector<int> v(clusters.size());
    std::iota(v.begin(), v.end(), 0);
    priority_queue<int,vector<int>,CompareClusters>
        clusterQ(v.begin(), v.end(), cmp);
    while(! clusterQ.empty()) {
        int currIndex = clusterQ.top();
        clusterQ.pop();
        Cluster& c = clusters[currIndex];
        if( c.isMerged() ) continue;
        set<int> intersect = findIntersect(c, postLSD);
        if( intersect.empty() ) continue;

        Cluster megaCluster;
        if(postLSD) {
            set<int>::const_iterator it=intersect.begin(), end=intersect.end();
            for(; it!=end; ++it) {
                set<int> temp; temp.insert(*it);
                if(mergeCluster(c, temp, megaCluster))
                    break;
            }
        } else
            mergeCluster(c, intersect, megaCluster);

        if( c.isMerged() ) {
            clusters.push_back(megaCluster);
            clusterQ.push(clusters.size()-1);
        }
    }
}

/// Append valid clusters of the ROI to \a filtered.
/// A valid cluster has not been merged and has a sufficiently large -log(NFA).
bool ROI::filterClusters(vector<Cluster>& filtered, Cimage& used,
                         double log_eps) {
    bool adding=false;
    for (size_t i=0; i<clusters.size(); i++) {
        if(clusters[i].isMerged() || clusters[i].getNFA() <= log_eps)
            continue;
        adding = true;
        clusters[i].setUsed(used);
        clusters[i].setIndex(filtered.size());
        filtered.push_back(clusters[i]);
    }
    return adding;
}

// --- Free functions ---

/// Refine \a rawSegments from previous image scale. The output is a collection
/// of clusters. The rectangle of each segment is upscaled and connected
/// components of pixels within at current scale with compatible direction
/// may me merged.
vector<Cluster> refineRawSegments(const vector<Segment>& rawSegments,
                                  Dimage angles, Dimage modgrad,
                                  Cimage& used,
                                  double logNT, double log_eps) {
    vector<Cluster> clusters;
    for(size_t i=0; i<rawSegments.size(); i++) {
        Segment seg = rawSegments[i].upscaled();
        ROI roi(seg, angles, modgrad, logNT, used); // Constructor with segment
        if (roi.isVoid())
            continue;
        roi.mergeClusters(false);
        if(! roi.filterClusters(clusters, used, log_eps))
            roi.setUsed(used); // Tag pixels inside to prevent detection by LSD
    }
    return clusters;
}

/// Merge greedily aligned clusters that should be merged (in NFA meaning).
/// Some of the input \a clusters are merged, generating some larger clusters.
/// At output, the number of clusters is less or equal.
void mergeClusters(vector<Cluster>& clusters,
                   Dimage angles, Dimage modgrad, Cimage& used,
                   double logNT, double log_eps) {
    ROI roi(clusters, angles, modgrad, logNT); // ROI is the full image
    roi.mergeClusters(true);
    clusters.clear();
    roi.filterClusters(clusters, used, log_eps);
}
