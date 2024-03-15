// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file cluster.hpp
 * @brief muLSD: clusters (segment candidates)
 * @author Yohann Salaun <yohann.salaun@imagine.enpc.fr>
 *         Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2016, 2023-2024
 */

#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include "segment.hpp"
extern "C" {
#include "lsd.h"
}
#include <set>
#include <vector>

typedef image_char   Cimage;
typedef image_double Dimage;

/** \brief Cluster=segment candidate. */
/// A cluster stores information about a segment according to LSD: a rectangular
/// section of image with pixels of coinciding gradient direction. This is
/// associated with an NFA. Initially, the cluster is a connected set of pixels,
/// issued from LSD. Afterwards, clusters can be merged. Finally, it is
/// converted to a Segment.
class Cluster {
    std::vector<point> pixels; ///< Pixels with coinciding direction
    rect rec;                  ///< Rectangle (issued from LSD)
    double nfa;                ///< -log_10(NFA)
    int index;                 ///< Unique identifier for the cluster
    bool merged;               ///< Merged clusters are ignored

public:
    Cluster() {} ///< The only sane behavior with this is to use operator= after
    Cluster(Dimage angles, Dimage modgrad, double logNT,
            const std::vector<point>& d, double t, double p, int idx);
    Cluster(Dimage angles, double logNT,
            const point* d, int dsize, rect& r, int idx);

    Cluster united(const std::vector<Cluster>& clusters,
                   const std::set<int>& indexToMerge,
                   Dimage angles, Dimage modgrad,
                   double logNT) const;

    double length() const;
    Segment toSegment() const;
    void setUsed(CImage& used) const;

    const std::vector<point>& getPixels() const { return pixels; }
    const rect& rectangle() const { return rec; }
    double getTheta() const { return rec.theta; }
    double getNFA() const { return nfa; }
    int getIndex() const { return index; }
    void setIndex(int i) { index=i; }
    bool isMerged() const { return merged; }
    void setMerged() { merged=true; }
};

/// Compute the segments at current scale with information from previous scale.
std::vector<Cluster> refineRawSegments(const std::vector<Segment>& rawSegments,
                                       Dimage angles,
                                       Dimage modgrad, Cimage& used,
                                       double logNT, double log_eps);

/// Merge clusters at same scale that belong to the same line.
void mergeClusters(std::vector<Cluster>& clusters,
                   Dimage angles, Dimage modgrad, Cimage& used,
                   double logNT, double log_eps);

#endif
