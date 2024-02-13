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

#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include "segment.hpp"
extern "C" {
#include "lsd.h"
}
#include <set>
#include <vector>

/// A Cluster stores information about a segment according to LSD: a rectangular
/// section of image with pixels of coinciding gradient direction. This is
/// associated with an NFA. Initially, the cluster is a connected set of pixels,
/// issued from LSD. Afterwards, clusters can be merged. Finally, it is
/// converted to a Segment.
class Cluster {
    std::vector<point> pixels; ///< Pixels with coinciding direction
    rect rec;                  ///< Rectangle (issued from LSD)
    double nfa;                ///< -log_10(NFA)
    int index;                 ///< Unique identifier for the cluster
    int scale;                 ///< Image scale of detection
    bool merged;               ///< Merged clusters are ignored

public:
    Cluster() {} ///< The only sane behavior with this is to use operator= after
    Cluster(image_double angles, image_double modgrad, double logNT,
            const std::vector<point>& d, double t, double p, int idx, int s);
    Cluster(image_double angles, double logNT,
            const point* d, int dsize, rect& r, int idx, int s);

    Cluster united(const std::vector<Cluster>& clusters,
                   const std::set<int>& indexToMerge,
                   image_double angles, image_double modgrad,
                   double logNT) const;

    double length() const;
    Segment toSegment() const;
    void setUsed(image_char& used) const;

    const std::vector<point>& getPixels() const { return pixels; }
    const rect& rectangle() const { return rec; }
    double getTheta() const { return rec.theta; }
    double getNFA() const { return nfa; }
    int getIndex() const { return index; }
    void setIndex(int i) { index=i; }
    int getScale() const { return scale; }
    bool isMerged() const { return merged; }
    void setMerged() { merged=true; }
};

/// Compute the segments at current scale with information from previous scale.
std::vector<Cluster> refineRawSegments(const std::vector<Segment>& rawSegments,
                                       std::vector<Segment>& finalLines,
                                       int i_scale,
                                       image_double angles,
                                       image_double modgrad, image_char& used,
                                       double logNT, double log_eps);

/// Merge clusters at same scale that belong to the same line.
void mergeClusters(std::vector<Cluster>& clusters,
                   double minLength, int i_scale,
                   image_double angles, image_double modgrad, image_char& used,
                   double logNT, double log_eps);

#endif
