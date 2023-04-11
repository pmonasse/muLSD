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

#ifndef MLSD_HPP
#define MLSD_HPP

#include "segment.hpp"
extern "C" {
#include "lsd.h"
}
#include <set>
#include <vector>

class Cluster {
    std::vector<point> pixels;
    rect rec;

    // NFA parameters
    double theta, prec, nfa;

    // fusion parameters
    int index;
    int scale;
    bool merged;

public:
    Cluster() {}
    Cluster(image_double angles, image_double modgrad, double logNT,
            const std::vector<point>& d, double t, double p, int i, int s);
    Cluster(image_double angles, double logNT,
            const point* d, const int dsize, rect &r, int i, int s);

    Cluster united(const std::vector<Cluster>& clusters,
                   const std::set<int>& indexToMerge,
                   image_double angles, image_double modgrad,
                   double logNT) const;

    double length() const;
    Segment toSegment() const;
    void setUsed(image_char& used) const;

    const std::vector<point>& getPixels() const { return pixels; }
    const rect& rectangle() const { return rec; }
    double getTheta() const { return theta; }
    double getPrec() const { return prec; }
    double getNFA() const { return nfa; }
    int getIndex() const { return index; }
    void setIndex(int i) { index=i; }
    int getScale() const { return scale; }
    bool isMerged() const { return merged; }
    void setMerged() { merged=true; }
};

/// Compute the segments at current scale with information from previous scale
std::vector<Cluster> refineRawSegments(const std::vector<Segment>& rawSegments,
                                       std::vector<Segment>& finalLines,
                                       int i_scale,
                                       image_double angles,
                                       image_double modgrad, image_char& used,
                                       double logNT, double log_eps);

/// Merge clusters at same scale that belong to the same line
void mergeClusters(std::vector<Cluster>& clusters,
                   double minLength, int i_scale,
                   image_double angles, image_double modgrad, image_char& used,
                   double logNT, double log_eps);

// filter some area in image for better SfM results and faster line detection
void denseGradientFilter(std::vector<int>& noisyTexture, int w, int h, 
                         image_double angles, image_char& used,
                         int xsize, int ysize, int N);

#endif
