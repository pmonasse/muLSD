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
#include <cmath>

struct Point2d{
    double x;
    double y;
    Point2d(double X,double Y):x(X),y(Y){}
    Point2d(){}
};

class NFA_params{
    double value, theta, prec;
    bool computed;
public:
    NFA_params();
    NFA_params(const double t, const double p);

    void computeNFA(rect &rec, image_double angles, const double logNT);

    double getValue() const;
    double getTheta() const;
    double getPrec() const;
};

class Cluster{
    // pixels parameters
    std::vector<point> data;
    rect rec;

    // NFA parameters
    NFA_params nfa;
    double nfa_separated_clusters;

    // fusion parameters
    int index;
    int scale;
    bool merged;

public:
    Cluster();
    Cluster(image_double angles, image_double modgrad, double logNT,
            const std::vector<point>& d, double t, double p, int i, int s,
            double nfaSeparate);
    Cluster(image_double angles, double logNT,
            const point* d, const int dsize, rect &r, int i, int s);

    Cluster mergedCluster(const std::vector<Cluster>& clusters, const std::set<int>& indexToMerge,
                          image_double angles, image_double modgrad, double logNT) const;
    bool isToMerge(image_double angles, double logNT);

    double length() const;
    Segment toSegment() const;

    bool isMerged() const;
    void setMerged();
    double getNFA() const;
    int getIndex() const;
    int getScale() const;
    double getTheta() const;
    double getPrec() const;
    Point2d getCenter() const;
    Point2d getSlope() const;
    const std::vector<point>& getData() const;
    void setUsed(image_char used) const;
    void setIndex(int i);
};

// filter some area in the image for better SfM results and a faster line detection
void denseGradientFilter(std::vector<int> &noisyTexture, int w, int h, 
                         const image_double &angles, const image_char &used,
                         int xsize, int ysize, int N);

// compute the segments at current scale with information from previous scale
std::vector<Cluster> refineRawSegments(const std::vector<Segment>& rawSegments,
                                       std::vector<Segment> &finalLines,
                                       int i_scale,
                                       image_double angles,
                                       image_double modgrad, image_char used,
                                       double logNT, double log_eps);

// merge segments at same scale that belong to the same line
void mergeSegments(std::vector<Cluster>& refinedLines,
                   double segment_length_threshold, int i_scale,
                   image_double angles, image_double modgrad, image_char& used,
                   double logNT, double log_eps);

#endif /* !MULTISCALE_LSD_HEADER */
/*----------------------------------------------------------------------------*/
