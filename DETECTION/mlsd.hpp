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
#include "interface.hpp"
#include "lsd.hpp"
#include <set>

// parameters for pyramid of images
const double prec = 3.0;
const double sigma_scale = 0.8;
const double scale_step = 2.0;
const int h_kernel = (unsigned int)ceil(sigma_scale * sqrt(2.0 * prec * log(10.0)));

struct Point2d{
    double x;
    double y;
    Point2d(double X,double Y):x(X),y(Y){}
    Point2d(){}
};

class NFA_params{
    double value, theta, prec_divided_by_pi, prec;
    bool computed;
public:
    NFA_params();
    NFA_params(const double t, const double p);
    NFA_params(const NFA_params &n);

    rect regionToRect(std::vector<point> &data, image_double modgrad);
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
    bool merged;
    int scale;

public:

    Cluster();
    Cluster(image_double angles, image_double modgrad, const double logNT,
            std::vector<point> &d, const double t, const double p, const int i, const int s, const double n);
    Cluster(image_double angles, image_double modgrad, const double logNT,
            point* d, const int dsize, rect &r, const int i, const int s);

    Cluster mergedCluster(const std::vector<Cluster> &clusters, const std::set<int> &indexToMerge,
                          image_double angles, image_double modgrad, const double logNT) const;
    bool isToMerge(image_double angles, const double logNT);

    double length() const;
    Segment toSegment();

    bool isMerged() const;
    void setMerged();
    double getNFA() const;
    int getIndex() const;
    int getScale() const;
    double getTheta() const;
    double getPrec() const;
  Point2d getCenter() const;
    Point2d getSlope() const;
    const std::vector<point>* getData() const;
    void setUsed(image_char used) const;
    void setIndex(int i);
};

// generate the number of scales needed for a given image
int scaleNb(std::shared_ptr<Image>  im, const bool multiscale);

// filter some area in the image for better SfM results and a faster line detection
void denseGradientFilter(std::vector<int> &noisyTexture, std::shared_ptr<Image>  im, 
                         const image_double &angles, const image_char &used,
                         const int xsize, const int ysize, const int N);

// compute the segments at current scale with information from previous scale
std::vector<Cluster> refineRawSegments(const std::vector<Segment> &rawSegments, std::vector<Segment> &finalLines, const int i_scale,
                                  image_double angles, image_double modgrad, image_char used,
                                  const double logNT, const double log_eps);

// merge segments at same scale that belong to the same line
void mergeSegments(std::vector<Cluster> &refinedLines, const double segment_length_threshold, const int i_scale,
                   image_double angles, image_double modgrad, image_char &used,
                   const double logNT, const double log_eps);

#endif /* !MULTISCALE_LSD_HEADER */
/*----------------------------------------------------------------------------*/
