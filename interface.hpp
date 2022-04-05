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

#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <iostream>

// strucuture to keep track of segment detected
struct Segment{
    // segment coordinates
    double x1, y1, x2, y2;

    // segment geometric attributes
    double width, length, angle;

    // NFA related arguments
    double log_nfa, prec;

    // scale of last detection from 0 (rawest) to n (finest)
    int scale;

    Segment(){}
    Segment(const double X1, const double Y1, const double X2, const double Y2,
            const double w, const double p, const double nfa, const double s);

    // For multiscale LSD
    void upscale(const double k);

    double qlength();
};

std::istream& operator>>(std::istream &str, Segment& s);
std::ostream& operator<<(std::ostream &str, const Segment& s);

#endif
