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

#include "segment.hpp"
#include <cmath>  

/// Constructor
Segment::Segment(double X1, double Y1, double X2, double Y2,
                 double w, double p, double nfa, double s) {
    x1 = X1; y1 = Y1;
    x2 = X2; y2 = Y2;
    width = w;
    prec = p;
    log_nfa = nfa;
    scale = s;
}

// For multiscale LSD
Segment Segment::upscaled() const {
    Segment s = *this;
    s.x1 *= 2; s.y1 *= 2;
    s.x2 *= 2; s.y2 *= 2;
    s.width *= 2;
    return s;
}

double Segment::length() const {
    return hypot(x1-x2,y1-y2);
}

double Segment::angle() const {
    return atan2(y2-y1, x2-x1);
}

std::ostream& operator<<(std::ostream& str, const Segment& s) {
    str << s.x1 << "  " << s.y1 << "  " << s.x2 << "  " << s.y2 << " "
        << s.width << " " << s.prec << " " << s.log_nfa << " " << s.scale << '\n';
    return str;
}
