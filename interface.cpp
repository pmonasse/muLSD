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

#include "interface.hpp"
#include <cmath>  

/// Constructor
Segment::Segment(double X1, double Y1, double X2, double Y2,
                 double w, double p, double nfa, double s) {
    x1 = X1; y1 = Y1;
    x2 = X2; y2 = Y2;
    width = w;
    length = std::sqrt(qlength());
    prec = p;
    log_nfa = nfa;
    angle = atan2(y2 - y1, x2 - x1);
    scale = s;
}

// For multiscale LSD
void Segment::upscale(double k){
    x1 *= k; y1 *= k;
    x2 *= k; y2 *= k;
    width *= k;
    length *= k;
}

// I/O METHODS for segments
std::istream& operator>>(std::istream &str, Segment& s){
    str >> s.x1 >> s.y1 >> s.x2 >> s.y2
        >> s.width >> s.prec >> s.log_nfa >> s.scale;
    s.scale = 0;
    s.angle = std::atan2(s.y2-s.y1, s.x2-s.x1);
    s.length = std::sqrt(s.qlength());
    return str;
}

std::ostream& operator<<(std::ostream &str, const Segment& s) {
    str << s.x1 << "  " << s.y1 << "  " << s.x2 << "  " << s.y2 << " "
        << s.width << " " << s.prec << " " << s.log_nfa << " " << s.scale << '\n';
    return str;
}

double Segment::qlength() {
    double dX = x1 - x2;
    double dY = y1 - y2;
    length = dX*dX + dY*dY;
    return length;
}
