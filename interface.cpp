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
#include <math.h>  
#include "interface.hpp"
using namespace std;


// for color randomization
/*vector<Scalar> rgb = {Scalar(255,0,0), Scalar(0,255,0), Scalar(0,0,255), Scalar(255,255,0), Scalar(0,255,255), Scalar(255,0,255), 
                      Scalar(0,150,0), Scalar(150,0,0), Scalar(0,0,150), Scalar(150,150,0), Scalar(0,150,150), Scalar(150,0,150)};
int countRandomRGB = 0;
*/
/*=================== SEGMENT ===================*/
Segment::Segment(const double X1, const double Y1, const double X2, const double Y2,
                 const double w, const double p, const double nfa, const double s){
    x1 = X1; y1 = Y1;
    x2 = X2; y2 = Y2;
    width = w;
    length = std::sqrt(qlength());
    prec = p;
    log_nfa = nfa;
    angle = atan2(y2 - y1, x2 - x1);
    scale = s;
}

// FOR MULTISCALE LSD
void Segment::upscale(const double k){
    x1 *= k; y1 *= k;
    x2 *= k; y2 *= k;
    width *= k;
    length *= k;
}

// I/O METHODS for segments
void Segment::readSegment(std::ifstream &file){
    file >> x1 >> y1 >> x2 >> y2 >> width >> prec >> log_nfa >> scale;
    scale = 0;
    angle = atan2(y2 - y1, x2 - x1);
    length = sqrt(qlength());
}

void Segment::saveSegment(std::ofstream &file) const{
    file << x1 << "  " << y1 << "  " << x2 << "  " << y2 << " " << width << " " << prec << " " << log_nfa << " " << scale << std::endl;
}

double Segment::qlength(){
    double dX = x1 - x2;
    double dY = y1 - y2;
    length = dX*dX + dY*dY;
    return length;
}
