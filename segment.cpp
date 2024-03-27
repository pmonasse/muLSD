// SPDX-License-Identifier: MPL-2.0
/**
 * @file segment.cpp
 * @brief muLSD: segment
 * @author Yohann Salaun <yohann.salaun@imagine.enpc.fr>
 * @date 2016, 2023-2024
 */

#include "segment.hpp"
#include <cmath>  

/// Constructor.
Segment::Segment(double X1, double Y1, double X2, double Y2,
                 double w, double p, double nfa) {
    x1 = X1; y1 = Y1;
    x2 = X2; y2 = Y2;
    width = w;
    prec = p;
    log_nfa = nfa;
}

/// Upscale by a factor 2.
Segment Segment::upscaled() const {
    Segment s = *this;
    s.x1 *= 2; s.y1 *= 2;
    s.x2 *= 2; s.y2 *= 2;
    s.width *= 2;
    return s;
}

/// Length (in pixels) of the segment
double Segment::length() const {
    return hypot(x1-x2,y1-y2);
}

/// Angle of segment with horizontal.
double Segment::angle() const {
    return atan2(y2-y1, x2-x1);
}

/// Output segment.
std::ostream& operator<<(std::ostream& str, const Segment& s) {
    str << s.x1 << "  " << s.y1 << "  " << s.x2 << "  " << s.y2 << " "
        << s.width << " " << s.prec << " " << s.log_nfa << '\n';
    return str;
}
