// SPDX-License-Identifier: MPL-2.0
/**
 * @file segment.hpp
 * @brief muLSD: segment
 * @author Yohann Salaun <yohann.salaun@imagine.enpc.fr>
 * @date 2016, 2023-2024
 */

#ifndef SEGMENT_HPP
#define SEGMENT_HPP

#include <iostream>

/// Strucuture to keep track of detected segment.
struct Segment{
    // segment coordinates
    double x1, y1, x2, y2;
    double width;

    // NFA related arguments
    double log_nfa, prec;

    Segment(){}
    Segment(double X1, double Y1, double X2, double Y2,
            double w, double p, double nfa);

    // For multiscale LSD
    Segment upscaled() const;

    double length() const;
    double angle() const;
};

std::ostream& operator<<(std::ostream& str, const Segment& s);

#endif
