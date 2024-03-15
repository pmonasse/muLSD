// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file mulsd.hpp
 * @brief muLSD: main interface
 * @author Yohann Salaun <yohann.salaun@imagine.enpc.fr>
 *         Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2016, 2023-2024
 */

#ifndef MLSD_HPP
#define MLSD_HPP

#include "segment.hpp"
#include "image.h"
#include <vector>

// muLSD interface: detect segments from image pyramid.
// \param imgPyr pyramid of scaled images.
// \param grad minimum gradient magniture (0=automatic).
std::vector<Segment> lsd_multiscale(const std::vector<Image<float>*>& imgPyr,
                                    float grad);

#endif
