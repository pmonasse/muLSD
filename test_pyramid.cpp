// SPDX-License-Identifier: GPL-3.0-or-later
/**
 * @file test_pyramid.cpp
 * @brief Compute Gaussian pyramid of image
 * @author Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2023-2024
 */

#include "image.h"
#include "cmdLine.h"
#include <iostream>

int main(int argc, char* argv[]) {
    int nScales=3;
    std::string prefix = "pyramid";
    CmdLine cmd;
    cmd.add( make_option('n', nScales).doc("number of scales") );
    cmd.add( make_option('p', prefix).doc("prefix for output image files") );
    try {
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << "Error: " << s << std::endl;
        return 1;
    }
    if(argc!=2) {
        std::cerr << "Usage: " << argv[0] << " image.png\nOptions:\n" << cmd;
        return 1;
    }
    Image<float> im(argv[1]);
    if(!im.data) {
        std::cerr << "Error reading image file " << argv[1] << std::endl;
        return 1;
    }
    std::vector<Image<float>*> P = gaussPyramid(im, nScales);
    for(int i=0; i<nScales; i++) {
        std::string fileName = prefix + std::to_string(i) + ".png";
        io_png_write_f32(fileName.c_str(), P[i]->data, P[i]->w, P[i]->h, 1);
        delete P[i];
        std::cout << "Writing image " << fileName << std::endl;
    }

    return 0;
}
