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

#include "mlsd.hpp"
#include "cmdLine.h"
#include <fstream>

using namespace std;

/// Adapt number of scales to the resolution of the picture.
/// With pictures of less than 1Mpx, the original LSD detector performs well.
int nbScales(int w, int h){
    int n = 1;
    int maxWH = max(w,h);
    while(maxWH > 500){
        maxWH /= 2;
        n++;
    }
    return n;
}

void saveLines(const std::vector<Segment> &lines, const char* name) {
    ofstream linesTxt(name, ofstream::out);
    linesTxt << lines.size() << endl;
    for (size_t i=0; i<lines.size(); i++)
        linesTxt << lines[i];
}

int main(int argc, char* argv[]) { 
    // parse arguments
    CmdLine cmd;

    int nScales=0;
    float grad=0;
    
    // options
    cmd.add( make_option('s', nScales, "scales")
             .doc("nb scales (0=automatic)") );
    cmd.add( make_option('g', grad, "gradient")
             .doc("Min gradient norm (0=automatic)") );

    try {
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << "Error: " << s << std::endl;
        return 1;
    }
    if(argc != 3) {
        cerr << "Usage: " << argv[0] << " [options] imgIn.png out.txt\n" << cmd;
        return 1;
    }

    Image<float> im(argv[1]);
    if(! im.data) {
        cerr << "Unable to load image " << argv[1] << endl;
        return 1;
    }

    clock_t t0 = clock();

    if(nScales==0)
        nScales = nbScales(im.w, im.h);
    vector<Image<float>*> imagePyramid = gaussPyramid(im, nScales);

    vector<Segment> segments = lsd_multiscale(imagePyramid, grad);
    saveLines(segments, argv[2]);

    cout << "Runtime: " << (clock()-t0)/float(CLOCKS_PER_SEC) << endl;

    for(size_t i=0; i<imagePyramid.size(); i++)
        delete imagePyramid[i];

    return 0;
}
