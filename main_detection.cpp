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

#include "detection.hpp"
#include "interface.hpp"
#include "cmdLine.h"

using namespace std;

void saveLines(const std::vector<Segment> &lines, const char* name) {
    ofstream linesTxt(name, ofstream::out);
    linesTxt << lines.size() << endl;
    for (size_t i=0; i<lines.size(); i++)
        lines[i].saveSegment(linesTxt);
}

int main(int argc, char* argv[]) { 
    // parse arguments
    CmdLine cmd;

    string picPath;

    double segment_length_threshold=0;//0.0009 ;
    bool multiscale = true;

    // options
    cmd.add( make_option('m', multiscale, "multiscale")
             .doc("multiscale option") );
    cmd.add( make_option('t', segment_length_threshold, "threshold")
             .doc("threshold for segment length") );

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

    std::shared_ptr<Image> im = Charger(argv[1]);
    if(! im) {
        cerr << "Unable to load image " << argv[1] << endl;
        return 1;
    }

    clock_t t0 = clock();

    std::shared_ptr<Image> ImGray=rgbtogray(im);
    vector<std::shared_ptr<Image>> imagePyramid =
        computeImagePyramid(ImGray, multiscale);

    vector<Segment> segments =
        lsd_multiscale(imagePyramid, segment_length_threshold, multiscale);
    saveLines(segments, argv[2]);

    cout << "Runtime: " << (clock()-t0)/float(CLOCKS_PER_SEC) << endl;

    return 0;
}
