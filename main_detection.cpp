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
#include "cmdLine/cmdLine.h"

using namespace std;

int main(int argc, char* argv[]){ 
    // Seed random function
    srand((unsigned int)(time(0)));

    // parse arguments
    CmdLine cmd;

    string dirPath;
    string picList;

    string name_image;

    string str="png";

    string picPath;

    bool consecutive = true;
    bool withDetection = false;
    double segment_length_threshold=0;//0.0009 ;
    bool multiscale = true;

    // required
    cmd.add( make_option('d', dirPath, "dirPath") );
    cmd.add( make_option('i', name_image, "inputPic") );
    
    // optional
    cmd.add( make_option('m', multiscale, "multiscale") );
    cmd.add( make_option('t', segment_length_threshold, "threshold") );

    try {
        if (argc == 1) throw std::string("Invalid command line parameter.");
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << "Usage: " << argv[0] << '\n'
                  << "[-d|--dirPath] feature path]\n"
                  << "[-i|--name_image] onepicture] \n"
                  << "\n[Optional]\n"
                  << "[-m|--multiscale] multiscale option (default = TRUE)\n"
                  << "[-t|--threshold] threshold for segment length (default = 0.05% of image size)\n"
                  << endl;
        return EXIT_FAILURE;
    }
    dirPath += "/";


    const string ext = (multiscale)? "" : "_no_mlsd";
    bool exist=name_image.find(str) != string::npos;
    if(exist==true){
        clock_t processing_time = clock();
        const char* parsName = name_image.c_str();
        std::shared_ptr<Image> im = Charger(parsName);

        //std::cout<<"parsName="<<parsName<<std::endl;


        std::shared_ptr<Image> ImGray=rgbtogray(im);
        //int k=Sauver(ImGray, "Image22.png");
        vector<std::shared_ptr<Image>> imagePyramid = computeImagePyramid(ImGray, multiscale);

        vector<Segment> segments = lsd_multiscale(imagePyramid, segment_length_threshold, multiscale);
        saveLines(segments, dirPath,name_image+ext);

        cout << "PROCESSED IN " << (clock() - processing_time) / float(CLOCKS_PER_SEC) << endl;

        /*Image* im1 = Charger(parsName);
    vector<Segment> segments1 = readLines(dirPath, name_image+ext);
    saveLinesPicture(segments1, *im1, dirPath, name_image+ext, false);*/
    }
    else std::cout<<"give a PNG image"<<std::endl;

    return 0;
}
