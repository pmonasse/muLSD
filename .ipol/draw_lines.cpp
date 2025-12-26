// SPDX-License-Identifier: MPL-2.0
/**
 * @file draw_lines.cpp
 * @brief Draw line segments with own random color on transparent image.
 * @author Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2025
 */

#include "io_png.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

struct Seg {
    float x1, y1, x2, y2;
};

std::vector<Seg> read_seg(const char* file) {
    std::vector<Seg> V;
    std::ifstream F(file);
    if(! F.is_open()) {
        std::cerr << "Error reading file of lines " << file << std::endl;
        return V;
    }
    int ignored=0;
    while(! F.eof()) {
        std::string str;
        std::getline(F, str);
        std::istringstream is(str);
        Seg s;
        is >> s.x1 >> s.y1 >> s.x2 >> s.y2;
        if(is.fail()) {
            if(!str.empty() && // Don't report blank lines
               str.find_first_not_of(" \t\n\v\f\r")!=std::string::npos)
                ++ignored;
        } else
            V.push_back(s);
    }
    if(ignored)
        std::cerr << "Ignored " << ignored
                  << " lines in file " << file << std::endl;
    return V;
}

void find_bb(const std::vector<Seg>& V, int& w, int& h) {
    float wf=0, hf=0;
    for(std::vector<Seg>::const_iterator i=V.begin(); i!=V.end(); ++i) {
        if(wf<i->x1) wf=i->x1;
        if(wf<i->x2) wf=i->x2;
        if(hf<i->y1) hf=i->y1;
        if(hf<i->y2) hf=i->y2;
    }
    w = (int)std::ceil(wf)+1;
    h = (int)std::ceil(hf)+1;
}

inline unsigned char randu8() { return (unsigned char)(std::rand()%RAND_MAX); }

void draw_line(const Seg& s, unsigned char* im, int w, int) {
    unsigned char col[4]={randu8(), randu8(), randu8(), 255};
    int p[] = {(int)std::round(s.x1), (int)std::round(s.y1)};
    int q[] = {(int)std::round(s.x2), (int)std::round(s.y2)};
    int delta[] = {q[0]-p[0], q[1]-p[1]};
    int d[]={1,1};
    if(delta[0]<0) { delta[0]=-delta[0]; d[0]=-1; }
    if(delta[1]<0) { delta[1]=-delta[1]; d[1]=-1; }
    int o = (delta[0]>=delta[1])? 0: 1;
    int plus=2*delta[1-o], sub=plus-delta[o], minus=sub-delta[o];
    while(p[o]!=q[o]) {
        for(int j=0; j<4; j++)
            im[4*(p[0]+p[1]*w)+j] = col[j];
        p[o] += d[o];
        if(sub<=0)
            sub += plus;
        else {
            sub += minus;
            p[1-o] += d[1-o];
        }
    }
    for(int j=0; j<4; j++)
        im[4*(p[0]+p[1]*w)+j] = col[j];
}

int main(int argc, char* argv[]) {
    if(argc!=3) {
        std::cerr << "Usage: " << argv[0] << " lines.txt img.png" << std::endl;
        return 1;
    }
    std::vector<Seg> V = read_seg(argv[1]);
    int w, h;
    find_bb(V, w,h);
    unsigned char* im= new unsigned char[4*w*h];
    for(int i=4*w*h-1; i>=0; i--)
        im[i]=0;
    for(std::vector<Seg>::const_iterator i=V.begin(); i!=V.end(); ++i)
        draw_line(*i, im, w, h);
    if(io_png_write_u8(argv[2], im, w, h, 4)) {
        std::cerr << "Error writing image file " << argv[2] << std::endl;
        return 1;
    }
    delete [] im;
    return 0;
}
