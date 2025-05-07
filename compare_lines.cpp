// SPDX-License-Identifier: MPL-2.0
/**
 * @file compare_lines.cpp
 * @brief Metric evaluation of the detections
 * @author Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2025
 */

#include "cmdLine.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

struct Seg {
    float x1, y1, x2, y2;
    float cx, cy, l, angle;
};

struct Param {
    float dist;    ///< Max orthogonal distance
    float angle;   ///< Max angle difference
    float overlap; /// < Min overlap ratio
    Param(): dist(1), angle(5), overlap(0.5) {}
};

std::istream& operator>>(std::istream& str, Seg& s) {
    str >> s.x1 >> s.y1 >> s.x2 >> s.y2;
    s.cx = 0.5f*(s.x1+s.x2);
    s.cy = 0.5f*(s.y1+s.y2);
    float dx=s.x1-s.x2, dy=s.y1-s.y2;
    s.l = std::hypot(dx,dy);
    s.angle = std::atan2(dy,dx);
    return str;
}

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
        is >> s;
        if(is.fail() && !F.eof()) {
            ++ignored;
        } else
            V.push_back(s);
    }
    if(ignored)
        std::cerr << "Ignored " << ignored
                  << " lines in file " << file << std::endl;
    return V;
}

float mod_pi(float angle) {
    if(angle<0)
        angle += M_PI;
    else if(angle>M_PI)
        angle -= M_PI;
    return angle;
}

float perp_dist(const Seg& s1, const Seg& s2, float c, float s) {
    return std::abs(-s*(s1.cx-s2.cx)+c*(s1.cy-s2.cy));
}

float angle_diff(float angle1, float angle2) {
    angle1 = mod_pi(angle1), angle2=mod_pi(angle2);
    float d = std::abs(angle1-angle2);
    if(d>M_PI/2) d-=M_PI/2;
    return d;
}

float intersect(const Seg& sGT, const Seg& sD) {
    float dx=sGT.x2-sGT.x1, dy=sGT.y2-sGT.y1;
    float d=hypot(dx,dy);
    dx/=d; dy/=d;
    float x0=0, x1=d;
    float x2=(sD.x1-sGT.x1)*dx+(sD.y1-sGT.y1)*dy;
    float x3=(sD.x2-sGT.x1)*dx+(sD.y2-sGT.y1)*dy;
    if(x2>x3) std::swap(x2,x3);
    return std::max(0.0f,std::min(x1,x3)-std::max(x0,x2));
}

bool eval_seg(const Seg& sGT, const std::vector<Seg>& D, const Param& param,
              float& tp, float& fp, float& tpIOU, float& fpIOU) {
    const float c=std::cos(sGT.angle), s=std::sin(sGT.angle);
    bool hit=false;
    for(std::vector<Seg>::const_iterator it=D.begin(); it!=D.end(); ++it) {
        if(! (perp_dist(sGT,*it,c,s)          <= param.dist &&
              angle_diff(sGT.angle,it->angle) <= param.angle))
            continue;
        float inter=intersect(sGT,*it);
        if(inter > 0) {
            hit=true;
            if(inter >= param.overlap*sGT.l)
                tp += it->l;
            else
                fp += it->l;
            tpIOU += inter;
            fpIOU += sGT.l+it->l-inter;
        }
    }
    return hit;
}

int main(int argc, char* argv[]) {
    Param param;
    CmdLine cmd;
    cmd.add( make_option('d',param.dist).doc("Max orthogonal dist") );
    cmd.add( make_option('a',param.angle).doc("Max angle diff in degree") );
    cmd.add( make_option('o',param.overlap).doc("Min overlap ratio") );
    try {
        cmd.process(argc, argv);
    } catch(const std::string& s) {
        std::cerr << "Error: " << s << std::endl;
        return 1;
    }
    if(argc!=3) {
        std::cerr<<"Usage: "<<argv[0]<<" [options] detect.txt gt.txt\n" << cmd;
        return 1;
    }
    param.angle *= M_PI/180; // Convert to radian
    std::vector<Seg> D, GT;
    D  = read_seg(argv[1]);
    GT = read_seg(argv[2]);
    float tp=0, fp=0, fn=0, tpIOU=0, fpIOU=0;
    for(std::vector<Seg>::const_iterator it=GT.begin(); it!=GT.end(); ++it) {
        if(! eval_seg(*it, D, param, tp, fp, tpIOU, fpIOU))
            fn += it->l;
    }
    float p=tp; if(tp+fp>0) p /= (tp+fp);
    float r=tp; if(tp+fn>0) r /= (tp+fn);
    float iou = fpIOU>0? tpIOU/fpIOU: 0;
    std::cout << "Prec,Rec,F1,IOU= "
              << p << ' ' << r << ' ' << 2*(p*r)/(p+r) << ' '
              << iou << std::endl;
    return 0;
}
