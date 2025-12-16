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
#include <list>
#include <iterator>
#include <cmath>

const std::string COL_TP = "#00ff00";
const std::string COL_FP = "#ff0000";
const std::string COL_FN = "#ffff00";

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

std::ostream& operator<<(std::ostream& str, const Seg& s) {
    return str << s.x1 << ' ' << s.y1 << ' ' << s.x2 << ' ' << s.y2;
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

class UnionIntervals {
    std::list<float> L;
public:
    void insert(float x0, float x1);
    float length() const;
};

void UnionIntervals::insert(float x0, float x1) {
    std::list<float>::iterator i,j;
    bool in0=false, in1=false;
    for(i=L.begin(); i!=L.end() && x0<*i; ++i)
        in0=!in0;
    for(in1=in0,j=i; j!=L.end() && x1<=*j; ++j)
        in1=!in1;
    if(i==j) {
        if(!in0)
            L.insert(i, {x0,x1});
        return;
    }
    if(!in0)
        *i++ = x0;
    if(!in1)
        *--j = x1;
    L.erase(i,j);
}

float UnionIntervals::length() const {
    float l=0;
    for(std::list<float>::const_iterator it=L.begin(); it!=L.end(); ++it) {
        int x0=*it++;
        l += *it-x0;
    }
    return l;
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
    if(d>M_PI/2) d = M_PI-d;
    return d;
}

float get_dir(const Seg& s, float& dx, float& dy) {
    dx=s.x2-s.x1, dy=s.y2-s.y1;
    float d=hypot(dx,dy);
    dx/=d; dy/=d;
    return d;
}

void project(const Seg& sGT, const Seg& sD,
             float dx, float dy, float& x0, float& x1) {
    x0=(sD.x1-sGT.x1)*dx+(sD.y1-sGT.y1)*dy;
    x1=(sD.x2-sGT.x1)*dx+(sD.y2-sGT.y1)*dy;
    if(x0>x1) std::swap(x0,x1);
}

// Length of [x0,x1] inter [x2,x3]
float overlap(float x0, float x1, float x2, float x3) {
    return std::max(0.0f,std::min(x1,x3)-std::max(x0,x2));
}

bool eval_seg(const Seg& sGT, const std::vector<Seg>& D, const Param& param,
              float& tp, float& fp, float& inter, float& uni,
              std::vector<bool>& match) {
    const float c=std::cos(sGT.angle), s=std::sin(sGT.angle);
    float x0=0, dx, dy, x1 = get_dir(sGT, dx, dy);
    float m=x0, M=x1; // min/max of overlapping segments
    bool hit=false;
    UnionIntervals U;
    int i=0;
    for(std::vector<Seg>::const_iterator it=D.begin(); it!=D.end(); ++it, ++i) {
        if(! (perp_dist(sGT,*it,c,s)          <= param.dist &&
              angle_diff(sGT.angle,it->angle) <= param.angle))
            continue;
        float x2, x3;
        project(sGT, *it, dx, dy, x2, x3);
        float over = overlap(x0, x1, x2, x3);
        if(over <= 0)
            continue;
        float l=sGT.l, L=it->l; if(l>L) std::swap(l,L); // min/max lengths
        if(over >= param.overlap * L && l >= param.overlap * over) {
            hit=true;
            match[i] = true;
            if(x2<m) m=x2;
            if(x3>M) M=x3;
            if(x2<x0) x2=x0;
            if(x3>x1) x3=x1;
            U.insert(x2,x3);
            tp += it->l;
        } else
            fp += it->l;
    }
    inter += U.length();
    uni += M-m;
    return hit;
}

int main(int argc, char* argv[]) {
    Param param;
    std::string fileName;
    CmdLine cmd;
    cmd.add( make_option('d',param.dist).doc("Max orthogonal dist") );
    cmd.add( make_option('a',param.angle).doc("Max angle diff in degree") );
    cmd.add( make_option('o',param.overlap).doc("Min overlap ratio") );
    cmd.add( make_option('f',fileName).doc("File for output TP, FP and FN") );
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
    if(param.overlap<0 || param.overlap>1) {
        std::cerr<<"Parameter overlap (option -o) must be in [0,1]"<<std::endl;
        return 1;
    }
    param.angle *= M_PI/180; // Convert to radian
    std::vector<Seg> D, GT;
    D  = read_seg(argv[1]);
    GT = read_seg(argv[2]);
    std::ofstream* file=0;
    if(!fileName.empty()) {
        file = new std::ofstream(fileName);
        if(!file->is_open()) {
            std::cerr << "Unable to create file " << fileName << std::endl;
            return 1;
        }
    }
    float tp=0, fp=0, fn=0, inter=0, uni=0;
    std::vector<bool> match(D.size(), false);
    for(std::vector<Seg>::const_iterator it=GT.begin(); it!=GT.end(); ++it) {
        std::string col=COL_TP;
        if(! eval_seg(*it, D, param, tp, fp, inter, uni, match)) {
            fn += it->l;
            col = COL_FN;
        }
        if(file)
            *file << *it << ' ' << col << std::endl;
    }
    if(file) {
        std::vector<Seg>::const_iterator i=D.begin();
        std::vector<bool>::const_iterator j=match.begin();
        for(; j!=match.end(); ++j, ++i)
            if(!*j)
                *file << *i << ' ' << COL_FP << std::endl;
    }
    delete file;
    float p=tp; if(tp+fp>0) p /= (tp+fp);
    float r=tp; if(tp+fn>0) r /= (tp+fn);
    float f1=(p*r==0)? 0: 2*(p*r)/(p+r);
    float iou = uni>0? inter/uni: 0;
    std::cout << "Prec,Rec,F1,IoU= " << p<<' '<<r<<' '<<f1<<' '<<iou<<std::endl;
    return 0;
}
