#include "image.h"
#include <cmath>
#include <string>
#include <cassert>

static const float SIGMA=1.6f;

/// Read PNG file as float image. Check field data to see if file was correct.
template<> Image<float>::Image(const char* filename) {
    size_t nx, ny;
    float* d = io_png_read_f32_gray(filename, &nx, &ny);
    if(! d) {
        w=h=0;
        data = 0;
        return;
    }
    w = (int)nx; h = (int)ny;
    data = new float[nx*ny];
    int n = w*h;
    for(int i=0; i<n; i++)
        data[i] = d[i];
    free(d);
}

/// Return Gaussian kernel of size 2*radius+1. \a radius is set to ceil(3*sigma)
/// with sigma a global constant.
float* gaussKernel(int& radius) {
    radius = (int)ceil(3*SIGMA);
    int sz = 2*radius+1;
    float sum=0, *ker = new float[sz];
    for(int i=0; i<=radius; i++) {
        ker[sz-i-1] = ker[i] = exp(-(i-radius)*(i-radius)/(2*SIGMA*SIGMA));
    }
    for(int i=0; i<sz; i++)
        sum += ker[i];
    for(int i=0; i<sz; i++)
        ker[i] /= sum;
    return ker;
}

/// x-convolution and 1/2 x-subsampling. The result is transposed.
void convx_sub2(const Image<float>& im, const float* ker, int radius,
                Image<float>& out) {
    assert(radius<=2*im.w);
    out.reset(im.h, (im.w+1)/2); // Beware of odd width
    for(int y=0; y<im.h; y++) {
        const float* in=im.data+y*im.w;
        for(int x=0; x<im.w; x+=2) {
            float c=0;
            const float* k = ker;
            for(int i=-radius; i<=+radius; i++) {
                int j = x+i;
                if(j<0) j = -j-1;
                if(j>=im.w) j = 2*im.w-j-1;
                c += in[j] * *k++;
            }
            out(y,x/2) = c;
        }
    }
}

/// Generate Gaussian pyramid. \a n is the number of scales. The first scale
/// is always the initial image.
std::vector<Image<float>*> gaussPyramid(const Image<float>& in, int n) {
    std::vector<Image<float>*> P;
    Image<float>* im = new Image<float>(in.w, in.h);
    for(int i=0; i< in.w*in.h; i++)
        im->data[i] = in.data[i];
    P.push_back(im);
    int radius;
    float* ker = gaussKernel(radius);
    Image<float> xconv;
    for(int i=1; i<n; i++) {
        convx_sub2(*P.back(), ker, 5, xconv);
        im = new Image<float>;
        convx_sub2(xconv, ker, radius, *im);
        P.push_back(im);
    }
    delete [] ker;
    return P;
}

/// Compute standard deviation of gradient norm.
float stdGradNorm(const Image<float>& im) {
    if(im.w<=2 || im.h<=2) // Unable to compute
        return 0;
    const int dx=1, dy=im.w;
    const int n = (im.w-2)*(im.h-2);
    float* grad = new float[n];
    float* g = grad;
    const float* p=im.data+dx+dy; // position (1,1)

    // Computing mean
    float m=0;
    for(int y=1; y+1<im.h; y++) {
        for(int x=1; x+1<im.w; x++, p++) {
            float gx = (p[+dx]-p[-dx])*.5f;
            float gy = (p[+dy]-p[-dy])*.5f;
            m += (*g++ = hypot(gx,gy));
        }
        p += 2;
    }
    m /= n;

    // std
    float s=0;
    for(int i=0; i<n; i++)
        s += (grad[i]-m)*(grad[i]-m);
    s /= n;
    delete [] grad;
    return sqrt(s);
}
