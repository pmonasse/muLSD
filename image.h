/**
 * SPDX-License-Identifier: GPL-3.0-or-later
 * @file image.h
 * @brief Elementary image class
 *
 * Copyright (c) 2021-2022 Robin Gay, Pascal Monasse
 */

#ifndef IMAGE_H
#define IMAGE_H

#include "io_png.h"
#include <vector>
#include <cstring>
#include <cstdlib>

/// A simple (simplistic?) image class
template <typename T>
struct Image {
    int w,h;
    T* data;
public:
    Image();
    Image(int w0, int h0);
    Image(int w0, int h0, unsigned char* rawdata);
    Image(const char* filename);
    ~Image();

    void reset(int w0, int h0);

    bool operator!() const { return (w==0 || h==0); }
    void fill(T value);

    T operator()(int i, int j) const;
    T& operator()(int i, int j);

    bool inside(int x, int y) const {
        return (0<=x && x<w && 0<=y && y<h);
    }

private:
    Image(const Image&);
    void operator=(const Image&);
};

std::vector<Image<float>*> gaussPyramid(const Image<float>& im, int n);

/// Constructor
template <typename T>
Image<T>::Image()
: w(0), h(0), data(0) {}

/// Constructor
template <typename T>
Image<T>::Image(int w0, int h0)
: w(w0), h(h0) {
    data = new T[w*h];
}

/// Constructor copying data
template <typename T>
Image<T>::Image(int w0, int h0, unsigned char* rawdata)
: w(w0), h(h0) {
    data = new T[w*h];
    std::memcpy(data, rawdata, w*h*sizeof(T));
}

/// Not implemented for generic type, specialized for float.
template <typename T>
Image<T>::Image(const char*) {
    w = h = 0;
    data = 0;
}
/// Specialization for float, see image.cpp
template<> Image<float>::Image(const char* filename);

/// Destructor
template <typename T>
Image<T>::~Image() {
    delete [] data;
}

/// Reset dimensions of image
template <typename T>
void Image<T>::reset(int w0, int h0) {
    delete [] data;
    w = w0;
    h = h0;
    data = new T[w*h];
}

/// Fill with constant value
template <typename T>
void Image<T>::fill(T value) {
    T* p=data;
    for(int j=0; j<h; j++)
        for(int i=0; i<w; i++)
            *p++ = value;
}

/// Pixel access (read-only)
template <typename T>
T Image<T>::operator()(int i, int j) const {
    return data[i+j*w];
}

/// Pixel access (read/write)
template <typename T>
T& Image<T>::operator()(int i, int j) {
    return data[i+j*w];
}

#endif
