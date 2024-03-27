// SPDX-License-Identifier: MPL-2.0
/**
 * @file mulsd.cpp
 * @brief muLSD: main interface
 * @author Yohann Salaun <yohann.salaun@imagine.enpc.fr>
 *         Pascal Monasse <pascal.monasse@enpc.fr>
 * @date 2016, 2023-2024
 */

#include "mulsd.hpp"
#include "cluster.hpp"
#include <cmath>
using namespace std;

/// Classical main entry point of LSD with two additions:
/// 1. Pre-processing: refine segments detected at previous scale;
/// 2. Post-processing: Merge close aligned segments.
/// \param im image at current scale.
/// \param[in,out] segments (in) at previous scale, (out) at current scale.
/// \param grad minimum gradient norm, still scaled by sine of angle precision.
/// \param ang_th angle tolerance (degree), normally 22.5 by default in LSD,
/// set at 45 in muLSD.
/// \param log_eps threshold for detection: -log10(NFA).
/// \param density_th min density of aligned pixels inside a segment.
/// \param n_bins number of bins in quantization of grad norm for linear sort.
void LineSegmentDetection(const Image<float>& im,
                          vector<Segment>& segments,
                          double grad, double ang_th=45.0,
                          double log_eps=0,
                          double density_th=0.70, int n_bins=1024) {
    const double p = ang_th / 180.0;      // probability
    const double prec = M_PI * p;         // angle tolerance
    const double rho = grad / sin(prec); // gradient magnitude threshold

    unsigned int xsize=im.w, ysize=im.h;
    const unsigned int N = im.w*im.h;

    // Compute angle of gradient direction at each pixel
    double* data = new double[N];
    for(unsigned int i=0; i<N; i++)
        data[i]=(double)im.data[i];
    image_double image = new_image_double_ptr(xsize, ysize, data);
    coorlist* list_p; void* mem_p; image_double modgrad;
    image_double angles = ll_angle(image, rho, &list_p, &mem_p, &modgrad,
                                   (unsigned int)n_bins);
    free((void*)image);
    delete[] data;

    const double logNT = 5*log10((double)N)/2+log10(11.0); // Number of Tests
    // Min number of points in region that can give a meaningful event
    const int min_reg_size = (int)(-logNT / log10(p));

    image_char used = new_image_char_ini(xsize, ysize, NOTUSED);
    point* reg = (point*)calloc((size_t)(xsize*ysize), sizeof(point));

    // 1. Pre-processing: Refine line segments from previous scale
    vector<Cluster> clusters =
        refineRawSegments(segments, angles, modgrad, used, logNT, log_eps);

    // 1.5 Regular search outside existing segments' area
    for(; list_p != NULL; list_p = list_p->next)
        if(used->data[list_p->x + list_p->y * xsize] == NOTUSED &&
           angles->data[list_p->x + list_p->y * xsize] != NOTDEF) {
            int reg_size; double reg_angle;
            region_grow(list_p->x, list_p->y, angles, reg, &reg_size,
                        &reg_angle, used, prec); // Find CC
            if(reg_size < min_reg_size) // reject too small regions
                continue;
            rect rec; // rectangular approximation for the region
            region2rect(reg, reg_size, modgrad, reg_angle, prec, p, &rec);
            if(! refine(reg, &reg_size, modgrad, reg_angle,
                        prec, p, &rec, used, angles, density_th))
                continue;
            double log_nfa = rect_improve(&rec, angles, logNT, log_eps);
            if(log_nfa <= log_eps)
                continue;
            clusters.push_back(Cluster(angles, logNT, reg, reg_size, rec,
                                       clusters.size()));
        }

    // 2. Post-processing: merge aligned clusters
    mergeClusters(clusters, angles,modgrad,used, logNT, log_eps);

    // Convert clusters into segments
    segments.clear();
    for(size_t i=0; i<clusters.size(); i++)
        if(! clusters[i].isMerged())
            segments.push_back(clusters[i].toSegment());

    free_image_double(angles);
    free_image_double(modgrad);
    free_image_char(used);
    free((void*)mem_p);
    free(reg);
}

/// In the image pyramid, the first one is the original image and subsequent
/// are scaled-down versions of factor 2.
/// When grad is zero, set it at standard deviation of gradient magnitude. 
std::vector<Segment> lsd_multiscale(const std::vector<Image<float>*>& imgs,
                                    float grad) {
    vector<Segment> segments;
    for(int i=(int)imgs.size()-1; i>=0; i--) {
        cout <<"scale:" <<i <<" (" <<imgs[i]->w <<'x'<<imgs[i]->h <<")" <<flush;
        float minGrad = (grad<=0? stdGradNorm(*imgs[i]): grad);
        if(grad<=0) cout << " grad:" << minGrad << flush;
        LineSegmentDetection(*imgs[i], segments, minGrad);
        cout << " #lines:" << segments.size() << endl;
    }
    return segments;
}
