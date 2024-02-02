/*----------------------------------------------------------------------------

  LSD - Line Segment Detector on digital images

  This code is part of the following publication and was subject
  to peer review:

  "LSD: a Line Segment Detector" by Rafael Grompone von Gioi,
  Jeremie Jakubowicz, Jean-Michel Morel, and Gregory Randall,
  Image Processing On Line, 2012. DOI:10.5201/ipol.2012.gjmr-lsd
  http://dx.doi.org/10.5201/ipol.2012.gjmr-lsd

  Copyright (c) 2007-2011 rafael grompone von gioi <grompone@gmail.com>

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program. If not, see <http://www.gnu.org/licenses/>.

  ----------------------------------------------------------------------------*/
#ifndef LSD_HPP
#define LSD_HPP

/** PI */
#ifndef M_PI
#define M_PI   3.14159265358979323846
#endif /* !M_PI */

/** Label for pixels with undefined gradient. */
#define NOTDEF -1024.0

/** Label for pixels not used in yet. */
#define NOTUSED 0

/** Label for pixels already used in detection. */
#define USED    1

/** Chained list of coordinates. */
struct coorlist
{
    int x, y;
    struct coorlist * next;
};

/** A point (or pixel). */
struct point {int x, y;};

/** Rectangle structure: line segment with width. */
struct rect
{
    double x1, y1, x2, y2; /* first and second point of the line segment */
    double width;          /* rectangle width */
    double x, y;           /* center of the rectangle */
    double theta;          /* angle */
    double dx, dy;         /* (dx,dy) is vector oriented as the line segment */
    double prec;           /* tolerance angle */
    double p;              /* probability of a point with angle within 'prec' */
};

/*----------------------------------------------------------------------------*/
/*----------------------------- Image Data Types -----------------------------*/
/*----------------------------------------------------------------------------*/

/** double image data type.

  The pixel value at (x,y) is accessed by:

  image->data[ x + y * image->xsize ]

  with x and y integer.
  */
typedef struct image_double_s
{
    double * data;
    unsigned int xsize, ysize;
} *image_double;

/** Create a new image_double of size 'xsize' times 'ysize'
  with the data pointed by 'data'.
  */
image_double new_image_double_ptr(unsigned int xsize,
                                  unsigned int ysize, double * data);

/** Free memory used in image_double 'i'. */
void free_image_double(image_double i);

/** char image data type.

  The pixel value at (x,y) is accessed by:

  image->data[ x + y * image->xsize ]

  with x and y integer.
  */
typedef struct image_char_s
{
    unsigned char * data;
    unsigned int xsize, ysize;
} *image_char;

/** Create a new image_char of size 'xsize' times 'ysize',
  initialized to the value 'fill_value'.
  */
image_char new_image_char_ini(unsigned int xsize, unsigned int ysize,
                              unsigned char fill_value);

/** Free memory used in image_char 'i'. */
void free_image_char(image_char i);

/** Computes the direction of the level line of 'in' at each point. */
image_double ll_angle(image_double in, double threshold,
                      struct coorlist ** list_p, void ** mem_p,
                      image_double * modgrad, unsigned int n_bins);

/*----------------------------------------------------------------------------*/
/*------------------------------------ NFA -----------------------------------*/
/*----------------------------------------------------------------------------*/

/** Compute a rectangle's NFA value. */
double rect_nfa(struct rect * rec, image_double angles, double logNT);

/** Build a region of pixels that share the same angle, up to a
    tolerance 'prec', starting at point (x,y).
 */
void region_grow( int x, int y, image_double angles, struct point * reg,
                  int * reg_size, double * reg_angle, image_char used,
                  double prec );

/** Computes a rectangle that covers a region of points. */
void region2rect(struct point * reg, int reg_size,
                 image_double modgrad, double reg_angle,
                 double prec, double p, struct rect * rec);

/** Refine a rectangle. */
int refine(struct point * reg, int * reg_size, image_double modgrad,
           double reg_angle, double prec, double p, struct rect * rec,
           image_char used, image_double angles, double density_th);

/** Try some rectangles variations to improve NFA value. Only if the
  rectangle is not meaningful (i.e., log_nfa <= log_eps).
  */
double rect_improve(struct rect * rec, image_double angles,
                    double logNT, double log_eps);

/*----------------------------------------------------------------------------*/
/*------------------------- Miscellaneous functions --------------------------*/
/*----------------------------------------------------------------------------*/

/** Absolute value angle difference. */
double angle_diff(double a, double b);

#endif
