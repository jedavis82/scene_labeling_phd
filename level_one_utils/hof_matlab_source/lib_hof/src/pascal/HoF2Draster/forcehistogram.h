#ifndef FORCEHISTOGRAM_H
#define FORCEHISTOGRAM_H

#include "objfunc.h"

#define A	0
#define B	1
#define AB	2

/*Line segment*/

typedef struct segment
{
	double start;       /*Start position of the segment along a straight line*/
	double end;         /*End position of the segment*/
	double length;      /*Length of the segment*/

}SEGMENT;
  
int getSegments(SEGMENT * p_seg, double ** p_obj, int * p_line, int p_offset, int p_start, int p_end, int p_n, int p_d);
int reverseSegments(SEGMENT * p_seg, int p_n);

double segmentForce(SEGMENT * p_seg_a, SEGMENT * p_seg_b, int p_a, int p_b, double (*F)(double x, double y, double z));

void force_histogram(double * histogram, int p_d, double (*F)(double x, double y, double z), int * line, int rows, int cols, double p_f,
					RAS_OBJ * p_obj_a, BOUND_BOX * box_a, SEGMENT * segment_a,
					RAS_OBJ * p_obj_b, BOUND_BOX * box_b, SEGMENT * segment_b);

double * forceHistogram(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b, double (*F)(double x, double y, double z));
double * forceHistogram_double(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b, double (*F)(double x, double y, double z));
double * forceHistogram_single(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b, double (*F)(double x, double y, double z), int type);


double segmentForce_new(SEGMENT * p_seg_a, SEGMENT * p_seg_b, int p_a, int p_b, double r_t, double (*F)(double x, double y, double z, double r));

void force_histogram_new(double * histogram, int p_d, double r_t, double (*F)(double x, double y, double z, double r), int * line, int rows, int cols, double p_f,
						 RAS_OBJ * p_obj_a, BOUND_BOX * box_a, SEGMENT * segment_a,
						 RAS_OBJ * p_obj_b, BOUND_BOX * box_b, SEGMENT * segment_b);

double * forceHistogram_new(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b, double r_t, double (*F)(double x, double y, double z, double r));
double * forceHistogram_double_new(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b, double r_t, double (*F)(double x, double y, double z, double r));
double * forceHistogram_single_new(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b, double r_t, double (*F)(double x, double y, double z, double r), int type);


#endif
