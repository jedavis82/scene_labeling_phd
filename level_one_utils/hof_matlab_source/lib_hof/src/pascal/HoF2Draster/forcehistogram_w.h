#ifndef W_FORCE_HISTOGRAM_H
#define W_FORCE_HISTOGRAM_H

#include "objfunc.h"

/*record of a node in the tree pyramid*/
typedef struct s_node
{
	double count_a;       /*number of pixels belonging to object A on a line*/
	double count_b;       /*number of pixels belonging to object B on  a line*/
	double force;         /*the force appears on the line*/
} NODE;

/*record a line segment*/
typedef struct seg
{
	int end_left;        /*the end of the left half of the segment*/
	int begin_right;     /*the start of the right half of the segment*/
}SEG;

int line_force_w(double * forward_f, double * backward_f, double ** obj_a, double ** obj_b, int * p_line, int p_offset, int p_start, int p_m1, int p_m2, int p_end, int p_n, int p_d);
void force_histogram_w(double * histogram, int p_d, int * line, int rows, int cols, RAS_OBJ * p_obj_a, BOUND_BOX * box_a, RAS_OBJ * p_obj_b, BOUND_BOX * box_b);
double * forceHistogram_w(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);

NODE * init_pyramid(int n, int factor); 

SEG ** create_seg(int factor);
void destroy_seg(SEG ** segment,int factor);
SEG ** build_seg(int factor);


void assign_pyramid_1(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n); // 0 - pi/4
void assign_pyramid_2(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n); // pi/4 - pi/2
void assign_pyramid_3(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n); // pi/2 - 3pi/4
void assign_pyramid_4(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n); // 3pi/4 - pi
void assign_pyramid_5(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n); // pi - 5pi/4
void assign_pyramid_6(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n); // 5pi/4 - 3pi/2
void assign_pyramid_7(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n); // 3pi/2 - 7pi/4
void assign_pyramid_8(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n); // 7pi/4 - 2pi


NODE * calculate_pyramid(NODE * py_low, NODE * py_high, SEG ** segment, int n, int factor);

void clean_pyramid(NODE * pyramid, int n);
void destroy_pyramid(NODE * pyramid);


double * forceHistogram_w_2(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);

#endif
