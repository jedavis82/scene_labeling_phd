#ifndef FRAMEWORK_H
#define FRAMEWORK_H

#include "objfunc.h"

RAS_OBJ * createLandscape_slow(RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);
RAS_OBJ * createLandscape_conv(RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);
RAS_OBJ * createLandscape_conv_opt(RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);

int landscape_angle_histogram(int p_d, double * p_histogram, RAS_OBJ * p_scape);
double * landscapeAngleHistogram(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);

double adjF0(double angle);
double adjF2(double angle);


int landscape_force_histogram_0(int p_d, double * p_histogram, RAS_OBJ * p_scape);
double * landscapeForceHistogram_0(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);

int landscape_force_histogram_2(int p_d, double * p_histogram, RAS_OBJ * p_scape);
double * landscapeForceHistogram_2(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);

int landscape_force_histogram_r(double r_t, int p_d, double * p_histogram, RAS_OBJ * p_scape);
double * landscapeForceHistogram_r(double r_t, int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);

int landscape_force_histogram_r_precise(double r_t, int p_d, double * p_histogram, RAS_OBJ * p_scape);
double * landscapeForceHistogram_r_precise(double r_t, int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b);

double F_less_1_opt(double y, double r);
double F_1_opt(double y, double r);
double F_1_2_opt(double y, double r);
double F_2_opt(double y, double r);
double F_bigger_2_opt(double y, double r);

#endif
