#ifndef FOURIER_H
#define FOURIER_H

#include "objfunc.h"

int * shuffle(int p_value, int * p_factor, int * p_length);
int * shuffle_reg(int p_factor, int p_length);
COMP_VEC * fourierDiagnal(int p_n, int p_is_reverse);
int fourierDiagnal_load(COMP_VEC * p_diagnal, int p_n, int p_is_reverse);
COMP_VEC * FFT1D(COMP_VEC * p_source, COMP_VEC * p_aid, COMP_VEC * p_diagnal, int * p_pos, int p_n, int p_factor);
RAS_OBJ * FFT2D_col(RAS_OBJ * p_source, RAS_OBJ * p_aid, COMP_VEC * p_diagnal, int * p_pos, int p_factor,
					int p_start, int p_end, int is_im);
RAS_OBJ * FFT2D_row(RAS_OBJ * p_source, RAS_OBJ * p_aid, COMP_VEC * p_diagnal, int * p_pos, int p_factor,
					int p_start, int p_end, int is_im);

RAS_OBJ * FFT2D(RAS_OBJ * p_source, RAS_OBJ * p_aid, 
				COMP_VEC * p_diagnal_vertical, int * p_pos_vertical, int p_factor_vertical,
				COMP_VEC * p_diagnal_horizontal, int * p_pos_horizontal, int p_factor_horizontal,
				int p_start, int p_end, int is_h, int is_im);

RAS_OBJ * fourierTran2D(RAS_OBJ * p_object, RAS_OBJ * p_aid, int p_is_reverse);
void setRange(BOUND_BOX * p_box, int rows, int cols, int * p_start, int * p_end, int * is_h);

RAS_OBJ * convolution(RAS_OBJ * p_obj_a, BOUND_BOX * p_box_a, RAS_OBJ * p_obj_b, BOUND_BOX * p_box_b, RAS_OBJ * p_aid_a, RAS_OBJ * p_aid_b, int p_factor_vertical, int p_factor_horizontal);

void fourierDiagnal2D(COMP_VEC * p_diagnal_vertical, int p_rows, int p_factor_vertical,
					  COMP_VEC * p_diagnal_horizontal, int p_cols, int p_factor_horizontal, int p_is_reverse);
int reverseFourierDiagnal(COMP_VEC * p_diagnal, int p_n);
#endif
