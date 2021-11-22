/*
This is the implementation dealing with Fourier and Fast Fourier transformation,
and digital convolution.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fourier.h"

/*
This is the old shuffle functions and not be used for the current
application.
*/
int * shuffle(int p_value, int * p_factor, int * p_length)
{
	int i, j, length, factor;
	int * pos = NULL;
	int count = 1;
	int step;
	factor = factorBy2(p_value);
	step = length = (int)pow(2,factor);
	pos = (int *) calloc (length, sizeof(int));
	
	pos[0] = 0;

	for(i = 0; i < factor; i++)
	{
		step = step >> 1;
		for(j = 0; j < count; j++)
		{
			pos[j+count] = pos[j] + step; 
		}
		count = count << 1;
	}
	* p_factor = factor;
	* p_length = length;
	return pos;
}

/*
For Fast Fourier transformaiton of a 1D vector, we have to first rearrange the positions of
every elements in the vector. The parameter p_length indicates the length of the 1D array.
Here, p_length = pow(2,p_factor). For example, if a 1D array is of length 8=pow(2,3), then
p_length = 8 and p_factor = 3. For the second (001) element in the array, its position in the shuffled
array will be 001 --> shuffle --> 100 = become the seventh element in the shuffled array.
*/

int * shuffle_reg(int p_factor, int p_length)
{
	int i, j;
	int * pos = NULL;
	int count = 1;
	int factor = p_factor;
	int length = p_length;

	pos = (int *) calloc (length, sizeof(int));
	
	pos[0] = 0;

	for(i = 0; i < factor; i++)
	{
		length = length >> 1;
		for(j = 0; j < count; j++)
		{
			pos[j+count] = pos[j] + length; 
		}
		count = count << 1;
	}

	return pos;
}

/*
No be used by the current application
*/

COMP_VEC * fourierDiagnal(int p_n, int p_is_reverse)
{
	int i;
	double factor = -1.0;
	double pi = acos(-1.0);
	double b_a = 2.0 * pi / (double) p_n;
	double angle;
	COMP_VEC * vector = createCompVec(p_n);
	if(p_is_reverse)
		factor = 1.0;
	for(i = 0; i < p_n; i++)
	{
		angle = (double)i * b_a;
		vector->re[i] = cos(angle);
		vector->im[i] = factor * sin(angle);
	}
	return vector;
}

/*
No be used by the current application
*/

int fourierDiagnal_load(COMP_VEC * p_diagnal, int p_n, int p_is_reverse)
{
	int i;
	double factor = -1.0;
	double pi = acos(-1.0);
	double b_a = 2.0 * pi / (double) p_n;
	double * re = p_diagnal->re;
	double * im = p_diagnal->im;
	double angle;
	if(p_is_reverse)
		factor = 1.0;
	for(i = 0; i < p_n; i++)
	{
		angle = (double)i * b_a;
		re[i] = cos(angle);
		im[i] = factor * sin(angle);
	}
	return 1;
}

/*
FFT for 1D vector.
p_source is the source vector.
p_aid and p_diagnal are the vectors for storing calculation results. They switch at
each iteration such that one stores the middle calculation results and the other stores
the final calculation results.
p_pos: stores the positions after shuffle.
p_n: is the size of the 1D vector and p_factor are the corresponding factor that p_n = pow(2,p_factor);

*/

COMP_VEC * FFT1D(COMP_VEC * p_source, COMP_VEC * p_aid, COMP_VEC * p_diagnal, int * p_pos, int p_n, int p_factor)
{
	int i, j, k0, k1, k, p, r;
	int n = p_n;
	int factor = 1;
	int dia_pos = 0;
	COMP_VEC * source = p_source;
	COMP_VEC * aid = p_aid;
	COMP_VEC * temp = NULL;
	double * source_re = source->re;
	double * source_im = source->im;
	double * aid_re = aid->re;
	double * aid_im = aid->im;
	double * dia_re = p_diagnal->re;
	double * dia_im = p_diagnal->im;
	double re_value, im_value;
	for(i = 0; i < p_n; i ++)
	{
		p = p_pos[i];
		aid_re[i] = source_re[p];
		aid_im[i] = source_im[p];
	}

	for(i = p_factor - 1; i >= 0; i--)
	{
		r = p_factor - i;
		n = n >> 1;
		temp = source;
		source = aid;
		aid = temp;
		source_re = source->re;
		source_im = source->im;
		aid_re = aid->re;
		aid_im = aid->im;
		for (k = 0; k < factor; k++)
		{
			dia_pos = k << i;
			for (j = 0; j < n; j++)
			{	
				k0 = k + (j << r);
				k1 = k0 + factor;
				re_value = source_re[k1]*dia_re[dia_pos] - source_im[k1]*dia_im[dia_pos];
				im_value = source_re[k1]*dia_im[dia_pos] + source_im[k1]*dia_re[dia_pos];
				aid_re[k0] = source_re[k0] + re_value;
				aid_im[k0] = source_im[k0] + im_value;
				aid_re[k1] = source_re[k0] - re_value;
				aid_im[k1] = source_im[k0] - im_value;
			}
		}
		factor = factor << 1;
		
	}
	return aid;
}

/*
For a given 2D matrix, in order to do the 2D FFT, we first perform the 1D FFT on each
column of the matrix. Some additional optimization is adopted, for example, for those columns
outside of the object's bounding rectangle, because all the values in them are 0s, they can be
skipped from computation.

p_source: is the source raster object
Samely, we use p_aid to ack a aid matrix storing middle calculation results.
According to FFT, during each iteration, we have to compute the values at the diagnal of the Fourier
Matrix according to some formula, the p_diagnal vector records these values. Essentially a Fourier transformation
is a matrix-vector product between the Fourier Matrix and the source vector (which is to be transformed).
Again p_pos records the shuffled positions.
p_factor means the length of each column = pow(2,p_factor).
p_start and p_end are the starting and ending position of the object's bounding rectangle along the 
horizontal direction.
is_im: indicates do we have to consider the values in the source matrix as complex numbers.
is_im = 1: the values are complex numbers.
is_im = 0: the values are real numbers.
*/

RAS_OBJ * FFT2D_col(RAS_OBJ * p_source, RAS_OBJ * p_aid, COMP_VEC * p_diagnal, int * p_pos, int p_factor,
					int p_start, int p_end, int is_im)
{
	int rows, cols, start, end;
	int i, j, t, k0, k1, k, p, r;
	int factor = 1;
	int dia_pos = 0;
	RAS_OBJ * source = p_source;
	RAS_OBJ * aid = p_aid;
	RAS_OBJ * temp = NULL;
	double ** source_re = source->re;
	double ** source_im = source->im;
	double ** aid_re = aid->re;
	double ** aid_im = aid->im;
	double * dia_re = p_diagnal->re;
	double * dia_im = p_diagnal->im;
	double re_value, im_value;

	rows = p_source->rows;
	cols = p_source->cols;

	if(p_start < 0)
	{
		start = 0;
		end = cols - 1;
	}
	else
	{
		start = p_start;
		end = p_end;
	}

	if(is_im)
	{
		for(i = 0; i < rows; i ++)
		{
			p = p_pos[i];
			for(t = start; t <= end; t++)
			{
				aid_re[i][t] = source_re[p][t];
				aid_im[i][t] = source_im[p][t];
			}
		}
	}
	else
	{
		for(i = 0; i < rows; i ++)
		{
			p = p_pos[i];
			for(t = start; t <= end; t++)
				aid_re[i][t] = source_re[p][t];
		}
	}
	for(i = p_factor - 1; i >= 0; i--)
	{
		r = p_factor - i;
		rows = rows >> 1;
		temp = source;
		source = aid;
		aid = temp;
		source_re = source->re;
		source_im = source->im;
		aid_re = aid->re;
		aid_im = aid->im;
		for (k = 0; k < factor; k++)
		{
			dia_pos = k << i;
			for (j = 0; j < rows; j++)
			{	
				k0 = k + (j << r);
				k1 = k0 + factor;
				for(t = start; t <= end; t++)
				{
					re_value = source_re[k1][t]*dia_re[dia_pos] - source_im[k1][t]*dia_im[dia_pos];
					im_value = source_re[k1][t]*dia_im[dia_pos] + source_im[k1][t]*dia_re[dia_pos];
					aid_re[k0][t] = source_re[k0][t] + re_value;
					aid_im[k0][t] = source_im[k0][t] + im_value;
					aid_re[k1][t] = source_re[k0][t] - re_value;
					aid_im[k1][t] = source_im[k0][t] - im_value;
				}
			}
		}
		factor = factor << 1;
		
	}
	return aid;
}

/*
For a given matrix (raster object), perform the FFT on each row.
*/

RAS_OBJ * FFT2D_row(RAS_OBJ * p_source, RAS_OBJ * p_aid, COMP_VEC * p_diagnal, int * p_pos, int p_factor,
					int p_start, int p_end, int is_im)
{
	int i, j, t, k0, k1, k, p, r, start, end;
	int cols, rows;
	int factor = 1;
	int dia_pos = 0;
	RAS_OBJ * source = p_source;
	RAS_OBJ * aid = p_aid;
	RAS_OBJ * temp = NULL;
	double ** source_re = source->re;
	double ** source_im = source->im;
	double ** aid_re = aid->re;
	double ** aid_im = aid->im;
	double * dia_re = p_diagnal->re;
	double * dia_im = p_diagnal->im;
	double re_value, im_value;

	cols = p_source->cols;
	rows = p_source->rows;

	if(p_start < 0)
	{
		start = 0;
		end = rows - 1;
	}
	else
	{
		start = p_start;
		end = p_end;
	}

	if(is_im)
	{
		for(i = 0; i < cols; i ++)
		{
			p = p_pos[i];
			for(t = start; t <= end; t++)
			{
				aid_re[t][i] = source_re[t][p];
				aid_im[t][i] = source_im[t][p];
			}
		}
	}
	else
	{
		for(i = 0; i < cols; i ++)
		{
			p = p_pos[i];
			for(t = start; t <= end; t++)
				aid_re[t][i] = source_re[t][p];
		}
	}

	for(i = p_factor - 1; i >= 0; i--)
	{
		r = p_factor - i;
		cols = cols >> 1;
		temp = source;
		source = aid;
		aid = temp;
		source_re = source->re;
		source_im = source->im;
		aid_re = aid->re;
		aid_im = aid->im;
		for (k = 0; k < factor; k++)
		{
			dia_pos = k << i;
			for (j = 0; j < cols; j++)
			{	
				k0 = k + (j << r);
				k1 = k0 + factor;
				for(t = start; t <= end; t++)
				{
					re_value = source_re[t][k1]*dia_re[dia_pos] - source_im[t][k1]*dia_im[dia_pos];
					im_value = source_re[t][k1]*dia_im[dia_pos] + source_im[t][k1]*dia_re[dia_pos];
					aid_re[t][k0] = source_re[t][k0] + re_value;
					aid_im[t][k0] = source_im[t][k0] + im_value;
					aid_re[t][k1] = source_re[t][k0] - re_value;
					aid_im[t][k1] = source_im[t][k0] - im_value;
				}
			}
		}
		factor = factor << 1;
		
	}
	return aid;
}

/*
2D FFT. For a given 2D matrix, we first perform 1D transform on each column (or row) and then
perform 1D transform on each row(column). This is the basic principle. 

p_source: is the source matrix (raster object)
p_aid: is the matrix for storing the middle calculatin results (results after the 1st round transfromation).
p_diagnal_vertical: The diagnal setting of the Fourier Matrix according to the length of each column.
p_pos_vertical: the shuffled positions according to the length of each column.
p_factor_vertical: means the length of each column = pow(2,p_factor_vertical).

p_diagnal_horizontal, p_pos_horizontal and p_factor_horizontal are the corresponding parameters regarding to each row.

is_h: indicate we do first transform on columns or on rows.
is_h = 1: we first perform 1D transform on columns
          p_start and p_end indicate the ends of bounding rectangle along horizontal direction.
          If p_start = -1, then we have to do the computations on all the columns
          
is_h = 0: we first perform 1D transform on rows
          p_start and p_end indicate the ends of bounding rectangle along vertical direction.
          if p_start = -1, we have to do the computations on all the rows.

is_im indicate whether the source object is complex and real.
is_im = 1: the object is complex (all the values have to be considered as complex numbers).
is_im = 0: the object is real (all the values are real numbers).

between the two transforms, we switch p_aid and p_source to let one be the aid matrix and one stores
the final results. This implementation will save the memory usage.
*/

RAS_OBJ * FFT2D(RAS_OBJ * p_source, RAS_OBJ * p_aid, 
				COMP_VEC * p_diagnal_vertical, int * p_pos_vertical, int p_factor_vertical,
				COMP_VEC * p_diagnal_horizontal, int * p_pos_horizontal, int p_factor_horizontal,
				int p_start, int p_end, int is_h, int is_im)
{
	RAS_OBJ * source = p_source;
	RAS_OBJ * aid = p_aid;
	RAS_OBJ * temp = NULL;
	if(is_h)
	{
		temp = FFT2D_col(source,aid,p_diagnal_vertical,p_pos_vertical,p_factor_vertical, p_start, p_end, is_im); /*1D transform on each column*/

		if(temp == aid)
				aid = source;                                                                                    /*switch p_aid an p_source*/
		return (FFT2D_row(temp,aid,p_diagnal_horizontal,p_pos_horizontal,p_factor_horizontal, -1, 0, 1));        /*1D transform on each row, 
                                                                                                                 this time we have to consider all the rows*/


	}
	else
	{
		temp = FFT2D_row(source,aid,p_diagnal_horizontal,p_pos_horizontal,p_factor_horizontal, p_start, p_end, is_im); /*1D transform on each row*/

		if(temp == aid)
				aid = source;                                                                                          /*switch the p_aid and p_source*/
		return (FFT2D_col(temp,aid,p_diagnal_vertical,p_pos_vertical,p_factor_vertical, -1, 0, 1));                    /*1D transform on each column
                                                                                                                       this time we have to consider all the columns*/


	}

	
}

/*
The application of 2D FFT.
p_object is the source raster object (either complex or real)
p_aid is again the aid matrix for storing the middle computation results.
p_is_reverse indicates whether the current operatin is reversed FFT or FFT:
1: reversed FFT
0: FFT.
*/

RAS_OBJ * fourierTran2D(RAS_OBJ * p_object, RAS_OBJ * p_aid, int p_is_reverse)
{
	RAS_OBJ * result;
	BOUND_BOX * box;
	
	int rows = p_object->rows;
	int cols = p_object->cols;
	int factor_vertical, factor_horizontal;

	COMP_VEC * diagnal_vertical, * diagnal_horizontal;
	int * pos_vertical, * pos_horizontal;
	
	box = loadBoundBox(p_object);
	factor_vertical = factorBy2(rows);
	factor_horizontal = factorBy2(cols);

	pos_vertical = shuffle_reg(factor_vertical, rows);
	pos_horizontal = shuffle_reg(factor_horizontal, cols);
	
	diagnal_vertical = createCompVec(rows);
	diagnal_horizontal = createCompVec(cols);
	
	fourierDiagnal2D(diagnal_vertical, rows, factor_vertical, diagnal_horizontal, cols, factor_horizontal, p_is_reverse);

	result = FFT2D(p_object, p_aid, diagnal_vertical, pos_vertical, factor_vertical, diagnal_horizontal, pos_horizontal, factor_horizontal, -1,0,1,1);

	destroyCompVec(diagnal_vertical);
	destroyCompVec(diagnal_horizontal);
	free(pos_vertical);
	free(pos_horizontal);

	return result;
	
}

/*
For a given raster object, calculate the length of its bounding rectangle along the vertical
and the horizontal directions.
And use certain formula to decide whether we should first do the column transform or the row transform in order to
save as much computational time as we can.

The results will be passed out by parameters p_start and p_end (the ends of the bounding rectangel), and is_h
indicating whether we do column transform first or row transform first.
*/

void setRange(BOUND_BOX * p_box, int rows, int cols, int * p_start, int * p_end, int * is_h)
{
	int range_vertical, range_horizontal;
	double value_vertical, value_horizontal;

	range_vertical = p_box->v1[0] - p_box->v3[0] + 1;
	range_horizontal = p_box->v3[1] - p_box->v1[1] + 1; 

	value_vertical = (double)range_vertical * log(cols);
	value_horizontal = (double)range_horizontal * log(rows);

	if(value_vertical > value_horizontal)
	{
		*is_h = 1;
		*p_start = p_box->v1[1];
		*p_end = p_box->v3[1];
	}
	else
	{
		*is_h = 0;
		*p_start = p_box->v3[0];
		*p_end = p_box->v1[0];
	}

	
}

/*
2D convolution between two raster objects (real objects).
p_obj_a: is the raster object A (reference object).
P_box_a: is the minimum bounding rectangle of A.
p_obj_b: is the raster object B (argument object).
p_box_b: is the minimum bounding rectangle of B.

p_aid_a and p_aid_b are the two empty objects for storing the middle computation results.

p_factor_vertical: means the height of the image = pow(2,p_factor_vertical);
p_factor_horizontal: means the width of the image = pow(2,p_factor_horizontal).
*/

RAS_OBJ * convolution(RAS_OBJ * p_obj_a, BOUND_BOX * p_box_a, RAS_OBJ * p_obj_b, BOUND_BOX * p_box_b, RAS_OBJ * p_aid_a, RAS_OBJ * p_aid_b, int p_factor_vertical, int p_factor_horizontal)
{
	int i,j, start, end, is_h;
	int rows, cols;
	double size;
	RAS_OBJ * obj_a = p_obj_a;
	RAS_OBJ * obj_b = p_obj_b;
	RAS_OBJ * aid_a = p_aid_a;
	RAS_OBJ * aid_b = p_aid_b;
	RAS_OBJ * fourier_obj_a = NULL;
	RAS_OBJ * fourier_obj_b = NULL;
	RAS_OBJ * result = NULL;
	int * pos_vertical = NULL;
	int * pos_horizontal = NULL;
	COMP_VEC * diagnal_vertical = NULL;
	COMP_VEC * diagnal_horizontal = NULL;

	double re_value, im_value;

	double ** obj_a_re, ** obj_a_im;
	double ** obj_b_re, ** obj_b_im;

	rows = p_obj_a->rows;
	cols = p_obj_a->cols;

	pos_vertical = shuffle_reg(p_factor_vertical, rows);              /*position shuffle according to the length of row*/
	pos_horizontal = shuffle_reg(p_factor_horizontal, cols);          /*position shuffle according to the length of column*/

	diagnal_vertical = createCompVec(rows);
	diagnal_horizontal = createCompVec(cols);

	/*
	Calculate the values at diagnal of Fourier Matrices (one for the row transform and the other for the column transform).
    */
    fourierDiagnal2D(diagnal_vertical, rows, p_factor_vertical, diagnal_horizontal, cols, p_factor_horizontal, 0);

	/*Compute the ends of bounding rectangel of object A and decide do column transform first or row transform first*/
    setRange(p_box_a, rows, cols, &start, &end, &is_h);                

	/* FFT transformation of object A*/
    fourier_obj_a = FFT2D(obj_a, aid_a, diagnal_vertical, pos_vertical, p_factor_vertical, diagnal_horizontal, pos_horizontal, p_factor_horizontal,
		                  start, end, is_h, 0);


    /*Compute the ends of bounding rectangel of object B and decide do column transform first or row transform first*/
	setRange(p_box_b, rows, cols, &start, &end, &is_h);

	/* FFT transformation of object B*/
    fourier_obj_b = FFT2D(obj_b, aid_b, diagnal_vertical, pos_vertical, p_factor_vertical, diagnal_horizontal, pos_horizontal, p_factor_horizontal,
		                  start, end, is_h, 0);


	obj_a_re = fourier_obj_a->re;
	obj_a_im = fourier_obj_a->im;
	obj_b_re = fourier_obj_b->re;
	obj_b_im = fourier_obj_b->im;
	
	size = (double) rows * (double) cols;

	
    /*Product of the transformed object A snd the transformed object B*/
    for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			re_value = obj_a_re[i][j] * obj_b_re[i][j] - obj_a_im[i][j] * obj_b_im[i][j];
			im_value = obj_a_re[i][j] * obj_b_im[i][j] + obj_a_im[i][j] * obj_b_re[i][j];
			obj_a_re[i][j] = re_value;
			obj_a_im[i][j] = im_value;
		}
	}
	
	reverseFourierDiagnal(diagnal_vertical, rows);       /*Compute the values at diagnal of the reversed Fourier Matrix for rows*/
	reverseFourierDiagnal(diagnal_horizontal, cols);     /*Compute the values at diagnal of the reversed Fourier Matrix for columns*/

	
	if(fourier_obj_a == aid_a)
		aid_a = obj_a;

	
    /*Compute the reversed FFT */
    result = FFT2D(fourier_obj_a, aid_a, diagnal_vertical, pos_vertical, p_factor_vertical, diagnal_horizontal, pos_horizontal, p_factor_horizontal,
					-1, 0, 1, 1);
    
	/*Normalize the results*/
    for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			result->re[i][j] = result->re[i][j] / size;

	free(pos_vertical);
	free(pos_horizontal);
	free(diagnal_vertical);
	free(diagnal_horizontal);

	return result;
}

/*
	Calculate the values at diagnal of Fourier Matrices (one for the row transform and the other for the column transform).
*/

void fourierDiagnal2D(COMP_VEC * p_diagnal_vertical, int p_rows, int p_factor_vertical,
					  COMP_VEC * p_diagnal_horizontal, int p_cols, int p_factor_horizontal, int p_is_reverse)
{
	int i;
	double * source_re, * source_im;
	double * target_re, * target_im;

	int factor, pos, min_l;

	

	if(p_factor_vertical > p_factor_horizontal)
	{
		factor = p_factor_vertical - p_factor_horizontal;
		min_l = p_cols; 
		fourierDiagnal_load(p_diagnal_vertical, p_rows, p_is_reverse);
		
		source_re = p_diagnal_vertical->re;
		source_im = p_diagnal_vertical->im;
		target_re = p_diagnal_horizontal->re;
		target_im = p_diagnal_horizontal->im;
	}
	else
	{
		factor = p_factor_horizontal - p_factor_vertical;
		min_l = p_rows;
		fourierDiagnal_load(p_diagnal_horizontal, p_cols, p_is_reverse);
		target_re = p_diagnal_vertical->re;
		target_im = p_diagnal_vertical->im;
		source_re = p_diagnal_horizontal->re;
		source_im = p_diagnal_horizontal->im;
		
	}

	for(i = 0; i < min_l; i++)
	{
		pos = i << factor; 
		target_re[i] = source_re[pos];
		target_im[i] = source_im[pos];
	}
}

int reverseFourierDiagnal(COMP_VEC * p_diagnal, int p_n)
{
	int i;
	double * im = p_diagnal->im;
	for(i = 0; i < p_n; i++)
	{
		im[i] = -1.0 * im[i];
	}
	return 1;
}
