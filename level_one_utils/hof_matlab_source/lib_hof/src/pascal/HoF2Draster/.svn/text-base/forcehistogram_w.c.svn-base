/*
Computing constant (r=0) force histogram using Wang's algorithms
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "util.h"
#include "forcehistogram_w.h"


/* 
Function for Wang's first algorithm
It calculates the force along a rasterized line in certain direction.
forward_f: records the force along the direction
backward_f: records the force along the opposite direction.
obj_a: the 2D array of raster object A
obj_b: the 2D array of raster object B
p_line and p_offset define a rasterized straight line.
p_start p_m1, p_m2 and p_end define the search range along the line.
p_n: the length of the line (along the horizontal or vertical direction).
p_d: whether the line is close to the horizontal or the vertical direction:
p_d = 1: close to the horizontal directin.
p_d = 0: close to the vertical direction.
*/
int line_force_w(double * forward_f, double * backward_f, double ** obj_a, double ** obj_b, int * p_line, int p_offset, int p_start, int p_m1, int p_m2, int p_end, int p_n, int p_d)
{
	double forward_force = 0.0;
	double backward_force = 0.0;
	double c_a = 0.0;
	double c_b = 0.0;
	double temp;
	int i,j;

	if(!p_d)
	{
		for(j = p_start; j <= p_m1; j++)
		{
			i = p_line[j] + p_offset;
			if( i < 0 || i >= p_n)
				continue;
			c_a += obj_a[i][j];
			c_b += obj_b[i][j];
			temp = obj_b[i][j] * obj_a[i][j] / 2.0;
			forward_force += (c_a * obj_b[i][j] - temp);
			backward_force +=(c_b * obj_a[i][j] - temp);
		}
		for(j = p_m2; j <= p_end; j++)
		{
			i = p_line[j] + p_offset;
			if( i < 0 || i >= p_n)
				continue;
			c_a += obj_a[i][j];
			c_b += obj_b[i][j];
			temp = obj_b[i][j] * obj_a[i][j] / 2.0;
			forward_force += (c_a * obj_b[i][j] - temp);
			backward_force +=(c_b * obj_a[i][j] - temp);
		}
		
	}
	else
	{
		for(j = p_start; j <= p_m1; j++)
		{
			i = p_line[j] + p_offset;
			if( i < 0 || i >= p_n)
				continue;
			c_a += obj_a[j][i];
			c_b += obj_b[j][i];
			temp = obj_b[j][i] * obj_a[j][i] / 2.0;
			forward_force += (c_a * obj_b[j][i] - temp);
			backward_force += (c_b * obj_a[j][i] - temp);
		}
		for(j = p_m2; j <= p_end; j++)
		{
			i = p_line[j] + p_offset;
			if( i < 0 || i >= p_n)
				continue;
			c_a += obj_a[j][i];
			c_b += obj_b[j][i];
			temp = obj_b[j][i] * obj_a[j][i] / 2.0;
			forward_force += (c_a * obj_b[j][i] - temp);
			backward_force += (c_b * obj_a[j][i] - temp);
		}
	}
	(*forward_f) = forward_force;
	(*backward_f) = backward_force;
	return 1;
}

/*
Compute the constant force histogram between raster (fuzzy) objects using Wang's first algorithm
histogram: recording the result constant force histogram.
p_d: the number of reference directions considered.
rows: the height of the image.
cols: the width of the image.
p_obj_a: the raster object A.
box_a: A's minimum bounding rectangle.
p_obj_b: the raster object B.
box_b: B's minimum bounding rectangle.
*/

void force_histogram_w(double * histogram, int p_d, int * line, int rows, int cols, RAS_OBJ * p_obj_a, BOUND_BOX * box_a, RAS_OBJ * p_obj_b, BOUND_BOX * box_b)
{
	double step_angle = 360.0 / (double) p_d;
	int d_8 = p_d / 8;
	int d_2 = p_d / 2;
	int i,j, pos, rev_pos;
	double ** a_re = p_obj_a->re;
	double ** b_re = p_obj_b->re;
	double tan_value;
	double angle, l_angle;
	double pi = acos(-1.0);
	int start_a, end_a, start_b, end_b, start, end;
	int from, middle1, middle2, to;
	double forward_value, backward_value;
	double forward_force, backward_force;
	double factor = 0.0;

	
	if(box_a->v1[1] < box_b->v1[1])
	{
		from = box_a->v1[1];
		middle1 = box_a->v3[1];
		middle2 = middle1 >= box_b->v1[1] ? middle1 + 1: box_b->v1[1];
	}
	else
	{
		from = box_b->v1[1];
		middle1 = box_b->v3[1];
		middle2 = middle1 >= box_a->v1[1] ? middle1 + 1: box_a->v1[1];
	}
	to = box_a->v3[1] < box_b->v3[1] ? box_b->v3[1] : box_a->v3[1];

	for(i = -1 * d_8; i<= d_8; i++)
	{
		forward_force = 0.0;
		backward_force = 0.0;
		
		angle = (double)i * step_angle;
		l_angle = angle * pi / 180.0;

		if(angle == -45)
			tan_value = -1.0;
		else
		{
			if(angle == 45)
				tan_value = 1.0;
			else
				tan_value = tan(l_angle);
		}
		factor = 1.0 / cos(l_angle);
		
		populateLine(line,tan_value,cols);
		if(angle>=0)
		{
			projection(&start_a, &end_a,box_a->v3, box_a->v1, tan_value, 0);
			projection(&start_b, &end_b,box_b->v3, box_b->v1, tan_value, 0);
			pos = i;
		}
		else
		{
			projection(&start_a, &end_a,box_a->v2, box_a->v4, tan_value, 0);
			projection(&start_b, &end_b,box_b->v2, box_b->v4, tan_value, 0);
			pos = i + p_d;
		}
		rev_pos = i + d_2;

		start = start_a > start_b ? start_a : start_b;
		end = end_a < end_b ? end_a : end_b;
		for( j = start; j <= end; j++)
		{
			line_force_w(&forward_value,&backward_value,p_obj_a->re,p_obj_b->re,line,j,from,middle1,middle2,to,rows,0);
			forward_force += forward_value;
			backward_force += backward_value;
		}
		histogram[pos] += (forward_force * factor);
		histogram[rev_pos] += (backward_force * factor);
	}

	if(box_a->v3[0] < box_b->v3[0])
	{
		from = box_a->v3[0];
		middle1 = box_a->v1[0];
		middle2 = middle1 >= box_b->v3[0] ? middle1 + 1 : box_b->v3[0];
	}
	else
	{
		from = box_b->v3[0];
		middle1 = box_b->v1[0];
		middle2 = middle1 >= box_a->v3[0] ? middle1 + 1 : box_a->v3[0];
	}
	to = box_a->v1[0] < box_b->v1[0] ? box_b->v1[0] : box_a->v1[0];
	for( i = d_8 + 1; i< (d_2 - d_8); i++)
	{
		forward_force = 0.0;
		backward_force = 0.0;

		angle = 90.0 - (double)i * step_angle;
		l_angle = angle * pi / 180.0;
		tan_value = tan( l_angle);
		factor = 1.0 / cos(l_angle);
		populateLine(line,tan_value,rows);
		if(angle>=0)
		{
			projection(&start_a, &end_a,box_a->v1, box_a->v3, tan_value, 1);
			projection(&start_b, &end_b,box_b->v1, box_b->v3, tan_value, 1);
		}
		else
		{
			projection(&start_a, &end_a,box_a->v2, box_a->v4, tan_value, 1);
			projection(&start_b, &end_b,box_b->v2, box_b->v4, tan_value, 1);
		}
		pos = i;
		rev_pos = i + d_2;

		start = start_a > start_b ? start_a : start_b;
		end = end_a < end_b ? end_a : end_b;
		for( j = start; j <= end; j++)
		{
			line_force_w(&forward_value,&backward_value,p_obj_a->re,p_obj_b->re,line,j,from,middle1,middle2,to,cols,1);
			forward_force += forward_value;
			backward_force += backward_value;
		}
		histogram[pos] += (forward_force * factor);
		histogram[rev_pos] += (backward_force * factor);
	}
}

/*
Implementation of Wang's first algorithm with complexity O(KN).
It calls function force_histogram_w.
*/

double * forceHistogram_w(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	int rows, cols, max;
	double * histogram;
	int * line;
	BOUND_BOX * box_a, * box_b;

	rows = p_obj_a->rows;
	cols = p_obj_a->cols;
	
	histogram = (double *) calloc (p_d, sizeof(double));

	if((box_a = loadBoundBox(p_obj_a)) == NULL)
		return histogram;
	if((box_b = loadBoundBox(p_obj_b)) == NULL)
	{
		free(box_b);
		return histogram;
	}

	max = rows > cols ? rows:cols;

	line = (int *) calloc (max,sizeof(int));
	
	force_histogram_w(histogram, p_d, line, rows, cols, p_obj_a, box_a,p_obj_b,box_b);

	free(line);
	free(box_a);
	free(box_b);
	return histogram;
}

/*
Initiate an empry pyramid.
It contains (3n - 2) * n nodes (segments).
n: is the length of the side of an image.
factor is not used in the function.
*/
NODE * init_pyramid(int n, int factor)
{
	NODE * pyramid = NULL;
	pyramid = (NODE *) calloc((3*n - 2) * n, sizeof(NODE));
	return pyramid;
}
/*
clean all the values in a pyramid.
*/
void clean_pyramid(NODE * pyramid, int n)
{
	int length = (3*n-2) * n;
	int i;
	for(i=0; i<length; i++)
	{
		pyramid[i].count_a = pyramid[i].count_b = pyramid[i].force = 0.0;
	}
}

/*
free the memory of a node.
*/

void destroy_pyramid(NODE * pyramid)
{
	free(pyramid);
}
/*
generates an array contains all the decompositions of an line segment.
The decomposition is conducted recrusively.
factor means : the length of the line = pow(2,factor).
*/

SEG ** create_seg(int factor)
{
	int i,t2;
	SEG ** segment = (SEG**) malloc(factor * sizeof(SEG*));
	for(i=0;i<factor;i++)
	{
		t2 = 1<<(i+1);
		segment[i] = (SEG*) malloc(t2 * sizeof(SEG));
	}
	return segment;
}

/*
free the memory charged by the decomposition segment array.
*/

void destroy_seg(SEG ** segment, int factor)
{
	int i;
	for(i=0;i<factor;i++)
	{
		free(segment[i]);
	}
	free(segment);
}

/*
Decomposing all the line segments whose directions are in the range [0,pi/4]
the value factor = log(the width of the image).
*/

SEG ** build_seg(int factor)
{
	int h,t2,k;
	double real_right, real_left;
	double middle,tan_value;
	SEG ** segment = create_seg(factor);
	SEG * s_1;
	/*For each iteration, the line segment we are dealing with have length pow(2,h+1)*/
    for(h = 0; h < factor; h++)
	{
		t2 = 1 << (h+1);
		middle = (double)t2 / 2.0;                        /*Cut the line in half*/
		s_1 = segment[h];                                 /*Record the 1D array entry*/
		/*The number of segments with length = pow(2,h+1) can have k different directions (in [0, pi/4])*/
        for(k = 0; k < t2; k++)
		{
			tan_value = (double)k / (double)(t2);         /*Calculate the tan value for each possible direction.*/
			real_left = tan_value * (middle-1);           /*Calculate the end of the left half of the line*/
			real_right = tan_value * (middle);            /*Calculate the begin of the right half of the line*/

			s_1[k].end_left = (int)(real_left + 0.5);     /*Some adjustment of half pixel*/
			s_1[k].begin_right = (int)(real_right + 0.5);
		}
	}
	return segment;
}

/*
calculate the pyramid, each pyramid node records the number of pixels belonging to object A
the number of pixels belonging to object B and the force values on a segment.

py_low and py_high are the two empty arrays of nodes for computation.
segment is the array of decomposion segments generated through build_seg
n is the length of the image side.
factor is log(n).
*/

NODE * calculate_pyramid(NODE * py_low, NODE * py_high, SEG ** segment, int n, int factor)
{
	NODE * temp;
	int h,i,j,k;
	int seg, width_low, width_high, step_low, height_high, height_low;
	int index_left, index_right, begin_right;
	int count,height;
	int low_left, low_right;
	SEG * s_1;

	height = 3*n-2;
	seg = n;
	width_low = 1;

	for(h = 1; h <= factor; h++)                           /*for each level of line lenght = pow(2,h)*/
	{
		width_high = width_low << 1;                       /*width_hight is the length of lines in horizontal direction*/
		height_low = height - width_low + 1;               /*width_low is the length of half line*/
		height_high = height - width_high + 1;             /*height_low is the number of starting pixels along one side of the image considering the length = width_low
                                                           height_high is the number of starting pixels along one side of the image considering the length = width_high*/

		step_low = width_low * height_low;                 /*Considering half line length width_low, step_low indicates the total number of line segments
                                                           whose dirctions are in range [0, pi/4]*/

		seg = seg >> 1;                                    /*Number of cuts of the images along the horizontal direction*/

		s_1 = segment[h-1];                                /*Entry of 1D array of nodes of the left most cut of the image*/

		count = 0;

		for(i = 0; i < seg; i++)                           /*For each horizontal cuts of the image*/
		{
			low_left =  2 * i * step_low;                  /*the location of the nodes of the current cuts in the pyramid*/

			for(j = 0; j < height_high; j++)
			{
				low_right = low_left + step_low;           /*The location of the nodes containing the second half of segments*/

				/*For each possible segments on the cuts of image and with this length, compute the force value
                and update the number of pixels of A and the number of pixels of B*/
                for(k = 0; k < width_high; k++)
				{
					
					index_left = low_left + s_1[k].end_left;
					begin_right = s_1[k].begin_right;
					index_right = low_right + begin_right * (width_low - 1) + k  ;

					py_high[count].count_a = py_low[index_left].count_a +
												  py_low[index_right].count_a;

					py_high[count].count_b = py_low[index_left].count_b +
												  py_low[index_right].count_b;

					py_high[count].force = py_low[index_left].force +
												py_low[index_right].force + 
												py_low[index_right].count_b *
												py_low[index_left].count_a;
					count++;

				}
				low_left += width_low;

			}
		}
		width_low = width_high;
		temp = py_high;
		py_high = py_low;
		py_low = temp;
	}
	return py_low;
}


/*Assign the image into an empry pyramid initially*/
void assign_pyramid_1(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n) // 0 - pi/4
{
	int i,j;
	int p1 = n-1;
	int p2;
	int height = 3*n - 2;
	int rows, cols;
	rows = obj_a->rows;
	cols = obj_a->cols;

	for(i = 0; i < cols; i++)
	{
		p2 = p1;
		for(j=0; j < rows; j++)
		{	
			pyramid[p2].count_a = obj_a->re[j][i];
			pyramid[p2].count_b = obj_b->re[j][i];
			pyramid[p2].force = obj_a->re[j][i] * obj_b->re[j][i] / 2.0;
			p2++;
		}
		p1 += height;
	}
}
/*Rotate the image pi/4 anti-clockwise and assign the rotated image into an empry pyramid*/
void assign_pyramid_2(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n) // pi/4 - pi/2
{
	int i,j;
	int p1 = n-1;
	int p2;
	int height = 3*n - 2;
	int rows, cols;
	rows = obj_a->rows;
	cols = obj_a->cols;

	for(i = 0; i < rows; i++)
	{
		p2 = p1;
		for(j=0; j < cols; j++)
		{

			pyramid[p2].count_a = obj_a->re[i][j];
			pyramid[p2].count_b = obj_b->re[i][j];
			pyramid[p2].force = obj_a->re[i][j] * obj_b->re[i][j] / 2.0;
			p2++;
		}
		p1 += height;
	}
}
/*Rotate the image pi/2 anti-clockwise and assign the rotated image into an empry pyramid*/
void assign_pyramid_3(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n) // pi/2 - 3pi/4
{
	int i,j,t;
	int p1 = n-1;
	int p2;
	int height = 3*n - 2;
	int rows, cols;
	rows = obj_a->rows;
	cols = obj_a->cols;
	for(i = 0; i < rows; i++)
	{
		p2 = p1;
		for(j=0; j < cols; j++)
		{
			t = cols - j - 1;
			pyramid[p2].count_a = obj_a->re[i][t];
			pyramid[p2].count_b = obj_b->re[i][t];
			pyramid[p2].force = obj_a->re[i][t] * obj_b->re[i][t] / 2.0;
			p2++;
		}
		p1 += height;
	}
}

/*Rotate the image 3pi/4 anti-clockwise and assign the rotated image into an empry pyramid*/
void assign_pyramid_4(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n) // 3pi/4 - pi
{
	int i,j,t;
	int p1 = n-1;
	int p2;
	int height = 3*n - 2;
	int rows, cols;
	rows = obj_a->rows;
	cols = obj_a->cols;

	for(i = 0; i < cols; i++)
	{
		p2 = p1;
		for(j=0; j < rows; j++)
		{
			t = cols - i - 1;
			pyramid[p2].count_a = obj_a->re[j][t];
			pyramid[p2].count_b = obj_b->re[j][t];
			pyramid[p2].force = obj_a->re[j][t] * obj_b->re[j][t] / 2.0;
			p2++;
		}
		p1 += height;
	}
}

/*Rotate the image pi anti-clockwise and assign the rotated image into an empry pyramid*/
void assign_pyramid_5(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n) // pi - 5pi/4
{
	int i,j,t,k;
	int p1 = n-1;
	int p2;
	int height = 3*n - 2;
	int rows, cols;
	rows = obj_a->rows;
	cols = obj_a->cols;

	for(i = 0; i < cols; i++)
	{
		p2 = p1;
		for(j=0; j < rows; j++)
		{
			t = rows - j - 1;
			k = cols - i - 1;
			pyramid[p2].count_a = obj_a->re[t][k];
			pyramid[p2].count_b = obj_b->re[t][k];
			pyramid[p2].force = obj_a->re[t][k] * obj_b->re[t][k] / 2.0;
			p2++;
		}
		p1 += height;
	}
}

/*Rotate the image 5pi/4 anti-clockwise and assign the rotated image into an empry pyramid*/
void assign_pyramid_6(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n) // 5pi/4 - 3pi/2
{
	int i,j,t,k;
	int p1 = n-1;
	int p2;
	int height = 3*n - 2;
	int rows, cols;
	rows = obj_a->rows;
	cols = obj_a->cols;
	for(i = 0; i < rows; i++)
	{
		p2 = p1;
		for(j=0; j < cols; j++)
		{
			t = rows - i - 1;
			k = cols - j - 1;
			pyramid[p2].count_a = obj_a->re[t][k];
			pyramid[p2].count_b = obj_b->re[t][k];
			pyramid[p2].force = obj_a->re[t][k] * obj_b->re[t][k] / 2.0;
			p2++;
		}
		p1 += height;
	}
}

/*Rotate the image 3pi/2 anti-clockwise and assign the rotated image into an empry pyramid*/
void assign_pyramid_7(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n) // 3pi/2 - 7pi/4
{
	int i,j,t;
	int p1 = n-1;
	int p2;
	int height = 3*n - 2;
	int rows, cols;
	rows = obj_a->rows;
	cols = obj_a->cols;
	for(i = 0; i < rows; i++)
	{
		p2 = p1;
		for(j=0; j < cols; j++)
		{
			t = rows - i - 1;
			pyramid[p2].count_a = obj_a->re[t][j];
			pyramid[p2].count_b = obj_b->re[t][j];
			pyramid[p2].force = obj_a->re[t][j] * obj_b->re[t][j] / 2.0;
			p2++;
		}
		p1 += height;
	}
}

/*Rotate the image 7pi/4 anti-clockwise and assign the rotated image into an empry pyramid*/
void assign_pyramid_8(NODE * pyramid, RAS_OBJ * obj_a, RAS_OBJ * obj_b, int n) // 7pi/4 - 2pi
{
	int i,j,t;
	int p1 = n-1;
	int p2;
	int height = 3*n - 2;
	int rows, cols;
	rows = obj_a->rows;
	cols = obj_a->cols;
	for(i = 0; i < cols; i++)
	{
		p2 = p1;
		for(j=0; j < rows; j++)
		{
			t = rows - j - 1;
			pyramid[p2].count_a = obj_a->re[t][i];
			pyramid[p2].count_b = obj_b->re[t][i];
			pyramid[p2].force = obj_a->re[t][i] * obj_b->re[t][i] / 2.0;
			p2++;
		}
		p1 += height;
	}
}

/*
Computer the constant force between two raster object 
p_d is the number of reference directions considered
p_obj_a is the raster object A
p_obj_b is the raster object B
*/
double * forceHistogram_w_2(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	int i,j;
	int factor, size, length;
	int rows, cols, count;
	int d_1_8 = (int)((double)p_d / 8.0);
	int d_2_8 = (int)((double)p_d / 4.0);
	int d_3_8 = (int)(3.0 * (double)p_d / 8.0);
	int d_4_8 = (int)((double)p_d / 2.0);
	int d_5_8 = (int)(5.0 * (double)p_d / 8.0);
	int d_6_8 = (int)(3.0 * (double)p_d / 4.0);
	int d_7_8 = (int)(7.0 * (double)p_d / 8.0);
	int ext_size;
	double pi = acos(-1.0);
	double angle;
	double step;
	double force;
	int end, begin;
	NODE * py_high;
	NODE * py_low;
	NODE * py_result;
	SEG ** segment;

	double * his = (double *)calloc(p_d,sizeof(double));
	rows = p_obj_a->rows;                               /*height of the image */
	cols = p_obj_a->cols;                               /*width of the image */

	length = rows > cols? rows : cols;

	factor = factorBy2(length);
	size = 1 << factor;                                 /*make the image square and be the factor of 2*/
	ext_size = 2*size - 1;

	py_high = init_pyramid(size,factor);                /*initialize two pyramids for use*/
	py_low = init_pyramid(size,factor);
	segment = build_seg(factor);                        /*decompose of line segments*/


	step = 2.0 * pi / (double)p_d;
	count = 0;

	assign_pyramid_1(py_low,p_obj_a,p_obj_b,size);                     /*Assign the image into the lowest level of pyramid py_low*/
	py_result = calculate_pyramid(py_low,py_high,segment,size,factor); /*Compute the pyramid*/
	/*get the force values for directions in [0,pi/4] from the highest level of the result pyramid*/
    for(i = 0; i<= d_1_8; i++)
	{
		angle = (double) i * step;
		end = (int)(tan(angle) * (double)(size-1) +0.5);
		force = 0.0;

		begin = size - 1 - end;
		for(j = 0; j<ext_size; j++)
		{
			force += py_result[j * size + end].force;
		}

		his[i] = force / cos(angle);                               /*adjust the force value*/
	}
	clean_pyramid(py_low,size);

	assign_pyramid_2(py_low,p_obj_a,p_obj_b,size);                 /*Assign the rotated image (pi/4 rotation) into the lowest level of pyramid py_low*/
	py_result = calculate_pyramid(py_low,py_high,segment,size,factor);
	/*get the force values for directions in [pi/4,pi/2] from the highest level of the result pyramid*/
    for(i = d_1_8 + 1; i<= d_2_8; i++)
	{
		angle = (double) (d_2_8 - i) * step;
		end = (int)(tan(angle) * (double)(size-1) +0.5);
		force = 0.0;
		begin = size - 1 - end;
		for(j = 0; j<ext_size; j++)
		{
			force += py_result[j * size + end].force;
		}
		his[i] = force / cos(angle);                                /*adjust the force value*/
	}
	clean_pyramid(py_low,size);
	
	assign_pyramid_3(py_low,p_obj_a,p_obj_b,size);                           /*Assign the rotated image (pi/2 rotation) into the lowest level of pyramid py_low*/
	py_result = calculate_pyramid(py_low,py_high,segment,size,factor);
	/*get the force values for directions in [pi/2, 3pi/4] from the highest level of the result pyramid*/
    for(i = d_2_8 + 1; i<= d_3_8; i++)
	{
		angle = (double) (i - d_2_8) * step;
		end = (int)(tan(angle) * (double)(size-1) +0.5);
		force = 0.0;
		begin = size - 1 - end;
		for(j = 0; j<ext_size; j++)
		{
			force += py_result[j * size + end].force;
		}
		his[i] = force / cos(angle);                                /*adjust the force value*/
	}
	clean_pyramid(py_low,size);

	assign_pyramid_4(py_low,p_obj_a,p_obj_b,size);                           /*Assign the rotated image (3pi/4 rotation) into the lowest level of pyramid py_low*/
	py_result = calculate_pyramid(py_low,py_high,segment,size,factor);
	/*get the force values for directions in [3pi/4,pi] from the highest level of the result pyramid*/
    for(i = d_3_8 + 1; i<= d_4_8; i++)
	{
		angle = (double) (d_4_8 - i) * step;
		end = (int)(tan(angle) * (double)(size-1) +0.5);
		force = 0.0;
		begin = size - 1 - end;
		for(j = 0; j<ext_size; j++)
		{
			force += py_result[j * size + end].force;
		}
		his[i] = force / cos(angle);                                 /*adjust the force value*/
	}
	clean_pyramid(py_low,size);


	assign_pyramid_5(py_low,p_obj_a,p_obj_b,size);                          /*Assign the rotated image (pi rotation) into the lowest level of pyramid py_low*/
	py_result = calculate_pyramid(py_low,py_high,segment,size,factor);
	/*get the force values for directions in [pi, 5pi/4] from the highest level of the result pyramid*/
    for(i = d_4_8 + 1; i<= d_5_8; i++)
	{
		angle = (double) ( i - d_4_8) * step;
		end = (int)(tan(angle) * (double)(size-1) +0.5);
		force = 0.0;
		begin = size - 1 - end;
		for(j = 0; j<ext_size; j++)
		{
			force += py_result[j * size + end].force;
		}
		his[i] = force / cos(angle);                                      /*Adjust the force value*/
	}
	clean_pyramid(py_low,size);

	assign_pyramid_6(py_low,p_obj_a,p_obj_b,size);                        /*Assign the rotated image (5pi/4 rotation) into the lowest level of pyramid py_low*/
	py_result = calculate_pyramid(py_low,py_high,segment,size,factor);
	/*get the force values for directions in [5pi/4, 3pi/2] from the highest level of the result pyramid*/
    for(i = d_5_8 + 1; i<= d_6_8; i++)
	{
		angle = (double) (d_6_8 - i) * step;
		end = (int)(tan(angle) * (double)(size-1) +0.5);
		force = 0.0;
		begin = size - 1 - end;
		for(j = 0; j<ext_size; j++)
		{
			force += py_result[j * size + end].force;
		}
		his[i] = force / cos(angle);                                         /*adjust the force value*/
	}
	clean_pyramid(py_low,size);

	assign_pyramid_7(py_low,p_obj_a,p_obj_b,size);                           /*Assign the rotated image (3pi/2 rotation) into the lowest level of pyramid py_low*/
	py_result = calculate_pyramid(py_low,py_high,segment,size,factor);
	/*get the force values for directions in [3pi/2,7pi/4] from the highest level of the result pyramid*/
    for(i = d_6_8 + 1; i<= d_7_8; i++)
	{
		angle = (double) (i - d_6_8) * step;
		end = (int)(tan(angle) * (double)(size-1) +0.5);
		force = 0.0;
		begin = size - 1 - end;
		for(j = 0; j<ext_size; j++)
		{
			force += py_result[j * size + end].force;
		}
		his[i] = force / cos(angle);                                        /*adjust the force value*/
	}
	clean_pyramid(py_low,size);

	assign_pyramid_8(py_low,p_obj_a,p_obj_b,size);                          /*Assign the rotated image (7pi/4 rotation) into the lowest level of pyramid py_low*/
	py_result = calculate_pyramid(py_low,py_high,segment,size,factor);
	/*get the force values for directions in [7pi/4,2pi] from the highest level of the result pyramid*/
    for(i = d_7_8 + 1; i< p_d; i++)
	{
		angle = (double) (p_d - i) * step;
		end = (int)(tan(angle) * (double)(size-1) +0.5);
		force = 0.0;
		begin = size - 1 - end;
		for(j = 0; j<ext_size; j++)
		{
			force += py_result[j * size + end].force;
		}
		his[i] = force / cos(angle);                                      /*adjust the force value*/
	}
	destroy_pyramid(py_low);
	destroy_pyramid(py_high);
	destroy_seg(segment,factor);
	return his;

}



/*
int main(int argc, char * argv[])
{
	int p_d;
	int grey_A, grey_B, grey_A_B, m;
	FILE * fp = NULL;
	double * histogram;
	IMAGE * image;
	RAS_OBJ * obj_a;
	RAS_OBJ * obj_b;
	clock_t tv1, tv2;
	double time = 0.0;

	if(argc < 8)
	{
		printf("parameters: file_name A B A_B d m his_file\n");
		return 0;
	}

	image = loadImage(argv[1]);
	sscanf(argv[2],"%d", &grey_A);
	sscanf(argv[3],"%d", &grey_B);
	sscanf(argv[4],"%d", &grey_A_B);
	sscanf(argv[5],"%d", &p_d);
	sscanf(argv[6],"%d", &m);

	obj_a = image2RasObjExtract(image,grey_A,grey_A_B);
	obj_b = image2RasObjExtract(image,grey_B,grey_A_B);
	if(m>0)
	{
		randomFuzzify(obj_a, 100, m);
		randomFuzzify(obj_b, 200, m);
	}

	tv1 = clock();
	histogram = forceHistogram_w(p_d,obj_a,obj_b);
	//histogram = forceHistogram_w_2(p_d,obj_a,obj_b);

	tv2 = clock();
	time = (((double)(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0))/1000.0);

	printf("Wang and Makedon's algorithm_1 O(KN):\n");
	printf("Computational time (seconds) is: %f\n",time);
	
	fp = fopen(argv[7], "w");
	histogramWrite(fp,histogram,p_d);
	fclose(fp);
	destroyRasObj(obj_a);
	destroyRasObj(obj_b);
	destroyImage(image);
	free(histogram);
	return 1;
}*/

/*
The main function
*/
int main(int argc, char * argv[])
{
	int p_d;
	int  m;
	FILE * fp = NULL;
	double * histogram;
	IMAGE * image_a;
	IMAGE * image_b;
	RAS_OBJ * obj_a;
	RAS_OBJ * obj_b;
	clock_t tv1, tv2;
	double time = 0.0;

	if(argc < 6)
	{
		printf("\nWang's algorithms for computing constant force histograms (r=0) \n\n");
        printf("parameters: file_name_A file_name_B d first_or_second(1/2) histogram_file\n\n");
		printf("file_name_A: the pgm image of object A (reference object);\n");
		printf("file_name_B: the pgm image of object B (argument object);\n");
		printf("d: the number of reference directions considered (integer);\n");
		printf("first_or_second: the type of algorithms (1: Wang's first algorithm O(KN); 2: Wang's second algorithm O(NlogN));\n");
		printf("histogram_file: the txt file for storing the histogram values;\n");
		return 0;
	}

	image_a = loadImage(argv[1]);                                        /*Load PGM image containing object A*/
	image_b = loadImage(argv[2]);                                        /*Load PGM image containing object B*/

	sscanf(argv[3],"%d", &p_d);                                          /*The number of reference directions considered*/
	sscanf(argv[4],"%d", &m);                                            /*Which algorithm to use*/

	obj_a = image2RasObj(image_a);                                       /*Load raster object A from Image*/
	obj_b = image2RasObj(image_b);                                       /*Load raster object B from Image*/
	
	if(m == 1)                                                           /*Use Wang's first algorithm with complexity O(KN)*/
	{
       tv1 = clock();
	   histogram = forceHistogram_w(p_d,obj_a,obj_b);
       tv2 = clock();
	   time = (((double)(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0))/1000.0);
	   printf("Wang and Makedon's algorithm_1 O(KN):\n");
	}
	else                                                                 /*Use Wang's second algorithm with complexity O(NlogN)*/
	{
        tv1 = clock();
        histogram = forceHistogram_w_2(p_d,obj_a,obj_b);
	    tv2 = clock();
	    time = (((double)(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0))/1000.0);
	    printf("Wang and Makedon's algorithm_2 O(NlogN):\n");     
    }
	printf("Computational time (seconds) is: %f\n",time);
	
	fp = fopen(argv[5], "w");
	histogramWrite(fp,histogram,p_d);                                    /*Write the histogram into a txt file*/
	fclose(fp);
	destroyRasObj(obj_a);
	destroyRasObj(obj_b);
	destroyImage(image_a);
	destroyImage(image_b);
	free(histogram);
	return 1;
}



