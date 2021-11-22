/*
Define the data structure for raster objects and implement some operation functions
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "objfunc.h"

/*
Functions loadAngleObject, image2AngleObject and destroyAngleObject are used specificly for 
generating angle histogram, which is not included in the force histogram project.
*/

/*
A_OBJ * loadAngleObject(RAS_OBJ * p_obj)
{
	int i, j, rows, cols;
	int size;
	double ** re;
	A_OBJ * angle_obj;
	double * value;
	int * y, * x;
	int count = 0;

	rows = p_obj->rows;
	cols = p_obj->cols;
	re = p_obj->re;
	
	size = rows * cols;

	angle_obj = (A_OBJ *) calloc (1, sizeof(A_OBJ));
	value = (double *) calloc (size, sizeof(double));
	y = (int *) calloc (size, sizeof(int));
	x = (int *) calloc (size, sizeof(int));

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			if(re[i][j] > 0.0)
			{
				value[count] = re[i][j];
				y[count] = i;
				x[count] = j;
				count ++;
			}
		}
	}
	angle_obj->value = value;
	angle_obj->y = y;
	angle_obj->x = x;
	angle_obj->n = count;
	return angle_obj;
}

A_OBJ * image2AngleObject(IMAGE * image)
{
	int i, j, rows , cols;
	double maxval;
	int size;
	unsigned char ** pixels = NULL;
	double * value;
	int * y, * x;
	int count = 0;
	A_OBJ * object = NULL;

	cols = image->width;
	rows = image->height;
	size = cols * rows;
	maxval = (double)(image->max_grey);
	pixels = image->pixels;

	object = (A_OBJ *)malloc(sizeof(A_OBJ));
	value = (double *) calloc (size, sizeof(double));
	y = (int *) calloc (size, sizeof(int));
	x = (int *) calloc (size, sizeof(int));

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j< cols; j++)
		{
			if(pixels[i][j] > 0.0)
			{
				value[count] = ((double) (pixels[i][j])) / maxval;
				y[count] = rows - i - 1;
				x[count] = j;
				count++;
			}
		}
	}
	object->n = count;
	object->value = value;
	object->x = x;
	object->y = y;
	return object;	

}

void destroyAngleObject(A_OBJ * p_obj)
{
	free(p_obj->value);
	free(p_obj->x);
	free(p_obj->y);
	free(p_obj);
}
*/

/*
Create an empty raster object, only initiate the real array
*/

RAS_OBJ * createRasObj(int p_rows, int p_cols, int p_org_y, int p_org_x)
{
	RAS_OBJ * object = (RAS_OBJ *) calloc (1,sizeof(RAS_OBJ));

	object->rows = p_rows;
	object->cols = p_cols;
	object->org_y = p_org_y;
	object->org_x = p_org_x;

	object->re = createMatrix(p_rows, p_cols);
	object->im = NULL;
	return object;
}
/*
Create an empty raster object, initiate both the real and the imaginary arrays.
*/
RAS_OBJ * createRasObjComplex(int p_rows, int p_cols, int p_org_y, int p_org_x)
{
	RAS_OBJ * object = createRasObj(p_rows,p_cols, p_org_y, p_org_x);
	complexRasObj(object);
	return object;
}

/*
create a 2D array with size p_cols by p_rows
*/

double ** createMatrix(int p_rows, int p_cols)
{
	int i;
	double ** data = (double **) calloc (p_rows,sizeof(double*));
	for( i = 0; i < p_rows; i++)
		data[i] = (double *) calloc (p_cols, sizeof(double));
	return data;
}

/*
Load a (possibly) fuzzy object from p_image, each membership value is calculated by
dividing the grey level by the max grey level of the image (usually is 255).
*/

RAS_OBJ * image2RasObj(IMAGE * p_image)
{
	int i, j, rows , cols;
	double maxval;
	unsigned char ** pixels = NULL;
	double ** re = NULL;
	RAS_OBJ * object = NULL;

	cols = p_image->width;
	rows = p_image->height;
	maxval = (double)(p_image->max_grey);
	pixels = p_image->pixels;

	object = createRasObj(rows, cols, 0, 0);

	re = object->re;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j< cols; j++)
		{
			re[rows - i - 1][j] = ((double) (pixels[i][j])) / maxval;
		}
	}
	return object;
}

/*
Assume p_image contains two crisp objects. One object is painted in grey level
p_mark and the joint part is painted in grey level p_mark_joint. This function
extract one object (with grey level p_mark) from p_image.
*/

RAS_OBJ * image2RasObjExtract(IMAGE * p_image, int p_mark, int p_mark_joint)
{
	int i, j, rows , cols;
	double maxval;
	unsigned char ** pixels = NULL;
	double ** re = NULL;
	RAS_OBJ * object = NULL;

	cols = p_image->width;
	rows = p_image->height;
	maxval = (double)(p_image->max_grey);
	pixels = p_image->pixels;

	
	object = createRasObj(rows, cols, 0, 0);
	

	re = object->re;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			if(pixels[i][j] == p_mark || pixels[i][j] == p_mark_joint)
				re[rows - i - 1][j] = 1.0;
		}
	}
	return object;
}

/*
Copy a raster object.
If p_only_re=1: only the real array will be copied.
If p_only_re=0: both the real and the imaginary arrays will be copied.
*/

RAS_OBJ * copyRasObj(RAS_OBJ * p_object, int p_only_re)
{
	int rows, cols;
	RAS_OBJ * tar_object = NULL;
	double ** org_data = NULL;
	double ** tar_data = NULL;

	rows = p_object->rows;
	cols = p_object->cols;
	org_data = p_object->re;

	tar_object = createRasObj(rows,cols, p_object->org_y, p_object->org_x);

	tar_data = tar_object->re;

	copyMatrix(org_data, tar_data, rows, cols);

	org_data = p_object->im;

	if(org_data == NULL || p_only_re)
		return tar_object;
	
	complexRasObj(tar_object);
	tar_data = tar_object->im;

	copyMatrix(org_data, tar_data, rows, cols);

	return tar_object;
}
/*
Copy a matrix p_org to another matrix p_tar with the same dimension.
*/

int copyMatrix(double ** p_org, double ** p_tar, int p_rows, int p_cols)
{
	int i,j;
	for(i = 0; i < p_rows; i++)
	{
		for(j = 0; j < p_cols; j++)
		{
			p_tar[i][j] = p_org[i][j];
		}
	}
	return 1;
}
/*
Extend the region of a raster object four times bigger and set the origin to the middle of
the extended range.
*/
int extendRasObj(RAS_OBJ * p_object)
{
	int rows, cols, t_rows, t_cols, org_y, org_x;
	double ** org_data = NULL;
	double ** tar_data = NULL;

	rows = p_object->rows;
	cols = p_object->cols;
	org_data = p_object->re;

	t_rows = rows * 2 - 1;
	t_cols = cols * 2 - 1;
	org_y = rows - 1;
	org_x = cols - 1;

	p_object->rows = t_rows;
	p_object->cols = t_cols;
	p_object->org_y = org_y;
	p_object->org_x = org_x;
	p_object->re = extendMatrix(org_data,rows, cols, t_rows, t_cols, org_y, org_x);

	org_data = p_object->im;
	if(org_data == NULL)
		return 1;
	
	p_object->im = extendMatrix(org_data,rows, cols, t_rows, t_cols, org_y, org_x);;
	return 1;
}

/*
Extend a certain region (with original size = bo_rows by bo_cols with the position of low-left corner = (region_x,region_y)) 
to a certain size. After calling the function, the region will have size = p_t_rows by p_t_cols, 
and have the position of the origin at (p_org_x, p_org_y).
*/

RAS_OBJ * resizeRasObj(RAS_OBJ * p_object, int p_t_rows, int p_t_cols, int p_org_y, int p_org_x, int region_y, int region_x, int bo_rows, int bo_cols)
{
	int rows, cols, i, j;
	RAS_OBJ * tar_object = NULL;
	double ** tar_re = NULL;
	double ** org_re = NULL;
	
	rows = bo_rows;
	cols = bo_cols;
	org_re = p_object->re;

	if((p_t_rows - p_org_y) < rows || (p_t_cols - p_org_x) < cols )
		return NULL;

	tar_object = createRasObj(p_t_rows, p_t_cols, p_org_y, p_org_x);
	tar_re = tar_object->re;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			tar_re[p_org_y + i][p_org_x + j] = org_re[i + region_y][j + region_x];
		}
	}
	return tar_object;
}

/*
Extend the object region to a certain size. The raster object is directed loaded from a Image using 
function image2RasObj or image2RasObjExtract without being extened or resized yet.
After calling the function, the new size of the object region will be p_t_rows by p_t_cols and the new
position of the origin will be at (p_org_x, p_org_y).
During the extension, the function will automatically adjust the minimum bounding rectangle of the object 
in the new region and save the new bounding rectangle in p_box.
*/


RAS_OBJ * resizeRasObj_boundBox(RAS_OBJ * p_object, BOUND_BOX * p_box, int p_t_rows, int p_t_cols, int p_org_y, int p_org_x)
{
	RAS_OBJ * tar_object = resizeRasObj(p_object, p_t_rows, p_t_cols, p_org_y, p_org_x,0,0, 
		                                p_object->rows, p_object->cols);
	p_box->v1[0] += p_org_y; p_box->v1[1] += p_org_x;
	p_box->v3[0] += p_org_y; p_box->v3[1] += p_org_x;
	p_box->v2[0] = p_box->v3[0]; p_box->v2[1] = p_box->v1[1];
	p_box->v4[0] = p_box->v1[0]; p_box->v4[1] = p_box->v3[1];
	return tar_object;
}

/*
Extend the region of object's minimum bounding rectangle.
The low-left corner of the region is (p_org_x, p_org_y) and the size of the region is bound_rows by bound_cols
After extension, the target size will be p_t_rows by p_t_cols and the target position of the origin will be at
(p_org_x, p_org_y).
At the same time, the function adjusts the minimum bounding rectangle of the raster object in the extended region.
*/

RAS_OBJ * resizeRasObj_boundBox_region(RAS_OBJ * p_object, BOUND_BOX * p_box, int p_t_rows, int p_t_cols, int p_org_y, int p_org_x, 
									   int bound_rows, int bound_cols)
{
	RAS_OBJ * tar_object = resizeRasObj(p_object, p_t_rows, p_t_cols, p_org_y, p_org_x, p_box->v2[0], p_box->v2[1],
		                                bound_rows, bound_cols);
	p_box->v1[0] += (p_org_y - p_box->v2[0]); p_box->v1[1] += (p_org_x - p_box->v2[1]);
	p_box->v3[0] += (p_org_y - p_box->v2[0]); p_box->v3[1] += (p_org_x - p_box->v2[1]);
	p_box->v2[0] = p_box->v3[0]; p_box->v2[1] = p_box->v1[1];
	p_box->v4[0] = p_box->v1[0]; p_box->v4[1] = p_box->v3[1];
	return tar_object;
}

/*
Extend the size of a 2D array. The original size of the array is p_rows by p_cols.
After extension the target size is p_t_rows by p_t_cols and the positon of the origin will be
at (p_org_x, p_org_y).
*/
double ** extendMatrix(double ** p_org, int p_rows, int p_cols, int p_t_rows, int p_t_cols, int p_org_y, int p_org_x)
{
	int i, j;
	double ** tar = NULL;
	tar = (double **) calloc (p_t_rows,sizeof(double *));

	for(i = 0; i < p_t_rows; i++)
		tar[i] = (double *) calloc (p_t_cols,sizeof(double));	

	for(i = 0; i < p_rows; i++)
	{
		for(j = 0; j < p_cols; j++)
		{
			tar[i + p_org_y][j + p_org_x] = p_org[i][j];
		}
	}
	for(i = 0; i < p_rows; i++)
		free(p_org[i]);
	free(p_org);
	return tar;
}

/*
In the region, reflect the raster object about the origin.
*/

int reflectRasObj(RAS_OBJ * p_object)
{
	int rows, cols, org_y, org_x;
	double ** data = NULL;
	rows = p_object->rows;
	cols = p_object->cols;
	org_y = p_object->org_y;
	org_x = p_object->org_x;

	
	if(org_y == 0 || org_x == 0)
		return 0;


	data = p_object->re;
	reflectMatrix(data, rows, cols, org_y, org_x);
		
	data = p_object->im;

	if(data == NULL)
		return 1;
	reflectMatrix(data, rows, cols, org_y, org_x);
	return 1;
}

/*
In the region, reflect the raster object about the origin, at the same time,
adjust the minimum bounding rectangle.
*/

int reflectRasObj_boundBox(RAS_OBJ * p_object, BOUND_BOX * p_box)
{
	int rows, cols, org_y, org_x, s_i, e_i, s_j, e_j, i, j, r_i, r_j;
	int org_y_2, org_x_2;
	double temp;
	double ** data = NULL;
	rows = p_object->rows;
	cols = p_object->cols;
	org_y = p_object->org_y;
	org_x = p_object->org_x;
	org_y_2 = org_y * 2;
	org_x_2 = org_x * 2;

	s_i = p_box->v3[0];
	e_i = p_box->v1[0];
	s_j = p_box->v1[1];
	e_j = p_box->v3[1];

	if(org_y == 0 || org_x == 0)
		return 0;

	
	data = p_object->re;
	temp = data[org_y][org_x];

	for(i = s_i; i <= e_i; i++)
	{
		for(j = s_j; j<= e_j; j++)
		{
			r_i = org_y_2 - i;
			r_j = org_x_2 - j;
			data[r_i][r_j] = data[i][j];
			data[i][j] = 0;
		}
	}
	data[org_y][org_x] = temp;
	s_i = org_y_2 - s_i;
	e_i = org_y_2 - e_i;
	s_j = org_x_2 - s_j;
	e_j = org_x_2 - e_j;

	p_box->v1[0] = s_i; p_box->v1[1] = e_j;
	p_box->v3[0] = e_i; p_box->v3[1] = s_j;
	p_box->v2[0] = p_box->v3[0]; p_box->v2[1] = p_box->v1[1];
	p_box->v4[0] = p_box->v1[0]; p_box->v4[1] = p_box->v3[1];

	return 1;

}

/*
Reflect a 2D matrix with size p_rows by p_cols about the origin locating at (p_org_x, p_org_y).
*/

int reflectMatrix(double ** p_matrix, int p_rows, int p_cols, int p_org_y, int p_org_x)
{
	int b_low_y, b_high_y, b_low_x, b_high_x, real_cols;
	int i, j, r_i, r_j;
	double t;
	if(p_org_y <= p_rows/2)
	{
		b_low_y = 0;
		b_high_y = p_org_y;
	}
	else
	{
		b_low_y = p_org_y + 1;
		b_high_y = p_rows;
	}

	if(p_org_x <= p_cols/2)
	{
		b_low_x = 0;
		b_high_x = p_org_x;
	}
	else
	{
		b_low_x = p_org_x + 1;
		b_high_x = p_cols;
	}
	real_cols = (b_high_x - b_low_x) * 2 +1;

	for(i = b_low_y; i < b_high_y; i++)
	{
		for(j = 0; j < real_cols; j++)
		{
			t = p_matrix[i][j];
			r_i = 2*p_org_y - i;
			r_j = 2*p_org_x - j;
			p_matrix[i][j] = p_matrix[r_i][r_j];
			p_matrix[r_i][r_j] = t;
		}
	}
	for(j = b_low_x; j < b_high_x; j++)
	{
		t = p_matrix[p_org_y][j];
		r_j = 2*p_org_x - j;
		p_matrix[p_org_y][j] = p_matrix[p_org_y][r_j];
		p_matrix[p_org_y][r_j] = t;
	}
	return 1;
}

/*
Initiate a raster object's 2D imaginary array.
*/

int complexRasObj(RAS_OBJ * p_object)
{
	p_object->im = createMatrix(p_object->rows, p_object->cols);
	return 1;
}

/*
For each complex value recorded in the arrays re and im, compute the magnitude as
sqrt(re*re + im*im).
*/

int magnitudeRasObj(RAS_OBJ * p_object)
{
	int rows, cols, i,j;
	double ** re = NULL;
	double ** im = NULL;

	rows = p_object->rows;
	cols = p_object->cols;
	re = p_object->re;
	im = p_object->im;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			re[i][j] = sqrt(pow(re[i][j], 2) + pow(im[i][j], 2));
		}
	}
	return 1;
}

/*
Normalize the re array, bring all the value inside into the range [0,1].
*/

int normalizRasObj(RAS_OBJ * p_object)
{
	int rows, cols, i, j;
	double ** re = NULL;
	double max_val = 0;

	rows = p_object->rows;
	cols = p_object->cols;
	re = p_object->re;
	
	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			if(re[i][j] > max_val)
				max_val = re[i][j];
		}
	}
	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			re[i][j] = re[i][j] / max_val;
		}
	}
	return 1;
}

/*
Write a raster object back to a Image.
*/

IMAGE * rasObj2Image(RAS_OBJ * p_object, int p_maxval)
{
	IMAGE * image = NULL;
	double ** re = NULL;
	unsigned char ** pixels = NULL;
	int rows, cols, i, j;
	double t;

	rows = p_object->rows;
	cols = p_object->cols;
	re = p_object->re;

	image = createImage(cols,rows,p_maxval);

	image->max_grey = p_maxval;
	pixels = image->pixels;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			t = re[i][j] * (double)p_maxval;
			
			if(t > p_maxval)
				t = (double) p_maxval;
			if(t < 0)
				t = 0.0;
				
			pixels[rows - i - 1][j] = (unsigned char) t;
		}
	}
	return image;
}

/*
Release the memory charged by a raster object.
*/

void destroyRasObj(RAS_OBJ * p_object)
{
	int rows;
	double ** data;
	rows = p_object->rows;
	data = p_object->re;

	destroyMatrix(data, rows);

	data = p_object->im;

	if(data != NULL)
		destroyMatrix(data, rows);
	free(p_object);
}

/*
Release the memory of a 2D array.
*/

void destroyMatrix(double ** p_matrix, int p_rows)
{
	int i;
	for(i = 0; i < p_rows; i++)
		free(p_matrix[i]);
	free(p_matrix);
}

/*
Create an empry vector for storing complex values.
*/

COMP_VEC * createCompVec(int p_length)
{
	COMP_VEC * vector = (COMP_VEC *) calloc (1, sizeof(COMP_VEC));
	vector->length = p_length;
	vector->re = (double *) calloc (p_length, sizeof(double));
	vector->im = (double *) calloc (p_length, sizeof(double));
	return vector;
}

/*Copy a vector in with all the values are complex numbers.*/

COMP_VEC * copyCompeVec(COMP_VEC * p_vector)
{
	int length;
	COMP_VEC * tar_vector = NULL;
	double * org = NULL;
	double * tar = NULL;

	length = p_vector->length;
	tar_vector = createCompVec(length);
	org = p_vector->re;
	tar = tar_vector->re;
	copyVec(org, tar, length);

	org = p_vector->im;
	tar = tar_vector->im;
	copyVec(org, tar, length);
	return tar_vector;
}

/*
Copy a 1D array.
*/

int copyVec(double * p_org, double * p_tar, int p_length)
{
	int i;
	for(i = 0; i < p_length; i++)
	{
		p_tar[i] = p_org[i];
	}
	return 1;
}

/*
Copy a 1D array p_vector to a specific row of a 2D array p_matrix.
And the position of the 1D array in p_matrix will be at row=p_row, and
in that row, the starting and ending positions are defined by parameters
start and end.
*/

int vec2Row(double * p_vector, int p_row, double ** p_matrix, int start, int end)
{
	int i;
	for(i = start; i < end; i++)
	{
		p_matrix[p_row][i] = p_vector[i];
	}
	return 1;
}

/*
Copy a 1D array p_vector to a specific column of a 2D array p_matrix.
Also in the column, the position of p_vector is defined by parameters start and end.
*/

int vec2Col(double * p_vector, int p_col, double ** p_matrix, int start, int end)
{
	int j;
	for(j = start; j < end; j++)
	{
		p_matrix[j][p_col] = p_vector[j];
	}
	return 1;
}

/*
Copy a specific part of a row (defined by p_row, start and end) of a 2D matrix p_matrix to
a 1D array.
*/

int row2Vec(double ** p_matrix, int p_row, double * p_vector, int start, int end)
{
	int i;
	for(i = start; i < end; i++)
	{
		p_vector[i] = p_matrix[p_row][i];
	}
	return 1;
}

/*
Copy a specific part of a column (defined by p_row, start and end) of a 2D matrix p_matrix to
a 1D array.
*/

int col2Vec(double ** p_matrix, int p_col, double * p_vector, int start, int end)
{
	int j;
	for(j = start; j < end; j++)
	{
		p_vector[j] = p_matrix[j][p_col];
	}
	return 1;
}

/*
Release the memory charged by a 1D complex vector.
*/

void destroyCompVec(COMP_VEC * p_vector)
{
	free(p_vector->re);
	free(p_vector->im);
	free(p_vector);
}

/*
Calculate the minimum bounding rectangle of a raster object.
The function called boundBox.
*/

BOUND_BOX * loadBoundBox(RAS_OBJ * p_object)
{
	
	BOUND_BOX * box = (BOUND_BOX *) malloc (sizeof(BOUND_BOX));
	if(!boundBox(box,p_object))
	{
		free(box);
		return NULL;
	}
	return box;
}

/*
Calculate the minimum bounding rectangle of a raster object.
*/

int boundBox(BOUND_BOX * p_box, RAS_OBJ * p_object)
{
	int i, j, rows, cols;
	int min_i, min_j, max_i, max_j;
	double ** re = NULL;

	rows = p_object->rows;
	cols = p_object->cols;
	re = p_object->re;

	min_i = rows - 1;
	min_j = cols - 1;
	max_i = 0;
	max_j = 0;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			if(re[i][j] > 0.0)
			{
				if(i < min_i)
					min_i = i;
				if(j < min_j)
					min_j = j;
				if(i > max_i)
					max_i = i;
				if(j > max_j)
					max_j = j;
			}
		}
	}
	if(max_i >= min_i && max_j >= min_j)
	{
		p_box->v1[Y] = max_i; p_box->v1[X] = min_j;
		p_box->v2[Y] = min_i; p_box->v2[X] = min_j;
		p_box->v3[Y] = min_i; p_box->v3[X] = max_j;
		p_box->v4[Y] = max_i; p_box->v4[X] = max_j;
		return 1;
	}
	return 0;

}

/*
Given a crisp raster object, the function randomly assign each pixel belonging to the object
a membership degree in range[0,1].
seed: is the random seed set by user.
p_number: indicates how many possible membership degrees.
For example, if p_number=5, the function will assign a pixel (belonging to the object) a membership
degree randomly picked from {0.2, 0.4, 0.6, 0.8, 1.0}.
*/

int randomFuzzify(RAS_OBJ * p_object, int seed, int p_number)
{
	int i,j, t;
	double * values;
	int rows, cols;
	double ** re;
	double step = 1.0 / (double)p_number;
	srand(seed);
	values = (double *) calloc (p_number,sizeof(double));
	for (i = 0;i<p_number; i++)
	{
		values[i] = 1.0 - (double)i * step;
	}
	rows = p_object->rows;
	cols = p_object->cols;
	re = p_object->re;

	for( i = 0; i < rows; i++)
	{
		for( j = 0; j < cols; j++)
		{
			if(re[i][j] > 0)
			{
				t = rand() % p_number;
				re[i][j] = values[t];
			}

		}
	}
	free(values);
	return 1;
}

/*
Create an empty list of alpha cuts.
*/

A_LIST * createAList()
{
	A_LIST * list = (A_LIST *) malloc (sizeof(A_LIST));
	list->head = NULL;
	list->n = 0;
	return list;
}

/*
Insert a value into the list and make the list ordered.
Each value only appears in the list once.
*/

int insertAList(A_LIST * p_list, double p_f)
{
	ALPHA * a = p_list->head;
	ALPHA * n_alpha = NULL;
	if(a == NULL)
	{
		n_alpha = (ALPHA *) malloc (sizeof(ALPHA));
		n_alpha->f = p_f;
		n_alpha->next = NULL;
		p_list->head = n_alpha;
		p_list->n ++;
		return 1;
	}
	if(p_f == a->f)
		return 0;
	if(p_f < a->f)
	{
		n_alpha = (ALPHA *) malloc (sizeof(ALPHA));
		n_alpha->f = p_f;
		n_alpha->next = a;
		p_list->head = n_alpha;
		p_list->n ++;
		return 1;
	}
	while(a->next != NULL)
	{
		if(p_f == a->next->f)
			return 0;
		if(p_f < a->next->f)
		{
			n_alpha = (ALPHA *) malloc (sizeof(ALPHA));
			n_alpha->f = p_f;
			n_alpha->next = a->next;
			a->next = n_alpha;
			p_list->n ++;
			return 1;
		}
		a = a->next;
	}
	n_alpha = (ALPHA *) malloc (sizeof(ALPHA));
	n_alpha->f = p_f;
	n_alpha->next = NULL;
	a->next = n_alpha;
	p_list->n ++;
	return 1;
}

/*
Print the list.
*/

void printfAList(A_LIST * p_list)
{
	ALPHA * a = p_list->head;
	while(a!=NULL)
	{
		printf("%f -> ", a->f);
		a = a->next;
	}
	printf("\n");
}

/*
Free the memory of the list.
*/

void destroyAList(A_LIST * p_list)
{
	ALPHA * a = p_list->head;
	while(a!=NULL)
	{
		p_list->head = a->next;
		free(a);
		a = p_list->head;
	}
}

/*
Build a list from a raster objects. The list stores all the possible
membership degrees appearing in the raster object.
*/

int detectAlphas(RAS_OBJ * p_object, A_LIST * p_list)
{
	int i,j;
	int min_i, max_i, min_j, max_j;
	BOUND_BOX * box = loadBoundBox(p_object);
	double ** re;

	re = p_object->re;
	min_i = box->v3[0]; max_i = box->v1[0];
	min_j = box->v1[1]; max_j = box->v3[1];

	for(i = min_i; i<=max_i; i++)
	{
		for(j = min_j; j<=max_j; j++)
		{
			if(re[i][j] > 0.0)
				insertAList(p_list,re[i][j]);
		}
	}
	return p_list->n;
}

/*
Create the crisp object p_cut from the fuzzy raster object p_object.
p_cut is the alpha cut (with alpha value = p_f) of p_object.
*/

int getAlphaCut(RAS_OBJ * p_object, RAS_OBJ * p_cut, double p_f)
{
	int i,j,rows,cols;
	double ** re_obj, ** re_cut;

	re_cut = p_cut->re;
	re_obj = p_object->re;

	rows = p_cut->rows;
	cols = p_cut->cols;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
		{
			if(re_obj[i][j] >= p_f)
				re_cut[i][j] = 1.0;
			else
				re_cut[i][j] = 0.0;
		}
	}
	return 1;
}

/*
Write a histogram into a txt file.
p_d is the number of values in the histogram.
*/

int histogramWrite(FILE * p_fp, double * p_histogram, int p_d)
{
	int i;
	fprintf(p_fp,"%d\n",p_d);
	for(i = 0; i < p_d; i++)
	{
		fprintf(p_fp, "%.10f\n",p_histogram[i]);
	}
	return 1;
}

/*
Swap two values a and b.
*/

void swap(int * a, int *b)
{
	int temp;
	temp = *a;
	*a = *b;
	*b = temp;
}

/*
The official Bresenham's line drawing algorithm.
Draw a rastrized line from vertex start to vertex end.
*/

void bresenham(unsigned char ** pixels, VERTEX * start, VERTEX * end, unsigned char color)
{
	int x0, x1, y0, y1, x, y;
	int steep, deltax, deltay, ystep;
	double error, deltaerr;

	x0 = start->x;
	y0 = start->y;
	x1 = end->x;
	y1 = end->y;

	if(abs(y1-y0) > abs(x1-x0))
		steep = 1;
	else
		steep = 0;
	if(steep)
	{
		swap(&x0,&y0);
		swap(&x1,&y1);
	}
	if(x0 > x1)
	{
		swap(&x0, &x1);
		swap(&y0, &y1);
	}
	deltax = x1-x0;
	deltay = abs(y1-y0);
	error = 0.0;
	deltaerr = (double)deltay / (double)deltax;
	y = y0;
	if(y0 < y1)
		ystep = 1;
	else
		ystep = -1;
	for(x = x0; x <= x1; x++)
	{
		if(steep)
			pixels[x][y] = color;
		else
			pixels[y][x] = color;
		error = error + deltaerr;
		if(error >= 0.5)
		{
			y = y + ystep;
			error = error - 1.0;
		}
	}

}

/*
Initialization for creating a koch snow flake.
It actually draws a triangle with the lengths of all three sizes = n.
*/

VERTEX * init_koch_snowflake(int n)
{
	VERTEX * list = NULL;
	VERTEX * p;
	int side = (int)((double)n * 3.0 / (2.0 * sqrt(3.0)));
	int dis = (n-side) / 2;
	int height = (int)((double)side * sqrt(3.0) / 2.0);
	printf("height = %d\n",height);
	list = (VERTEX *) malloc (sizeof(VERTEX));
	p = list;
	p->x = dis + side-1;
	p->y = (n-1) - height + 1;
	p->next = (VERTEX *) malloc(sizeof(VERTEX));
	p = p->next;
	p->x = dis;
	p->y = (n-1) - height + 1;
	p->next = (VERTEX *) malloc(sizeof(VERTEX));
	p = p->next;
	p->x = dis + side/2 - 1;
	p->y = n-1;
	p->next = NULL;
	return list;
}

/*
Compute the angle of the line linking vertex start and vertex end.
*/

double angle_calc(VERTEX * start, VERTEX * end)
{
	int x0, x1, y0, y1;
	double pi = acos(-1);
	double ang = pi/6.0;
	x0 = start->x;
	y0 = start->y;
	x1 = end->x;
	y1 = end->y;
	if(x1 == x0)
	{
		if(y1 > y0)
			return (pi / 4.0 + ang);
		else
			return (3 * pi / 4.0 + ang);
	}
	else
	{
		if(x1 > x0)
			return ( atan((double)(y1 - y0) / (double)(x1 - x0)) + ang );
		else
			return ( atan((double)(y1 - y0) / (double)(x1 - x0)) + ang + pi );
	}
}

/*According to Koch snow flake, interpolate each side to make more vertexes*/

VERTEX * interpolation_koch(VERTEX * start, VERTEX * end)
{
	double angle;
	double length;
	VERTEX * v = (VERTEX *) malloc (sizeof(VERTEX));
	v->next = NULL;

	angle = angle_calc(start,end);
	length = sqrt( pow(start->x - end->x + 1, 2.0) + pow(start->y - end->y + 1, 2.0) ) / sqrt(3);
	v->x = (int)(length * cos(angle)) + start->x;
	v->y = (int)(length * sin(angle)) + start->y;
	return v;
}

/*
Generate a koch snow flake with iter = the number of iterations.
*/

void iteration_koch_snowflake(int iter, VERTEX * list)
{
	int i;
	VERTEX * start, * end;
	VERTEX * p = NULL;
	VERTEX * v = NULL;
	double step_x, step_y;
	for(i=0;i<iter;i++)
	{
		p = list;
		while(p->next != NULL)
		{
			start = p;
			end = p->next;
			
			step_x = (double)(end->x - start->x + 1) / 3.0;
			step_y = (double)(end->y - start->y + 1) / 3.0;
			
			
			v = (VERTEX *) malloc (sizeof(VERTEX));
			v->x = (int)step_x + start->x;
			v->y = (int)step_y + start->y;

			v->next = p->next;
			p->next = v;
			p = p->next;

			v = interpolation_koch(start,end);
			v->next = p->next;
			p->next = v;
			p = p->next;

			v = (VERTEX *) malloc (sizeof(VERTEX));
			v->x = (int)(2.0 * step_x) + start->x;
			v->y = (int)(2.0 * step_y) + start->y;
			
			v->next = p->next;
			p->next = v;
			p = p->next->next;	
		}
		start = p;
		end = list;
			
		step_x = (double)(end->x - start->x + 1) / 3.0;
		step_y = (double)(end->y - start->y + 1) / 3.0;
			
			
		v = (VERTEX *) malloc (sizeof(VERTEX));
		v->x = (int)step_x + start->x;
		v->y = (int)step_y + start->y;

		v->next = p->next;
		p->next = v;
		p = p->next;

		v = interpolation_koch(start,end);
		v->next = p->next;
		p->next = v;
		p = p->next;

		v = (VERTEX *) malloc (sizeof(VERTEX));
		v->x = (int)(2.0 * step_x) + start->x;
		v->y = (int)(2.0 * step_y) + start->y;
			
		v->next = p->next;
		p->next = v;	
	}
}

/*
draw the generated koch snow flake on a raw image.
*/

void draw_koch_snowflake(int rows, int cols, unsigned char ** pixels, VERTEX * list, unsigned b_color, unsigned char f_color)
{
	int i,j;
	VERTEX * p = list;

	for(i = 0; i < rows; i++)
		for(j = 0; j < cols; j++)
			pixels[i][j] = b_color;

	while(p->next != NULL)
	{
		bresenham(pixels,p,p->next,f_color);
		p = p->next;
	}
	bresenham(pixels,p,list,f_color);
}




