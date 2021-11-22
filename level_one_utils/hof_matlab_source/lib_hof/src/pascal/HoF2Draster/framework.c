/*Compute the force histogram between two (fuzzy) raster objects using the new framework */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "framework.h"
#include "fourier.h"
#include "util.h"



/*The formulas for computing the forces between two pixels 
according to different types of forces considered*/

double F_less_1_opt(double y, double r)
{
	double t = 2.0 - r;
	if(y>1)
	{
		return ( pow(y+1,t) - 2.0*pow(y,t) + pow(y-1,t) ) / (t*(1.0-r));
	}
	else
	{
		if(y==1)
			return (pow(2.0,t)-2.0) / (t*(1.0-r));
		else
			return 1.0 / (t*(1.0-r));
	}
}
double F_1_opt(double y, double r)
{
	if(y>1)
	{
		return (y+1) * log(y+1) - 2.0 * y * log(y) + (y-1)*log(y-1);
	}
	else
	{
		if(y==1)
			return 2.0 * log(2);
		else
			return 0.0;
	}
}
double F_1_2_opt(double y, double r)
{
	double t = 2.0 - r;
	if(y>1)
	{
		return (pow(y+1,t) - 2.0*pow(y,t) + pow(y-1,t)) / (t*(1.0-r));
	}
	else
	{
		if(y==1)
			return (pow(2.0,t)-2.0) / (t*(1.0-r));
		else
			return 0.0;
	}
}
double F_2_opt(double y, double r)
{
	if(y>1)
		return log( pow(y,2.0) / ( (y-1)*(y+1) ) );
	else
		return 0.0;
}
double F_bigger_2_opt(double y, double r)
{
	double t = 2.0 - r;
	if(y>1)
	{
		return (pow(y+1,t) - 2.0*pow(y,t) + pow(y-1,t)) / (t*(1.0-r));
	}
	else
		return 0.0;
}

/*
Compute the spatial correlation between two raster objects using standard algorithm
with complexity O(N*N).
p_obj_a: the object A;
p_obj_b: the object B.
*/

RAS_OBJ * createLandscape_slow(RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	int i_a,j_a, i_b, j_b;
	int rows, cols, t_rows, t_cols, org_y, org_x;
	double ** obj_a_re;
	double ** obj_b_re;
	double ** result_re;
	double value;
	RAS_OBJ * result;
	
	obj_a_re = p_obj_a->re;
	obj_b_re = p_obj_b->re;
	rows = p_obj_a->rows;
	cols = p_obj_b->cols;

	t_rows = 2 * rows - 1;
	t_cols = 2 * cols - 1;
	org_y = rows - 1;
	org_x = cols - 1;

	result = createRasObj(t_rows,t_cols,org_y,org_x);
	result_re = result->re;

	for(i_a = 0; i_a < rows; i_a++)
	{
		for(j_a = 0; j_a < cols; j_a++)
		{
			value = obj_a_re[i_a][j_a];
			for(i_b = 0; i_b < rows; i_b++)
			{
				for(j_b = 0; j_b < cols; j_b++)
				{
					result_re[org_y + i_b - i_a][org_x + j_b - j_a] += (obj_b_re[i_b][j_b] * value);
				}	
			}
		}	
	}
	return result;
}
/*
Old implementation of fast computation of spatial correlation
not used any more.
*/
RAS_OBJ * createLandscape_conv(RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	int i, j;
	int rows, cols, t_rows, t_cols, org_y, org_x;
	int n_cols, n_rows, factor_vertical, factor_horizontal;

	RAS_OBJ * extend_a, * extend_b, * aid_a, * aid_b;
	RAS_OBJ * con_res;
	RAS_OBJ * result;
	BOUND_BOX * box_a, * box_b;

	double ** sou_re;
	double ** tar_re;

	rows = p_obj_a->rows;
	cols = p_obj_a->cols;
	t_rows = 2 * rows - 1;
	t_cols = 2 * cols - 1;
	org_y = rows - 1;
	org_x = cols - 1;
	
	factor_vertical = factorBy2(t_rows);
	factor_horizontal = factorBy2(t_cols);

	n_rows = 1 << factor_vertical;
	n_cols = 1 << factor_horizontal;

	box_a = loadBoundBox(p_obj_a);
	box_b = loadBoundBox(p_obj_b);

	extend_a = resizeRasObj_boundBox(p_obj_a, box_a, n_rows,n_cols,org_y,org_x);
	extend_b = resizeRasObj_boundBox(p_obj_b, box_b, n_rows,n_cols,0,0);
	

	complexRasObj(extend_a);
	complexRasObj(extend_b);

	aid_a = createRasObjComplex(n_rows,n_cols,org_y,org_x);
	aid_b = createRasObjComplex(n_rows,n_cols,0,0);

	reflectRasObj_boundBox(extend_a, box_a);

	

	con_res = convolution(extend_a, box_a, extend_b, box_b, aid_a, aid_b, factor_vertical,factor_horizontal);

	result = createRasObj(t_rows,t_cols,org_y,org_x);

	
	sou_re = con_res->re;
	tar_re = result->re;
	
	for(i = 0; i < t_rows; i++)
	{
		for(j = 0; j < t_cols; j++)
				tar_re[i][j] = sou_re[i][j];
	}

	destroyRasObj(extend_a);
	destroyRasObj(extend_b);
	destroyRasObj(aid_a);
	destroyRasObj(aid_b);
	free(box_a);
	free(box_b);

	return result;

}

/*
Compute the spatial correlation using convolution and fast FFT.
*/

RAS_OBJ * createLandscape_conv_opt(RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	int i, j;
	int rows, cols, t_rows, t_cols, org_y, org_x, small_rows, small_cols, t_small_rows, t_small_cols, small_org_y, small_org_x;
	int n_cols, n_rows, factor_vertical, factor_horizontal, a_rows, a_cols, b_rows, b_cols;
	int dis_y, dis_x, start_i, end_i, start_j, end_j;
	RAS_OBJ * extend_a, * extend_b, * aid_a, * aid_b;
	RAS_OBJ * con_res;
	RAS_OBJ * result;
	BOUND_BOX * box_a, * box_b;

	double ** sou_re;
	double ** tar_re;


	rows = p_obj_a->rows;
	cols = p_obj_a->cols;
	t_rows = 2 * rows - 1;
	t_cols = 2 * cols - 1;
	org_y = rows - 1;
	org_x = cols - 1;
	

	box_a = loadBoundBox(p_obj_a);                   /*Compute A's minimum bounding rectangle*/
	box_b = loadBoundBox(p_obj_b);                   /*Compute B's minimum bounding rectangle*/

	a_rows = box_a->v1[0] - box_a->v3[0] + 1;
	b_rows = box_b->v1[0] - box_b->v3[0] + 1;
	a_cols = box_a->v3[1] - box_a->v1[1] + 1;
	b_cols = box_b->v3[1] - box_b->v1[1] + 1;

	small_rows = a_rows > b_rows ? a_rows : b_rows;
	small_cols = a_cols > b_cols ? a_cols : b_cols;

	t_small_rows = 2 * small_rows - 1;
	t_small_cols = 2 * small_cols - 1;

	small_org_y = small_rows - 1;
	small_org_x = small_cols - 1;

	dis_y = box_b->v2[0] - box_a->v2[0];
	dis_x = box_b->v2[1] - box_a->v2[1];

	factor_vertical = factorBy2(t_small_rows);
	factor_horizontal = factorBy2(t_small_cols);


	n_rows = 1 << factor_vertical;
	n_cols = 1 << factor_horizontal;


    /*Extend the image region of object A 4 times bigger and set the origin in the middle of the extended image.
    And shift the object A from the low left quarter region to the up right quarter region*/
	extend_a = resizeRasObj_boundBox_region(p_obj_a, box_a, n_rows, n_cols, small_org_y, small_org_x, a_rows,a_cols);
	/*Extend the image region of object B 4 times bigger and set the origin in the middle of the extended image.*/
	extend_b = resizeRasObj_boundBox_region(p_obj_b, box_b, n_rows, n_cols, 0 , 0 , b_rows, b_cols);

	complexRasObj(extend_a);
	complexRasObj(extend_b);

	aid_a = createRasObjComplex(n_rows,n_cols,small_org_y,small_org_x);
	aid_b = createRasObjComplex(n_rows,n_cols,0,0);

	/*reflect the object a about the origin*/
    reflectRasObj_boundBox(extend_a, box_a);
	
    /*Compute the convolution between object A and the reflected object A.*/
	con_res = convolution(extend_a, box_a, extend_b, box_b, aid_a, aid_b, factor_vertical,factor_horizontal);


	/*--------------------------------

	result = createRasObj(t_small_rows,t_small_cols,small_org_y,small_org_x);

	
	sou_re = con_res->re;
	tar_re = result->re;
	
	for(i = 0; i < t_small_rows; i++)
	{
		for(j = 0; j < t_small_cols; j++)
				tar_re[i][j] = sou_re[i][j];
	}

	--------------------------------*/

    
	result = createRasObj(t_rows, t_cols,org_y,org_x);

	
	sou_re = con_res->re;
	tar_re = result->re;

	org_y = org_y - small_org_y + dis_y;
	org_x = org_x - small_org_x + dis_x;

	start_i = 0 > (-1 * org_y) ? 0: -1 * org_y;
	start_j = 0 > (-1 * org_x) ? 0: -1 * org_x;

	end_i = t_small_rows < (t_rows - org_y) ? t_small_rows : t_rows - org_y;
	end_j = t_small_cols < (t_cols - org_x) ? t_small_cols : t_cols - org_x;
    
    /*Adjust the results*/
	for(i = start_i; i < end_i; i++)
	{
		for(j = start_j; j < end_j; j++)
			tar_re[i + org_y][j + org_x] = sou_re[i][j] > 0.0 ? sou_re[i][j] : 0.0;
	}

	
	destroyRasObj(extend_a);
	destroyRasObj(extend_b);
	destroyRasObj(aid_a);
	destroyRasObj(aid_b);
	free(box_a);
	free(box_b);

	return result;

}
/*
Compute the angle histogram according to the new framework.
*/

int landscape_angle_histogram(int p_d, double * p_histogram, RAS_OBJ * p_scape)
{
	int i,j,org_x,org_y,rows,cols,c_d, range_j, range_i, i_1, i_2, j_1, j_2;
	int p_d_2, p_d_4, p_d_3_4;
	double ** re;
	double step = 360 / (double) p_d;
	double pi = acos(-1.0);
	double rate = 180.0 / pi;
	double angle, value_0;

	p_d_2 = p_d / 2;
	p_d_4 = p_d / 4;
	p_d_3_4 = 3 * p_d_4;

	rows = p_scape->rows;
	cols = p_scape->cols;
	org_y = p_scape->org_y;
	org_x = p_scape->org_x;
	re = p_scape->re;

	range_j = cols - org_x;
	range_i = rows - org_y;

	if(range_j > (org_x + 1))
		range_j = org_x + 1;

	if(range_i > (org_y + 1))
		range_i = org_y + 1;

	for(i = 1; i < range_i; i++)
	{
		i_1 = org_y + i;
		i_2 = org_y - i;
		for(j = 1; j < range_j; j++)
		{
			j_1 = org_x + j;
			j_2 = org_x - j;

			angle = rate * atan((double) i / (double) j);
			c_d = (int)(angle / step + 0.5);
			p_histogram[c_d % p_d] += re[i_1][j_1];

			c_d = (int)((360.0 - angle) / step + 0.5);
			p_histogram[c_d % p_d] += re[i_2][j_1];
			
			c_d = (int)((180.0 - angle) / step + 0.5);
			p_histogram[c_d % p_d] += re[i_1][j_2];

			c_d = (int)((180.0 + angle) / step + 0.5);
			p_histogram[c_d % p_d] += re[i_2][j_2];
		}
	}

	for(j = 1; j < range_j; j++)
	{
		p_histogram[0] += re[org_y][org_x + j];
		p_histogram[p_d_2] += re[org_y][org_x - j];		
	}
	for(i = 1; i < range_i; i++)
	{
		p_histogram[p_d_4] += re[org_y + i][org_x];
		p_histogram[p_d_3_4] += re[org_y - i][org_x];
	}

	value_0 = re[org_y][org_x];
	if(value_0 > 0)
	{
		for(i = 0; i < p_d; i++)
		{
			p_histogram[i] += value_0;
		}
	}
	return 1;
}

/*
Application of angle histogram computation using the new framework.
It calls function landscape_angle_histogram.
*/

double * landscapeAngleHistogram(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	double * histogram = (double *) calloc(p_d, sizeof(double));
	RAS_OBJ * scape = createLandscape_conv(p_obj_a, p_obj_b);
	landscape_angle_histogram(p_d, histogram, scape);
	destroyRasObj(scape);
	return histogram;
}

double adjF0(double angle)
{
	return (1.0 / cos(angle) );
}
double adjF2(double angle)
{
	return cos(angle);
}

/*Not used any more*/

int landscape_force_histogram_0(int p_d, double * p_histogram, RAS_OBJ * p_scape)
{
	int rows, cols, org_y, org_x, range_i, range_j, loc;
	double ** re;
	double step_angle = 360.0 / (double) p_d;

	int i,j,d, pos, r_pos, comp, r_comp;
	double tan_value;
	double angle, l_angle;
	double pi = acos(-1.0);
	double factor = 1.0;

	int d_8 = p_d / 8;
	int d_2 = p_d / 2;
	int d_4 = p_d / 4;
	int d_3_4 = 3 * d_4;

	rows = p_scape->rows;
	cols = p_scape->cols;
	org_y = p_scape->org_y;
	org_x = p_scape->org_x;
	re = p_scape->re;

	range_j = cols - org_x;
	range_i = rows - org_y;

	if(range_j > (org_x + 1))
		range_j = org_x + 1;

	if(range_i > (org_y + 1))
		range_i = org_y + 1;

	for(j = 0; j < range_j; j++)
	{
		p_histogram[0] +=  re[org_y][org_x + j]; 
		p_histogram[d_2] += re[org_y][org_x - j]; 
	}

	for(i = 0; i< range_i; i++)
	{
		p_histogram[d_4] += re[org_y + i][org_x];
		p_histogram[d_3_4] += re[org_y - i][org_x];
	}

	for(d = 1; d <= d_8; d++)
	{
		pos = d;
		r_pos = d_2 + pos;
		comp = d_2 - pos;
		r_comp = d_2 + comp;


		angle = (double)d * step_angle;
		l_angle = angle * pi / 180.0;

		if(angle == 45)
			tan_value = 1.0;
		else
			tan_value = tan(l_angle);

		factor = 1.0 / cos(l_angle); 
		for(j = 0; j < range_j; j++)
		{
			
			loc = (int)((double)j * tan_value + 0.5);
			if(loc >= range_i)
				break;
			p_histogram[pos] += re[org_y + loc][org_x + j]; 
			p_histogram[r_pos] += re[org_y - loc][org_x - j]; 
			p_histogram[comp] += re[org_y + loc][org_x - j]; 
			p_histogram[r_comp] += re[org_y - loc][org_x + j]; 
		}
		p_histogram[pos] *= factor; 
		p_histogram[r_pos] *= factor; 
		p_histogram[comp] *= factor; 
		p_histogram[r_comp] *= factor; 
	}

	for(d = d_8 + 1; d < d_4; d++)
	{
		pos = d;
		r_pos = d_2 + pos;
		comp = d_2 - pos;
		r_comp = d_2 + comp;

		angle = 90.0 - (double)d * step_angle;
		l_angle = angle * pi / 180.0;

		tan_value = tan(l_angle);
		factor = 1.0 / cos(l_angle);  
		for(i = 0; i < range_i; i++)
		{
			loc = (int)((double)i * tan_value + 0.5);
			if(loc >= range_j)
				break;	
			p_histogram[pos] += re[org_y + i][org_x + loc]; 
			p_histogram[r_pos] += re[org_y - i][org_x - loc]; 
			p_histogram[comp] += re[org_y + i][org_x - loc]; 
			p_histogram[r_comp] += re[org_y - i][org_x + loc]; 
		}
		p_histogram[pos] *= factor; 
		p_histogram[r_pos] *= factor; 
		p_histogram[comp] *= factor; 
		p_histogram[r_comp] *= factor; 
	}
	return 1;

}

/*
Not used any more.
*/
double * landscapeForceHistogram_0(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	double * histogram = (double *) calloc(p_d, sizeof(double));
	RAS_OBJ * scape = createLandscape_conv_opt(p_obj_a, p_obj_b);
	landscape_force_histogram_0(p_d, histogram, scape);
	destroyRasObj(scape);
	return histogram;
}

/*
Not used any more.
*/

int landscape_force_histogram_2(int p_d, double * p_histogram, RAS_OBJ * p_scape)
{
	int rows, cols, org_y, org_x, range_i, range_j, loc, max;
	double ** re;
	double * dis;
	double step_angle = 360.0 / (double) p_d;

	int i,j,d, pos, r_pos, comp, r_comp;
	double tan_value;
	double angle, l_angle;
	double pi = acos(-1.0);
	double factor = 1.0;
	double r;

	int d_8 = p_d / 8;
	int d_2 = p_d / 2;
	int d_4 = p_d / 4;
	int d_3_4 = 3 * d_4;

	rows = p_scape->rows;
	cols = p_scape->cols;
	org_y = p_scape->org_y;
	org_x = p_scape->org_x;
	re = p_scape->re;

	range_j = cols - org_x;
	range_i = rows - org_y;

	if(range_j > (org_x + 1))
		range_j = org_x + 1;

	if(range_i > (org_y + 1))
		range_i = org_y + 1;

	max = rows > cols ? rows : cols;

	dis = (double *) malloc (max * sizeof(double));

	for(i = 0; i < max; i ++)
		dis[i] = i * i;

	for(j = 1; j < range_j; j++)
	{
		r = dis[j];
		p_histogram[0] +=  (re[org_y][org_x + j] / r); 
		p_histogram[d_2] += (re[org_y][org_x - j] / r); 
	}

	for(i = 1; i< range_i; i++)
	{
		r = dis[i];
		p_histogram[d_4] += (re[org_y + i][org_x] / r);
		p_histogram[d_3_4] += (re[org_y - i][org_x] / r);
		
	}

	for(d = 1; d <= d_8; d++)
	{
		pos = d;
		r_pos = d_2 + pos;
		comp = d_2 - pos;
		r_comp = d_2 + comp;


		angle = (double)d * step_angle;
		l_angle = angle * pi / 180.0;

		if(angle == 45)
			tan_value = 1.0;
		else
			tan_value = tan(l_angle);

		factor = cos(l_angle); 
		for(j = 1; j < range_j; j++)
		{
			r = dis[j];
			loc = (int)((double)j * tan_value + 0.5);
			if(loc >= range_i)
				break;
			p_histogram[pos] += (re[org_y + loc][org_x + j] / r); 
			p_histogram[r_pos] += (re[org_y - loc][org_x - j]  / r); 
			p_histogram[comp] += (re[org_y + loc][org_x - j]  / r); 
			p_histogram[r_comp] += (re[org_y - loc][org_x + j]  / r); 
		}
		p_histogram[pos] *= factor; 
		p_histogram[r_pos] *= factor; 
		p_histogram[comp] *= factor; 
		p_histogram[r_comp] *= factor; 
	}

	for(d = d_8 + 1; d < d_4; d++)
	{
		pos = d;
		r_pos = d_2 + pos;
		comp = d_2 - pos;
		r_comp = d_2 + comp;

		angle = 90.0 - (double)d * step_angle;
		l_angle = angle * pi / 180.0;

		tan_value = tan(l_angle);
		factor = cos(l_angle);  
		for(i = 1; i < range_i; i++)
		{
			r = dis[i];
			loc = (int)((double)i * tan_value + 0.5);
			if(loc >= range_j)
				break;	
			p_histogram[pos] += (re[org_y + i][org_x + loc] / r); 
			p_histogram[r_pos] += (re[org_y - i][org_x - loc] / r); 
			p_histogram[comp] += (re[org_y + i][org_x - loc] / r); 
			p_histogram[r_comp] += (re[org_y - i][org_x + loc] / r); 
		}
		p_histogram[pos] *= factor; 
		p_histogram[r_pos] *= factor; 
		p_histogram[comp] *= factor; 
		p_histogram[r_comp] *= factor; 
	}
	free(dis);
	return 1;

}

/*
Not used any more.
*/
double * landscapeForceHistogram_2(int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	double * histogram = (double *) calloc(p_d, sizeof(double));
	RAS_OBJ * scape = createLandscape_conv_opt(p_obj_a, p_obj_b);
	landscape_force_histogram_2(p_d, histogram, scape);
	destroyRasObj(scape);
	return histogram;
}

/*
Old implementatin, not used any more.
*/

int landscape_force_histogram_r(double r_t, int p_d, double * p_histogram, RAS_OBJ * p_scape)
{
	int rows, cols, org_y, org_x, range_i, range_j, loc, max;
	double ** re;
	double * dis;
	double step_angle = 360.0 / (double) p_d;

	int i,j,d, pos, r_pos, comp, r_comp;
	double tan_value;
	double angle, l_angle;
	double pi = acos(-1.0);
	double factor = 1.0;
	double r;

	int d_8 = p_d / 8;
	int d_2 = p_d / 2;
	int d_4 = p_d / 4;
	int d_3_4 = 3 * d_4;

	rows = p_scape->rows;
	cols = p_scape->cols;
	org_y = p_scape->org_y;
	org_x = p_scape->org_x;
	re = p_scape->re;

	range_j = cols - org_x;
	range_i = rows - org_y;

	if(range_j > (org_x + 1))
		range_j = org_x + 1;

	if(range_i > (org_y + 1))
		range_i = org_y + 1;

	max = rows > cols ? rows : cols;

	dis = (double *) malloc (max * sizeof(double));

	for(i = 0; i < max; i ++)
		dis[i] = pow((double)i,r_t);

	for(j = 1; j < range_j; j++)
	{
		r = dis[j];
		p_histogram[0] +=  (re[org_y][org_x + j] / r); 
		p_histogram[d_2] += (re[org_y][org_x - j] / r); 
	}

	for(i = 1; i< range_i; i++)
	{
		r = dis[i];
		p_histogram[d_4] += (re[org_y + i][org_x] / r);
		p_histogram[d_3_4] += (re[org_y - i][org_x] / r);
		
	}

	for(d = 1; d <= d_8; d++)
	{
		pos = d;
		r_pos = d_2 + pos;
		comp = d_2 - pos;
		r_comp = d_2 + comp;


		angle = (double)d * step_angle;
		l_angle = angle * pi / 180.0;

		if(angle == 45)
			tan_value = 1.0;
		else
			tan_value = tan(l_angle);

		factor = cos(l_angle); 
		factor = factor / pow(factor,2.0 - r_t);

		for(j = 1; j < range_j; j++)
		{
			r = dis[j];
			loc = (int)((double)j * tan_value + 0.5);
			if(loc >= range_i)
				break;
			p_histogram[pos] += (re[org_y + loc][org_x + j] / r); 
			p_histogram[r_pos] += (re[org_y - loc][org_x - j]  / r); 
			p_histogram[comp] += (re[org_y + loc][org_x - j]  / r); 
			p_histogram[r_comp] += (re[org_y - loc][org_x + j]  / r); 
		}
		p_histogram[pos] *= factor; 
		p_histogram[r_pos] *= factor; 
		p_histogram[comp] *= factor; 
		p_histogram[r_comp] *= factor; 
	}

	for(d = d_8 + 1; d < d_4; d++)
	{
		pos = d;
		r_pos = d_2 + pos;
		comp = d_2 - pos;
		r_comp = d_2 + comp;

		angle = 90.0 - (double)d * step_angle;
		l_angle = angle * pi / 180.0;

		tan_value = tan(l_angle);
		factor = cos(l_angle); 
		factor = factor / pow(factor,2.0 - r_t);
		for(i = 1; i < range_i; i++)
		{
			r = dis[i];
			loc = (int)((double)i * tan_value + 0.5);
			if(loc >= range_j)
				break;	
			p_histogram[pos] += (re[org_y + i][org_x + loc] / r); 
			p_histogram[r_pos] += (re[org_y - i][org_x - loc] / r); 
			p_histogram[comp] += (re[org_y + i][org_x - loc] / r); 
			p_histogram[r_comp] += (re[org_y - i][org_x + loc] / r); 
		}
		p_histogram[pos] *= factor; 
		p_histogram[r_pos] *= factor; 
		p_histogram[comp] *= factor; 
		p_histogram[r_comp] *= factor; 
	}
	free(dis);
	return 1;

}

/*
Old implementation not used any more.
*/

double * landscapeForceHistogram_r(double r_t, int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	double * histogram = (double *) calloc(p_d, sizeof(double));
	RAS_OBJ * scape = createLandscape_conv_opt(p_obj_a, p_obj_b);
	landscape_force_histogram_r(r_t, p_d, histogram, scape);
	destroyRasObj(scape);
	return histogram;
}

/*
Current implementation of force histogram computation using the new framework.
r_t: the type of force considered
p_d: the number of reference directions considered.
p_histogram: stores the result histogram.
p_scape: the spatial correlation between the object A and the object B.
*/

int landscape_force_histogram_r_precise(double r_t, int p_d, double * p_histogram, RAS_OBJ * p_scape)
{
	int rows, cols, org_y, org_x, range_i, range_j, loc, max;
	double ** re;
	double * dis;
	double step_angle = 360.0 / (double) p_d;

	int i,j,d, pos, r_pos, comp, r_comp;
	double tan_value;
	double angle, l_angle;
	double pi = acos(-1.0);
	double factor = 1.0;
	double r;

	int d_8 = p_d / 8;
	int d_2 = p_d / 2;
	int d_4 = p_d / 4;
	int d_3_4 = 3 * d_4;
	
	

	double (* F)(double x, double y, double z, double r);
    
    /*According the type of the force, decide the formula for computing the force between two pixels*/
	if(r_t < 1.0)
		F = F_less_1;
	else
	{
		if(r_t == 1.0)
			F = F_1;
		else
		{
			if(r_t < 2.0)
				F = F_1_2;
			else
			{
				if(r_t == 2.0)
					F = F_2;
				else
					F = F_bigger_2;
			}
		}
	}


	rows = p_scape->rows;
	cols = p_scape->cols;
	org_y = p_scape->org_y;
	org_x = p_scape->org_x;
	re = p_scape->re;

	range_j = cols - org_x;
	range_i = rows - org_y;

	if(range_j > (org_x + 1))
		range_j = org_x + 1;

	if(range_i > (org_y + 1))
		range_i = org_y + 1;

	max = rows > cols ? rows : cols;

	dis = (double *) malloc (max*sizeof(double));

    /*Initialize the force values between pixels with different distances */
	for(i = 0; i < max; i ++)
		dis[i] = F(1.0,(double)i-1.0,1.0,r_t);

    /*Calculate first the force values along the primitive dirctions: left and right*/
	//for(j = 1; j < range_j; j++)
	for(j = 0; j < range_j; j++)
	{
		r = dis[j];
		p_histogram[0] +=  (re[org_y][org_x + j] * r); 
		p_histogram[d_2] += (re[org_y][org_x - j] * r); 
	}
	
    /*Calculate first the force values along the primitive dirctions: up and below*/
	//for(i = 1; i< range_i; i++)
	for(i = 0; i< range_i; i++)
	{
		r = dis[i];
		p_histogram[d_4] += (re[org_y + i][org_x] * r);
		p_histogram[d_3_4] += (re[org_y - i][org_x] * r);
		
	}
    /*Compute the force values along the directions in (0, pi/4]*/
	for(d = 1; d <= d_8; d++)
	{
		pos = d;
		r_pos = d_2 + pos;
		comp = d_2 - pos;
		r_comp = d_2 + comp;


		angle = (double)d * step_angle;
		l_angle = angle * pi / 180.0;

		if(angle == 45)
			tan_value = 1.0;
		else
			tan_value = tan(l_angle);

		factor = pow(cos(l_angle),r_t - 1.0);
		//factor = cos(l_angle); 
		//factor = factor / pow(factor,2.0 - r_t);

		//for(j = 1; j < range_j; j++)
		for(j = 0; j < range_j; j++)
		{
			r = dis[j];
			loc = (int)((double)j * tan_value + 0.5);                 /*Rasterize the line along the current direction*/
			if(loc >= range_i)
				break;
			p_histogram[pos] += (re[org_y + loc][org_x + j] * r);     /*Compute the force value along the direction*/
			p_histogram[r_pos] += (re[org_y - loc][org_x - j]  * r);  /*Compute the force value along the opposite of the direction in range (pi, 5pi/4]*/
			p_histogram[comp] += (re[org_y + loc][org_x - j]  * r);   /*Compute the force value along the corresponding direction in range [3pi/4, pi)*/
			p_histogram[r_comp] += (re[org_y - loc][org_x + j]  * r); /*Compute the force value along the corresponding direction in range [7pi/4,2pi)*/
		}
		/*Adjust the force values*/
        p_histogram[pos] *= factor; 
		p_histogram[r_pos] *= factor; 
		p_histogram[comp] *= factor; 
		p_histogram[r_comp] *= factor; 
	}
    /*Compute the force values along the directions in (pi/4, pi/2)*/
	for(d = d_8 + 1; d < d_4; d++)
	{
		pos = d;
		r_pos = d_2 + pos;
		comp = d_2 - pos;
		r_comp = d_2 + comp;

		angle = 90.0 - (double)d * step_angle;
		l_angle = angle * pi / 180.0;

		tan_value = tan(l_angle);
		factor = pow(cos(l_angle),r_t - 1.0);
		//factor = cos(l_angle); 
		//factor = factor / pow(factor,2.0 - r_t);
		//for(i = 1; i < range_i; i++)
		for(i = 0; i < range_i; i++)
		{
			r = dis[i];
			loc = (int)((double)i * tan_value + 0.5);                /*Rasterize the line along the current direction*/
			if(loc >= range_j)
				break;	
			p_histogram[pos] += (re[org_y + i][org_x + loc] * r);    /*Compute the force value along the direction*/  
			p_histogram[r_pos] += (re[org_y - i][org_x - loc] * r);  /*Compute the force value along the opposite of the direction in range (5pi/4, 3pi/2)*/
			p_histogram[comp] += (re[org_y + i][org_x - loc] * r);   /*Compute the force value along the corresponding direction in range (pi/2, 3pi/4)*/
			p_histogram[r_comp] += (re[org_y - i][org_x + loc] * r); /*Compute the force value along the corresponding direction in range (3pi/2, 7pi/4)*/
		}
		/*Adjust force values*/
        p_histogram[pos] *= factor; 
		p_histogram[r_pos] *= factor; 
		p_histogram[comp] *= factor; 
		p_histogram[r_comp] *= factor; 
	}
	free(dis);
	return 1;

}

/*
The application for computing force histogram using the new framework.
It calls function landscape_force_histogram_r_precise.
*/

double * landscapeForceHistogram_r_precise(double r_t, int p_d, RAS_OBJ * p_obj_a, RAS_OBJ * p_obj_b)
{
	double * histogram = (double *) calloc(p_d, sizeof(double));
	RAS_OBJ * scape = createLandscape_conv_opt(p_obj_a, p_obj_b);
	landscape_force_histogram_r_precise(r_t, p_d, histogram, scape);
	destroyRasObj(scape);
	return histogram;
}

/*
int main(int argc, char * argv[])
{
	int p_d;
	int grey_A, grey_B, grey_A_B,m;
	float f;
	FILE * fp = NULL;
	double * histogram;
	IMAGE * image;
	RAS_OBJ * obj_a;
	RAS_OBJ * obj_b;
	clock_t tv1, tv2;
	double time = 0.0;


	if(argc < 9)
	{
		printf("parameters: file_name A B A_B d f m his_file\n");
		return 0;
	}
	image = loadImage(argv[1]);


	sscanf(argv[2],"%d", &grey_A);
	sscanf(argv[3],"%d", &grey_B);
	sscanf(argv[4],"%d", &grey_A_B);
	sscanf(argv[5],"%d", &p_d);
	sscanf(argv[6],"%f", &f);
	sscanf(argv[7],"%d", &m);

	printf("This is the optmized version of new framework.\n");

	obj_a = image2RasObjExtract(image,grey_A,grey_A_B);
	obj_b = image2RasObjExtract(image,grey_B,grey_A_B);
	if(m > 0)
	{
		randomFuzzify(obj_a, 100, m);
		randomFuzzify(obj_b, 200, m);
	}
	
	//if(f == 0)
	//	histogram = landscapeForceHistogram_0(p_d,obj_a,obj_b);
	//else
	//	histogram = landscapeForceHistogram_2(p_d,obj_a,obj_b);
	
	tv1 = clock();

	histogram = landscapeForceHistogram_r_precise((double)f, p_d,obj_a,obj_b);

	tv2 = clock();
	time = (((double)(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0))/1000.0);

	printf("New framework computing force histogram:\n");
	printf("Computation time (seconds) is: %f\n",time);

	//histogram = landscapeAngleHistogram(d,obj_a,obj_b);
	fp = fopen(argv[8], "w");
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
	float r;
	FILE * fp = NULL;
	double * histogram;
	IMAGE * image_A;
	IMAGE * image_B;
	RAS_OBJ * obj_a;
	RAS_OBJ * obj_b;
	clock_t tv1, tv2;
	double time = 0.0;


	if(argc < 6)
	{
		printf("\nComputing different types of force histograms using the New Framework\n\n");
        printf("parameters: file_name_A file_name_B d r histogram_file\n\n");
		printf("file_name_A: the pgm image of object A (reference object);\n");
		printf("file_name_B: the pgm image of object B (argument object);\n");
		printf("d: the number of reference directions considered (integer);\n");
		printf("r: the type of forces (0: constant forces; 2:gravitational forces; or other real numbers);\n");
		printf("histogram_file: the txt file for storing the histogram values;\n");
		return 0;
	}
	image_A = loadImage(argv[1]);       /*Load pgm image containing object A*/
	image_B = loadImage(argv[2]);       /*Load pgm iamge containing object B*/


	sscanf(argv[3],"%d", &p_d);         /*The number of reference directions*/
	sscanf(argv[4],"%f", &r);           /*The type of forces*/

	printf("This is the optmized version of new framework.\n");

	obj_a = image2RasObj(image_A);       /*Load raster object A from Image*/
	obj_b = image2RasObj(image_B);       /*Load raster object B from Image*/
	
	
    tv1 = clock();

	histogram = landscapeForceHistogram_r_precise((double)r, p_d,obj_a,obj_b);     /*Compute the force histogram*/

	tv2 = clock();
	time = (((double)(tv2 - tv1)/(CLOCKS_PER_SEC / (double) 1000.0))/1000.0);

	printf("New framework computing force histogram:\n");
	printf("Computation time (seconds) is: %f\n",time);

	fp = fopen(argv[5], "w");
	histogramWrite(fp,histogram,p_d);     /*Write the result histogram into a txt file.*/
	fclose(fp);
	destroyRasObj(obj_a);
	destroyRasObj(obj_b);
	destroyImage(image_A);
	destroyImage(image_B);
	free(histogram); 
	return 1;
}


