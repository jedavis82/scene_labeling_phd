/*
Implemented some utility functions 
*/

#include <math.h>
#include "util.h"

/*
Calculate the smallest factor m such that pow(2,m) > p_value;
*/
int factorBy2(int p_value)
{
	int factor = 0;
	int residue = 0;
	int value = p_value;
	
	if(p_value < 1)
		return -1;

	while(value > 1)
	{	
		if((value % 2)!=0)
			residue = 1;
		value = value >> 1;
		factor ++;
	}
	return(factor + residue);
}

/*
Rasterize a line using Bresenham algorithm:
*/
int populateLine(int * p_line, double p_tan, int p_n)
{
	int i;
	double base = 0.5;
	
	if(p_tan < 0.0)
	{
		base = -0.5;
	}
	for(i = 0; i < p_n; i++)
	{
		p_line[i] = (int)base;
		base += p_tan;
	}
	return p_n;
}

/*
Rasterize a line with any direction using the Bresenham algorithm:
*/

int rasterizeLine(int * p_line, double p_angle, int p_rows, int p_cols)
{
	double angle = p_angle;
	double pi = acos(-1.0);
	if(angle < 45 || angle > 135)
		return populateLine(p_line, tan(angle * pi / 180.0), p_cols);
	if(angle == 45)
		return populateLine(p_line, 1.0, p_cols);
	if(angle == 135)
		return populateLine(p_line, -1.0, p_cols);
	return populateLine(p_line, tan((90.0 - angle) * pi / 180.0), p_rows);
	
}
/*
Projecting a segment on to one of the horizontal and vertical axies:
*/

int projection(int * p_start, int * p_end, int * p_v_start, int * p_v_end, double p_tan, int p_d)
{
	double start, end;
	int x = 1 - p_d;
	start = (double) p_v_start[p_d] - p_tan * (double) p_v_start[x];
	end = (double) p_v_end[p_d] - p_tan * (double) p_v_end[x];
	* p_start = (int) (start - 1.0);
	* p_end = (int) (end + 1.0);
	return 1;
}

/*
Formula for calculating the constant force between segments.
*/
double F0(double x, double y, double z)
{
	double value;
	if(y > 0)
		value = x * z;
	else
	{
		if( x + y > 0)
		{
			if( y + z > 0)
				value = x * z - y * y / 2.0;
			else
				value = (x + y + z / 2.0) * z;
		}
		else
		{
			if( y + z > 0)
				value = (x / 2.0 + y + z) * x;
			else
			{
				if( x + y + z > 0)
					value = pow((x + y + z), 2);
				else
					value = 0.0;
			}
		}
	}
	return value;
}

/*
Formula for calculating the gravitational force between segments.
*/

double F2(double x, double y, double z)
{
	if( y > 0)
		return log( (y + x) * (y + z) / (y * (y + x + z)) );
	else
	{
		if( x == 0 || z == 0 || x + y + z <=0)
			return 0.0;
		else
			return 0.0;	
	}

}

/*
Calculating the force between segments with r<1;
*/

double F_less_1(double x, double y, double z, double r)
{
	double value;
	double t = 2-r;
	if(y > 0)
		value = (pow(y,t) - pow(x+y,t) - pow(y+z,t) + pow(x+y+z,t) ) / (t*(1-r));
	else
	{
		if( x + y > 0)
		{
			if( y + z > 0)
				value = (pow(x+y+z,t) - pow(x+y,t) - pow(y+z,t) ) / (t*(1-r));
			else
				value = (pow(x+y+z,t) - pow(x+y,t)) / (t*(1-r));
		}
		else
		{
			if( y + z > 0)
				value = (pow(x+y+z,t) - pow(y+z,t)) / (t*(1-r));
			else
			{
				if( x + y + z > 0)
					value = pow(x+y+z,t) / (t*(1-r));
				else
					value = 0.0;
			}
		}
	}
	return value;
}

/*
Calculating the force between segments with r=1;
*/

double F_1(double x, double y, double z, double r)
{
	if( y >= 0)
	{
		if(y==0)
			return (y+z)*log( (x+y+z)/(y+z) ) + x*log( (x+y+z)/(x+y) );
		else
			return (y+z)*log( (x+y+z)/(y+z) ) + x*log( (x+y+z)/(x+y) ) + y*log( y/(x+y) );
	}
	else
	{
		if( x == 0 || z == 0 || x + y + z <=0)
			return 0.0;
		else
			return 0.0;	
	}
}

/*
Calculating the force between segments with 1<r<2;
*/

double F_1_2(double x, double y, double z, double r)
{
	double t = 2-r;
	if( y >= 0)
		return ( (pow(y,t) - pow(x+y,t) - pow(y+z,t) + pow(x+y+z,t)) / (t*(1-r)) );
	else
	{
		if( x == 0 || z == 0 || x + y + z <=0)
			return 0.0;
		else
			return 0.0;	
	}
}

/*
Calculating the force between segments with r=2 (gravitational force);
*/

double F_2(double x, double y, double z, double r)
{
	if( y > 0)
		return log( (y + x) * (y + z) / (y * (y + x + z)) );
	else
	{
		if( x == 0 || z == 0 || x + y + z <=0)
			return 0.0;
		else
			return 0.0;	
	}
}

/*
Calculating the force between segments with r>2;
*/

double F_bigger_2(double x, double y, double z, double r)
{
	double t = 2-r;
	if( y > 0)
		return ( (pow(y,t) - pow(x+y,t) - pow(y+z,t) + pow(x+y+z,t)) / (t*(1-r)) );
	else
	{
		if( x == 0 || z == 0 || x + y + z <=0)
			return 0.0;
		else
			return 0.0;	
	}
}



