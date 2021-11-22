//-------------------------------------------------------------------------
// F02Histogram_CrispRaster.cpp
// 2/11/2011
// Author: Andrew Buck
//
// MATLAB MEX-file wrapper for the F-hybrid histogram in the crisp raster
// case.  Intended for use with the MATLAB wrapper function 'hof_raster.m'.
//
// WARNING!  ALL TYPE CHECKING IS PERFORMED IN 'hof_raster.m'  USE THIS
// FUNCTION ONLY WITH DATA FORMATTED IN THIS STYLE!
//-------------------------------------------------------------------------

#include "HoF_Raster.hpp"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//Get values from MATLAB
	unsigned char *imageA = (unsigned char *)mxGetData(prhs[0]);
	unsigned char *imageB = (unsigned char *)mxGetData(prhs[1]);
	int numberDirections = (int)mxGetScalar(prhs[2]);
	double p0 = (double)mxGetScalar(prhs[3]);
	double p1 = (double)mxGetScalar(prhs[4]);

	//Get matrix sizes
	int M = mxGetM(prhs[0]);
	int N = mxGetN(prhs[0]);

	//Prepare space for output histogram
	plhs[0] = mxCreateDoubleMatrix(1, numberDirections+1, mxREAL);
	double *histogram = mxGetPr(plhs[0]);

	//Calculate HoF
	hof::HoF_Raster raster_obj;
	raster_obj.F02Histogram_CrispRaster(histogram, numberDirections, imageA, imageB, M, N, p0, p1);
}