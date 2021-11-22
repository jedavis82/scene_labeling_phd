// HoF_Raster_mex.cpp
// 1/27/2011
// Author: Andrew Buck
//
// MATLAB mex wrapper class for raster histograms of forces.


#include "HoF_Raster.hpp"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	//Check for correct number of arguments
    if(nrhs != 2 && nrhs != 3)
        mexErrMsgTxt("Usage:\nhof_raster(imageA, imageB)\nhof_raster(imageA, imageB, numDir)");
    if(nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");

	//Get number of directions for evaluation
	int numDir = 180;
	if(nrhs == 3)
		numDir = (int)mxGetScalar(prhs[2]);

    //Get matrix sizes
	mwSize mA = mxGetM(prhs[0]);
	mwSize nA = mxGetN(prhs[0]);
	mwSize mB = mxGetM(prhs[1]);
	mwSize nB = mxGetN(prhs[1]);

	//Check to make sure matricies are the same size
	if(mA != mB || nA != nB)
		mexErrMsgTxt("Image dimensions must be equal.");

    //Construct image A
    unsigned char *imageA;
    imageA = (unsigned char *)mxCalloc(mA*nA, sizeof(unsigned char));
    if(mxIsUint8(prhs[0])) {
        //Already in uint8 format
        
        //Get pointer to input data as unsigned char array
        unsigned char *t_imgA = (unsigned char *)mxGetData(prhs[0]);
        
        //Transpose t_imgA to imageA
        int offset = 0;
        for(int i = 0; i < mA; i++) {
            for(int j = 0; j < nA; j++) {
                *(imageA + offset) = *(t_imgA + mA*j + i);
                offset++;
            }
        }
    }
    else if(mxIsDouble(prhs[0])) {
        //Convert to uint8 format
        
        //Get pointer to input data as double array
        double *t_imgA = mxGetPr(prhs[0]);

        //Transpose t_imgA to imageA
        int offset = 0;
        for(int i = 0; i < mA; i++) {
            for(int j = 0; j < nA; j++) {
                *(imageA + offset) = (unsigned char)(*(t_imgA + mA*j + i));
                offset++;
            }
        }
    }
    else
        mexErrMsgTxt("Images must be in double or uint8 format.");
    
    //Construct image B
    unsigned char *imageB;
    imageB = (unsigned char *)mxCalloc(mB*nB, sizeof(unsigned char));
    if(mxIsUint8(prhs[1])) {
        //Already in uint8 format
        
        //Get pointer to input data as unsigned char array
        unsigned char *t_imgB = (unsigned char *)mxGetData(prhs[1]);
        
        //Transpose t_imgB to imageB
        int offset = 0;
        for(int i = 0; i < mB; i++) {
            for(int j = 0; j < nB; j++) {
                *(imageB + offset) = *(t_imgB + mB*j + i);
                offset++;
            }
        }
    }
    else if(mxIsDouble(prhs[1])) {
        //Convert to uint8 format
        
        //Get pointer to input data as double array
        double *t_imgB = mxGetPr(prhs[1]);
        
        //Transpose t_imgB to imageB
        int offset = 0;
        for(int i = 0; i < mB; i++) {
            for(int j = 0; j < nB; j++) {
                *(imageB + offset) = (unsigned char)(*(t_imgB + mB*j + i));
                offset++;
            }
        }
    }
    else
        mexErrMsgTxt("Images must be in double or uint8 format.");
    
    //Prepare space for output histograms
	double *F0histogram = (double *)mxCalloc(numDir+1, sizeof(double));
    double *F2histogram = (double *)mxCalloc(numDir+1, sizeof(double));
    double *F02histogram = (double *)mxCalloc(numDir+1, sizeof(double));
    
    //Calculate HoF
	hof::HoF_Raster raster_obj;
    
// 	   raster_obj.FRHistogram_FuzzyRaster(F0histogram, numDir, 0, imageA, imageB, nA, mA, 1);
//     raster_obj.FRHistogram_FuzzyRaster(F2histogram, numDir, 2, imageA, imageB, nA, mA, 1);
//     raster_obj.F02Histogram_FuzzyRaster(F02histogram, numDir, imageA, imageB, nA, mA, 0.01, 3, 1);
    
    raster_obj.FRHistogram_CrispRaster(F0histogram, numDir, 0, imageA, imageB, nA, mA);
    raster_obj.FRHistogram_CrispRaster(F2histogram, numDir, 2, imageA, imageB, nA, mA);
    raster_obj.F02Histogram_CrispRaster(F02histogram, numDir, imageA, imageB, nA, mA, 0.01, 3);
    
    //Construct output histogram
	double *hof;
	plhs[0] = mxCreateDoubleMatrix(3, numDir+1, mxREAL);
    hof = mxGetPr(plhs[0]);
    int offset = 0;
    for(int i = 0; i < numDir+1; i++) {
        *(hof + offset++) = F0histogram[i];
        *(hof + offset++) = F2histogram[i];
        *(hof + offset++) = F02histogram[i];
    }
    
    //Free memory
    mxFree(imageA);
    mxFree(imageB);
    mxFree(F0histogram);
    mxFree(F2histogram);
    mxFree(F02histogram);
}