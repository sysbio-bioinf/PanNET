#ifndef THRESHOLDerrorMinimizer_H
#define THRESHOLDerrorMinimizer_H

#include<R.h>
#include<Rinternals.h>




void THRESHOLDerrorMinimizer(	double* data,
				double* weights, 
				double* sortedFeature, 
				double*	thresholds, 
				int* 	numFea, 
				int* 	last_pos, 
				int* 	last_neg, 
				int*	numThresh,
                int*    testGEQ,
                int*    testLEQ,
				double* output);



void THRESHOLDprediction(	double* data,
				int*    features, 
				int*    directions, 
				double*	thresholds, 
				int* 	numFea, 
				int* 	numSam, 
				int* 	numMod, 
				double* output);

#endif

