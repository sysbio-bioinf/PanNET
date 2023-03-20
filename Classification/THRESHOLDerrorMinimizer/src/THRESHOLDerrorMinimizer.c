#include "THRESHOLDerrorMinimizer.h"


void THRESHOLDerrorMinimizer(	double* data,
				double* weights,  
				double*	sortedFeature,
				double*	thresholds,  
				int* 	numFea,
				int* 	last_pos, 
				int* 	last_neg, 
				int*	numThresh,
                int*    testGEQ,
                int*    testLEQ,
                double* output)
{

//counter

int i;	//Iterator Thresholds
int j;	//Iterator Samples

int k;	//Iterator of feature



double* actWeight;	// Current weight

double* feature;	// Pointer to feature dimension

double* actSortedFea;

double* actFea;		// Pointer to current sample feature dimension

double* thresh;
double* actThreshold;

double* actOutput;	// Pointer to current position in output

double	actScore;	// Score of current classifier

double  bestScore;	 // Best score for feature combination
double  bestDir;	 // Begin of best interval (for feature combination)
double  bestThr;         // End of best interval (for feature combination)


//int 	verbose;


//verbose=0;

actOutput= output;

thresh = thresholds;

actSortedFea = sortedFeature;

	for (k =0; k < *numFea; k++)
	{
		feature      = &data[*last_neg*k];
		actSortedFea = sortedFeature;
		
		for(i=0; i < *last_neg; i++)
		{
			*actSortedFea++ = *feature++;
		}

		feature      = &data[*last_neg*k];
		actSortedFea = sortedFeature;

		R_rsort (actSortedFea,*last_neg);


		//printf("************************************\n");
		for(i=0; i < (*last_neg-1); i++)
		{
			thresh[i] = (actSortedFea[i] + actSortedFea[i+1])/2.0;
		//	printf("%f \n",thresh[i]);
		}
		
		thresh[(*last_neg-1)] = actSortedFea[0]-0.5;
		thresh[*last_neg]     = actSortedFea[(*last_neg-1)]+0.5;

		

		actThreshold = thresh;

		bestScore = 99999999999999999.9;
		bestDir   = 1.0;
		bestThr   = -1.0;    

      
        for(i=0; i < *numThresh; i++)
        {
            if(*testGEQ == 1)
            {
                //****************************************************
                //**
                //** Analyse positive direction

                actFea   = feature;
                actWeight = weights;

                actScore = 0;

                // Classify positive examples
                for (j =0; j < *last_pos; j++)
                {
                    //Misclassified Sample of class 1
                    if( *actFea++ <= *actThreshold)
                    {
                        actScore = actScore + (*actWeight);
                    }
		
                    actWeight++;
                }

                // Classify negative examples
                for (j =*last_pos; j < *last_neg; j++)
                {
                    if( *actFea++ > *actThreshold)
                    {
                        actScore = actScore + (*actWeight);
                    }

                    actWeight++;
                }

                // Updating currentComb and bestScoreComb
                if(actScore < bestScore)
                {
                    bestScore = actScore;
                    bestDir   = 1;
                    bestThr   = *actThreshold;
                    //printf("%f \n",bestThr);
                }
            }

            if(*testLEQ == 1)
            {
                //****************************************************
                //**
                //** Analyse negative direction

                actFea    = feature;
                actWeight = weights;

                actScore = 0;

                // Classify positive examples
                for (j =0; j < *last_pos; j++)
                {
                    //Misclassified Sample of class 1
                    if( *actFea++ >= *actThreshold)
                    {
                        actScore = actScore + (*actWeight);
                    }
		
                    actWeight++;
                }

                // Classify negative examples
                for (j =*last_pos; j < *last_neg; j++)
                {
                    if( *actFea++ < *actThreshold)
                    {
                        actScore = actScore + (*actWeight);
                    }

                    actWeight++;
                }

                // Updating currentComb and bestScoreComb
                if(actScore < bestScore)
                {
                    bestScore = actScore;
                    bestDir   = -1.0;
                    bestThr   = *actThreshold;
                    //printf("%f \n",bestThr);
                }
            }

			actThreshold++;
        }

		*actOutput++ = (double) k;//feature
		*actOutput++ = bestDir;//direction
		*actOutput++ = bestThr;//threshold
		*actOutput++ = bestScore;//score
	}//END for (k =0; k < (*numFea)-1; k++)

}

void THRESHOLDprediction(	double* data,
				int*    features, 
				int*    directions, 
				double*	thresholds, 
				int* 	numFea, 
				int* 	numSam, 
				int* 	numMod, 
				double* output)
{
	int i;
	int j;	

	double* actData;
	int* 	actFea;
	int*	actDir;
	double* actThr;
	double* actOut;

	actFea = features;
	actDir = directions;
	actThr = thresholds;
	actOut = output;

	for(i=0; i < *numMod; i++)
	{
		actData = &data[(*actFea++)*(*numSam)];

		if(*actDir++ == 1)
		{
			for(j =0; j < *numSam;j++)
			{
				if(*actData++ > *actThr)
				{
					*actOut = 1.0;
				}
				actOut++;
			}
		}else{
			for(j =0; j < *numSam;j++)
			{
				if(*actData++ < *actThr)
				{
					*actOut = 1.0;
				}
				actOut++;
			}
		}

		actThr++;
	}
}
