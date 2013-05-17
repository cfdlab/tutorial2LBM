#ifndef _MAIN_C_
#define _MAIN_C_


#include <stdio.h>
#include <stdlib.h>
#include "LBDefinitions.h"
#include "collision.h"
#include "streaming.h"
#include "initLB.h"
#include "helper.h"
#include "visualLB.h"
#include "boundary.h"
#include "math.h"


int main (int argc, char *argv[]){
	double *collideField=NULL;
	double *streamField=NULL;
	int *flagField=NULL;
	int xlength;
	double tau;
	double velocityWall[3];
	int timesteps;
	int timestepsPerPlotting;
	int t;

	if(readParameters(&xlength, &tau, velocityWall, &timesteps, &timestepsPerPlotting, argc , argv[1])==1){

		/* Allocate memory for the collide, stream and flag fields */
		collideField = (double *)  malloc((size_t)( Q *(xlength+2)*(xlength+2) *(xlength+2)* sizeof( double )));
		streamField = (double *)  malloc((size_t)( Q *(xlength+2)*(xlength+2) *(xlength+2)* sizeof( double )));
		flagField = (int *) malloc((size_t)(xlength+2)*(xlength+2) *(xlength+2)* sizeof( int ));


		/* Initialise the fields with lattice weights and with the corresponding flags and check that there was no errors*/

		initialiseFields(collideField,streamField,flagField,xlength);

		/* Run this cycle for the number of timesteps required */
		for(t = 0; t < timesteps; t++){
			/* Create a temporary pointer to store swap the stream and collide pointers */
			double *swap=NULL;
			/* Do the streaming step using the collide field as input */
			doStreaming(collideField,streamField,flagField,xlength);
			/* Swap the streaming field with the collide field */
			swap = collideField;
			collideField = streamField;
			streamField = swap;
			/* Do the collision step */
			doCollision(collideField,flagField,&tau,xlength);
			/* Do the boundary treatment */
			treatBoundary(collideField,flagField,velocityWall,xlength);
			/* Create the output file depending on how many timesteps are defined */
			if (t%timestepsPerPlotting==0){
				writeVtkOutput(collideField,flagField,argv[0],t,xlength);
			}
		}

		/*Kill the pointers*/
		free(collideField);
		free(streamField);
		free(flagField);
	}
return 0;
}

#endif

