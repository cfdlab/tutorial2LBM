#include "initLB.h"
#include "helper.h"
#include "LBDefinitions.h"

/* reads the parameters for the lid driven cavity scenario from a config file */
int readParameters(
		int *xlength,
		double *tau,
		double *velocityWall,
		int *timesteps,
		int *timestepsPerPlotting,
		int argc,
		char *argv
){
	double velocityWallx;
	double velocityWally;
	double velocityWallz;
	/* Check if there is one and only one input argument which should be the data file  */
	if(argc==2){
		/* Read the values */
		READ_INT( argv, *xlength );
		READ_INT ( argv, *timesteps);
		READ_INT( argv, *timestepsPerPlotting);
		READ_DOUBLE( argv, *tau );
		/* Since the velocity is a vector of 1 x 3, we read the three different values and then save them in one array */
		READ_DOUBLE( argv, velocityWallx);
		READ_DOUBLE( argv, velocityWally);
		READ_DOUBLE( argv, velocityWallz);
		velocityWall[0]=velocityWallx;
		velocityWall[1]=velocityWally;
		velocityWall[2]=velocityWallz;
	}
	else{
		/* In case there was only one argument print an error and return 0*/
		ERROR("One and only one argument (data file) should be passed to the function" );
		return 0;
	}
	return 1;
}

/* Initialises the particle distribution function fields collideField, flagField and streamField. */
void initialiseFields(double *collideField, double *streamField, int *flagField, int xlength){
	/*i-th distribution function in the cell (x, y, z) is (Q * (z * xlength * xlength + y * xlength + x)) + i; */
	int i, x, y, z;

	/* We initialize the flagField values directly here, checking in which boundary they are. */
	for (z = 0; z < xlength + 2; z++){
		for (y = 0; y < xlength + 2; y++){
			for (x = 0; x < xlength + 2; x++){
				flagField[((z * (xlength+2) * (xlength+2) + y * (xlength+2) + x))] = 0;
			}
		}
	}


	for (z = 0; z < xlength + 2; z++){
		for (y = 0; y < xlength + 2; y++){
			flagField[((z * (xlength+2) * (xlength+2) + y * (xlength+2) ))] = 1;
		}
	}
	for (z = 0; z < xlength + 2; z++){
		for (y = 0; y < xlength + 2; y++){
			flagField[((z * (xlength+2) * (xlength+2) + y * (xlength+2) + (xlength+1)))] = 1;
		}
	}
	for (z = 0; z < xlength + 2; z++){
		for (x = 0; x < xlength + 2; x++){
			flagField[((z * (xlength+2) * (xlength+2) + x))] = 1;
		}
	}

	for (z = 0; z < xlength + 2; z++){
		for (x = 0; x < xlength + 2; x++){
			flagField[((z * (xlength+2) * (xlength+2) + (xlength+1) * (xlength+2) + x))] = 1;
		}
	}

	for (y = 0; y < xlength + 2; y++){
		for (x = 0; x < xlength + 2; x++){
			flagField[((y * (xlength+2) + x))] = 1;
		}
	}

	for (y = 0; y < xlength + 2; y++){
		for (x = 0; x < xlength + 2; x++){
			flagField[(((xlength+1) * (xlength+2) * (xlength+2) + y * (xlength+2) + x))] = 2;
		}
	}

	for (z = 0; z < xlength + 2; z++){
		for (y = 0; y < xlength + 2; y++){
			for (x = 0; x < xlength + 2; x++){

				for (i = 0; i < Q; i++){
					/* Initialize the fields to the values of the Lattice weights */
					collideField[(Q * (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x)) + i] = LATTICEWEIGHTS[i];
					streamField[(Q * (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x)) + i] = LATTICEWEIGHTS[i];
				}
			}
		}
	}




}


