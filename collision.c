#include "collision.h"
#include "computeCellValues.h"
#include "LBDefinitions.h"

/** computes the post-collision distribution functions according to the BGK update rule and
 *  stores the results again at the same position.
 */
void computePostCollisionDistributions(double *currentCell, const double * const tau, const double *const feq){
	int i;
	for (i = 0; i < Q; i++){
		currentCell[i] = currentCell[i] - (1.0/ (*tau))*(currentCell[i]-feq[i]);
	}
}

/** carries out the whole local collision process. Computes density and velocity and
 *  equilibrium distributions. Carries out BGK update.
 */
void doCollision(double *collideField, int *flagField,const double * const tau,int xlength){

	int x,y,z ;
	int counter;
	double density;
	double velocity[3] ;
	double feq[Q] ;

    /* we loop over all the inner cells i.e. fluid cells */
	for (z = 1; z < xlength+1; z++) {
		for (y = 1; y < xlength+1; y++) {
			for (x = 1; x < xlength+1; x++) {
                /* get the current index */
				counter  = Q*(z*(xlength+2)*(xlength+2) + y * (xlength+2) + x );
				/* Compute density, velocity, f_eq to finally get the postcollision distribution */
				computeDensity (&collideField[counter], &density) ;
                computeVelocity(&collideField[counter], &density,velocity) ;
				computeFeq(&density,velocity,feq) ;
				computePostCollisionDistributions(&collideField[counter],tau,feq);
			}
		}
	}
}

