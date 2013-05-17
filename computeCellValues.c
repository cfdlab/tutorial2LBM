#include "computeCellValues.h"
#include "LBDefinitions.h"
#include <math.h>

/** computes the density from the particle distribution functions stored at currentCell.
 *  currentCell thus denotes the address of the first particle distribution function of the
 *  respective cell. The result is stored in density.
 */
void computeDensity(const double *const currentCell, double *density){
	double sum = 0.0 ;
	int i ;
    /* The density is the sum of the distribution functions across all directions in a particular cell  */
	for (i = 0; i < Q; i++) {
		sum = sum + currentCell[i] ;
	}
	*density = sum ;
}

/** computes the velocity according to Eq.(9) within currentCell and stores the result in velocity */
void computeVelocity(const double * const currentCell, const double * const density, double *velocity){

	int i, j;
    double sum[3];

	for (j = 0; j < 3; j++) {
		sum[j]  = 0.0;
		for (i = 0; i < Q; i++) {
			sum[j] = sum[j]  +  currentCell[i] * LATTICEVELOCITIES[i][j] ;
		}
	}
	for (j = 0; j < 3; j++) {
		velocity[j] = sum[j] /	(*density) ;
	}
}

/** computes the equilibrium distributions for all particle distribution functions of one
 *  cell from density and velocity and stores the results in feq.
 */
void computeFeq(const double * const density, const double * const velocity, double *feq){

	int i, j;
	double innerProductCU;
	double innerProductUU;

    /* First we calculate the innerproducts between c.u and u.u */
	
    /* Then we solve the equation for f_eq */
	for (i = 0; i < Q; i++) {
        innerProductCU = 0;
        innerProductUU = 0;
		for(j = 0; j < 3; j++){
            innerProductCU = innerProductCU + velocity[j]*(double)LATTICEVELOCITIES[i][j];
        }
        for(j = 0; j < 3; j++){
            innerProductUU = innerProductUU + velocity[j]*velocity[j];
        }
        feq[i] = LATTICEWEIGHTS[i] * (*density) * ( 1 + innerProductCU/(C_S*C_S)+
				(innerProductCU*innerProductCU)/(2*C_S*C_S*C_S*C_S)- innerProductUU/(2*C_S*C_S));
	}
}

