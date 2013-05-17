#include "boundary.h"
#include "LBDefinitions.h"
#include "computeCellValues.h"
#include <stdio.h>

/* treatBoundary
 Carries out the boundary treatment. Therefore, we loop over the outer boundary cells, check
 each cell for its state (NO SLIP or MOVING WALL), and set the respective distribution
 functions inside this cell according to Eq.(16) and Eq.(18).
 */
void treatBoundary(double *collideField, int* flagField, const double * const wallVelocity, int xlength){
	int x, y, z;
	int i;
	double density;
	int counter;
	int counter_wall;
    
	/*The following array are the values for i that are pointing in the direction of the inner cells
	 * for example in the border where y = 0, the c vectors pointing inside must be the ones which
	 * y component is = 1, in this case the ones with index 4, 11, 12, 13, 18, and the array is denoted
	 * as y_min. The same case is for when y = xlength+1 and the array is y_max, and similarly for the
	 * other 4 planes that form the border of the cube. */
	int y_min[5]={4, 11, 12, 13, 18};
	int y_max[5]={14, 7, 6, 5, 0};
	int x_min[5]={3, 7, 10, 13, 17};
	int x_max[5]={15, 11, 8, 5, 1};
	int z_min[5]={14, 15, 16, 17, 18};
	int z_max[5]={4, 3, 2, 1, 0};
    
	/*This is divided in 6 similar for loops for each of the planes, starting with the plane where y = 0*/
	y = 0;
	for(z = 0; z < xlength + 1; z++){
		for(x = 0; x < xlength + 2; x++){
			/* Store an index for the current cell */
			counter = (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x);
			/* Check if the boundary condition is NO SLIP*/
			if(flagField[counter]==1){
				for(i = 0; i<5; i++){
					/*Check that the updated directions are pointing to the inner cells*/
					if(((x+LATTICEVELOCITIES[y_min[i]][0]) > 0) && ((z+LATTICEVELOCITIES[y_min[i]][2]) > 0) &&
                       ((x+LATTICEVELOCITIES[y_min[i]][0])< xlength+1)&&((z+LATTICEVELOCITIES[y_min[i]][2])< xlength+1)){
						collideField[Q*counter+y_min[i]]= collideField[Q*((z+LATTICEVELOCITIES[y_min[i]][2])*
                                                                          (xlength+2) * (xlength+2) + (y+LATTICEVELOCITIES[y_min[i]][1]) * (xlength+2) +
                                                                          (x+LATTICEVELOCITIES[y_min[i]][0]))+(Q-y_min[i]-1)] ;
					}
				}
			}
			if(flagField[counter]==2){
				for(i = 0; i < 5; i++){
					/* Check that there exists an inner cell in that direction*/
					if(((x+LATTICEVELOCITIES[y_min[i]][0]) > 0) && ((z+LATTICEVELOCITIES[y_min[i]][2]) > 0) &&
                       ((x+LATTICEVELOCITIES[y_min[i]][0])< xlength+1)&&((z+LATTICEVELOCITIES[y_min[i]][2])< xlength+1)){
						/*treat the boundary as MOVING WALL according to Eq. (18)*/
						counter_wall  = Q*((z+LATTICEVELOCITIES[y_min[i]][2]) * (xlength+2) * (xlength+2) +
                                           (y+LATTICEVELOCITIES[y_min[i]][1]) * (xlength+2) + (x+LATTICEVELOCITIES[y_min[i]][0]));
						computeDensity (&collideField[counter_wall] , &density) ;
						collideField[Q*counter+y_min[i]]= collideField[counter_wall+(Q-y_min[i]-1)]+
                        2*LATTICEWEIGHTS[y_min[i]]*density/(C_S*C_S)*((LATTICEVELOCITIES[y_min[i]][0]*wallVelocity[0])+
                                                                      (LATTICEVELOCITIES[y_min[i]][1]*wallVelocity[1])+(LATTICEVELOCITIES[y_min[i]][2]*wallVelocity[2]));
                        
					}
				}
			}
		}
	}
    
	y = xlength+1;
	for(z = 0; z < xlength + 1; z++){
		for(x = 0; x < xlength + 2; x++){
			counter = (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x);
			if(flagField[counter]==1){
				for(i = 0; i<5; i++){
					if(((x+LATTICEVELOCITIES[y_max[i]][0]) > 0) && ((z+LATTICEVELOCITIES[y_max[i]][2]) > 0) &&
                       ((x+LATTICEVELOCITIES[y_max[i]][0])< xlength+1)&&((z+LATTICEVELOCITIES[y_max[i]][2])< xlength+1)){
						collideField[Q*counter+y_max[i]]= collideField[Q*((z+LATTICEVELOCITIES[y_max[i]][2])*
                                                                          (xlength+2) * (xlength+2) + (y+LATTICEVELOCITIES[y_max[i]][1]) * (xlength+2) +
                                                                          (x+LATTICEVELOCITIES[y_max[i]][0]))+(Q-y_max[i]-1)] ;
					}
				}
			}
			if(flagField[counter]==2){
				for(i = 0; i < 5; i++){
					if(((x+LATTICEVELOCITIES[y_max[i]][0]) > 0) && ((z+LATTICEVELOCITIES[y_max[i]][2]) > 0) &&
                       ((x+LATTICEVELOCITIES[y_max[i]][0])< xlength+1)&&((z+LATTICEVELOCITIES[y_max[i]][2])< xlength+1)){
						counter_wall  = Q*((z+LATTICEVELOCITIES[y_max[i]][2]) * (xlength+2) * (xlength+2) +
                                           (y+LATTICEVELOCITIES[y_max[i]][1]) * (xlength+2) + (x+LATTICEVELOCITIES[y_max[i]][0]));
						computeDensity (&collideField[counter_wall] , &density) ;
						collideField[Q*counter+y_max[i]]= collideField[counter_wall+(Q-y_max[i]-1)]+
                        2*LATTICEWEIGHTS[y_max[i]]*density/(C_S*C_S)*((LATTICEVELOCITIES[y_max[i]][0]*wallVelocity[0])+
                                                                      (LATTICEVELOCITIES[y_max[i]][1]*wallVelocity[1])+(LATTICEVELOCITIES[y_max[i]][2]*wallVelocity[2]));
                        
					}
				}
			}
		}
	}
    
	x = 0;
	for(z = 0; z < xlength + 1; z++){
		for(y = 0; y < xlength + 2; y++){
			counter = (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x);
			if(flagField[counter]==1){
				for(i = 0; i<5; i++){
					if(((y+LATTICEVELOCITIES[x_min[i]][1]) > 0) && ((z+LATTICEVELOCITIES[x_min[i]][2]) > 0) &&
                       ((y+LATTICEVELOCITIES[x_min[i]][1])< xlength+1)&&((z+LATTICEVELOCITIES[x_min[i]][2])< xlength+1)){
						collideField[Q*counter+x_min[i]]= collideField[Q*((z+LATTICEVELOCITIES[x_min[i]][2])*
                                                                          (xlength+2) * (xlength+2) + (y+LATTICEVELOCITIES[x_min[i]][1]) * (xlength+2) +
                                                                          (x+LATTICEVELOCITIES[x_min[i]][0]))+(Q-x_min[i]-1)] ;
					}
				}
			}
			if(flagField[counter]==2){
				for(i = 0; i < 5; i++){
					/* Check that there exists a neighbor cell in that direction*/
					if(((y+LATTICEVELOCITIES[x_min[i]][1]) > 0) && ((z+LATTICEVELOCITIES[x_min[i]][2]) > 0) &&
                       ((y+LATTICEVELOCITIES[x_min[i]][1])< xlength+1)&&((z+LATTICEVELOCITIES[x_min[i]][2])< xlength+1)){
						/*treat the boundary as MOVING WALL according to Eq. (18)*/
						counter_wall  = Q*((z+LATTICEVELOCITIES[x_min[i]][2]) * (xlength+2) * (xlength+2) +
                                           (y+LATTICEVELOCITIES[x_min[i]][1]) * (xlength+2) + (x+LATTICEVELOCITIES[x_min[i]][0]));
						computeDensity (&collideField[counter_wall] , &density) ;
						collideField[Q*counter+x_min[i]]= collideField[counter_wall+(Q-x_min[i]-1)]+
                        2*LATTICEWEIGHTS[x_min[i]]*density/(C_S*C_S)*((LATTICEVELOCITIES[x_min[i]][0]*wallVelocity[0])+
                                                                      (LATTICEVELOCITIES[x_min[i]][1]*wallVelocity[1])+(LATTICEVELOCITIES[x_min[i]][2]*wallVelocity[2]));
                        
					}
				}
			}
		}
	}
    
	x = xlength+1;
	for(z = 0; z < xlength + 1; z++){
		for(y = 0; y < xlength + 2; y++){
			counter = (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x);
			if(flagField[counter]==1){
				for(i = 0; i<5; i++){
					if(((y+LATTICEVELOCITIES[x_max[i]][1]) > 0) && ((z+LATTICEVELOCITIES[x_max[i]][2]) > 0) &&
                       ((y+LATTICEVELOCITIES[x_max[i]][1])< xlength+1)&&((z+LATTICEVELOCITIES[x_max[i]][2])< xlength+1)){
						collideField[Q*counter+x_max[i]]= collideField[Q*((z+LATTICEVELOCITIES[x_max[i]][2])*
                                                                          (xlength+2) * (xlength+2) + (y+LATTICEVELOCITIES[x_max[i]][1]) * (xlength+2) +
                                                                          (x+LATTICEVELOCITIES[x_max[i]][0]))+(Q-x_max[i]-1)] ;
					}
				}
			}
			if(flagField[counter]==2){
				for(i = 0; i < 5; i++){
					/* Check that there exists a neighbor cell in that direction*/
					if(((y+LATTICEVELOCITIES[x_max[i]][1]) > 0) && ((z+LATTICEVELOCITIES[x_max[i]][2]) > 0) &&
                       ((y+LATTICEVELOCITIES[x_max[i]][1])< xlength+1)&&((z+LATTICEVELOCITIES[x_max[i]][2])< xlength+1)){
						/*treat the boundary as MOVING WALL according to Eq. (18)*/
						counter_wall  = Q*((z+LATTICEVELOCITIES[x_max[i]][2]) * (xlength+2) * (xlength+2) +
                                           (y+LATTICEVELOCITIES[x_max[i]][1]) * (xlength+2) + (x+LATTICEVELOCITIES[x_max[i]][0]));
						computeDensity (&collideField[counter_wall] , &density) ;
						collideField[Q*counter+x_max[i]]= collideField[counter_wall+(Q-x_max[i]-1)]+
                        2*LATTICEWEIGHTS[x_max[i]]*density/(C_S*C_S)*((LATTICEVELOCITIES[x_max[i]][0]*wallVelocity[0])+
                                                                      (LATTICEVELOCITIES[x_max[i]][1]*wallVelocity[1])+(LATTICEVELOCITIES[x_max[i]][2]*wallVelocity[2]));
                        
					}
				}
			}
		}
	}
    
	z = 0;
	for(y = 0; y < xlength + 2; y++){
		for(x = 0; x < xlength + 2; x++){
			counter = (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x);
			if(flagField[counter]==1){
				for(i = 0; i<5; i++){
					if(((x+LATTICEVELOCITIES[z_min[i]][0]) > 0) && ((y+LATTICEVELOCITIES[z_min[i]][1]) > 0) &&
                       ((x+LATTICEVELOCITIES[z_min[i]][0])< xlength+1)&&((y+LATTICEVELOCITIES[z_min[i]][1])< xlength+1)){
						collideField[Q*counter+z_min[i]]= collideField[Q*((z+LATTICEVELOCITIES[z_min[i]][2])*
                                                                          (xlength+2) * (xlength+2) + (y+LATTICEVELOCITIES[z_min[i]][1]) * (xlength+2) +
                                                                          (x+LATTICEVELOCITIES[z_min[i]][0]))+(Q-z_min[i]-1)] ;
						if(collideField[Q*counter+z_max[i]]<0 ||collideField[Q*counter+z_max[i]] > 2 ){
						printf("x=%i y=%i z=%i value=%f\n",x,y,z,collideField[Q*counter+z_min[i]]);
						}
					}
				}
			}
			
			if(flagField[counter]==2){
				for(i = 0; i < 5; i++){
					/* Check that there exists a neighbor cell in that direction*/
					if(((x+LATTICEVELOCITIES[z_min[i]][0]) > 0) && ((y+LATTICEVELOCITIES[z_min[i]][1]) > 0) &&
                       ((x+LATTICEVELOCITIES[z_min[i]][0])< xlength+1)&&((y+LATTICEVELOCITIES[z_min[i]][1])< xlength+1)){
						/*treat the boundary as MOVING WALL according to Eq. (18)*/
						counter_wall  = Q*((z+LATTICEVELOCITIES[z_min[i]][2]) * (xlength+2) * (xlength+2) +
                                           (y+LATTICEVELOCITIES[z_min[i]][1]) * (xlength+2) + (x+LATTICEVELOCITIES[z_min[i]][0]));
						computeDensity (&collideField[counter_wall] , &density) ;
						collideField[Q*counter+z_min[i]]= collideField[counter_wall+(Q-z_min[i]-1)]+
                        2*LATTICEWEIGHTS[z_min[i]]*density/(C_S*C_S)*((LATTICEVELOCITIES[z_min[i]][0]*wallVelocity[0])+
                                                                      (LATTICEVELOCITIES[z_min[i]][1]*wallVelocity[1])+(LATTICEVELOCITIES[z_min[i]][2]*wallVelocity[2]));
                        if(collideField[Q*counter+z_max[i]]<0 ||collideField[Q*counter+z_max[i]] > 2 ){
						printf("x=%i y=%i z=%i value=%f\n",x,y,z,collideField[Q*counter+z_min[i]]);
						}
					}
				}
			}
		}
	}
    
	z = xlength+1;
	for(y = 0; y < xlength + 2; y++){
		for(x = 0; x < xlength + 2; x++){
			counter = (z * (xlength+2) * (xlength+2) + y * (xlength+2) + x);
			if(flagField[counter]==1){
				for(i = 0; i<5; i++){
					if(((x+LATTICEVELOCITIES[z_max[i]][0]) > 0) && ((y+LATTICEVELOCITIES[z_max[i]][1]) > 0) &&
                       ((x+LATTICEVELOCITIES[z_max[i]][0])< xlength+1)&&((y+LATTICEVELOCITIES[z_max[i]][1])< xlength+1)){
						collideField[Q*counter+z_max[i]]= collideField[Q*((z+LATTICEVELOCITIES[z_max[i]][2])*
                                                                          (xlength+2) * (xlength+2) + (y+LATTICEVELOCITIES[z_max[i]][1]) * (xlength+2) +
                                                                          (x+LATTICEVELOCITIES[z_max[i]][0]))+(Q-z_max[i]-1)] ;
					}
				}
			}
			if(flagField[counter]==2){
				for(i = 0; i < 5; i++){
					/* Check that there exists a neighbor cell in that direction*/
					if(((x+LATTICEVELOCITIES[z_max[i]][0]) > 0) && ((y+LATTICEVELOCITIES[z_max[i]][1]) > 0) &&
                       ((x+LATTICEVELOCITIES[z_max[i]][0]) < xlength+1)&&((y+LATTICEVELOCITIES[z_max[i]][1])< xlength+1)){
						/*treat the boundary as MOVING WALL according to Eq. (18)*/
						counter_wall  = Q*((z+LATTICEVELOCITIES[z_max[i]][2]) * (xlength+2) * (xlength+2) +
                                           (y+LATTICEVELOCITIES[z_max[i]][1]) * (xlength+2) + (x+LATTICEVELOCITIES[z_max[i]][0]));
						computeDensity (&collideField[counter_wall] , &density) ;						
						collideField[Q*counter+z_max[i]]= collideField[counter_wall+(Q-z_max[i]-1)]+
                        2*LATTICEWEIGHTS[z_max[i]]*density/(C_S*C_S)*((LATTICEVELOCITIES[z_max[i]][0]*wallVelocity[0])+
                                                                      (LATTICEVELOCITIES[z_max[i]][1]*wallVelocity[1])+(LATTICEVELOCITIES[z_max[i]][2]*wallVelocity[2]));					
					}
				}
			}
		}
	}
}
