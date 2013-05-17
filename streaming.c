#include "streaming.h"
#include "LBDefinitions.h"

/** carries out the streaming step and writes the respective distribution functions from
 *  collideField to streamField.
 */
void doStreaming(double *collideField, double *streamField,int *flagField,int xlength){
	int x ;
	int y;
	int z ;
	int i ;
	int currentCell;
    /* Loop through the inner cells */
	for (z = 1; z < xlength+1; z++ ) {
		for (y = 1; y < xlength+1; y++) {
			for (x = 1; x < xlength+1; x++) {
                /* Compute the index for the current cell */
				currentCell = Q*(z*(xlength+2)*(xlength+2) + y * (xlength+2) + x );
				for (i = 0; i < Q; i++) {
                    /* Carries out the streaming step. For each FLUID cell, the distributions fi
                     * from ALL neighbouring cells ~x + ~ci are copied from the collideField to the
                     * i-th position in the streamingField. */
					streamField[currentCell+i] = collideField[Q*(((z-LATTICEVELOCITIES[i][2])*(xlength+2)*(xlength+2)) +
							((y-LATTICEVELOCITIES[i][1]) * (xlength+2)) + (x-LATTICEVELOCITIES[i][0])) + i] ;
				}
			}
		}
	}
}

