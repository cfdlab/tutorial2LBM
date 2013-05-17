#include <stdio.h>
#include "visualLB.h"
#include "LBDefinitions.h"
#include "helper.h"
#include "computeCellValues.h"

/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. We re-used parts of the code
 *  from visual.c (VTK output for Navier-Stokes solver) and modified it for 3D datasets.
 */
void writeVtkOutput(const double * const collideField,
		const int * const flagField,
		const char* filename,
		unsigned int t, int xlength){
	int x, y, z;
	char szFileName[200];
	FILE *fp=NULL;
	int counter;
	double density;
	double velocity[3] ;

	/* Create the new vtk file */
	sprintf( szFileName, "%s.%i.vtk",filename, t );
	fp = fopen( szFileName, "w");
	if( fp == NULL )
	{
		char szBuff[80];
		sprintf( szBuff, "Failed to open %s", szFileName );
		ERROR( szBuff );
		return;
	}

	/* Write the VTK file header information and the geometry information */
	write_vtkHeader( fp, xlength, xlength, xlength);

	/* Write the velocity vectors to the VTK file*/
	fprintf(fp,"\nPOINT_DATA %i \n", (xlength+2)*(xlength+2)*(xlength+2) );
	fprintf(fp, "VECTORS velocity float\n");
	for(z = 0; z < xlength+2; z++) {
		for(y = 0; y < xlength+2; y++) {
			for(x = 0; x < xlength+2; x++) {
				if(x!=0 && x!=xlength+1 && y!=0 && y!=xlength+1 && z!=0 && z!=xlength+1){
					/* Get the index for current cell */
					counter  = Q*(z*(xlength+2)*(xlength+2) + y * (xlength+2) + x );
					/* Compute the velocity of the current cell */
					computeDensity (&collideField[counter] , &density) ;
					computeVelocity(&collideField[counter], &density,velocity) ;
					/* Print the values to the file */
					fprintf(fp, "%f %f %f\n", velocity[0], velocity[1] , velocity[2]);
				}
				/*else if(z == (xlength + 1)){
					fprintf(fp, "1 0 0\n");
				}*/
				else{
					fprintf(fp, "0 0 0\n");
				}
			}
		}
	}

	/* Write the density values to the vtk file */
	fprintf(fp,"\n");
	fprintf(fp, "SCALARS density double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(z = 0; z < xlength+2; z++) {
		for(y = 0; y < xlength+2; y++) {
			for(x = 0; x < xlength+2; x++) {
				if(x!=0 && x!=xlength+1 && y!=0 && y!=xlength+1 && z!=0 && z!=xlength+1){
					/* Get the index for current cell */
					counter  = Q*(z*(xlength+2)*(xlength+2) + y * (xlength+2) + x );
					/* Compute the density of the current cell */
					computeDensity (&collideField[counter] , &density) ;
					/* Print the value to the file */
					fprintf(fp, "%f\n", density);
				}
				else{
					fprintf(fp, "0\n");
				}
			}
		}
	}

	/* Write the flagfield for each cell*/
	fprintf(fp,"\n");
	fprintf(fp, "SCALARS flagfield int 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for(z = 0; z < xlength+2; z++) {
		for(y = 0; y < xlength+2; y++) {
			for(x = 0; x < xlength+2; x++) {
				/* Get the index for current cell */
				counter  = ((z * (xlength+2) * (xlength+2) + y * (xlength+2) + x));
				/* Print the value to the file */
				fprintf(fp, "%i\n", flagField[counter]);
			}
		}
	}

	/* Try to close file and show an error message if it fails.  */
	if( fclose(fp) )
	{
		char szBuff[80];
		sprintf( szBuff, "Failed to close %s", szFileName );
		ERROR( szBuff );
	}
}

/* auxiliary function to write the header and the geometry for the vtk file.*/
void write_vtkHeader( FILE *fp, int xlength, int ylength, int zlength ) {
	if( fp == NULL )
	{
		char szBuff[80];
		sprintf( szBuff, "Null pointer in write_vtkHeader" );
		ERROR( szBuff );
		return;
	}

	fprintf(fp,"# vtk DataFile Version 2.0\n");
	fprintf(fp,"generated for CFD-lab course output (Based in code by Tobias Neckel) \n");
	fprintf(fp,"ASCII\n");
	fprintf(fp,"\n");
	fprintf(fp,"DATASET STRUCTURED_POINTS\n");
	fprintf(fp,"DIMENSIONS  %i %i %i \n", xlength+2, ylength + 2, zlength + 2);
	fprintf(fp,"ORIGIN 0 0 0\n");
	fprintf(fp,"SPACING 1 1 1\n");
	fprintf(fp,"\n");
}


