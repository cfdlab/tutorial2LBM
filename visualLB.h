#ifndef _VISUALLB_H_
#define _VISUALLB_H_


/** writes the density and velocity field (derived from the distributions in collideField)
 *  to a file determined by 'filename' and timestep 't'. */
void writeVtkOutput(const double * const collideField,
		const int * const flagField,
		const char *filename,
		unsigned int t, int xlength);

/* auxiliary function to write the header and the geometry for the vtk file.*/
void write_vtkHeader( FILE *fp, int xlength, int ylength, int zlength);


#endif

