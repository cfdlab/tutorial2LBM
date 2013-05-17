#ifndef _LBDEFINITIONS_H_
#define _LBDEFINITIONS_H_
#define Q 19
#include <math.h>

/* Define constant values that will be used in the computations including Q and the Lattice weights and velocities */
  static const int LATTICEVELOCITIES[19][3]={{0,-1,-1},{-1,0,-1},{0,0,-1},{1,0,-1},
		  {0, 1, -1},{-1, -1, 0},{0, -1, 0},{1, -1, 0},{-1, 0, 0},{0, 0, 0},{1, 0, 0},
		  {-1, 1, 0},{0, 1, 0},{1, 1, 0},{0, -1, 1},{-1, 0, 1},{0,0,1},{1,0,1},{0,1,1}};
  static const double LATTICEWEIGHTS[19]={1/36.0, 1/36.0, 2/36.0, 1/36.0, 1/36.0, 1/36.0,
		  2/36.0, 1/36.0, 2/36.0, 12/36.0, 2/36.0, 1/36.0, 2/36.0, 1/36.0, 1/36.0, 1/36.0, 2/36.0,
		  1/36.0, 1/36.0};

  /* The following threw an error at compilation time so it was defined in the functions where C_S is used:*/
  static const double C_S = 0.57735026918963;

#endif

