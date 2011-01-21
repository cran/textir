/*  Linear algebra tools (wrappers for lapack and blas) 
 *  and some basic array tools.
 *  Everything is COLUMN MAJOR for communication with fortran.
 *  
 */

#ifndef __CUBIC_H__
#define __CUBIC_H__ 

#include <stdio.h>
#include <stdlib.h>

double cubed(double x);
double cubeRoot(double x);
double* solvecubic(double a, double b, double c);
void Rcubic(double *coef, double *num);

#endif

