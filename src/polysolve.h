/* simple tools for finding cubic and quadratic equation roots */

#ifndef __POLYSOLVE_H__
#define __POLYSOLVE_H__ 

#include <stdio.h>
#include <stdlib.h>

double cubed(double x);
double cubeRoot(double x);
double* solvecubic(double a, double b, double c);
double* solvequadratic(double b, double c);
void Rcubic(double *coef, double *num);
void Rquadratic(double *coef, double *num);

#endif

