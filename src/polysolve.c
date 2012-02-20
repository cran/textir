#include "latools.h"
#include <assert.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

double cubed(double x)
{ return x*x*x; }

double cubeRoot(double x)
{
    if (x < 0)
        return -pow(-x, 1.0/3.0);
    else
        return pow(x, 1.0/3.0);
}

/* x^3 + ax^2 + bx +c trig method solution */
double* solvecubic(double a, double b, double c)
{
  double *roots = new_dvec(4);
  roots[0] = 0;
  double j, k, l, m, n, o, p, q, r;

  p = (3.0*b  - a*a)/3.0;
  q = (2.0*a*a*a  - 9.0*a*b)/27.0  + c;
  r = q*q/4.0 + p*p*p/27.0;
  
  if (p == 0 && q == 0)
    { // one single root
      roots[0] = 1;
      roots[1] = roots[2] = roots[3] = -cubeRoot(c); }
  else if (r > 0)
    {
      // 1 real root and 2 complex roots     
      l = cubeRoot( -q/2.0 + sqrt(r) );
      m = cubeRoot( -q/2.0 - sqrt(r) );
      n = -a/3.0;
      
      roots[0] = 2;
      roots[1] = (l + m) + n;
      roots[2] = -(l + m)/2 + n; // real part
      roots[3] = (l - m)*sqrt(3)*0.5;  // plus or minus imaginary part 'i'
    }
  else if (r <= 0)
    {
      // 3 real roots
      j = sqrt(q*q/4.0 - r);
      k = cubeRoot(j);
      l = acos(-(q / (2.0 * j)));
      m = cos(l / 3.0);
      n = sqrt(3) * sin(l/3.0);
      o = -a/3.0;

      roots[0] = 3;
      roots[1] = 2.0*k*m + o;
      roots[2] = -k*(m + n) + o;
      roots[3] = -k*(m - n) + o;
    }
  return roots;
}

/* x^2 + bx + c root finder 
   roots[0] is the number of real roots */
double *solvequadratic(double b, double c)
{
  double *roots = new_dvec(3);
  double q = b*b - 4*c;
  if(q==0){ // one real
    roots[0] = 1.0;
    roots[1] = roots[2] = -0.5*b;
  }
  else if(q>0){ // two real
    roots[0] = 2.0;
    roots[1] = -0.5*(b + sqrt(q));
    roots[2] = -0.5*(b - sqrt(q));
  }
  else{ // two complex
    roots[0]  = 0.0;
    roots[1] = -0.5*b;
    roots[2] = -0.5*sqrt(-q);
  }
  return roots;
}

void Rcubic(double *coef, double *num)
{
  double *roots = solvecubic(coef[0], coef[1], coef[2]);
  *num = roots[0];
  coef[0] = roots[1];
  coef[1] = roots[2];
  coef[2] = roots[3];
  free(roots);
}
  

void Rquadratic(double *coef, double *num)
{
  double *roots = solvequadratic(coef[0], coef[1]);
  *num = roots[0];
  coef[0] = roots[1];
  coef[1] = roots[2];
  free(roots);
}
 
