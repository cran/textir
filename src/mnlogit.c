#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
#include "latools.h"
#include "rhelp.h"
#include "cubic.h"
# ifdef _OPENMP
#include <omp.h>
# endif


int n;
int p;
int d;
int *m;
double *X;
int **xi;
int N;
double **V;
double **H;
double **D;
double **eta;
double *denom;
double **B;
double **G;
double Bsum;
double nregpar;

double lambda;
double *lampar;
int maplam;

/* un-normalized negative log posterior */

double neglogpost(){
  
  double L = 0.0;
  int i;
  
  #pragma omp parallel for private(i) reduction(+: L) 
  for(i=0; i<N; i++) L += -X[i]*(eta[xi[1][i]][xi[0][i]] - log(denom[xi[0][i]]));

  if(maplam) L += (nregpar + lampar[0] -1.0)*Bsum/(lampar[1]+Bsum);
  else L += Bsum*lambda;

  return L;
}

/* update eta, Bsum, B, and return the new objective value */

void update(int j, int k, double bnew)
{
  assert(j!=0);
  int i;
  #pragma omp parallel for private(i) 
  for(i=0; i<n; i++)
    { 
      denom[i] += -exp(eta[j][i]);
      eta[j][i] = eta[j][i] + V[k][i]*(bnew-B[k][j]);
      denom[i] += exp(eta[j][i]);
    }

  if(k!=0) Bsum += fabs(bnew) - fabs(B[k][j]);
  B[k][j] = bnew;

  if(maplam)  lambda = (nregpar + lampar[0] -1.0)/(lampar[1] + Bsum);

}

/* Hessian bound updates  */

void calcH(int j, int k){
  int i;
  double Hkj = 0.0;
  assert(D[k][j] >= 0.0);

  #pragma omp parallel for private(i) reduction(+: Hkj)
  for(i=0; i<n; i++)
    {
      double F, E, erm, erp; 
      E = denom[i] - exp(eta[j][i]);      
      erm = exp(eta[j][i]-D[k][j]);
      erp = exp(eta[j][i]+D[k][j]);
      
      if(E < erm) F = erm/E + E/erm;
      else if(E > erp) F = erp/E + E/erp;
      else F = 2.0; 
         
      Hkj += V[k][i]*V[k][i]*((double) m[i])/(F+2.0);   
    }

  H[k][j] = Hkj;
}

/* the moveable part of our bounding function */

double bound(double z, int j, int k, double grad) 
{
  double c;
  if(maplam)
    c = (lampar[0] + nregpar - 1.0)*(Bsum - fabs(B[k][j]) + fabs(z))/(lampar[1] + Bsum - fabs(B[k][j]) + fabs(z));
  else c = lambda*(Bsum - fabs(B[k][j]) + fabs(z));

  return grad*(z-B[k][j]) + 0.5*(z-B[k][j])*(z-B[k][j])*H[k][j] + c;
}

/* the cubic rooter */

double findroot(int j, int k, double grad, double sgn)
{
  assert(maplam);

  double a, b, c, s, r, ghb;
  
  s = (lampar[0] + nregpar - 1.0)*lampar[1];
  r = lampar[1] + Bsum - fabs(B[k][j]);
  ghb = grad/H[k][j]-B[k][j];
    
  a = ghb + sgn*2.0*r;
  b = sgn*ghb*2.0*r + r*r;
  c = ghb*r*r + sgn*s/H[k][j];

  int ns = 0;
  double sln;

  double *roots = solvecubic(a,b,c);
  int nr = roots[0];
  if(nr==2) nr = 1; // ignore complex

  for(int h=1; h<=nr; h++)
    if(sgn == sign(roots[h]))
      if( ( H[k][j] - 2.0*s/cubed(r+fabs(roots[h])) ) > 0.0 )
	{ sln = roots[h]; ns++; } 
     
  assert(ns <= 1);  

  if(ns == 1)
    { if(  bound(0.0, j, k, grad) < bound(sln, j, k, grad) ) sln = 0.0; }
  else sln = 0.0;

  free(roots);
  return sln;
}

/* The gradient descent move for given direction */ 

double Bmove(int j, int k, double grad, double sgn, double pen)
{
  double dbet;

  if(H[k][j] == 0.0) return -B[k][j]; // happens only if you have all zero predictors

  if(pen == 0.0){ // unpenalized (intercept) parameters
    dbet = -grad/H[k][j]; 
  }
  else if(!maplam){ // update with fixed lambda 
    if(sgn != 0.0)
      { dbet = -(grad+sgn*pen)/H[k][j]; 
	if(sgn != sign(B[k][j]+dbet)) dbet = -B[k][j]; }
    else{ // check both sides
      double obj = bound(B[k][j], j, k, grad);
      dbet =  -(grad + pen)/H[k][j];
      if( !(bound(B[k][j] + dbet, j, k, grad) < obj) )
	{ dbet = -(grad - pen)/H[k][j];
	  if( !(bound(B[k][j] + dbet, j, k, grad) < obj) ) dbet = 0.0; } }
  }
  else{ // joint lambda/beta update
    if(sgn != 0.0) dbet = findroot(j, k, grad, sgn) - B[k][j];
    else{  // check both sides
      assert(B[k][j] == 0.0);
      double sln = findroot(j, k, grad, +1.0);
      if(sln == 0.0) sln = findroot(j, k, grad, -1.0);
      dbet = sln - B[k][j]; }
  }      

  if(fabs(dbet) > D[k][j]) dbet = sign(dbet)*D[k][j]; // trust region bounds 
  
  return dbet;
}

/* 
 * cgd_mnlogit:
 *
 * Cyclic coordinate descent for m.a.p. estimation of multinomial
 * logistic regression coefficients under a laplace prior.  
 *
 * Assumes that the first column of V is the (unpenalized) intercept
 *
 */

void cgd_mnlogit(int *n_in, int *p_in, int *d_in, int *m_in, double *tol_in, int *niter,
		 int *N_in, double *X_in, int *xi_in, double *V_in, double *beta_vec, double *Lout, 
		 int *maplam_in, double *lampar_in, double *lamout, double *dmin, double *delta, int *verbalize)
{
  int i, j, k, l, h, t, verb;
  double tol, tmax, grad, diff, penalty;

 
  /** Build everything **/
  verb = *verbalize;
  tol = *tol_in;
  tmax = *niter;

  n = *n_in;
  p = *p_in;
  d = *d_in;  
  nregpar = ((double) (d*p-p));

  m = new_dup_ivec(m_in, n); 

  N = *N_in;
  X = new_dup_dvec(X_in, N);
  xi = new_imat_fromv(N, 2, xi_in);
  V = new_mat_fromv(n, d, V_in);

  H = new_zero_mat(p+1, d);  
  D = new_mat_fromv(p+1, d, delta);

  maplam = *maplam_in;
  lambda = lampar_in[0];
  if(maplam) 
    { lampar = new_dup_dvec(&lampar_in[1], 2);
      lamout[0] = lambda;  }
  
  B = new_mat_fromv(p+1, d, beta_vec);
  Bsum = 0.0;
  for(j=1; j<=p; j++) for(k=1; k<d; k++) Bsum += fabs(B[k][j]);

  eta = new_zero_mat(n, p+1);
  la_dgemm( 0, 1, n, d, p+1, d, n, p+1, V, B, eta, 1.0, 0.0 ); 
  denom = new_dzero(n);
#pragma omp parallel for private(i,j)
  for(i=0; i<n; i++)
    { for(j=0; j<=p; j++) denom[i] += exp(eta[j][i]); }
  if(eta[0][0]!=0.0) myprintf(stdout, "You've input nonzero null category betas; these are not updated.\n");

  G = new_zero_mat(p+1, d);
  for(i=0; i<N; i++) for(k=0; k<d; k++) G[k][xi[1][i]] += -X[i]*V[k][xi[0][i]]; 

  Lout[0] = neglogpost();
  
  /* introductory print statements */
  if(verb)
    { myprintf(stdout, "\n ** Logit-Lasso m.a.p. estimation for a %d x %d beta matrix ** \n\n", p, d);
# ifdef _OPENMP
#pragma omp parallel
      if(omp_get_thread_num() == 0) myprintf(stdout, "OpenMP enabled with %d threads. \n", omp_get_num_threads()); 
# endif
      if(maplam) myprintf(stdout, "Joint estimation for lambda under a gamma(%g,%g) prior\n", lampar[0], lampar[1]);
      myprintf(stdout, "Trust region starts at +/- %g for upper bound on info.\n\n", delta[1]);
      myprintf(stdout, "Objective L initialized at %g\n", Lout[0]); }
  
  /* optimize until objective stops improving or t>tmax */

  diff = tol*100.0;
  t = 0;
  int dozero = 1; // switch for full set updating
  int numzero;
  while(t<tmax && diff > tol){
    numzero = 0;
    for(j=1; j<=p; j++)
      for(k=0; k<d; k++)
	{ 
	  if(k==0) penalty = 0.0;
	  else penalty = lambda;	  

	  grad = G[k][j];
	  #pragma omp parallel for private(i) reduction(+: grad)
	  for(i=0; i<n; i++) grad += ((double) m[i])*exp(eta[j][i] - log(denom[i]))*V[k][i];

	  calcH(j, k); 

	  double bnew = B[k][j];
	  if(B[k][j] != 0.0 || penalty==0) bnew = B[k][j] + Bmove(j, k, grad, sign(B[k][j]), penalty);
	  else if( (runif(0.0,1.0) < 0.25) || dozero) bnew = B[k][j] + Bmove(j, k, grad, 0.0, penalty);	
     
	  if(bnew != B[k][j])
	    { D[k][j] = fmax(*dmin,fmax(0.5*D[k][j], 2.0*fabs(B[k][j] - bnew)));
	      update(j, k,  bnew); assert(D[k][j] > 0.0); }
	  
	  if(B[k][j] == 0.0) numzero++;
	  
	}
    
    t++;
    
    if(maplam) lamout[t] = lambda; 
    
    Lout[t] = neglogpost();
    diff = Lout[t-1] - Lout[t];
    
    if(Lout[t]!=Lout[t]){ diff = 0.0;  myprintf(stdout, "\nL is nan! Probably due to tiny lambda. \n"); }
    if(verb)
      { myprintf(stdout, "t = %d: L = %g (diff of %g) with lambda = %g & %d zero betas.\n", 
		 t, Lout[t], diff, lambda, numzero);
	if(t==tmax) myprintf(stdout, "Terminating optimization; exhausted max number of iterations.\n"); }
   
    if(dozero == 1) dozero = 0;
    else if(diff < tol)
      { diff += tol+1.0; 
	dozero = 1; } 

  }
  
  /* return new beta (and finishing deltas)*/

  for(j=0; j<=p; j++) 
    for(k=0; k<d; k++)
      { beta_vec[k*(p+1) + j] = B[k][j];
	delta[k*(p+1) + j] = D[k][j]; }

  *niter = t+1;

  /* destroy everything */

  free(m);
  free(X);
  delete_imat(xi);
  delete_mat(V);
  delete_mat(eta);
  free(denom);
  delete_mat(B);
  delete_mat(G);
  delete_mat(H);
  delete_mat(D); 
  if(maplam) free(lampar);
}




