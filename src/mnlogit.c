// need to make sure both lams are k indexed (not k-1)
#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
#include "latools.h"
#include "rhelp.h"
#include "cubic.h"
#include <time.h>

/* global variables */

int n, p, d, N;

int dirty = 0;
int *m = NULL;
double *X = NULL;
int **xi = NULL;
double **V = NULL;
double *bmax = NULL;
double **eta = NULL;
double *denom = NULL;
double **B = NULL;
double **G = NULL;
double **H = NULL;
double **D = NULL;
double **lam = NULL;
int *maplam = NULL;

/* un-normalized negative log posterior */

double neglogpost(){  
  
  double L = 0.0;
  int i, j, k;
  
  for(i=0; i<N; i++) L += -X[i]*(eta[xi[1][i]][xi[0][i]] - log(denom[xi[0][i]]));

  for(k=1; k<d; k++){ 
    if(maplam[k]==1) 
      for(j=0;j<=p;j++) L += lam[k][0]*fabs(B[k][j])/(lam[k][1]+fabs(B[k][j]));
    else if(maplam[k]==0) 
      for(j=0;j<=p;j++) L += lam[k][0]*fabs(B[k][j]);
  }

  return L;
}


/* calculates and returns the pearson residuals */
void residuals(double *xhat, double *r){
  int i;
  // x.hat fitted counts
  for(i=0; i<N; i++) xhat[i] = m[xi[0][i]]*exp(eta[xi[1][i]][xi[0][i]] - log(denom[xi[0][i]]));
  // residuals
  for(i=0; i<N; i++) r[i] = (X[i]-xhat[i])/sqrt(xhat[i]*(m[xi[0][i]]-xhat[i])/m[xi[0][i]]);
}


/* update eta + B, and return the new objective value */

void update(int j, int k, double bnew)
{
  assert(j!=0);
  int i;
  for(i=0; i<n; i++)
    { 
      denom[i] += -exp(eta[j][i]);
      eta[j][i] = eta[j][i] + V[k][i]*(bnew-B[k][j]);
      denom[i] += exp(eta[j][i]);
    }

  B[k][j] = bnew;
}

/* Hessian bound updates  */

void calcH(int j, int k){
  int i;
  double Hkj = 0.0;
  assert(D[k][j] >= 0.0);

  // #pragma omp parallel for private(i) reduction(+: Hkj)
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
  double c = 0.0;
  if(maplam[k]==1)
    c = lam[k][0]*fabs(z)/(lam[k][1] + fabs(z));
  else if(maplam[k]==0) 
    c = lam[k][0]*fabs(z);

  return grad*(z-B[k][j]) + 0.5*(z-B[k][j])*(z-B[k][j])*H[k][j] + c;
}

/* the cubic rooter */

double findroot(int j, int k, double grad, double sgn)
{
  assert(maplam[k]>=0);

  double a, b, c, s, r, ghb;
  
  s = lam[k][0]*lam[k][1];
  r = lam[k][1];

  ghb = grad/H[k][j]-B[k][j];
    
  a = ghb + sgn*2.0*r;
  b = sgn*ghb*2.0*r + r*r;
  c = ghb*r*r + sgn*s/H[k][j];

  int ns = 0;
  double sln = 0.0;

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

double Bmove(int j, int k, double grad, double sgn)
{
  double dbet;
  assert(maplam[k] >= 0);
  if(H[k][j] == 0.0) return -B[k][j]; // happens only if you have all zero predictors

  if(lam[k][0]==0){ // unpenalized (e.g. intercept) parameters
    dbet = -grad/H[k][j]; 
  }
  else if(maplam[k]==0){ // update with fixed lambda 
    if(sgn != 0.0)
      { dbet = -(grad+sgn*lam[k][0])/H[k][j]; 
	if(sgn != sign(B[k][j]+dbet)) dbet = -B[k][j]; }
    else{ // check both sides
      double obj = bound(B[k][j], j, k, grad);
      dbet =  -(grad + lam[k][0])/H[k][j];
      if( !(bound(B[k][j] + dbet, j, k, grad) < obj) )
	{ dbet = -(grad - lam[k][0])/H[k][j];
	  if( !(bound(B[k][j] + dbet, j, k, grad) < obj) ) dbet = 0.0; } }
  }
  else{ // lambda/beta update
    if(sgn != 0.0) dbet = findroot(j, k, grad, sgn) - B[k][j];
    else{  // check both sides
      assert(B[k][j] == 0.0);
      double sln = findroot(j, k, grad, +1.0);
      if(sln == 0.0) sln = findroot(j, k, grad, -1.0);
      dbet = sln - B[k][j]; }
  }      

  if(fabs(dbet) > D[k][j]) dbet = sign(dbet)*D[k][j]; // trust region bounds 
  if(fabs(B[k][j]+dbet) > bmax[k]) dbet = sign(B[k][j]+dbet)*bmax[k] - B[k][j]; // numerical overload bounds
  return dbet;
}


/* cleanup function */
void mnlm_cleanup(){
  if(!dirty) return;

  if(m){ free(m); m = NULL; }
  if(X){ free(X); X = NULL; }
  if(xi){ delete_imat(xi); xi = NULL; }
  if(V){ delete_mat(V); V = NULL; }
  if(bmax){ free(bmax); bmax = NULL; }
  if(eta){ delete_mat(eta); eta = NULL; }
  if(denom){ free(denom); denom = NULL; }
  if(B){ delete_mat(B); B = NULL; }
  if(G){ delete_mat(G); G = NULL; }
  if(H){ delete_mat(H); H = NULL; }
  if(D){ delete_mat(D); D = NULL; }
  if(maplam){ free(maplam); maplam = NULL; }
  if(lam){ delete_mat(lam); lam = NULL; }
}

/* 
 * Main Function: Rmnlogit
 *
 * Cyclic coordinate descent for m.a.p. estimation of multinomial
 * logistic regression coefficients under a laplace prior.  
 *
 * Assumes that the first column of V is the (unpenalized) intercept
 *
 */

void Rmnlogit(int *n_in, int *p_in, int *d_in, int *m_in, double *tol_in, int *niter,
	      int *N_in, double *X_in, int *xi_in, double *V_in, 
	      double *beta_vec, double *bmax_in, double *Lout, double *resids, double *fitted,
	      int *maplam_in, double *lam_in, double *dmin, double *delta, int *verbalize)
{
  dirty = 1; // flag to say the function has been called
  time_t itime = time(NULL);  // time stamp for periodic R interaction 

  int i, j, k, t, verb;
  double tol, tmax, grad, diff;

  /** Build everything **/
  verb = *verbalize;
  tol = *tol_in;
  tmax = *niter;
  n = *n_in;
  p = *p_in;
  d = *d_in;  

  m = new_dup_ivec(m_in, n); 
  N = *N_in;
  X = new_dup_dvec(X_in, N);
  xi = new_imat_fromv(N, 2, xi_in);
  V = new_mat_fromv(n, d, V_in);
  bmax = new_dup_dvec(bmax_in, d);
  H = new_zero_mat(p+1, d);  
  D = new_mat(p+1, d);
  for(j=0; j<=p; j++) for(k=0; k<d; k++) D[k][j] = delta[0]; 

  maplam = new_dup_ivec(maplam_in, d);
  lam = new_mat_fromv(2, d, lam_in);
  
  B = new_mat_fromv(p+1, d, beta_vec);

  eta = new_zero_mat(n, p+1);
  la_dgemm( 0, 1, n, d, p+1, d, n, p+1, *V, *B, *eta, 1.0, 0.0 ); 
  denom = new_dzero(n);
  for(i=0; i<n; i++) for(j=0; j<=p; j++) denom[i] += exp(eta[j][i]); 
  if(eta[0][0]!=0.0) myprintf(stdout, "You've input nonzero null category betas; these are not updated.\n");

  G = new_zero_mat(p+1, d);
  for(i=0; i<N; i++) for(k=0; k<d; k++) G[k][xi[1][i]] += -X[i]*V[k][xi[0][i]]; 

  Lout[0] = neglogpost();
  if(isinf(Lout[0])){  
    myprintf(stdout, "\nInfinite initial fit; starting at zero instead.  Perhaps try `normalize=TRUE'.\n");
    for(j=0; j<=p; j++) for(k=0; k<d; k++) B[k][j] = 0.0;
    for(i=0; i<n; i++)
      { for(j=0; j<=p; j++) eta[j][i] = 0.0;
	denom[i] = ((double) p) + 1.0; }
    Lout[0]  = neglogpost(); }
    
  diff = tol*100.0;
  t = 0;
  int dozero = 1; 
  double numzero, nregpar;
  double bnew;

  /* introductory print statements */
  if(verb)
    { myprintf(stdout, "*** Logistic Regression with a %d x %d coefficient matrix ***\n", p, d);
      for(k=0; k<d; k++){
	if(maplam[k] == 1) myprintf(stdout, "V[%d]: Joint MAP penalty estimation under a gamma(%g,%g) prior\n", k, lam[k][0], lam[k][1]);
	else if (maplam[k] == 0) myprintf(stdout, "V[%d]: Estimation with L1 penalty of %g\n", k, lam[k][0]);
	else if (maplam[k] == -1) myprintf(stdout, "V[%d]: Fixed coefficients\n", k);
      }
      myprintf(stdout, "Objective L initialized at %g\n", Lout[0]); }
  
  /* optimize until objective stops improving or t>tmax */
  while(t<tmax && diff > tol){
    numzero = 0.0;
    nregpar = 0.0;

    // loop through coefficient
    for(j=1; j<=p; j++)
      for(k=0; k<d; k++)
	if(maplam[k]>=0){ 
	  if(B[k][j] != 0.0 || dozero || lam[k][0]==0 || t < 3){
	    // gradient
	    grad = G[k][j];
	    for(i=0; i<n; i++) grad += ((double) m[i])*exp(eta[j][i] - log(denom[i]))*V[k][i];
	    // curvature 
	    calcH(j, k); 
	    // conditional newton update
	    bnew = B[k][j] + Bmove(j, k, grad, sign(B[k][j]));
	    // check and update dependencies
	    if(bnew != B[k][j])
	      { D[k][j] = fmax(*dmin,fmax(0.5*D[k][j], 2.0*fabs(B[k][j] - bnew)));
		update(j, k,  bnew); assert(D[k][j] > 0.0); }
	  }
	  // sum the zeros and check for escape from R
	  if(k>0)
	    { nregpar++;
	      if(B[k][j] == 0.0) numzero++; }
	  itime = my_r_process_events(itime); }
    
    // iterate
    t++;
    Lout[t] = neglogpost();
    diff = Lout[t-1] - Lout[t];
    
    // print 
    if(Lout[t]!=Lout[t] || !isfinite(Lout[t]) || Lout[t] < 0){ 
      diff = 0.0;  dozero=1; 
      myprintf(stdout, "L is NaN!  Try a larger `penalty' or use normalize=TRUE. \n"); }
    else if(verb)
      { myprintf(stdout, "t = %d: L = %g (diff of %g) with %g%% zero loadings.\n", 
		 t, Lout[t], diff, 100.0*(numzero/nregpar));
	if(t==tmax) myprintf(stdout, "Terminating optimization; exhausted max number of iterations.\n"); }
   
    if(diff < 0.0){
      myprintf(stdout, "WARNING: non-monotonic convergence. \nThis means you have a non-concave posterior, leading to an unstable MAP. \nTry increasing the hyperprior rate on your penalty or using fixed penalties.\n"); }

    // check for active set update
    if(dozero == 1) dozero = 0;
    else if(diff < tol)
      { diff += tol+1.0; 
	dozero = 1; } 
  }
  
  /* clean, collect, and exit */
  *niter = t+1; // total iterations
  for(j=0; j<=p; j++) for(k=0; k<d; k++) beta_vec[k*(p+1) + j] = B[k][j];  // beta fit
  residuals(fitted, resids);
  mnlm_cleanup(); // clean 
  dirty = 0;  // declare normal exit
}

/* Fast binning of variables */
void Rbin(int *nobs, int *dim, int *nbin, double *B, double *V, int *O){

  int n = *nobs;
  int d = *dim;
  int b = *nbin;

  double **vals = new_mat_fromv(n, d, V);
  double **bins = new_mat_fromv(b, d, B);
  int **ords = new_imat_fromv(n, d, O);

  int i,j,k;

  for(j=0; j<d; j++){
    k = 0;
    for(i=0; i<n; i++){
      while( vals[j][ords[j][i]] > bins[j][k] ) k++;
      O[j*n + ords[j][i]] = k+1;
    }
  }

  delete_imat(ords);
  delete_mat(bins);
  delete_mat(vals);
}


