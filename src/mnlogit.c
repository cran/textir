// need to make sure both lams are k indexed (not k-1)
#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
#include "latools.h"
#include "rhelp.h"
#include "polysolve.h"
#include <time.h>
# ifdef _OPENMP
#include <omp.h>
# endif

/* global variables */

int n, p, d, Nx, Nv, RE;

int dirty = 0;
double *m = NULL;
double *X = NULL;
int **xind = NULL;
int *xi = NULL;
double *V = NULL;
int **vind = NULL;
int *vk = NULL;
double **eta = NULL;
double *denom = NULL;
double **B = NULL;
double **G = NULL;
double **H = NULL;
double **D = NULL;
double **lam = NULL;
int *maplam = NULL;
double *revar = NULL;
double *remean = NULL;
double **U = NULL;
double *nvec;

/* un-normalized negative log posterior */

double neglogpost(){  
  
  double L = 0.0;
  int i, j, k;
  
  for(i=0; i<Nx; i++) L += -X[i]*(eta[xind[1][i]][xind[0][i]] - log(denom[xind[0][i]]));

  for(k=0; k<d; k++){ 
    if(maplam[k]==1) 
      for(j=0;j<=p;j++) L += lam[k][0]*log( 1.0 + fabs(B[k][j])/lam[k][1] );
    else if(maplam[k]==0) 
      for(j=0;j<=p;j++) L += lam[k][0]*fabs(B[k][j]);
  }
  if(RE) for(i=0;i<n;i++) for(j=1;j<=p;j++) L += (U[j][i]-remean[i])*(U[j][i]-remean[i])*0.5/revar[i];

  return L;

}

/* update eta + B, and return the new objective value */

void update(int j, int k, double bnew)
{
  assert(j!=0);
  int i;
  copy_dvec(nvec, eta[j], n);
  for(i=vk[k]; i<vk[k+1]; i++) nvec[vind[0][i]] += V[i]*(bnew-B[k][j]);
  for(i=0; i<n; i++)
    { denom[i] += exp(nvec[i]) - exp(eta[j][i]);
      eta[j][i] = nvec[i]; }
  B[k][j] = bnew;
}

/* lub factor */
double calcF(int i, int j, double delta){
  double E = denom[i] - exp(eta[j][i]);      
  double erm = exp(eta[j][i]-delta);
  double erp = exp(eta[j][i]+delta);
      
  if(E < erm) return erm/E + E/erm;
  else if(E > erp) return erp/E + E/erp;
  else return 2.0;
} 

/* draw random effects and update eta+denom */

void update_rei(int i){
  assert(RE);
  int j, l;
  double g, h, dbet, delta;
  delta = 0.1; // fixed trust to avoid storage

  for(j=1; j<=p; j++){ 
    g = 0.0;
    for(l=xi[i]; l<xi[i+1]; l++) if(xind[1][l]==j){ g += -X[l]; break; }
    g += m[i]*exp(eta[j][i] - log(denom[i])) + (U[j][i] - remean[i])/revar[i];
    if(g == 0.0) continue;
    h = m[i]/(calcF(i, j, delta)+2.0) + 1.0/revar[i];
    dbet = -g/h;
    if(fabs(dbet) < 0.0001) continue; // numerical overload otherwise
    if(fabs(dbet) > delta) dbet = sign(dbet)*delta; 
    U[j][i] += dbet;
    eta[j][i] += dbet;
    denom[i] += exp(eta[j][i]) - exp(eta[j][i]-dbet);
  }
}


/* Hessian bound updates  */

void calcH(int j, int k){
  int i;
  double Hkj = 0.0;
  assert(D[k][j] >= 0.0);

  for(i=0; i<n; i++) nvec[i] = calcF(i, j, D[k][j]);
  for(i=vk[k]; i<vk[k+1]; i++) Hkj += V[i]*V[i]*m[vind[0][i]]/(nvec[vind[0][i]]+2.0);   
  H[k][j] = Hkj;
}

/* the moveable part of our bounding function */

double bound(double z, int j, int k, double grad) 
{
  double c = 0.0;
  if(maplam[k]==1)
    c = lam[k][0]*log( 1.0 + fabs(z)/lam[k][1] );
  else if(maplam[k]==0) 
    c = lam[k][0]*fabs(z);

  return grad*(z-B[k][j]) + 0.5*(z-B[k][j])*(z-B[k][j])*H[k][j] + c;
}

/* the quadratic rooter */

double findroot(int j, int k, double grad, double sgn)
{
  if(maplam[k]<0) error("mnlogit -> findroot is trying to solve for fixed coeficients.  Please report this bug.");

  double b, c, ghb;
  ghb = grad/H[k][j]-B[k][j];
  b = ghb + sgn*lam[k][1];
  c = ghb*sgn*lam[k][1] + lam[k][0]/H[k][j];

  double *roots = solvequadratic(b,c);
  int ns = 0;

  double sln = 0.0;
  for(int h=1; h<=((int) roots[0]); h++) // check 1 or 2 roots
    if(sgn == sign(roots[h])) // same sign as current phi
      if( ( H[k][j] - lam[k][0]/((lam[k][1]+fabs(roots[h]))*(lam[k][1]+fabs(roots[h]))) ) > 0.0 ) // local min
	if(  bound(0.0, j, k, grad) >= bound(roots[h], j, k, grad) ) // better than zero
	  { sln = roots[h]; ns++; } 
  
  if(ns > 1) warning("multiple candidate solutions in findroot.");  
  free(roots);
  return sln;
}

/* The gradient descent move for given direction */ 

double Bmove(int j, int k, double grad, double sgn)
{
  double dbet;
  if(maplam[k]<0) error("mnlogit is trying to solve for fixed coeficients.  Please report this bug.");
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
      if(B[k][j] != 0.0) error("accounting error in mnlogit -> Bmove.  Please report this bug.");
      double sln = findroot(j, k, grad, +1.0);
      if(sln == 0.0) sln = findroot(j, k, grad, -1.0);
      dbet = sln - B[k][j]; }
  }      

  if(fabs(dbet) > D[k][j]) dbet = sign(dbet)*D[k][j]; // trust region bounds 
  return dbet;
}


/* cleanup function */
void mnlm_cleanup(){
  if(!dirty) return;

  if(m){ free(m); m = NULL; }
  if(nvec){ free(nvec); nvec = NULL; }
  if(X){ free(X); X = NULL; }
  if(xind){ delete_imat(xind); xind = NULL; }
  if(xi){ free(xi); xi = NULL; }
  if(V){ free(V); V = NULL; }
  if(vind){ delete_imat(vind); vind = NULL; }
  if(vk){ free(vk); vk = NULL; }
  if(eta){ delete_mat(eta); eta = NULL; }
  if(denom){ free(denom); denom = NULL; }
  if(B){ delete_mat(B); B = NULL; }
  if(G){ delete_mat(G); G = NULL; }
  if(H){ delete_mat(H); H = NULL; }
  if(D){ delete_mat(D); D = NULL; }
  if(maplam){ free(maplam); maplam = NULL; }
  if(lam){ delete_mat(lam); lam = NULL; }
  if(remean){ free(remean); remean = NULL; }
  if(revar){ free(revar); revar = NULL; }
  if(U){ delete_mat(U); U = NULL; }
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

void Rmnlogit(int *n_in, int *p_in, int *d_in, double *m_in, double *tol_in, 
	      int *Nx_in, double *X_in, int *xind_in, int *xi_in,
	      int *Nv_in, double *V_in, int *vind_in, int *vk_in,
	      double *beta_vec, double *fitted,
	      int *maplam_in, double *lam_in, 
	      double *dmin, double *dinit, 
	      double *Gout, int *RE_in, double *randeff,
	      int *verbalize)
{
  dirty = 1; // flag to say the function has been called
  time_t itime = time(NULL);  // time stamp for periodic R interaction 

  int i, j, k, t, verb;
  double tol, grad, diff;
  double Lold, Lnew;

  /** Build everything **/
  verb = *verbalize;
  tol = *tol_in;
  n = *n_in;
  p = *p_in;
  d = *d_in;  

  nvec = new_dzero(n);
  m = new_dup_dvec(m_in, n); 

  Nx = *Nx_in;
  X = new_dup_dvec(X_in, Nx);
  xind = new_imat_fromv(Nx, 2, xind_in);

  Nv = *Nv_in;
  V = new_dup_dvec(V_in, Nv);
  vind = new_imat_fromv(Nv, 2, vind_in);
  vk = new_dup_ivec(vk_in, d+1); 

  H = new_zero_mat(p+1, d);  
  D = new_mat(p+1, d);
  for(j=0; j<=p; j++) for(k=0; k<d; k++) D[k][j] = dinit[0]; 

  maplam = new_dup_ivec(maplam_in, d);
  lam = new_mat_fromv(2, d, lam_in);
  
  B = new_mat_fromv(p+1, d, beta_vec);

  RE = *RE_in;
  if(RE){
    xi = new_dup_ivec(xi_in, n+1); 
    remean = new_dup_dvec(randeff, n);
    revar = new_dup_dvec(&randeff[n], n);
    U = new_mat(n,p+1);
    zero_dvec(U[0], n);
    for(j=1; j<=p; j++) copy_dvec(U[j], remean, n);
    eta = new_dup_mat(n, p+1, U);
  } else{ eta = new_zero_mat(n, p+1); }
  
  for(j=0; j<=p; j++) for(i=0; i<Nv; i++) eta[j][vind[0][i]] += V[i]*B[vind[1][i]][j];
  denom = new_dzero(n);
  for(i=0; i<n; i++) for(j=0; j<=p; j++) denom[i] += exp(eta[j][i]); // includes 1.0 for null category
  if(eta[0][0]!=0.0) myprintf(mystdout, "You've input nonzero null category betas; these are not updated.\n");

  G = new_mat_fromv(p+1, d, Gout);

  Lnew = neglogpost();
  if(isinf(Lnew) || isnan(Lnew)){  
    warning(" Infinite or NaN initial fit; starting at zero instead.  \n  Try `normalize=TRUE' or a larger penalty shape or rate.\n");
    for(j=0; j<=p; j++) for(k=0; k<d; k++) B[k][j] = 0.0;
    if(RE){
      copy_mat(n, p+1, eta, U);
      zero_dvec(denom, n);
      for(i=0; i<n; i++) for(j=0; j<=p; j++) denom[i] += exp(eta[j][i]); }
    else{
      zero_mat(eta, n, p+1);
      for(i=0; i<n; i++) denom[i] = ((double) p) + 1.0; }
    Lnew  = neglogpost(); 
    if(isinf(Lnew) || isnan(Lnew)){
      warning(" Probabilities of exactly zero in your likelihood.  \n   You need to use a larger penalty.\n");
      Lnew = 100000000.0; 
    }
  }
    

  diff = tol*100.0;
  t = 0;
  int dozero = 1; 
  double numzero, nregpar;
  double bnew;
  double rateincrease = 0.0; 


  /* introductory print statements */
  if(verb)
    { myprintf(mystdout, "*** Logistic Regression with a %d x %d coefficient matrix ***\n", p, d);
      myprintf(mystdout, "Objective L initialized at %g\n", Lnew); }
  
  /* optimize until objective stops improving */
  while(diff > tol | diff < 0){
    numzero = 0.0;
    nregpar = 0.0;

    // loop through coefficient
    for(j=1; j<=p; j++)
      for(k=0; k<d; k++)
	if(maplam[k]>=0){ 
	  if(B[k][j] != 0.0 || dozero || lam[k][0]==0 || t < 3 || runif(0,1) < 0.1){
	    // gradient
	    grad = G[k][j];
	    for(i=vk[k]; i<vk[k+1]; i++) grad += m[vind[0][i]]*exp(eta[j][vind[0][i]] - log(denom[vind[0][i]]))*V[i]; 
	    // curvature 
	    calcH(j, k); 
	    // conditional newton update
	    bnew = B[k][j] + Bmove(j, k, grad, sign(B[k][j]));
	    // check and update dependencies
	    if(bnew != B[k][j])
	      { D[k][j] = fmax(*dmin,fmax(0.5*D[k][j], 2.0*fabs(B[k][j] - bnew)));
		update(j, k,  bnew); 
		if(D[k][j] <= 0.0) error("negative bound window in mnlogit. Please report this bug."); }
	  }
	  // sum the zeros and check for escape from R
	  if(k>0)
	    { nregpar++;
	      if(B[k][j] == 0.0) numzero++; }
	  itime = my_r_process_events(itime); }

    if(RE){
#pragma omp parallel for private(i) 
      for(i=0; i<n; i++) update_rei(i); 
    }

    // iterate
    t++;
    Lold = Lnew;
    Lnew = neglogpost();
    diff = Lold - Lnew;
    
    // print 
    if(Lnew!=Lnew || !isfinite(Lnew) || Lnew < 0){ 
      warning("The algorithm did not converge (non-finite likelihood).  \n   Try `normalize=TRUE' or a larger penalty shape or rate.\n"); 
      dozero = 1;
      diff = 0.0;
    }
    else if(verb)
      { myprintf(mystdout, "t = %d: L = %g (diff of %g) with %g%% zero loadings.\n", 
		 t, Lnew, diff, 100.0*(numzero/nregpar)); }
   
    if(diff < 0.0){
      i = 0;
      for(k=0; k<d; k++)
	if(maplam[k]==1)
	  { lam[k][0] *= 2.0;
	    lam[k][1] *= 2.0; 
	    i++; }
      if(i>0 && rateincrease < 8){
	if(verb) myprintf(mystdout, "WARNING: non-monotonic convergence, probably due to a non-concave posterior.  \n");
	rateincrease += 1.0;
	Lnew = neglogpost();}
      else{
	 warning("The algorithm did not converge (non-decreasing likelihood).  \n  Try `normalize=TRUE' or a larger penalty shape or rate.\n");
	 dozero = 1;
	 diff = 0.0;
      }
    }
    lam_in[0] = rateincrease;

    // check for active set update
    if(dozero == 1) dozero = 0;
    else if(diff < tol)
      { diff += tol+1.0; 
	dozero = 1; } 
  }
  
  /* clean, collect, and exit */
  for(j=0; j<=p; j++) for(k=0; k<d; k++) beta_vec[k*(p+1) + j] = B[k][j];  // beta fit
  for(i=0; i<Nx; i++) fitted[i] = exp(eta[xind[1][i]][xind[0][i]] - log(denom[xind[0][i]] - 1.0*((double) p > 1))); // exclude null for > 2 cat
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


