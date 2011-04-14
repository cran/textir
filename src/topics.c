#include <stdlib.h>
#include <assert.h>
#include <Rmath.h>
#include "latools.h"
#include "rhelp.h"
# ifdef _OPENMP
#include <omp.h>
# endif


double wllhd(int nwrd, double *x, double *q){
  double l = 0.0;
  int j;
  for(j=0; j<nwrd; j++) l += x[j]*log(q[j]);
  return l;
}

void wrdprob(double *q, int nwrd, int K, int *wrd, double **theta, double *w){ 
  int j, k;
  zero_dvec(q,nwrd);
  for(j=0; j<nwrd; j++)
    for(k=0; k<K; k++)
      q[j] += w[k]*theta[k][wrd[j]];
}

void wgrad(double *grad, int nwrd, int K, int *wrd, double *x, double *q, double **theta, double *w, int nef){
  int j, k;
  zero_dvec(grad,K);
  for(k=0; k<K; k++)
    { for(j=0; j<nwrd; j++) grad[k] += theta[k][wrd[j]]*x[j]/q[j];
      if(nef) grad[k] += 1.0/(w[k]*((double) K)); } 
}


void wneghess_lowtri(double *nH, int nwrd, int K, int *wrd, double *x, double *q, double **theta, double *w, int nef){
  int j, h, k;
  zero_dvec(nH, K*K);
  for(k=0; k<K; k++){
    if(nef) nH[k*(K+1)] += 1.0/(w[k]*w[k]*((double) K));
    for(h=k; h<K; h++)  
	for(j=0; j<nwrd; j++)
	  nH[k*K+h] += x[j]*theta[k][wrd[j]]*theta[h][wrd[j]]/(q[j]*q[j]); }

}

int wtighten(int K, double *w, int *active){
  int h, k, nact;
  // impose practical bounds > 0 for efficiency 
  nact = sum_ivec(active, K);

  for(h=0; h<K; h++)
    if(!active[h]){
      if(w[h] <= .00000000001)
	{ w[h]=0.0; active[h] = 1; nact++; }
      if(w[h] >= 1.0-.00000000001)
	{ zero_dvec(w,K); w[h] = 1.0; 
	  for(k=0; k<K; k++) active[k] = 1;
	  return 0; } }

  return K - nact;
}

int wupdate(int K, double *w, double *B, int *active){
    double delmin = 1.0;
    double delta;
    int h, d;

    d = 0;
    for(h=0; h<K; h++)
      if(!active[h])
	{ assert(w[h]>0);
	  delta = 1.0;
	  if(B[d] < -w[h]) delta = -w[h]/B[d];
	  if(B[d] > 1.0-w[h]) delta = (1.0-w[h])/B[d];
	  if(delta < delmin) delmin = delta;
	  d++; }
    
    d = 0;
    for(h=0; h<K; h++)
      if(!active[h]) 
	{ w[h] += delmin*B[d]; d++; }

    return  wtighten(K, w, active);
}

/* main program: sequential quadratic programming for topic weights */
int sqpw(int p, int nwrd, int K,  int *wrd, double *x, 
	 double **theta, double *w, int nef,
	 double tol, int tmax, int verb){

  int h,k,c,d;

  // small-sample solution
  if(nwrd==1){
    zero_dvec(w,K);
    k = 0;
    for(h=1;h<K;h++) 
      if(theta[h][*wrd] > theta[k][*wrd]) k = h;
    w[k] = 1.0;
    return 1; }

  double diff = tol+1.0;

  double *q = new_dvec(nwrd);
  wrdprob(q, nwrd, K, wrd, theta, w);
  double L = wllhd(nwrd, x, q);
  double Lnew;

  double *A = new_dzero( (K+1)*(K+1) );
  double *B = new_dzero(K+1);
  double *nhess = new_dvec(K*K);
  double *grad = new_dvec(K);
  double *wold = new_dup_dvec(w, K);

  int *active = new_izero(K);
  int nw = K;
  int info;
  int dosysv = 1;
  int iter = 1;
   
  while(diff > tol && iter < tmax && nw > 0){

    /* update gradient, and negative hessian */
    wgrad(grad, nwrd, K, wrd, x, q, theta, w, nef);
    wneghess_lowtri(nhess, nwrd, K, wrd, x, q, theta, w, nef);
    // NB: assumed prior concentration of 1/K if in NEF parametrization

    /* build  the linear system */
    d = 0;
    for(k=0; k<K; k++)
      if(!active[k]){ 
	B[d] = grad[k];
	c = d;
	for(h=k; h<K; h++) 
	  if(!active[h]){ A[d*(nw+1) + c] = nhess[k*K + h]; c++; }
	d++;
      }
    assert(d==nw);

    for(c=0; c<nw; c++) A[c*(nw+1) + nw] = 1.0;
    A[(nw+1)*(nw+1)-1] = 0.0;
    B[nw] = 0.0;

	  
    /* solve it  (sysv can be unstable for small samples; switch to gesv if necessary) */
    if(dosysv){ info = la_dsysv(nw+1, 1, A, B); 
      if(info!=0){ dosysv = 0; continue; } }
    else{
      for(d=1; d<=nw; d++) for(c=0; c<d; c++) A[d*(nw+1) + c] = A[c*(nw+1) +d];
      info = la_dgesv(nw+1, 1, A, B);
    }

    /* check and update */
    if(fabs(sum_dvec(B, nw)) > .001 || (info !=0 && dosysv == 0)) nw = wtighten(K, w, active);
    else nw = wupdate(K, w, B, active); 


    /* iterate */
    wrdprob(q, nwrd, K, wrd, theta, w);
    Lnew =  wllhd(nwrd, x, q);
    diff = 0.0;
    for(h=0; h<K; h++) diff += fabs(wold[h]-w[h]);
    copy_dvec(wold, w, K);

    if(verb==1){ 
      myprintf(stdout, 
	       "Omega Fit: iter = %d, unique words = %d, active=%d,  log llhd = %g, diff = %g \n",
	       iter, nwrd, K-nw, Lnew, diff);
      myprintf(stdout, "          W = "); print_dvec(w, K, stdout); }
    
    iter++;
    L = Lnew;
  }

  free(wold);
  free(active);
  free(q);
  free(B);
  free(grad);
  free(A);
  free(nhess);

  return 1;
}

void Romega( int *n, int *p, int *K, 
	     int *doc, int *wrd, double *X, 
	     double *theta_vec, double *W, int *nef,
	     double *tol, int *tmax, int *verb){
 
  int speakup = *verb;

# ifdef _OPENMP
  speakup=0; // otherwise it'll crash R trying to print all at once
# endif

 double **theta = new_mat_fromv(*p, *K, theta_vec);
 int i;

#pragma omp parallel for private(i) 
 for(i=0; i<(*n); i++) 
   if(sqpw(*p, doc[i+1]-doc[i], *K, &wrd[doc[i]], &X[doc[i]], theta,  &W[i*(*K)], *nef, *tol, *tmax, speakup) == 0)
     myprintf(stdout, "Failed to converge for omega at i = %d\n", i+1); 

 delete_mat(theta);
}


void RcalcQ(int *n_in, int *p_in, int *K_in, int *doc, int *wrd, 
	    int *N, double *omega, double *theta, double *q ){

  int n, p, K, l, h;

  n = *n_in;
  p = *p_in;
  K = *K_in;

  for(l=0; l<(*N); l++){
    q[l] = 0.0;
    for(h=0; h<K; h++) q[l] += omega[h*n + doc[l]]*theta[h*p + wrd[l]];
  }
}

void RcalcHW(int *n_in, int *p_in, int *K_in, double *omeg, double *thet,
	     int *doc, int *wrd, double *cnt, double *q, int *N, double *H){

  int n, p, K, K2, l, h, k;

  n = *n_in;
  p = *p_in;
  K = *K_in;
  K2 = K*K;
  double ttbyq2;

  zero_dvec(H, n*K2);

  for(l=0; l<(*N); l++)
    for(k=0; k<K; k++)
      {
	for(h=k; h<K; h++)
	  { ttbyq2 = exp( log(thet[k*p+wrd[l]])+log(thet[h*p+wrd[l]]) - 2*log(q[l]) );
	    H[doc[l]*K2 + K*k +h] += -cnt[l]*omeg[k*n+doc[l]]*omeg[h*n+doc[l]]*(1.0 - ttbyq2); 
	    H[doc[l]*K2 + K*k +h] += ((double) (K+1))*omeg[k*n+doc[l]]*omeg[h*n+doc[l]]; }
	
	H[doc[l]*K2 + (K+1)*k] += -cnt[l]*omeg[k*n+doc[l]]*( thet[k*p+wrd[l]]-q[l] )/q[l]
	  + ((double) (K+1))*omeg[k*n+doc[l]]; 	
      }

  for(l=0; l<n; l++)
    for(k=0; k<K; k++)
      for(h=0; h<k; h++) H[l*K2 + K*k +h] = H[l*K2 + K*h +k]; 
}


void logit(int d, double *eta, double *prob){
  int j;
  double lp0 = log(prob[0]);
  for(j=1; j<d; j++) eta[j-1] = log(prob[j]) - lp0; }

void expit(int d, double *prob, double *eta){
  int j;
  prob[0] = 1; 
  for(j=1; j<d; j++) prob[j] = exp(eta[j-1]);
  double p0 = sum_dvec(prob, d);
  for(j=0; j<d; j++) prob[j] *= 1.0/p0;  }


void RtoNEF(int *n_in, int *p_in, int *K_in, double *Y, double *theta, double *tomega){
  int l;
  int n = *n_in;
  int p = *p_in;
  int K = *K_in;

  for(l=0; l<K; l++) logit(p, &Y[l*(p-1)], &theta[l*p]);
  for(l=0; l<n; l++) logit(K, &Y[K*(p-1) + l*(K-1)], &tomega[l*K]); }


void RfromNEF(int *n_in, int *p_in, int *K_in, double *Y, double *theta, double *tomega){
  int l;
  int n = *n_in;
  int p = *p_in;
  int K = *K_in;

  for(l=0; l<K; l++) expit(p, &theta[l*p], &Y[l*(p-1)]);
  for(l=0; l<n; l++) expit(K, &tomega[l*K], &Y[K*(p-1) + l*(K-1)]); }


void Rzhat(int *n_in, int *p_in, int *K_in, int *N_in, double *xhat, 
	   int *doc, int *wrd, double *zhatj, double *zhati){
  int h,l;
  int n = *n_in;
  int p = *p_in;
  int K = *K_in;
  int N = *N_in;

  // zhat's must be pre-zero'd
  for(l=0; l<N; l++)
    for(h=0; h<K; h++)
      { zhatj[h*p + wrd[l]] += xhat[N*h +l];
	zhati[h*n + doc[l]] += xhat[N*h +l]; }
}