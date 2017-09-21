#include <math.h>
#include <string.h>
#include "Rinternals.h"
#include "R_ext/Rdynload.h"
#include <R.h>
#include <R_ext/Applic.h>

double crossprod(double *X, double *y, int n, int j);
double sum(double *x, int n);
int checkConvergence(double *beta, double *beta_old, double eps, int l, int J);
double SSL(double z, double beta, double lambda0, double lambda1, double theta, double v, int n, double delta, double sigma2);
double pstar(double x, double theta, double lambda1, double lambda0);
double lambdastar(double x, double theta, double lambda1, double lambda0);
double expectation_approx(double *beta, double a, double b, int p, int l);
double update_sigma(double *r, int n, double nu, double xi, int l);
double g(double x, double theta, double sigma2, double lambda1, double lambda0, double n);
double threshold(double theta, double sigma2, double lambda1, double lambda0, int n);

// Memory handling, output formatting (Gaussian)

SEXP cleanupG(double *a, double *r, int *e1, int *e2, double *z, double *thetas, double *sigmas, SEXP beta, SEXP loss, SEXP iter, SEXP thetas_export, SEXP sigmas_export) {

  Free(a);
  Free(r);
  Free(e1);
  Free(e2);
  Free(z);
  
  SEXP res;
  PROTECT(res = allocVector(VECSXP, 5));

  SET_VECTOR_ELT(res, 0, beta);
  SET_VECTOR_ELT(res, 1, loss);
  SET_VECTOR_ELT(res, 2, iter);
  SET_VECTOR_ELT(res, 3, thetas_export);
  SET_VECTOR_ELT(res, 4, sigmas_export);

  UNPROTECT(6);

  return(res);
}

// Coordinate descent for gaussian models

SEXP SSL_gaussian(SEXP X_, SEXP y_, SEXP penalty_, SEXP variance_, SEXP lambda0s_, SEXP eps_, SEXP max_iter_, 
                  SEXP lambda1_, SEXP theta_, SEXP sigma_, SEXP counter_, SEXP a_,
		              SEXP b_, SEXP nu_, SEXP xi_) {

  // Declarations
  
  double *X = REAL(X_);

  double *y = REAL(y_);

  int n = length(y_);

  int p = length(X_)/n;

  int L = length(lambda0s_);

  SEXP res, beta, loss, iter, thetas_export, sigmas_export;

  PROTECT(beta = allocVector(REALSXP, L*p));

  PROTECT(thetas_export = allocVector(REALSXP, L));

  PROTECT(sigmas_export = allocVector(REALSXP, L));
  
  double *b = REAL(beta);

  double *thetas = REAL(thetas_export);
  
  double *sigmas = REAL(sigmas_export);

  for (int j=0; j<(L*p); j++) b[j] = 0;

  PROTECT(loss = allocVector(REALSXP, L));

  PROTECT(iter = allocVector(INTSXP, L));

  double theta =REAL(theta_)[0];
  
  double sigma2 = pow(REAL(sigma_)[0], 2);
  
  int count_max =INTEGER(counter_)[0];

  double delta= 0;
  
  for (int i=0; i<L; i++){ 
    
    INTEGER(iter)[i] = 0;

  }

  double *a = Calloc(p, double); // Beta from previous iteration

  for (int j = 0; j < p; j++){
    a[j] = 0;
  }

  const char *penalty = CHAR(STRING_ELT(penalty_, 0));
  
  const char *variance = CHAR(STRING_ELT(variance_, 0));

  double *lambda0s = REAL(lambda0s_);

  double eps = REAL(eps_)[0];
  
  int max_iter = INTEGER(max_iter_)[0];
  
  double *r = Calloc(n, double);
  
  for (int i=0; i<n; i++) r[i] = y[i];
  
  for (int j = 0; j < p; j++){
    for(int i = 0; i < n; i++) r[i] -= a[j]*X[j*n + i];
  }
  
  double *z = Calloc(p, double);
  
  for (int j=0; j<p; j++) z[j] = crossprod(X, r, n, j);
  
  int *e1 = Calloc(p, int); // Index of an active set
  
  for (int j=0; j<p; j++) e1[j] = 1;
  
  int *e2 = Calloc(p, int); // Index of an eligible set from the strong rule
  
  for (int j=0; j<p; j++) e2[j] = 0;
  

  int converged, counter=0;
 
  double aa= REAL(a_)[0];;

  double bb= REAL(b_)[0];;
  
  double nu = REAL(nu_)[0];
  
  double xi = REAL(xi_)[0];

  double lambda1 =REAL(lambda1_)[0];

  double lambda0;
 
  // Regularization Path

  for (int l=0;l<L;l++) {

    R_CheckUserInterrupt();

    lambda0=lambda0s[l];



	while (INTEGER(iter)[l] < max_iter) {

	  // Solve over the active set

	  INTEGER(iter)[l]++;
	  
	  

	  for (int j=0; j<p; j++) {


	       delta = threshold(theta, sigma2, lambda1, lambda0, n);

	      z[j] = crossprod(X, r, n, j) + n*a[j];

	      // Update beta_j
	   
	      b[l*p+j] = SSL(z[j], a[j],lambda0,lambda1,theta,1,n,delta, sigma2);

	      // Update r
	    
	     double shift = b[l*p+j] - a[j];
	      
	     if (shift !=0) for (int i=0;i<n;i++) r[i] -= shift*X[j*n+i];

	     counter++;
	    

	    if (counter==count_max){
	      if(strcmp(penalty, "adaptive")==0){
	        theta=expectation_approx(b, aa, bb,p,l);
	      }
	      if(strcmp(variance, "unknown")==0){
	        sigma2=update_sigma(r, n, nu, xi, l);
	        
	      }

	      
	      counter=0;
	    }

	  } 

	  // Check for convergence
	  
	  converged = checkConvergence(b, a, eps, l, p);
	  
	  if ( n * sigma2 < 1) {
	    for (int j=0; j<p; j++) a[j] = 0;
	    for (int i=0; i<n; i++) r[i] = y[i];
	  } else{ 
	    for (int j=0; j<p; j++) a[j] = b[l*p+j];
	  } 
	  

	  if (converged) {
	    thetas[l]=theta;
      sigmas[l] = sqrt(sigma2);
	    break;}
	}


	if (INTEGER(iter)[l] == max_iter){
	  thetas[l] = NAN;
	  sigmas[l] = NAN;
	}

  
  }

  res = cleanupG(a, r, e1, e2, z, thetas, sigmas, beta, loss, iter, thetas_export, sigmas_export);

  return(res);
}
