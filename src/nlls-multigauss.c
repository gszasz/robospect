/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "robo.h"

/* #define DEBUG 0 */
static inline double f(double x, double m, double s, double A) {
  return(A * exp(-0.5 * pow((x - m)/s,2)));
}
static inline double dfdA2(double xms) {
  return(exp(-0.5 * pow(xms,2)));
}
static inline double dfdA(double x, double m, double s, double A) {
  return(dfdA2((x - m)/s));
}
static inline double dfdm2(double xms,double Am) {
  return(Am * pow(xms,2) * exp(-0.5 * pow(xms,2)));
}
static inline double dfdm(double x, double m, double s, double A) {
  if ((x - m) == 0) {
    return(0.0);
  }
  else {
    return(dfdm2((x-m)/s, A / (x-m)));
  }
}
static inline double dfds2(double xms,double Am) {
  return(Am * pow(xms,3) * exp(-0.5 * pow(xms,2)));
}
static inline double dfds(double x, double m, double s, double A) {
  if ((x - m) == 0) {
    return(0);
  }
  else {
    return(dfds2((x - m)/s,A / (x - m)));
  }
}

typedef struct {
  double m;
  double s;
  double A;

  double dm;
  double ds;
  double dA;
  
  double **a;
  double L[3];
  double **C;
  double **ata;
  double atb[3];
} multigauss_solver;

multigauss_solver *alloc_multigauss_solver(double m, double s, double A, int N) {
  multigauss_solver *S = malloc(sizeof(multigauss_solver));
  int i;
  
  S->a = malloc(N * sizeof(double *));
  S->a[0] = malloc(N * 3 * sizeof(double));
  for (i = 1; i < N; i++) {
    S->a[i] = S->a[0] + i * 3;
  }

  S->C = malloc(3 * sizeof(double *));
  S->C[0] = malloc(9 * sizeof(double));
  S->ata = malloc(3 * sizeof(double *));
  S->ata[0] = malloc(9 * sizeof(double));
  for (i = 1; i < 3; i++) {
    S->C[i] = S->C[0] + i * 3;
    S->ata[i] = S->ata[0] + i * 3;
  }

  S->m = m;
  S->s = s;
  S->A = A;

  S->L[0] = m * 0.01;
  S->L[1] = s * 0.01;
  S->L[2] = A * 0.01;

  return(S);
}

void free_multigauss_solver(multigauss_solver *S) {
  free(S->a[0]);
  free(S->a);
  free(S->C[0]);
  free(S->C);
  free(S->ata[0]);
  free(S->ata);
  free(S);
}

double eval_multigauss_solver(multigauss_solver *S, double x) {
  return(f(x,S->m,S->s,S->A));
}



/* Linear least squares fitting code.                                    */
int multigauss(double *X, double *Y, double *E, int N,
	       double *m, double *s, double *A, int M,
	       double *dm, double *ds, double *dA,
	       int *vm, int *vs, int *vA,
	       double relax,
	       double *chi, double tolerance, int max_iter) {
  int iteration = 0;
  gsl_vector *dB = gsl_vector_alloc(N);        /* Array of deviations */
  gsl_matrix *a  = gsl_matrix_alloc(N,M*3);    /* Matrix of derivatives at each position. */
  gsl_vector *L  = gsl_vector_alloc(M*3);      /* Proposed update values */
  gsl_matrix *C  = gsl_matrix_alloc(M*3,M*3);
  gsl_matrix *W  = gsl_matrix_alloc(M*3,M*3);
  gsl_matrix *ata = gsl_matrix_alloc(M*3,M*3);
  gsl_vector *atb = gsl_vector_alloc(M*3);
  gsl_permutation *P = gsl_permutation_alloc(M*3);

  gsl_matrix *SV_V = gsl_matrix_alloc(M*3,M*3);
  gsl_vector *SV_S = gsl_vector_alloc(M*3);
  gsl_vector *SV_WORK = gsl_vector_alloc(M*3);
  gsl_error_handler_t *old_handler;
  int is_converged = 0;
  double val = 0.0;
  double F,dF,xms;  /* Calculation speed tools. */
  int i; /* Data iterator */
  int j; /* Component iterator */

  double *mi, *si, *Ai;
  double chi_pixel = 0.0;
  double old_chi = 99e99;
  /* Save initial values in case we need to reset. */
  mi = malloc(sizeof(double) * M);
  si = malloc(sizeof(double) * M);
  Ai = malloc(sizeof(double) * M);
  for (j = 0; j < M; j++) {
    mi[j] = m[j];
    si[j] = s[j];
    A[j]  = A[j] * pow(s[j] * sqrt(2 * M_PI),-1);  /* Convert from flux */
    Ai[j] = A[j];
  }

  /* Iterate until either the fractional proposed update is smaller than */
  /* the tolerance, or we reach 1000 iterations (which suggests no fit). */
  *chi = 99e99;
  while ((is_converged == 0)&&
	 (iteration <= max_iter)) {
    is_converged = 1;
    old_chi = *chi;
    *chi = 0;
    chi_pixel = 0.0;

    /* Calculate current deviations, and the current derivative matrix. */
    gsl_vector_set_zero(dB);
    gsl_matrix_set_zero(a);
    gsl_matrix_set_zero(ata);
    gsl_vector_set_zero(atb);
    gsl_vector_set_zero(L);
    gsl_matrix_set_zero(W);
    gsl_matrix_set_zero(C);
    for (j = 0; j < M; j++) {
#ifdef DEBUG
	fprintf(stderr,"gaussfit: current values: iter: %d (m,s,A) = (%g %g %g) (delta) = (%g %g %g)\n",
		iteration,m[j],s[j],A[j],dm[j],ds[j],dA[j]);
#endif
    }
    for (i = 0; i < N; i++) {
      val = 0.0;
      for (j = 0; j < M; j++) { /* Accumulate all components */
	F = f(X[i],m[j],s[j],A[j]);
	dF = dfdA(X[i],m[j],s[j],A[j]) / E[i];
	gsl_matrix_set(a,i,2 + j * 3,dF);
	
	if (dF == 0.0) {
	  gsl_matrix_set(a,i,0 + j * 3,0.0);
	  gsl_matrix_set(a,i,1 + j * 3,0.0);
	}
	else {
	  xms = (X[i] - m[j]) / s[j];
	  gsl_matrix_set(a,i,0 + j * 3,dF * A[j] * xms / s[j]);
	  gsl_matrix_set(a,i,1 + j * 3,xms * gsl_matrix_get(a,i,0 + j * 3));
	}
	val += F;
      }	
      gsl_vector_set(dB,i,(Y[i] - val) / E[i]);
#ifdef DEBUG
      fprintf(stderr,"%d %g %g %g %g %g\n",i,gsl_vector_get(dB,i),val,X[i],Y[i],E[i]);
#endif

      if (fabs(val) > 9.8659e-10) { /* 6 sigma */
	*chi += pow(gsl_vector_get(dB,i),2);
	chi_pixel += 1.0;
      }
      else { /* If we're not within the 6 sigma limit, we don't want to include anything. */
	gsl_vector_set(dB,i,0.0);
	for (j = 0; j < M; j++) {
	  gsl_matrix_set(a,i,0+j*3,0.0);
	  gsl_matrix_set(a,i,1+j*3,0.0);
	  gsl_matrix_set(a,i,2+j*3,0.0);
	}
      }
    }
    if (chi_pixel > 0.0) {
      *chi /= chi_pixel;
#ifdef DEBUG
      fprintf(stderr,"gaussfit: iteration %d chi: %g %g\n",iteration,chi_pixel,*chi);
#endif
    }
    else {
      iteration = 200;
      continue;
    }

    /* Convert system of Ndata equations into system of Nparm equations*/
    /* {dB}       = {A} {dL}                                           */
    /* {A^T} {dB} = {A^T} {A} {dL}                                     */
    /* This is done, because {A} is in general non-square, and the     */
    /* inverse is only defined for a square matrix.  However, {A^T} {A}*/
    /* is always square.                                               */
    /* A^T A */
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,a,a,0.0,ata);
    /* A^T B */
    gsl_blas_dgemv(CblasTrans,1.0,a,dB,0.0,atb);
    /* dL = (A^T B) \ (A^T A) */
#ifdef DEBUG
    fprintf(stderr,"%g %g %g %g\n",gsl_matrix_min(ata),gsl_matrix_max(ata),gsl_vector_min(atb),gsl_vector_max(atb));
#endif

    old_handler = gsl_set_error_handler_off();
    j = 0;
    j = gsl_linalg_SV_decomp(ata,SV_V,SV_S,SV_WORK);
    if (j != 0) {
      iteration = 100000;
      continue;
    }
    j = gsl_linalg_SV_solve(ata,SV_V,SV_S,atb,L);
    if (j != 0) {
      iteration = 100000;
      continue;
    }
    gsl_set_error_handler(old_handler);

    /* Calculate the covariance matrix, which is just the inverse of ata */

    for (j = 0; j < M; j++) {
      gsl_matrix_set(W,0 + 3 * j,0 + 3 * j,pow(gsl_vector_get(SV_S,0 + 3 * j),-2.0));
      gsl_matrix_set(W,1 + 3 * j,1 + 3 * j,pow(gsl_vector_get(SV_S,1 + 3 * j),-2.0));
      gsl_matrix_set(W,2 + 3 * j,2 + 3 * j,pow(gsl_vector_get(SV_S,2 + 3 * j),-2.0));
    }
    gsl_matrix_set_zero(ata);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,SV_V,W,0.0,ata);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,ata,SV_V,0.0,C);
       
    /* Update the values with the proposed changes, and recalculate errors*/
    for (j = 0; j < M; j++) {
      if (vm[j]) {
	m[j] += gsl_vector_get(L,0 + 3 * j);
      }
      if (vs[j]) {
	s[j] += gsl_vector_get(L,1 + 3 * j);
      }
      if (vA[j]) { 
	A[j] += gsl_vector_get(L,2 + 3 * j);
      }
      
      dm[j] = sqrt(*chi * gsl_matrix_get(C,0 + 3 * j,0 + 3 * j));
      ds[j] = sqrt(*chi * gsl_matrix_get(C,1 + 3 * j,1 + 3 * j));
      dA[j] = sqrt(*chi * gsl_matrix_get(C,2 + 3 * j,2 + 3 * j));
      
      if (!isfinite( m[j])) {  m[j] = mi[j]; }
      if (!isfinite( s[j])) {  s[j] = si[j]; }
      if (!isfinite( A[j])) {  A[j] = Ai[j]; }
      if (!isfinite(dm[j])) { dm[j] = 0.0; }
      if (!isfinite(ds[j])) { ds[j] = 0.0; }
      if (!isfinite(dA[j])) { dA[j] = 0.0; }

      if (fabs(m[j] - mi[j]) > fabs(mi[j])) { m[j] = mi[j]; }
      if (fabs(s[j]) > 10 * fabs(si[j])) { s[j] = si[j]; }
      if (fabs(A[j]) > 10 * fabs(Ai[j])) { A[j] = Ai[j]; }

      if (fabs(m[j] - mi[j]) > LARGE_SHIFT * si[j]) {
	m[j] = mi[j];
	vm[j] = 0;
      }
      if (s[j] < 0.0) { s[j] = si[j]; }
      
      
      if ((fabs(gsl_vector_get(L,0 + 3 * j) / m[j]) > tolerance)&&
	  (fabs(gsl_vector_get(L,0 + 3 * j)) > tolerance))  { is_converged = 0; }
      if ((fabs(gsl_vector_get(L,1 + 3 * j) / s[j]) > tolerance)&&
	  (fabs(gsl_vector_get(L,1 + 3 * j)) > tolerance))  { is_converged = 0; }
      if ((fabs(gsl_vector_get(L,2 + 3 * j) / A[j]) > tolerance)&&
	  (fabs(gsl_vector_get(L,2 + 3 * j)) > tolerance))  { is_converged = 0; }

      if (fabs((old_chi - *chi)/old_chi) < tolerance) { is_converged = 1; }
      
      if ((!isfinite( m[j]))||(!isfinite( s[j]))||(!isfinite( A[j]))||
	  (!isfinite(dm[j]))||(!isfinite(ds[j]))||(!isfinite(dA[j]))) {
	iteration = 100000;
      }
#ifdef DEBUG
      fprintf(stderr,"gaussfit: current values: iter: %d (m,s,A) = (%g %g %g) (delta) = (%g %g %g) (%g %g %g @ %g) \n",
	      iteration,m[j],s[j],A[j],dm[j],ds[j],dA[j],
	      fabs(gsl_vector_get(L,0 + 3 * j) / m[j]),
	      fabs(gsl_vector_get(L,1 + 3 * j) / s[j]),
	      fabs(gsl_vector_get(L,2 + 3 * j) / A[j]),tolerance
	      );
#endif
    }
    iteration++;
  }
  /* Once the fit has converged:                                        */
  /* Free allocated memory.                                             */
  gsl_vector_free(dB);
  gsl_matrix_free(a);
  gsl_vector_free(L);
  gsl_matrix_free(C);
  gsl_matrix_free(ata);
  gsl_vector_free(atb);
  gsl_permutation_free(P);
  free(mi);
  free(si);
  free(Ai);
  for (j = 0; j < M; j++) {
     A[j] =  A[j] / pow(s[j] * sqrt(2 * M_PI),-1);
    dA[j] = dA[j] / pow(s[j] * sqrt(2 * M_PI),-1);
  }
  
  return(iteration);
}
