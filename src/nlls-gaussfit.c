/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/* Panda has a terrible libc install that doesn't have NAN? */
#ifndef NAN
# define NAN (__builtin_nanf (""))
#endif


#include "robo.h"
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

/* Invert a 3x3 matrix directly.                                       */
int invert3x3(double **in, double **out) {
  int i,j;
  double det= (in[0][0] * (in[1][1] * in[2][2] -
			   in[1][2] * in[2][1]) -
	       in[1][0] * (in[0][1] * in[2][2] -
			   in[2][1] * in[0][2]) +
	       in[2][0] * (in[0][1] * in[2][1] -
			   in[1][1] * in[0][2]));
  if (det == 0) {
    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++) {
	out[i][j] = 0;
	return(0);
      }
    }
  }
  else {
    out[0][0] = (in[1][1] * in[2][2] - in[1][2] * in[2][1]) / det;
    out[0][1] = (in[0][2] * in[2][1] - in[0][1] * in[2][2]) / det;
    out[0][2] = (in[0][1] * in[1][2] - in[0][2] * in[1][1]) / det;

    out[1][0] = (in[1][2] * in[2][0] - in[1][0] * in[2][2]) / det;
    out[1][1] = (in[0][0] * in[2][2] - in[0][2] * in[2][0]) / det;
    out[1][2] = (in[0][2] * in[1][0] - in[0][0] * in[1][2]) / det;

    out[2][0] = (in[1][0] * in[2][1] - in[1][1] * in[2][0]) / det;
    out[2][1] = (in[0][1] * in[2][0] - in[0][0] * in[2][1]) / det;
    out[2][2] = (in[0][0] * in[1][1] - in[0][1] * in[1][0]) / det;
    return(1);
  }
  return(-1);
}

/* #define DEBUG 1 */

/* Linear least squares fitting code.                                    */
int gaussfit(double *X, double *Y, double *E, int N,
	     double *m, double *s, double *A,
	     double *dm, double *ds, double *dA,
	     int vm, int vs, int vA, double relax,
	     double *chi, double tolerance, int max_iter) {
  int iteration = 0;
  double *dB;                  /* Array of deviations                     */
  double **a;                  /* Matrix of derivatives at each position. */
  double L[3];                 /* Proposed update values                  */
  double **C;                  /* Covariance Matrix                       */
  double **ata;
  double atb[3];

  int j,k;
  int i;
  double chi_pixel = 0.0;
  double xms;

#ifdef DEBUG
    double mi,si,Ai;
#endif

  /* Modify *A to from the F that's passed: */
  *A = *A * pow(*s * sqrt(2 * M_PI),-1);
  /* Allocate space for matrices. */
  dB = malloc(N * sizeof(double));
  a  = malloc(N * sizeof(double *));
  a[0] = malloc(N * 3 * sizeof(double));
  for (i = 1; i < N; i++) {
    a[i] = a[0] + i * 3;
  }

  C = malloc(3 * sizeof(double *));
  C[0] = malloc(9 * sizeof(double));
  ata = malloc(3 * sizeof(double *));
  ata[0] = malloc(9 * sizeof(double));
  for (i = 1; i < 3; i++) {
    C[i] = C[0] + i * 3;
    ata[i] = ata[0] + i * 3;
  }

  /* Assume the first guess isn't too bad, and set first updates to 1%. */
  L[0] = *m * 0.01;
  L[1] = *s * 0.01;
  L[2] = *A * 0.01;


  /* Backup input values in case we need to reset. */
#ifdef DEBUG
  mi = *m;
  si = *s;
  Ai = *A;
  fprintf(stderr,"gaussfit: initial values: %g %g %g %d\n",
	  *m,*s,*A,N);
#endif

  /* Iterate until either the fractional proposed update is smaller than */
  /* the tolerance, or we reach 1000 iterations (which suggests no fit). */
  while (((fabs(L[0] / *m) > tolerance)||
	  (fabs(L[1] / *s) > tolerance)||
	  (fabs(L[2] / *A) > tolerance))&&
	 (iteration <= max_iter)) {
    *chi = 0;
    chi_pixel = 0.0;

#ifdef DEBUG
    fprintf(stderr,"gaussfit: iteration %d: %f %f (%f %f %f)\n",iteration,*chi,chi_pixel,L[0],L[1],L[2]);
#endif
    /* Calculate current deviations, and the current derivative matrix. */
    for (i = 0; i < N; i++) {
      if (!isfinite(E[i])||(E[i] == 0.0)) {
	continue;
      }
      xms = (X[i] - *m) / *s;
      dB[i] = (Y[i] - f(X[i],*m,*s,*A)) / E[i];
      /* a[i][0] = dfdm(X[i],*m,*s,*A) / E[i]; */
      /* a[i][1] = dfds(X[i],*m,*s,*A) / E[i]; */
      /* a[i][2] = dfdA(X[i],*m,*s,*A) / E[i]; */
      a[i][2] = dfdA(X[i],*m,*s,*A) / E[i];
      if (a[i][2] == 0) {
	a[i][0] = 0;
	a[i][1] = 0;
      }
      else {
	a[i][0] = a[i][2] * *A * xms / *s;
	a[i][1] = a[i][0] * xms;
      }
      /* fprintf(stderr,"FF: %d %g %g %g %g :: %g %g %g\n",i,xms,a[i][0],a[i][1],a[i][2], */
      /* 	      dfdm(X[i],*m,*s,*A) / E[i], */
      /* 	      dfds(X[i],*m,*s,*A) / E[i], */
      /* 	      dfdA(X[i],*m,*s,*A) / E[i]); */

      if (fabs(f(X[i],*m,*s,*A)) > 9.8659e-10) { /* 6 sigma */
	*chi += pow(dB[i],2);
	chi_pixel += 1.0;
      }
    }
    if (chi_pixel > 0.0) {
      *chi /= chi_pixel;
#ifdef DEBUG
      fprintf(stderr,"gaussfit: iteration %d chi: %g %g\n",iteration,chi_pixel,*chi);
#endif
    }
    /* Convert system of Ndata equations into system of Nparm equations*/
    /* {dB}       = {A} {dL}                                           */
    /* {A^T} {dB} = {A^T} {A} {dL}                                     */
    /* This is done, because {A} is in general non-square, and the     */
    /* inverse is only defined for a square matrix.  However, {A^T} {A}*/
    /* is always square.                                               */
    for (j = 0; j < 3; j++) {
      atb[j] = 0;
      for (k = 0; k < 3; k++) {
	ata[j][k] = 0;
	for (i = 0; i < N; i++) {
	  ata[j][k] += a[i][j] * a[i][k];
	  if (k == 0) {
	    atb[j] += a[i][j] * dB[i];
	  }
	}
      }
    }
    /* Calculate the covariance matrix, which is just the inverse of ata*/
    /* The diagonals of this matrix hold the estimates of the           */
    /* uncertainties in the parameters.                                 */
    j = invert3x3(ata,C);
#ifdef DEBUG
    fprintf(stderr,"gaussfit: iteration %d inverse: %d\n",iteration,j);
#endif

    /* Use Gaussian Elimination to solve for the proposed updates:      */
    atb[1] -= atb[0] * ata[1][0] / ata[0][0];
    atb[2] -= atb[0] * ata[2][0] / ata[0][0];
    for (k = 0; k < 3; k++) {
      ata[1][k] -= ata[0][k] * ata[1][0] / ata[0][0];
      ata[2][k] -= ata[0][k] * ata[2][0] / ata[0][0];
    }
    atb[2] -= atb[1] * ata[2][1] / ata[1][1];
    for (k = 0; k < 3; k++) {
      ata[2][k] -= ata[1][k] * ata[2][1] / ata[1][1];
    }

    L[2] = atb[2] / ata[2][2];
    L[1] = (atb[1] - (L[2] * ata[1][2])) / ata[1][1];
    L[0] = (atb[0] - (L[2] * ata[0][2]) - (L[1] * ata[0][1])) / ata[0][0];

    /* Update the values with the proposed changes, and recalculate errors*/
    iteration++;

    if (vm) {
      *m += L[0] / relax;
    }
    if (vs) {
      *s += L[1] / relax;
    }
    if (vA) { 
      *A += L[2] / relax;
    }

    *dm = sqrt(C[0][0] / chi_pixel);
    *ds = sqrt(C[1][1] / chi_pixel);
    *dA = sqrt(C[2][2] / chi_pixel);

    if ((!isfinite(*m))||(!isfinite(*s))||(!isfinite(*A))||
	(!isfinite(*dm))||(!isfinite(*ds))||(!isfinite(*dA))) {
      iteration = 100000;
    }
    if (fabs(*m) > 1e8) {      *m = NAN;    }
    if (fabs(*s) > 1e8) {      *s = NAN;    }
    if (fabs(*A) > 1e8) {      *A = NAN;    }
#ifdef DEBUG
    fprintf(stderr,"gaussfit: current values: iter: %d (m,s,A) = (%g %g %g) (delta) = (%g %g %g) on (%g %g %g)\n",
	    iteration,*m,*s,*A,*dm,*ds,*dA,L[0],L[1],L[2]);
#endif
    
  }
  /* Once the fit has converged:                                        */
  /* Free allocated memory.                                             */
  free(dB);
  free(a[0]);
  free(a);
  *A = *A / pow(*s * sqrt(2 * M_PI),-1);
  *dA = *dA / pow(*s * sqrt(2 * M_PI),-1);
#ifdef DEBUG
    fprintf(stderr,"gaussfit: final values: iter: %d (m,s,A) = (%g %g %g) (delta) = (%g %g %g)\n",
	    iteration,*m,*s,*A,*dm,*ds,*dA);
#endif

  return(iteration);
}
