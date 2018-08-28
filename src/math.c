/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "robo.h"

#ifndef __MEDIAN_H__
#define __MEDIAN_H__

double gaussian(double x, double m, double s, double A) {
  return(A * exp(-0.5 * pow((x - m)/s,2)));
}

double median(double *values, long number) {
  int index = number / 2;
  int i,j,l,m;

  double result = 0;
  double temp;
  l = 0;
  m = number - 1;
  while (l < m) {
    result = values[index];
    i = l;
    j = m;
    while (i <= j) {
      while (values[i] < result) {
	i++;
      }
      while (values[j] > result) {
	j--;
      }
      if (i <= j) {
	temp = values[j];
	values[j] = values[i];
	values[i] = temp;
	i++;
	j--;
      }
    }
    if (j < index) {
      l = i;
    }
    if (index < i) {
      m = j;
    }
  }
  return(result);
}

#endif

int comp (const void *a, const void *b) {
  const double *da = (const double *) a;
  const double *db = (const double *) b;
  return (*da > *db) - (*da < *db);
}

stats *array_stats_safe(double *x, int N) {
  int i;
  stats *S = malloc(sizeof(stats));

  S->N = N;
  S->mean = 0;
  S->sigma = 0;
  S->min = 99e99;
  S->max = -99e99;
  S->med = 0;
  
  for (i = 0; i < N; i++) {
    S->mean += x[i];
    if (x[i] < S->min) {
      S->min = x[i];
    }
    if (x[i] > S->max) {
      S->max = x[i];
    }
  }
  S->mean /= N;
  
  for (i = 0; i < N; i++) {
    S->sigma += pow( x[i] - S->mean,2);
  }
  S->sigma = sqrt(S->sigma / (N - 1));

  /* S->med = median(x,(long) N); */

  return(S);
}
stats *array_stats(double *x, int N) {
  int i;
  stats *S = malloc(sizeof(stats));
  double *MADarr = malloc(sizeof(double) * N);
  S->N = N;
  S->mean = 0;
  S->sigma = 0;
  S->min = 99e99;
  S->max = -99e99;
  S->med = 0;
  S->MAD = 0;
  
  for (i = 0; i < N; i++) {
    S->mean += x[i];
    if (x[i] < S->min) {
      S->min = x[i];
    }
    if (x[i] > S->max) {
      S->max = x[i];
    }
  }
  S->mean /= N;
  
  for (i = 0; i < N; i++) {
    S->sigma += pow( x[i] - S->mean,2);
  }
  S->sigma = sqrt(S->sigma / (N - 1));

  S->med = median(x,(long) N);
  for (i = 0; i < N; i++) {
    MADarr[i] = fabs(x[i] - S->med);
  }
  S->MAD = median(MADarr,(long) N);
  free(MADarr);
  return(S);
}

stats *histogram_and_gaussfit(double *x, int N) {
  int i,j;
  int M = (int) sqrt(N);
  stats *S = malloc(sizeof(stats));
  double *X = malloc(sizeof(double) * M);
  double *Y = malloc(sizeof(double) * M);
  double *E = malloc(sizeof(double) * M);
  double F,dm,ds,dA,chi;
  S->N = N;
  S->mean = 0;
  S->sigma = 0;
  S->min = 99e99;
  S->max = -99e99;
  S->med = 0;
  
  for (i = 0; i < N; i++) {
    S->mean += x[i];
    if (x[i] < S->min) {
      S->min = x[i];
    }
    if (x[i] > S->max) {
      S->max = x[i];
    }
  }
  S->mean /= N;
  
  for (i = 0; i < N; i++) {
    S->sigma += pow( x[i] - S->mean,2);
  }
  S->sigma = sqrt(S->sigma / (N - 1));

  for (i = 0; i < M; i++) {
    X[i] = (S->mean -3.0 * S->sigma) + (1.0 * i)/M * 6.0 * S->sigma;
    Y[i] = 0;
  }
  j = 0;
  for (i = 0; i < N; i++) {
    j = (x[i] - X[0]) / (X[M-1] - X[0]) * M;
    if ((j >= 0)&&(j < M)) {
      Y[j] ++;
    }
  }
  F = 0;
  for (i = 0; i < M; i++) {
    if (Y[i] > 0) {
      E[i] = sqrt(Y[i]);
      if (Y[i] > F) {
	F = Y[i];
      }
    }
    else {
      E[i] = 1.0;
    }
  }

  gaussfit(X,Y,E,M,
	   &(S->mean),&(S->sigma),&F,
	   &dm,&ds,&dA,
	   1,1,1,2.0,
	   &chi,1e-3,20);
  free(X);
  free(Y);
  free(E);
  return(S);
}

/* double gaussian(double x, double m, double s, double F) { */
/*   return(pow(s * sqrt(2 * M_PI),-1) * exp(-0.5 * pow( (x - m) / s,2))); */
/* } */


double equivalent_width(lines *L, int i) {
  return(-1.0 * L->F[i] * 1000); /* Return mA like Julie uses, apparently. */
}
double equivalent_width_alt(lines *L, int i) {
  return(-1.0 * L->Fp[i] * 1000); /* Return mA like Julie uses, apparently. */
}
double equivalent_width_error(lines *L, int i) {
  return(L->dF[i] * 1000); /* Return mA like Julie uses, apparently. */
}

double *make_psf(double s, int N) {
  double *v = malloc(N * sizeof(double));
  int i;

  for (i = 0; i < N; i++) {
    if (i < N / 2) {
      v[i] = gaussian(i,0.0,s,1.0);
    }
    else {
      v[i] = gaussian(i,1.0 * N,s,1.0);
    }
  }
  return(v);
}

void robust_linear_fit (double *x, double *y, int N, double *m, double *b, double tolerance) {
  double m_in = *m;
  double b_in = *b;

  double m_old = m_in;
  double b_old = b_in;

  double delta_m = 99e99;
  double delta_b = 99e99;

  double S,Sx,Sy,Sxx,Syy,Sxy,D;
  int iterations = 0;
  /*  double *e = malloc(sizeof(double) * N); */
  double *w = malloc(sizeof(double) * N);
  int i;

  for (i = 0; i < N; i++) {
    /*    e[i] = fabs((m_in * x[i] + b_in) - y[i]); */
    w[i] = pow((m_in * x[i] + b_in) - y[i],-2);
    if (!isfinite(w[i])) {
      w[i] = 1e6;
    }
  }

  while ((delta_m > tolerance)&&(delta_b > tolerance)) {
    S = 0; Sx = 0; Sy = 0;
    Sxx = 0; Sxy = 0; Syy = 0;

    for (i = 0; i < N; i++) {
      S   += w[i];
      Sx  += x[i] * w[i];
      Sy  += y[i] * w[i];
      Sxx += pow(x[i],2) * w[i];
      Syy += pow(y[i],2) * w[i];
      Sxy += x[i] * y[i] * w[i];
    }

    D = (S * Sxx - Sx * Sx);
    if (D == 0) {
      *m = m_old;
      *b = b_old;
      return;
    }

    m_old = *m;
    b_old = *b;

    *m = (S * Sxy - Sx * Sy) / D;
    *b = (Sy * Sxx - Sx * Sxy) / D;

    delta_m = fabs(m_old - *m);
    delta_b = fabs(b_old - *b);

    for (i = 0; i < N; i++) {
      w[i] = pow((*m * x[i] + *b) - y[i],-2);
      if (!isfinite(w[i])) {
	w[i] = 1e-10;
      }
    }
    iterations++;
  }
  return;
}
    

  
double interp(double x0, double *X, double *Y, int N) {
  double v = 0.0;
  int i;
 
  for (i = 0; i < N - 1; i++) {
    if ((X[i] < x0)&&(X[i+1] > x0)) {
      v = Y[i] + (x0 - X[i]) * (Y[i+1] - Y[i]) / (X[i+1] - X[i]);
      break;
    }
  }
  return(v);
}
