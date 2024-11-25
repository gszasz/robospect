/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>



#include "robo.h"

struct lmmultigauss_parm {
  double *m;
  double *s;
  double *A;

  double *dm;
  double *ds;
  double *dA;

  double chi;
  int iteration;
  
  int M;

  double *xdata;
  double *ydata;
  double *error;
  int N;
};

double multi_f(double x, double *m, double *s, double *A, int M) {
  double V = 0;
  int i;
  
  for (i = 0; i < M; i++) {
    V += A[i] * exp(-0.5 * pow((x - m[i])/s[i],2));
  }
  return(V);
}
double multi_dfdA(double x, double *m, double *s, double *A, int j) {
  return(exp(-0.5 * pow((x - m[j])/s[j],2)));
}
double multi_dfdm(double x, double *m, double *s, double *A, int j) {
  return(A[j] * (x - m[j])/pow(s[j],2) * exp(-0.5 * pow((x - m[j])/s[j],2)));
}
double multi_dfds(double x, double *m, double *s, double *A, int j) {
  return(A[j] * pow(x - m[j],2)/pow(s[j],3) * exp(-0.5 * pow((x - m[j])/s[j],2)));
}


int lm_multigaussfit_f(const gsl_vector *x, void *params, gsl_vector *f) {
  struct lmmultigauss_parm data = *(struct lmmultigauss_parm *) params;
  int i;
  double V;
  double R = 0;
  for (i = 0; i < data.M; i++) {
    data.m[i] = gsl_vector_get(x,3 * i);
    data.s[i] = gsl_vector_get(x,3 * i + 1);
    data.A[i] = gsl_vector_get(x,3 * i + 2);
  }
  for (i = 0; i < data.N; i++) {
    V = multi_f(data.xdata[i],data.m,data.s,data.A,data.M);
    if (fabs(V) > 9.8659e-10) {
      gsl_vector_set(f,i,(data.ydata[i] - V)/data.error[i]);
      /* if (V != 0) { */
      /* 	fprintf(stderr,"%g -> %g %g\n",data.xdata[i],V,data.ydata[i]); */
      R += data.ydata[i] - V;
    }
    else {
      gsl_vector_set(f,i,0.0);
    }
  }
  fprintf(stderr,"    %g\n",R);
  return(GSL_SUCCESS);
}

int lm_multigaussfit_df(const gsl_vector *x, void *params, gsl_matrix *J) {
  struct lmmultigauss_parm data = *(struct lmmultigauss_parm *) params;
  int i,j;
  double V;
  
  for (i = 0; i < data.M; i++) {
    data.m[i] = gsl_vector_get(x,3 * i);
    data.s[i] = gsl_vector_get(x,3 * i + 1);
    data.A[i] = gsl_vector_get(x,3 * i + 2);
  }

  for (i = 0; i < data.N; i++) {
    V = multi_f(data.xdata[i],data.m,data.s,data.A,data.M);
    for (j = 0; j < data.M; j++) {
      if (fabs(V) > 9.8659e-10) {
	gsl_matrix_set(J,i,3 * j,multi_dfdm(data.xdata[i],data.m,data.s,data.A,j) / data.error[i]);
	gsl_matrix_set(J,i,3 * j + 1,multi_dfds(data.xdata[i],data.m,data.s,data.A,j) / data.error[i]);
	gsl_matrix_set(J,i,3 * j + 2,multi_dfdA(data.xdata[i],data.m,data.s,data.A,j) / data.error[i]);
      }
      else {
      	gsl_matrix_set(J,i,3 * j,0.0);
      	gsl_matrix_set(J,i,3 * j + 1,0.0);
      	gsl_matrix_set(J,i,3 * j + 2,0.0);
      }
    }
  }
  return(GSL_SUCCESS);
}

int lm_multigaussfit_fdf(const gsl_vector *x, void *params, gsl_vector *f, gsl_matrix *J) {
  lm_multigaussfit_f(x,params,f);
  lm_multigaussfit_df(x,params,J);
  return(GSL_SUCCESS);
}

int lm_multigauss(double *X, double *Y, double *E, int N,
	       double *m, double *s, double *A, int M,
	       double *dm, double *ds, double *dA,
	       int *vm, int *vs, int *vA,
	       double relax,
	       double *chi, double tolerance, int max_iter) {
  struct lmmultigauss_parm data;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *S;
  int status = GSL_CONTINUE;

  gsl_matrix *covar = gsl_matrix_alloc(3 * M,3 * M);
  gsl_multifit_function_fdf f;
  double *x_init;
  gsl_vector_view x;
  int i,j;
  
  x_init = malloc(3 * M * sizeof(double));
  if (x_init == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  data.m = m;
  data.s = s;
  data.A = A;
  data.dm = dm;
  data.ds = ds;
  data.dA = dA;
  data.M = M;
  
  data.xdata = X;
  data.ydata = Y;
  data.error = E;
  data.N = N;

  for (i = 0; i < M; i++) {
    x_init[3 * i]     = m[i];
    x_init[3 * i + 1] = s[i];
    x_init[3 * i + 2] = A[i];
    fprintf(stderr,"LOAD: %d %f %f %f\n",i,x_init[3 * i],x_init[3 * i + 1],x_init[3 * i + 2]);
  }
  
  f.f = &lm_multigaussfit_f;
  f.df = &lm_multigaussfit_df;
  f.fdf = &lm_multigaussfit_fdf;

  f.n = data.N;
  f.p = 3 * data.M;
  f.params = &data;

  x = gsl_vector_view_array(x_init,3 * M);
  T = gsl_multifit_fdfsolver_lmder;
  S = gsl_multifit_fdfsolver_alloc(T,f.n,f.p);
  data.iteration = 0;
  gsl_multifit_fdfsolver_set(S,&f,&x.vector);

  status = GSL_CONTINUE;

  for (j = 0; j < M; j++) {
    fprintf(stderr,"mml: (%d/%d) %d %f %f (%f,%f,%f) +- (%f,%f,%f)\n",
	    j,M,data.iteration,-1.0,0.0,data.m[j],data.s[j],data.A[j],data.dm[j],data.ds[j],data.dA[j]);
  }
  
  while ((status == GSL_CONTINUE)&&(data.iteration < max_iter)) {
    data.iteration++;
    status = gsl_multifit_fdfsolver_iterate(S);
    status = gsl_multifit_test_delta(S->dx,S->x,1e-35,1e-35);
    gsl_matrix *J = gsl_matrix_alloc(f.n, f.p);
    gsl_multifit_fdfsolver_jac(S, J);
    gsl_multifit_covar(J, 0.0, covar);
    for (j = 0; j < M; j++) {
      fprintf(stderr,"mml: (%d/%d) %d %f %f (%f,%f,%f) +- (%f,%f,%f)\n",
	      j,M,data.iteration,-1.0,0.0,data.m[j],data.s[j],data.A[j],data.dm[j],data.ds[j],data.dA[j]);
    }
  }

  for (i = 0; i < M; i++) {
    m[i] = x_init[3 * i];
    s[i] = x_init[3 * i + 1];
    A[i] = x_init[3 * i + 2];
    dm[i] = sqrt(gsl_matrix_get(covar,3 * i,3 * i));
    ds[i] = sqrt(gsl_matrix_get(covar,3 * i + 1,3 * i + 1));
    dA[i] = sqrt(gsl_matrix_get(covar,3 * i + 2,3 * i + 2));
  }
  
  *chi = pow(gsl_blas_dnrm2(S->f),2.0) / (f.n - (f.p));

  for (j = 0; j < M; j++) {
    fprintf(stderr,"mmlF: (%d/%d) %d %f %f (%f,%f,%f) +- (%f,%f,%f)\n",
	    j,M,data.iteration,-1.0,0.0,data.m[j],data.s[j],data.A[j],data.dm[j],data.ds[j],data.dA[j]);
    }

  
  gsl_multifit_fdfsolver_free(S);
  return(data.iteration);
}
