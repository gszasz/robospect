/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

double gmcps_UR(double x1, double x2) {
  double Rsq = pow(x1 - x2,2);
  if (Rsq <= 0) {
    return(0);
  }
  else {
    return(Rsq * 0.5 * log(Rsq));
  }
}

void generate_model_continuum_penalized_spline(opts *options, data *D, model *M) {
  double lambda = options->continuum_box;
  int i,j;
  gsl_matrix *P = gsl_matrix_alloc(D->N,2);
  gsl_vector *b = gsl_vector_alloc(D->N + 2);
  gsl_vector *w = gsl_vector_alloc(D->N);
  gsl_vector *a = gsl_vector_alloc(2);

  gsl_matrix *L = gsl_matrix_alloc(D->N + 2, D->N + 2);
  gsl_permutation *p = gsl_permutation_alloc(D->N + 2);
  gsl_vector *x = gsl_vector_alloc(D->N + 2);

  gsl_matrix *K = gsl_matrix_alloc(D->N,D->N);
  gsl_vector *y = gsl_vector_alloc(D->N);
    
  double alpha = 0.0;
  double E = 0.0;
  double z;
  /* tps load */
  fprintf(stderr,"Loading model.");
  for (i = 0; i < D->N; i++) {
    gsl_matrix_set(P,i,0,1);
    gsl_matrix_set(P,i,1,D->x[i]);
    gsl_vector_set(b,i,D->y[i]);    
  }
  fprintf(stderr,".");
  for (i = D->N; i < D->N + 2; i++) {
    gsl_vector_set(b,i,0);
  }
  fprintf(stderr,".");
  for (i = 0; i < D->N; i++) {
    for (j = 0; j < D->N; j++) {
      alpha += abs(D->x[i] - D->x[j]);
    }
  }
  alpha /= pow(D->N,2);
  fprintf(stderr,".\n");
  /* tps solve */
  fprintf(stderr,"Solving model.");
  for (i = 0; i < D->N + 2; i++) {
    for (j = 0; j < D->N + 2; j++) {
      if ((i < D->N)&&(j < D->N)) {
	if (i == j) {
	  gsl_matrix_set(L,i,j,gmcps_UR(D->x[i],D->x[j]) + pow(alpha,2) * lambda);
	  gsl_matrix_set(K,i,j,gsl_matrix_get(L,i,j));
	}
	else {
	  gsl_matrix_set(L,i,j,gmcps_UR(D->x[i],D->x[j]));
	}
      }
      else if ((i < D->N)&&(j >= D->N)) {
	gsl_matrix_set(L,i,j,gsl_matrix_get(P,i,j - D->N));
      }
      else if ((j < D->N)&&(j >= D->N)) {
	gsl_matrix_set(L,i,j,gsl_matrix_get(P,j,i - D->N));
      }
      else {
	gsl_matrix_set(L,i,j,0);
      }
    }
  }
  fprintf(stderr,".(dcomp)");
  gsl_linalg_LU_decomp(L,p,&i);
  fprintf(stderr,".(solve)");
  gsl_linalg_LU_solve(L,p,b,x);
  fprintf(stderr,".(transform)");

  for (i = 0; i < D->N; i++) {
    gsl_vector_set(w,i,gsl_vector_get(x,i));
    gsl_vector_set(y,i,0);
  }
  for (i = D->N; i < D->N + 2; i++) {
    gsl_vector_set(a,i - D->N, gsl_vector_get(x,i));
  }
  fprintf(stderr,".(solveE)");
  E = 0.0;
  gsl_blas_dgemv(CblasNoTrans,1.0,K,w,0.0,y);
  gsl_blas_ddot(w,y,&E);
  fprintf(stderr,".(done)\n");
  gsl_matrix_free(K);
  gsl_vector_free(y);
  gsl_matrix_free(L);
  gsl_permutation_free(p);
  gsl_vector_free(x);

  /* tps apply */
  fprintf(stderr,"Apply.");
  for (i = 0; i < D->N; i++) {
    z = gsl_vector_get(a,0) + gsl_vector_get(a,1) * D->x[i];
    for (j = 0; j < D->N; j++) {
      z += gsl_vector_get(w,i) * gmcps_UR(D->x[i],gsl_matrix_get(P,i,1));
    }
    M->continuum[i] = z;
  }
  fprintf(stderr,".(return)\n");
}
