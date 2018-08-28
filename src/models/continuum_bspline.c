/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void generate_model_continuum_bspline(opts *options,data *D, model *M) {
  int i,j;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  gsl_vector *x;
  gsl_vector *y;

  gsl_vector *c,*w;
  gsl_matrix *X, *cov;
  gsl_multifit_linear_workspace *mw;

  int nbreak = (D->N / options->continuum_box) + 1;
  int ncoeff = nbreak + 2;

  double chi,v,dv;

  bw = gsl_bspline_alloc(4, nbreak);
  B  = gsl_vector_alloc(ncoeff);

  x = gsl_vector_alloc(D->N);
  y = gsl_vector_alloc(D->N);
  w = gsl_vector_alloc(D->N);

  for (i = 0; i < D->N; i++) {
    gsl_vector_set(x,i,D->x[i]);
    gsl_vector_set(y,i,D->y[i]);
    gsl_vector_set(w,i,pow(0.05,-2));
  }

  X = gsl_matrix_alloc(D->N,ncoeff);
  c = gsl_vector_alloc(ncoeff);
  cov = gsl_matrix_alloc(ncoeff,ncoeff);
  mw = gsl_multifit_linear_alloc(D->N, ncoeff);


  fprintf(stderr,"in bspline: %ld %g %d %d\n",D->N,options->continuum_box,nbreak,ncoeff);
  
  gsl_bspline_knots_uniform(D->x[0],D->x[D->N - 1],bw);

  for (i = 0; i < D->N; i++) {
    gsl_bspline_eval(D->x[i],B,bw);
    for (j = 0; j < ncoeff; j++) {
      gsl_matrix_set(X,i,j,gsl_vector_get(B,j));
    }
  }
  gsl_multifit_wlinear(X,w,y,c,cov,&chi,mw);
  fprintf(stderr,"solved: %g %g %ld\n",
	  chi,
	  chi,
	  D->N - ncoeff);
  
  for (i = 0; i < D->N; i++) {
    gsl_bspline_eval(D->x[i],B,bw);
    gsl_multifit_linear_est(B,c,cov,&v,&dv);
    M->continuum[i] = v;
    if (i % 10 == 0) {
      fprintf(stderr,"%d %g %g %g\n",i,D->x[i],v,dv);
    }
  }

  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);  
  
}
