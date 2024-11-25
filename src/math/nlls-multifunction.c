/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


#define DEBUG 0
#define USE_GSL_SOLVER 0
#define USE_NEW_CODE 0


#if (! USE_GSL_SOLVER)

#if (USE_NEW_CODE)
int eval_chi(double *X, double *Y, double *E, int N,
	     int fit_function, int fit_parameters,
	     double *m, double *s, double *A, double *eta, int M,
	     gsl_vector *dB,
	     gsl_matrix *a,
	     double *chi) {
  int i;
  int j;
  double val,F,d0,d1,d2,d3;
  int keep;
  double chi_pixel = 0.0;
  for (i = 0; i < N; i++) {
    val = 0.0;
    for (j = 0; j < M; j++) { /* Accumulate all components */
      /* Set matrix values based on function and derivatives */
      if (fit_function == FUNCTION_GAUSSIAN) {
	gaussian_model(X[i],m[j],s[j],A[j],eta[j],
		       &F,&d0,&d1,&d2,&d3);
	d3 = 0.0;
      }
      else if (fit_function == FUNCTION_LORENTZIAN) {
	lorentzian_model(X[i],m[j],s[j],A[j],eta[j],
			 &F,&d0,&d1,&d2,&d3);
	d3 = 0.0;
      }
      else if (fit_function == FUNCTION_PSEUDOVOIGT) {
	pseudovoigt_model(X[i],m[j],s[j],A[j],eta[j],
			  &F,&d0,&d1,&d2,&d3);
      }
      else if (fit_function == FUNCTION_HJERTING) {
	hjerting_model(X[i],m[j],s[j],A[j],eta[j],
		       &F,&d0,&d1,&d2,&d3);
      }
      else if (fit_function == FUNCTION_HUMLICEK) {
	humlicek_model(X[i],m[j],s[j],A[j],eta[j],
		       &F,&d0,&d1,&d2,&d3);
      }
      else if (fit_function == FUNCTION_SKEWGAUSS) {
	skewgauss_model(X[i],m[j],s[j],A[j],eta[j],
		       &F,&d0,&d1,&d2,&d3);
      }
      else if (fit_function == FUNCTION_PLANCK) {
	planck_model(X[i],m[j],s[j],A[j],eta[j],
		     &F,&d0,&d1,&d2,&d3);
      }
      else {
	F = 0.0;
	d0 = 0.0;
	d1 = 0.0;
	d2 = 0.0;
	d3 = 0.0;
      }
      gsl_matrix_set(a,i,0 + j * fit_parameters,d0 / E[i]);
      gsl_matrix_set(a,i,1 + j * fit_parameters,d1 / E[i]);
      gsl_matrix_set(a,i,2 + j * fit_parameters,d2 / E[i]);
      gsl_matrix_set(a,i,3 + j * fit_parameters,d3 / E[i]);
      
      /*      gsl_matrix_set(ZZ,i,j,F); */
      val += F;
    } /* End loop over components. */
    gsl_vector_set(dB,i,(Y[i] - val) / E[i]);
    
    /* Should we use this point in the calculations? */
    keep = 0;
    for (j = 0; j < M; j++) {
      if (fabs((X[i] - m[j]) / (sqrt(2.0) * s[j])) < 7) {
	keep = 1;
      }
    }
    if (fit_function == FUNCTION_PLANCK) {
      keep = 1;
    }
    if ((fit_function == FUNCTION_HUMLICEK)&&(keep == 0)) {
      double F_core;
      double F_test;
      humlicek_model(X[i],m[j],s[j],A[j],eta[j],
		     &F_test,&d0,&d1,&d2,&d3);
      humlicek_model(m[j],m[j],s[j],A[j],eta[j],
		     &F_core,&d0,&d1,&d2,&d3);
      if (F_test / F_core > 1e-2) {
	keep = 1;
      }
    }
    
    /* Calculate chi */
    if (keep == 1) {
      *chi += pow(gsl_vector_get(dB,i),2);
      chi_pixel += 1.0;
    }
    else { /* If we're not within the 6 sigma limit, we don't want to include anything. */
      gsl_vector_set(dB,i,0.0);
      for (j = 0; j < M; j++) {
	gsl_matrix_set(a,i,0+j*fit_parameters,0.0);
	gsl_matrix_set(a,i,1+j*fit_parameters,0.0);
	gsl_matrix_set(a,i,2+j*fit_parameters,0.0);
	gsl_matrix_set(a,i,3+j*fit_parameters,0.0);
      }
    }
  } /* End loop over pixels */
  
  if (chi_pixel > 0.0) {
    *chi /= chi_pixel;
  }
  else {
    *chi = NAN;
  }
  return(0);
}

int multifunction(double *X, double *Y, double *E, int N,
		  int fit_function,
		  double *m, double *s, double *A, double *eta, int M,
		  double *dm, double *ds, double *dA, double *deta, 
		  int *vm, int *vs, int *vA, int *veta,
		  double relax,
		  double *chi, double tolerance, int max_iter) {
  int iteration = 0;
  int fit_parameters = 4;
  gsl_vector *dB;        /* Array of deviations */
  gsl_matrix *a ;        /* Matrix of derivatives at each position. */
  gsl_vector *L1 ;        /* Proposed update values */
  gsl_vector *L2 ;        /* Proposed update values */
  gsl_vector *L = NULL;
  gsl_matrix *C ; 
  gsl_matrix *W ; 
  gsl_matrix *ata; 
  gsl_vector *atb; 
  gsl_permutation *P;

  gsl_matrix *SV_V;
  gsl_vector *SV_S;
  gsl_vector *SV_WORK;

  gsl_error_handler_t *old_handler;

  gsl_matrix *ZZ;
  gsl_vector *dchi;
  gsl_vector *dchiN;
  
  int is_converged = 0;
  double val = 0.0;
  double lambda_counter = 0;
  /*  double F,dF,xms;  // Calculation speed tools. */
  double F,d0,d1,d2,d3;
  int i; /* Data iterator */
  int j; /* Component iterator */
  int keep;
  double *mi, *si, *Ai, *etai;
  double *m1, *s1, *A1, *eta1;
  double *m2, *s2, *A2, *eta2;
  double chi_pixel = 0.0;
  double old_chi = 99e99;

  double lambda = 1.0;
  double lambda0 = lambda;
  double nu = 10;
  double chi1,chi2;

  /* Futz the tolerance */
  if (fit_function == FUNCTION_HUMLICEK) {
    tolerance /= 5.0;
  }
  
  /* Allocate vectors and matrices */
  dB = gsl_vector_alloc(N);        /* Array of deviations */
  a  = gsl_matrix_alloc(N,M*fit_parameters);    /* Matrix of derivatives at each position. */
  L1  = gsl_vector_alloc(M*fit_parameters);      /* Proposed update values */
  L2  = gsl_vector_alloc(M*fit_parameters);      /* Proposed update values */
  C  = gsl_matrix_alloc(M*fit_parameters,M*fit_parameters);
  W  = gsl_matrix_alloc(M*fit_parameters,M*fit_parameters);
  ata = gsl_matrix_alloc(M*fit_parameters,M*fit_parameters);
  atb = gsl_vector_alloc(M*fit_parameters);
  P = gsl_permutation_alloc(M*fit_parameters);

  SV_V = gsl_matrix_alloc(M*fit_parameters,M*fit_parameters);
  SV_S = gsl_vector_alloc(M*fit_parameters);
  SV_WORK = gsl_vector_alloc(M*fit_parameters);

  ZZ = gsl_matrix_alloc(N,M);
  dchi = gsl_vector_alloc(M);
  dchiN = gsl_vector_alloc(M);
  
  /* Save initial values in case we need to reset. */
  mi = malloc(sizeof(double) * M);
  if (mi == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  si = malloc(sizeof(double) * M);
  if (si == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  Ai = malloc(sizeof(double) * M);
  if (Ai == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  etai = malloc(sizeof(double) * M);
  if (etai == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  for (j = 0; j < M; j++) {
    mi[j] = m[j];
    si[j] = s[j];
    A[j]  = A[j] * pow(s[j] * sqrt(2 * M_PI),-1);  /* Convert from flux */
    Ai[j] = A[j];
    etai[j] = eta[j];
  }

  /* Prepare test values */
  m1 = malloc(sizeof(double) * M);
  if (m1 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  s1 = malloc(sizeof(double) * M);
  if (s1 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  A1 = malloc(sizeof(double) * M);
  if (A1 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  eta1 = malloc(sizeof(double) * M);
  if (eta1 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  m2 = malloc(sizeof(double) * M);
  if (m2 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  s2 = malloc(sizeof(double) * M);
  if (s2 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  A2 = malloc(sizeof(double) * M);
  if (A2 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  eta2 = malloc(sizeof(double) * M);
  if (eta2 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
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
    double chi0 = 0.0;
    double lambda0 = lambda;
    double chiP = 0.0;
    double nu = 10;
    double lambdaP = lambda0 / nu;
    double chi1,chi2;
    /* Initial chi position */
    eval_chi(X,Y,E,N,fit_function,fit_parameters,
	     m,s,A,eta,M,
	     dB,a,&chi0);

    /*  */
    /* Take a step using lambda0 */
    /* A^T A + lambda*diag(A^T A) */
    gsl_matrix_set_zero(ata);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,a,a,0.0,ata);
    for (j = 0; j < M; j++) {
      for (k = 0; k < M; k++) {
	if (j != k) {
	  gsl_matrix_set(ata,j,k,0.0);
	}
      }
    }
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,a,a,lambda0,ata);
    /* A^T B  */
    gsl_blas_dgemv(CblasTrans,1.0,a,dB,0.0,atb);
    /* dL = (A^T B) \ (A^T A) */
    old_handler = gsl_set_error_handler_off();
    j = 0;
    j = gsl_linalg_SV_decomp(ata,SV_V,SV_S,SV_WORK);
    if (j != 0) {
      iteration = 100000;
      continue;
    }
    j = gsl_linalg_SV_solve(ata,SV_V,SV_S,atb,L1);
    if (j != 0) {
      iteration = 100000;
      continue;
    }
    gsl_set_error_handler(old_handler);

    /*  */
    /* Take a step using lambdaP */
    /* A^T A + lambda*diag(A^T A) */
    gsl_matrix_set_zero(ata); /* Needed for ata to be set up correctly for dgemm  */
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,a,a,0.0,ata);
    for (j = 0; j < M; j++) {
      for (k = 0; k < M; k++) {
	if (j != k) {
	  gsl_matrix_set(ata,j,k,0.0);
	}
      }
    }

    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,a,a,lambdaP,ata);
    /* A^T B  */
    gsl_blas_dgemv(CblasTrans,1.0,a,dB,0.0,atb);
    /* dL = (A^T B) \ (A^T A)  */
    old_handler = gsl_set_error_handler_off();
    j = 0;
    j = gsl_linalg_SV_decomp(ata,SV_V,SV_S,SV_WORK);
    if (j != 0) {
      iteration = 100000;
      continue;
    }
    j = gsl_linalg_SV_solve(ata,SV_V,SV_S,atb,L2);
    if (j != 0) {
      iteration = 100000;
      continue;
    }
    gsl_set_error_handler(old_handler);
    
    for (j = 0; j < M; j++) {
      m1[j] = m[j] + vm[j] * gsl_vector_get(L1,0 + fit_parameters * j) * pow(0.5,lambda_counter > 2);
      s1[j] = s[j] + vs[j] * gsl_vector_get(L1,1 + fit_parameters * j)  * pow(0.5,lambda_counter > 2);
      A1[j] = A[j] + vA[j] * gsl_vector_get(L1,2 + fit_parameters * j)  * pow(0.5,lambda_counter > 2);
      eta1[j] = eta[j] + veta[j] * gsl_vector_get(L1,3 + fit_parameters * j)  * pow(0.5,lambda_counter > 2);
      if (eta1[j] < 0) { eta1[j] = 0.0; }
      
      m2[j] = m[j] + vm[j] * gsl_vector_get(L2,0 + fit_parameters * j)  * pow(0.5,lambda_counter > 2);
      s2[j] = s[j] + vs[j] * gsl_vector_get(L2,1 + fit_parameters * j)  * pow(0.5,lambda_counter > 2);
      A2[j] = A[j] + vA[j] * gsl_vector_get(L2,2 + fit_parameters * j)  * pow(0.5,lambda_counter > 2);
      eta2[j] = eta[j] + veta[j] * gsl_vector_get(L2,3 + fit_parameters * j)  * pow(0.5,lambda_counter > 2);
      if (eta2[j] < 0) { eta2[j] = 0.0; }
    }

    eval_chi(X,Y,E,N,fit_function,fit_parameters,
	     m1,s1,A1,eta1,M,
	     dB,a,&chi1);
    eval_chi(X,Y,E,N,fit_function,fit_parameters,
	     m2,s2,A2,eta2,M,
	     dB,a,&chi2);

    if (!isfinite(chi1)) { chi1 = chi0 * 100; }
    if (!isfinite(chi2)) { chi2 = chi0 * 100; }
    
    printf("%d\t%f %f %f (%g %g)\n",
    	   iteration,chi0,chi1,chi2,lambda0,lambdaP);
    for (j = 0; j < M; j++) {
      printf("    %f %f %f %f\n",m[j],s[j],A[j],eta[j]);
      printf("    %f %f %f %f\n",m1[j],s1[j],A1[j],eta1[j]);
      printf("    %f %f %f %f\n",m2[j],s2[j],A2[j],eta2[j]);
    }

    
      if ((chi1 > chi0)&&
	  (chi2 > chi0)) { /* Bad choice  */
	lambda = lambda0 * nu;
	val = 1000; /* Something big to denote non-convergence. */
	L = NULL;
	iteration--;
	lambda_counter++;
	*chi = chi0;
      }
      else {
	if (chi1 < chi2) { /* lambda0 case was better */
	  for (j = 0;j < M; j++) {
	    m[j] = m1[j];
	    s[j] = s1[j];
	    A[j] = A1[j];
	    eta[j] = eta1[j];
	  }
	  *chi = chi1;
	  L = L1;
	  gsl_blas_ddot(L1,L1,&val);
	}
	else {             /* lambdaP case was better, promote lambdaP */
	  lambda = lambdaP;
	  for (j = 0;j < M; j++) {
	    m[j] = m2[j];
	    s[j] = s2[j];
	    A[j] = A2[j];
	    eta[j] = eta2[j];
	  }
	}
	*chi = chi2;
	L = L2;
	gsl_blas_ddot(L2,L2,&val);
      }

    for (j = 0; j < M; j++) {
      /* Failsafe checks that we haven't gone crazy. */
      if (!isfinite( m[j])) {  m[j] = mi[j]; }
      if (!isfinite( s[j])) {  s[j] = si[j]; }
      if (!isfinite( A[j])) {  A[j] = Ai[j]; }
      if (!isfinite(eta[j])){ eta[j] = etai[j]; }
      if (!isfinite(dm[j])) { dm[j] = 0.0; }
      if (!isfinite(ds[j])) { ds[j] = 0.0; }
      if (!isfinite(dA[j])) { dA[j] = 0.0; }
      if (!isfinite(deta[j])) { deta[j] = 0.0; }
      
      if (fabs(m[j] - mi[j]) > fabs(mi[j])) { m[j] = mi[j]; }
      if (fabs(s[j]) > 10 * fabs(si[j])) { s[j] = si[j]; }
      if (fabs(A[j]) > 10 * fabs(Ai[j])) { A[j] = Ai[j]; }
      if ((fit_function != FUNCTION_HUMLICEK)&&
	  (fit_function != FUNCTION_SKEWGAUSS)) {
	if (eta[j] > 1.0) { eta[j] = 1.0; }
      }
      if (fit_function != FUNCTION_SKEWGAUSS) {
	if (eta[j] < 0.0) { eta[j] = 0.0; }
      }
      
      if (fabs(m[j] - mi[j]) > LARGE_SHIFT * si[j]) {
	m[j] = mi[j];
	vm[j] = 0;
      }
      if (s[j] < 0.0) { s[j] = si[j]; }
      if (s[j] < 0.0) { s[j] = tolerance; }
      
      /* Check to see if we look converged */
      if (sqrt(val) < tolerance) {
	is_converged = 1;
	
	if (fabs(gsl_vector_get(L,0 + fit_parameters * j) / m[j]) > tolerance) { is_converged = 0; }
	if (fabs(gsl_vector_get(L,1 + fit_parameters * j) / s[j]) > tolerance) { is_converged = 0; }
	if (fabs(gsl_vector_get(L,2 + fit_parameters * j) / A[j]) > tolerance) { is_converged = 0; }
      }
      else {
	is_converged = 0;
      }
    
      if ((!isfinite( m[j]))||(!isfinite( s[j]))||(!isfinite( A[j]))||
	  (!isfinite(dm[j]))||(!isfinite(ds[j]))||(!isfinite(dA[j]))) {
	iteration = 100000;
      }
    }
    if (lambda_counter > 10) { is_converged = 1; }
    
    iteration++;
  } /* End iteration loop  */

    /* Construct the covariance matrix */
  /* for (j = 0; j < M; j++) { */
  /*   gsl_matrix_set(W,0 + fit_parameters * j,0 + fit_parameters * j,pow(gsl_vector_get(SV_S,0 + fit_parameters * j),-2.0)); */
  /*   gsl_matrix_set(W,1 + fit_parameters * j,1 + fit_parameters * j,pow(gsl_vector_get(SV_S,1 + fit_parameters * j),-2.0)); */
  /*   gsl_matrix_set(W,2 + fit_parameters * j,2 + fit_parameters * j,pow(gsl_vector_get(SV_S,2 + fit_parameters * j),-2.0)); */
  /*   gsl_matrix_set(W,3 + fit_parameters * j,3 + fit_parameters * j,pow(gsl_vector_get(SV_S,3 + fit_parameters * j),-2.0)); */
  /* } */
  /* gsl_matrix_set_zero(ata); */
  /* gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,SV_V,W,0.0,ata); */
  /* gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,ata,SV_V,0.0,C); */
  
  /* // Calculate errors */
  /* dm[j] = sqrt(*chi * gsl_matrix_get(C,0 + fit_parameters * j,0 + fit_parameters * j)); */
  /* ds[j] = sqrt(*chi * gsl_matrix_get(C,1 + fit_parameters * j,1 + fit_parameters * j)); */
  /* dA[j] = sqrt(*chi * gsl_matrix_get(C,2 + fit_parameters * j,2 + fit_parameters * j)); */
  /* deta[j] = sqrt(*chi * gsl_matrix_get(C,3 + fit_parameters * j,3 + fit_parameters * j)); */
  
#define CHICONTOUR 0
#if CHICONTOUR

  if ((lambda_counter > 5)&&
      (fabs(m[0] - 4863) < 1.0)) {
    const gsl_rng_type * T;
    gsl_rng * r;
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    for (k = 0; k <= 10000; k++) {
      for (j = 0; j < M; j++) {
	m1[j] = m[j] + gsl_ran_flat(r,-1.0,1.0);
	m1[j] = m[j];
	s1[j] = s[j] * pow(10.0,gsl_ran_flat(r,-2,-1));
	/*	s1[j] = s[j]; */
	A1[j] = A[j] * gsl_ran_flat(r,0.5,2.01);

	eta1[j] = eta[j] + gsl_ran_flat(r,0,55);
	A1[j] = A[j] * eta1[j];
      }
      eval_chi(X,Y,E,N,fit_function,fit_parameters,
	       m1,s1,A1,eta1,M,
	       dB,a,&chi0);
      gsl_matrix_set_zero(ata); /* Needed for ata to be set up correctly for dgemm */
    gsl_blas_dgemv(CblasTrans,1.0,a,dB,0.0,atb);
    /*    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,a,dB,0.0,ata); */

      printf("CONTOUR: %f %f \t",
	     *chi,chi0);
      printf("DERIV:   %f %f %f %f\t",
	     gsl_vector_get(atb,0),
	     gsl_vector_get(atb,1),
	     gsl_vector_get(atb,2),
	     gsl_vector_get(atb,3));

      for (j = 0; j < M; j++) {
	printf(" %f %f %f %f ",m1[j],s1[j],A1[j] / pow(s[j] * sqrt(2 * M_PI),-1),eta1[j]);
      }
      printf("\n");
      
      /* if (chi0 < *chi) { */
      /* 	for (j = 0; j < M; j++) { */
      /* 	  m[j] = m[j]; */
      /* 	  s[j] = s1[j]; */
      /* 	  A[j] = A1[j]; */
      /* 	  eta[j] = eta1[j]; */
      /* 	} */
      /* } */
    }
  }

#endif


  /* Once the fit has converged:                                        */
  /* Free allocated memory.                                             */
  gsl_vector_free(dB);
  gsl_matrix_free(a);
  gsl_vector_free(L1);
  gsl_vector_free(L2);
  gsl_matrix_free(C);
  gsl_matrix_free(W);
  gsl_matrix_free(ata);
  gsl_vector_free(atb);
  gsl_permutation_free(P);

  gsl_matrix_free(SV_V);
  gsl_vector_free(SV_S);
  gsl_vector_free(SV_WORK);
  
  gsl_matrix_free(ZZ);
  gsl_vector_free(dchi);
  gsl_vector_free(dchiN);
  
  free(mi);
  free(si);
  free(Ai);
  free(etai);
  for (j = 0; j < M; j++) {
    A[j] =  A[j] / pow(s[j] * sqrt(2 * M_PI),-1);
    dA[j] = dA[j] / pow(s[j] * sqrt(2 * M_PI),-1);
  }
  
  return(iteration);
}



#endif  

#if (! USE_NEW_CODE) 
/* Linear least squares fitting code.                                    */
int multifunction(double *X, double *Y, double *E, int N,
		  int fit_function,
		  double *m, double *s, double *A, double *eta, int M,
		  double *dm, double *ds, double *dA, double *deta, 
		  int *vm, int *vs, int *vA, int *veta,
		  double relax,
		  double *chi, double tolerance, int max_iter) {
  int iteration = 0;
  int fit_parameters = 4;
  gsl_vector *dB;        /* Array of deviations */
  gsl_matrix *a ;        /* Matrix of derivatives at each position. */
  gsl_vector *L ;        /* Proposed update values */
  gsl_matrix *C ; 
  gsl_matrix *W ; 
  gsl_matrix *ata; 
  gsl_vector *atb; 
  gsl_permutation *P;

  gsl_matrix *SV_V;
  gsl_vector *SV_S;
  gsl_vector *SV_WORK;

  gsl_error_handler_t *old_handler;

  gsl_matrix *ZZ;
  gsl_vector *dchi;
  gsl_vector *dchiN;
  
  int is_converged = 0;
  double val = 0.0;
  /*  double F,dF,xms;  // Calculation speed tools.  */
  double F,d0,d1,d2,d3;
  int i; /* Data iterator */
  int j; /* Component iterator */
  int keep;
  double *mi, *si, *Ai, *etai;
  double chi_pixel = 0.0;
  double old_chi = 99e99;

  double lambda = 1.0;
  double lambda0 = lambda;

  /* Futz the tolerance */
  if (fit_function == FUNCTION_HUMLICEK) {
    tolerance /= 5.0;
  }
  
  /* Allocate vectors and matrices */
  dB = gsl_vector_alloc(N);        /* Array of deviations */
  a  = gsl_matrix_alloc(N,M*fit_parameters);    /* Matrix of derivatives at each position. */
  L  = gsl_vector_alloc(M*fit_parameters);      /* Proposed update values */
  C  = gsl_matrix_alloc(M*fit_parameters,M*fit_parameters);
  W  = gsl_matrix_alloc(M*fit_parameters,M*fit_parameters);
  ata = gsl_matrix_alloc(M*fit_parameters,M*fit_parameters);
  atb = gsl_vector_alloc(M*fit_parameters);
  P = gsl_permutation_alloc(M*fit_parameters);

  SV_V = gsl_matrix_alloc(M*fit_parameters,M*fit_parameters);
  SV_S = gsl_vector_alloc(M*fit_parameters);
  SV_WORK = gsl_vector_alloc(M*fit_parameters);

  ZZ = gsl_matrix_alloc(N,M);
  dchi = gsl_vector_alloc(M);
  dchiN = gsl_vector_alloc(M);
  
  /* Save initial values in case we need to reset. */
  mi = malloc(sizeof(double) * M);
  if (mi == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  si = malloc(sizeof(double) * M);
  if (si == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  Ai = malloc(sizeof(double) * M);
  if (Ai == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  etai = malloc(sizeof(double) * M);
  if (etai == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  for (j = 0; j < M; j++) {
    mi[j] = m[j];
    si[j] = s[j];
    A[j]  = A[j] * pow(s[j] * sqrt(2 * M_PI),-1);  /* Convert from flux */
    Ai[j] = A[j];
    etai[j] = eta[j];
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
    {
      /* Calculate current deviations, and the current derivative matrix. */
      gsl_vector_set_zero(dB);
      gsl_matrix_set_zero(a);
      gsl_matrix_set_zero(ata);
      gsl_vector_set_zero(atb);
      gsl_vector_set_zero(L);
      gsl_matrix_set_zero(W);
      gsl_matrix_set_zero(C);
      
      for (i = 0; i < N; i++) {
	val = 0.0;
	for (j = 0; j < M; j++) { /* Accumulate all components */
	  /* Set matrix values based on function and derivatives */
	  if (fit_function == FUNCTION_GAUSSIAN) {
	    gaussian_model(X[i],m[j],s[j],A[j],eta[j],
			   &F,&d0,&d1,&d2,&d3);
	    d3 = 0.0;
	  }
	  else if (fit_function == FUNCTION_LORENTZIAN) {
	    lorentzian_model(X[i],m[j],s[j],A[j],eta[j],
			     &F,&d0,&d1,&d2,&d3);
	    d3 = 0.0;
	  }
	  else if (fit_function == FUNCTION_PSEUDOVOIGT) {
	    pseudovoigt_model(X[i],m[j],s[j],A[j],eta[j],
			      &F,&d0,&d1,&d2,&d3);
	  }
	  else if (fit_function == FUNCTION_HJERTING) {
	    hjerting_model(X[i],m[j],s[j],A[j],eta[j],
			   &F,&d0,&d1,&d2,&d3);
	  }
	  else if (fit_function == FUNCTION_HUMLICEK) {
	    humlicek_model(X[i],m[j],s[j],A[j],eta[j],
			   &F,&d0,&d1,&d2,&d3);
	  }
	  else if (fit_function == FUNCTION_SKEWGAUSS) {
	    skewgauss_model(X[i],m[j],s[j],A[j],eta[j],
			    &F,&d0,&d1,&d2,&d3);
	  }
	  else if (fit_function == FUNCTION_PLANCK) {
	    planck_model(X[i],m[j],s[j],A[j],eta[j],
			 &F,&d0,&d1,&d2,&d3);
	  }
	  else {
	    F = 0.0;
	    d0 = 0.0;
	    d1 = 0.0;
	    d2 = 0.0;
	    d3 = 0.0;
	  }
	  gsl_matrix_set(a,i,0 + j * fit_parameters,d0 / E[i]);
	  gsl_matrix_set(a,i,1 + j * fit_parameters,d1 / E[i]);
	  gsl_matrix_set(a,i,2 + j * fit_parameters,d2 / E[i]);
	  gsl_matrix_set(a,i,3 + j * fit_parameters,d3 / E[i]);
	  
	  gsl_matrix_set(ZZ,i,j,F);
	  val += F;
	} /* End loop over components. */
	gsl_vector_set(dB,i,(Y[i] - val) / E[i]);
	
	/* Should we use this point in the calculations? */
	keep = 0;
	for (j = 0; j < M; j++) {
	  if (fabs((X[i] - m[j]) / (sqrt(2.0) * s[j])) < 7) {
	    keep = 1;
	  }
	}
	if (fit_function == FUNCTION_PLANCK) {
	  keep = 1;
	}
	if ((fit_function == FUNCTION_HUMLICEK)&&(keep == 0)) {
	  if (fabs((X[i] - m[j]) / (sqrt(2.0) * ((1 + eta[j]) * s[j]))) < 10) {
	    keep = 1;
	  }
	}
	
	/* Calculate chi */
	if (keep == 1) {
	  *chi += pow(gsl_vector_get(dB,i),2);
	  chi_pixel += 1.0;
	}
	else { /* If we're not within the 6 sigma limit, we don't want to include anything. */
	  gsl_vector_set(dB,i,0.0);
	  for (j = 0; j < M; j++) {
	    gsl_matrix_set(a,i,0+j*fit_parameters,0.0);
	    gsl_matrix_set(a,i,1+j*fit_parameters,0.0);
	    gsl_matrix_set(a,i,2+j*fit_parameters,0.0);
	    gsl_matrix_set(a,i,3+j*fit_parameters,0.0);
	  }
	}
      } /* End loop over pixels */
      
      if (chi_pixel > 0.0) {
	*chi /= chi_pixel;
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
      /* A^T A  */
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,a,a,0.0,ata);
      
      /* LM addition. */
      if ((fit_function != FUNCTION_PLANCK))
      	{
	  if (*chi > old_chi) {
	    lambda = 1.0;
	  }
	  else {
	    if ((old_chi - *chi)/old_chi < 0.05) {
	      lambda /= 2.0;
	    }
	    else {
	      lambda *= 2.0;
	    }
	  }
      	  /*      if (lambda > 10) { lambda = 10; } */
      	}
      
      for (i = 0; i < M * fit_parameters; i++) {
	for (j = 0; j < M * fit_parameters; j++) {
	  if (i != j) {
	    gsl_matrix_set(ata,i,j,0.0);
	  }
	}
      }
      /* A^T A + lambda*diag(A^T A) */
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,a,a,lambda0,ata);
      
      /* A^T B  */
      gsl_blas_dgemv(CblasTrans,1.0,a,dB,0.0,atb);
      
      /* dL = (A^T B) \ (A^T A) */
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

    }

    /* Construct the covariance matrix */
    for (j = 0; j < M; j++) {
      gsl_matrix_set(W,0 + fit_parameters * j,0 + fit_parameters * j,pow(gsl_vector_get(SV_S,0 + fit_parameters * j),-2.0));
      gsl_matrix_set(W,1 + fit_parameters * j,1 + fit_parameters * j,pow(gsl_vector_get(SV_S,1 + fit_parameters * j),-2.0));
      gsl_matrix_set(W,2 + fit_parameters * j,2 + fit_parameters * j,pow(gsl_vector_get(SV_S,2 + fit_parameters * j),-2.0));
      gsl_matrix_set(W,3 + fit_parameters * j,3 + fit_parameters * j,pow(gsl_vector_get(SV_S,3 + fit_parameters * j),-2.0));
    }
    gsl_matrix_set_zero(ata);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,SV_V,W,0.0,ata);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,ata,SV_V,0.0,C);
       
    /* Update the values with the proposed changes, and recalculate errors*/
    for (j = 0; j < M; j++) {
      /* Update parameters */
      if (vm[j]) {
	m[j] += gsl_vector_get(L,0 + fit_parameters * j);
      }
      if (vs[j]) {
	s[j] += gsl_vector_get(L,1 + fit_parameters * j);
      }
      if (vA[j]) { 
	A[j] += gsl_vector_get(L,2 + fit_parameters * j);
      }
      if (veta[j])  {
	eta[j] += gsl_vector_get(L,3 + fit_parameters * j);
      }

      /* Calculate errors */
      dm[j] = sqrt(*chi * gsl_matrix_get(C,0 + fit_parameters * j,0 + fit_parameters * j));
      ds[j] = sqrt(*chi * gsl_matrix_get(C,1 + fit_parameters * j,1 + fit_parameters * j));
      dA[j] = sqrt(*chi * gsl_matrix_get(C,2 + fit_parameters * j,2 + fit_parameters * j));
      deta[j] = sqrt(*chi * gsl_matrix_get(C,3 + fit_parameters * j,3 + fit_parameters * j));

      /* Failsafe checks that we haven't gone crazy. */
      if (!isfinite( m[j])) {  m[j] = mi[j]; }
      if (!isfinite( s[j])) {  s[j] = si[j]; }
      if (!isfinite( A[j])) {  A[j] = Ai[j]; }
      if (!isfinite(eta[j])){ eta[j] = etai[j]; }
      if (!isfinite(dm[j])) { dm[j] = 0.0; }
      if (!isfinite(ds[j])) { ds[j] = 0.0; }
      if (!isfinite(dA[j])) { dA[j] = 0.0; }
      if (!isfinite(deta[j])) { deta[j] = 0.0; }

      if (fabs(m[j] - mi[j]) > fabs(mi[j])) { m[j] = mi[j]; }
      if (fabs(s[j]) > 10 * fabs(si[j])) { s[j] = si[j]; }
      if (fabs(A[j]) > 10 * fabs(Ai[j])) { A[j] = Ai[j]; }
      if ((fit_function != FUNCTION_HUMLICEK)&&
	  (fit_function != FUNCTION_SKEWGAUSS)) {
	if (eta[j] > 1.0) { eta[j] = 1.0; }
      }
      if (fit_function != FUNCTION_SKEWGAUSS) {
	if (eta[j] < 0.0) { eta[j] = 0.0; }
      }
      
      if (fabs(m[j] - mi[j]) > LARGE_SHIFT * si[j]) {
      	m[j] = mi[j];
      	vm[j] = 0;
      }
      if (s[j] < 0.0) { s[j] = si[j]; }
      if (s[j] < 0.0) { s[j] = tolerance; }

      /* Check to see if we look converged */
      if (fabs(gsl_vector_get(L,0 + fit_parameters * j) / m[j]) > tolerance) { is_converged = 0; }
      if (fabs(gsl_vector_get(L,1 + fit_parameters * j) / s[j]) > tolerance) { is_converged = 0; }
      if (fabs(gsl_vector_get(L,2 + fit_parameters * j) / A[j]) > tolerance) { is_converged = 0; }

      if ((!isfinite( m[j]))||(!isfinite( s[j]))||(!isfinite( A[j]))||
	  (!isfinite(dm[j]))||(!isfinite(ds[j]))||(!isfinite(dA[j]))) {
	iteration = 100000;
      }
#define CHICONTOUR 0
#if CHICONTOUR
      gsl_blas_ddot(L,L,&val);
      fprintf(stderr,"%d %d %g (%g %g %g %g) (%g %g %g %g) %g %g %g %g %g\n",
	      iteration,j,mi[j],m[j],s[j],A[j],eta[j],
	      gsl_vector_get(L,0 + fit_parameters * j),
	      gsl_vector_get(L,1 + fit_parameters * j),
	      gsl_vector_get(L,2 + fit_parameters * j),
	      gsl_vector_get(L,3 + fit_parameters * j),
	      tolerance,lambda,sqrt(val),*chi,old_chi);
#endif
    }
    iteration++;
    /*    if (fabs((old_chi - *chi)/old_chi) < tolerance / lambda) { is_converged = 1; } */
    gsl_blas_ddot(L,L,&val);
    if (sqrt(val) < tolerance) { is_converged = 1; }
    /* fprintf(stderr,"%g %g %g %d\n", */
    /* 	    sqrt(val),sqrt(val) * lambda,tolerance,is_converged); */

  }
  /* Once the fit has converged:                                        */
  /* Free allocated memory.                                             */
  gsl_vector_free(dB);
  gsl_matrix_free(a);
  gsl_vector_free(L);
  gsl_matrix_free(C);
  gsl_matrix_free(W);
  gsl_matrix_free(ata);
  gsl_vector_free(atb);
  gsl_permutation_free(P);

  gsl_matrix_free(SV_V);
  gsl_vector_free(SV_S);
  gsl_vector_free(SV_WORK);
  
  gsl_matrix_free(ZZ);
  gsl_vector_free(dchi);
  gsl_vector_free(dchiN);
  
  free(mi);
  free(si);
  free(Ai);
  free(etai);
  for (j = 0; j < M; j++) {
    A[j] =  A[j] / pow(s[j] * sqrt(2 * M_PI),-1);
    dA[j] = dA[j] / pow(s[j] * sqrt(2 * M_PI),-1);
  }
  
  return(iteration);
}
#endif

#endif

#if USE_GSL_SOLVER

struct fdata {
  int N;
  double *X;
  double *Y;
  double *E;
};

int ff(const gsl_vector *x, void *data, gsl_vector *f) {
  int i;
  double m = gsl_vector_get(x,0);
  double s = gsl_vector_get(x,1);
  double A = gsl_vector_get(x,2);
  double eta=gsl_vector_get(x,3);
  double F,d0,d1,d2,d3;

  int N = ((struct fdata *) data)->N;
  double *X = ((struct fdata *) data)->X;
  double *Y = ((struct fdata *) data)->Y;
  double *E = ((struct fdata *) data)->E;

  for (i = 0; i < N; i++) {  
    humlicek_model(X[i],m,s,A,eta,
		   &F,&d0,&d1,&d2,&d3);
    gsl_vector_set(f,i,(F - Y[i])/E[i]);
  }
  return(GSL_SUCCESS);
}
int df(const gsl_vector *x, void *data, gsl_matrix *J) {
  int i;
  double m = gsl_vector_get(x,0);
  double s = gsl_vector_get(x,1);
  double A = gsl_vector_get(x,2);
  double eta=gsl_vector_get(x,3);
  double F,d0,d1,d2,d3;

  int N = ((struct fdata *) data)->N;
  double *X = ((struct fdata *) data)->X;
  double *Y = ((struct fdata *) data)->Y;
  double *E = ((struct fdata *) data)->E;

  for (i = 0; i < N; i++) {  
    humlicek_model(X[i],m,s,A,eta,
		   &F,&d0,&d1,&d2,&d3);
    gsl_matrix_set(J,i,0,d0/E[i]);
    gsl_matrix_set(J,i,1,d1/E[i]);
    gsl_matrix_set(J,i,2,d2/E[i]);
    gsl_matrix_set(J,i,3,d3/E[i]);
  }
  return(GSL_SUCCESS);
}
int fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J) {
  ff(x,data,f);
  df(x,data,J);
  return(GSL_SUCCESS);
}

int multifunction(double *X, double *Y, double *E, int N,
		  int fit_function,
		  double *m, double *s, double *A, double *eta, int M,
		  double *dm, double *ds, double *dA, double *deta, 
		  int *vm, int *vs, int *vA, int *veta,
		  double relax,
		  double *chi, double tolerance, int max_iter) {
  int i,j,iter;
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *solver;
  int status;

  /*  gsl_matrix *covar = gsl_matrix_alloc(4,4); */
  struct fdata data;
  data.N = N;
  data.X = X;
  data.Y = Y;
  data.E = E;
  gsl_multifit_function_fdf f;

  f.f = &ff;
  f.df = &df;
  f.fdf = &fdf;
  f.n = N;
  f.p = 4;
  f.params = &data;

  T = gsl_multifit_fdfsolver_lmsder;
  solver = gsl_multifit_fdfsolver_alloc(T,N,4);

  for (j = 0; j < M; j++) {
    double x_init[4] = {m[j], s[j], A[j] * pow(s[j] * sqrt(2 * M_PI),-1) , eta[j]};
    gsl_vector_view x = gsl_vector_view_array(x_init, 4);

    gsl_multifit_fdfsolver_set(solver,&f,&x.vector);
    iter = 0;
    do {
      iter++;
      status = gsl_multifit_fdfsolver_iterate(solver);

      if (status) { break; }
      status = gsl_multifit_test_delta(solver->dx,solver->x,tolerance,tolerance);
    }
    while (status == GSL_CONTINUE && iter < max_iter);

    m[j] = gsl_vector_get(solver->x,0);
    s[j] = gsl_vector_get(solver->x,1);
    A[j] = gsl_vector_get(solver->x,2) / pow(s[j] * sqrt(2 * M_PI),-1);
    eta[j] = gsl_vector_get(solver->x,3);
    
  }
  gsl_multifit_fdfsolver_free(solver);
  
  return(0);
}

#endif
