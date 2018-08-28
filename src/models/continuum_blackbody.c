/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

double planck(double lambda, double tau, double A) {
  return(A * pow(lambda,-5) * pow(exp(pow(tau * lambda,-1)) - 1,-1));
}

double dPdA(double lambda, double tau, double A) {
  return(planck(lambda,tau,A) / A);
}
double dPdT(double lambda, double tau, double A) {
  double E = exp(pow(tau * lambda,-1));
  return(planck(lambda,tau,A) * pow(lambda,-1) * pow(tau,-2) * E / (E - 1));
}


#define HCKB 143877506.0
void generate_model_continuum_blackbody(opts *options, data *D, model *M) {
  int i;

  double A,tau,lT;
  double *EE = malloc(sizeof(double) * D->N);
  double l0,l1;
  double f0,f1;
  double K;
  double Z;
  double minZ = 99e99;
  double T;

  double s,eta;
  double dtau,ds,dA,deta;
  int vtau,vs,vA,veta;
  int iterations = 0;
  double chi;
  /* Calculate initial guesses. */
  l0 = D->x[0];
  l1 = D->x[D->N - 1];
  f0 = D->y[0];
  f1 = D->y[D->N - 1];

  K = (f0 * pow(l0,5)) / (f1 * pow(l1,5));
  T = 5778;
  for (lT = 2; lT <= 8; lT += 0.1) {
    tau = pow(lT,10) / HCKB;
    Z = fabs(K -
	     (exp(pow(l1 * tau,-1)) - 1) /
	     (exp(pow(l0 * tau,-1)) - 1));
    if (Z < minZ) {
      minZ = Z;
      T = tau;
    }
    else {
      if (minZ < 99e99) {
	continue;
      }
    }
  }
  tau = T;
  A = 0.5 * ((f0 * pow(l0,5) * (exp(pow(l0 * tau,-1)) - 1)) +
	     (f1 * pow(l1,5) * (exp(pow(l1 * tau,-1)) - 1)));
  s = 1 / sqrt(2.0 * M_PI);
  eta = 0.0;
  dtau = 0.0;
  ds = 0.0;
  dA = 0.0;
  deta = 0.0;
  vtau = 1;
  vs = 0;
  vA = 1;
  veta = 0;

  for (i = 0; i < D->N; i++) {
    EE[i] = 0.1 * D->y[i];
  }
  /* Do real fit */
  log_comment(options,ROBO_VERBOSE_CONTINUUM,
	      "gmcbb: (%d/%d) %g (%g %g %g) (PREFIT)",
	      i,D->N,
	      tau * HCKB,s,A,eta,
	      chi,iterations
	      );

  iterations = multifunction(D->x,D->y,EE,D->N,
			     FUNCTION_PLANCK,
			     &tau,&s,&A,&eta,1,
			     &dtau,&ds,&dA,&deta,
			     &vtau,&vs,&vA,&veta,
			     D->O->relax,
			     &chi,
			     1e-50,100);

  free(EE);
  for (i = 0; i < D->N; i++) {
    M->continuum[i] = planck(D->x[i],tau,A);
  }
  
  log_comment(options,ROBO_VERBOSE_CONTINUUM,
	      "gmcbb: (%d/%d) %g (%g %g %g) (%g  %ld) ",
	      i,D->N,
	      tau * HCKB,s,A,eta,
	      chi,iterations
	      );
}
