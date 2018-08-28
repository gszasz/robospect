/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

double skewness_from_eta(double eta) {
  double delta = eta / sqrt(1 + pow(eta,2));
  double gamma = (4.0 - M_PI)/2.0 * pow(delta * sqrt(2 / M_PI),3) / pow(1 - (2 / M_PI) * pow(delta,2),1.5);
  return(gamma);
}

double skewgaussian(double x, double m, double s, double A, double eta) {
  return(gaussian(x,m,s,A) * (1.0 + erf(eta * (x - m)/s * M_SQRT1_2)));
}


void skewgauss_model(double x, double m, double s, double A, double eta,
		     double *f,
		     double *dfdm, double *dfds, double *dfdA, double *dfdeta) {
  double xms = (x - m)/s;
  double G   = exp(-0.5 * pow(xms,2));
  
  *dfdA = G * (1.0 * erf(eta * xms * M_SQRT1_2));
  *f    = *dfdA * A;

  *dfdeta = A * 1/sqrt(0.5 * M_PI) * xms * G * exp(-0.5 * pow(eta * xms,2));
  *dfdm = *f * xms / s        - *dfdeta * (eta / (x - m));
  *dfds = *f * pow(xms,2) / s - *dfdeta * (eta / s);

}

double skewgauss_line(double x, lines *L, int i) {
  return(skewgaussian(x,L->m[i],L->s[i],L->F[i],L->eta[i]) *
	 pow(L->s[i] * sqrt(2 * M_PI),-1));
}

