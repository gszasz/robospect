/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

/* double gaussian(double x, double m, double s, double A) { */
/*   return(A * exp(-0.5 * pow( (x - m)/s,2))); */
/* } */

/* double dGdA(double x, double m, double s, double A, double eta) { */
/*   return(exp(-0.5 * pow( (x - m)/s,2))); */
/* } */

/* double dGdm(double x, double m, double s, double A, double eta) { */
/*   return(A * (x - m) / pow(s,2) * exp(-0.5 * pow( (x - m) / s,2))); */
/* } */

/* double dGds(double x, double m, double s, double A, double eta) { */
/*   return(A * pow(x - m,2) * pow(s,-3) * exp(-0.5 * pow( (x - m) / s,2))); */
/* } */

void gaussian_model(double x, double m, double s, double A, double eta,
		      double *f,
		      double *dfdm, double *dfds, double *dfdA, double *dfdeta) {
  double xms = (x - m)/s;

  *dfdA = exp(-0.5 * pow(xms,2));
  *f    = A * *dfdA;
  *dfdm = *f * xms / s;
  *dfds = *dfdm * xms;
  *dfdeta = 0.0;
}

double gaussian_line(double x, lines *L, int i) {
  return(gaussian(x,L->m[i],L->s[i],L->F[i]) * pow(L->s[i] * sqrt(2 * M_PI),-1));
}

double gaussian_line_alt(double x, lines *L, int i) {
  return(gaussian(x,L->mp[i],L->sp[i],L->Fp[i]) * pow(L->sp[i] * sqrt(2 * M_PI),-1));
}
