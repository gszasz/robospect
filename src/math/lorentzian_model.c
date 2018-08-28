/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

double lorentzian(double x, double m, double s, double A) {
  return(A * s / (2 * M_PI * ( pow(x - m,2) + pow(s / 2.0,2))));
}

/* double dLdA(double x, double m, double s, double A) { */
/*   return(s / (2 * M_PI * ( pow(x - m, 2) + pow(s / 2.0,2)))); */
/* } */

/* double dLdm(double x, double m, double s, double A) { */
/*   return(A * s * (x - m) / (M_PI * ( pow(x - m, 2) + pow(s / 2.0,2)))); */
/* } */

/* double dLds(double x, double m, double s, double A) { */
/*   return(A / (2 * M_PI) * ( pow( pow(x - m, 2) + pow(s / 2.0,2),-1) - */
/* 			    0.5 * pow(s,2) * pow( pow(x - m, 2) + pow(s / 2.0,2),-2) )); */
/* } */

void lorentzian_model(double x, double m, double s, double A, double eta,
		      double *f,
		      double *dfdm, double *dfds, double *dfdA, double *dfdeta) {
  double xms = pow( (pow(x - m,2) + pow(s / 2,2)), -1);

  *dfdA = (s / (2 * M_PI)) * xms;
  *f    = A * *dfdA;
  *dfdm = A * (s / M_PI) * xms * (x - m);
  *dfds = (A / (2 * M_PI)) * (xms - 0.5 * pow(s * xms,2));
  *dfdeta = 0.0;
}

double lorentzian_line(double x, lines *L, int i) {
  return(lorentzian(x,L->m[i],L->s[i],L->F[i]) * pow(L->s[i] * sqrt(2 * M_PI),-1));
}
