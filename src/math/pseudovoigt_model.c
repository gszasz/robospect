/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

double pseudovoigt(double x, double m, double s, double A, double eta) {
  return(A * (eta * s / (2 * M_PI * ( pow(x - m,2) + pow(s / 2,2))) +
	      (1 - eta) * exp(-0.5 * pow( (x - m)/s,2))));
}

/* double dPVdA(double x, double m, double s, double A, double eta) { */
/*   return((eta * s / (2 * M_PI * ( pow(x - m,2) + pow(s / 2,2))) + */
/* 	  (1 - eta) * exp(-0.5 * pow( (x - m)/s,2)))); */
/* } */

/* double dPVdm(double x, double m, double s, double A, double eta) { */
/*   return(A * (eta *  s * (x - m) / (M_PI * ( pow(x - m, 2) + pow(s / 2.0,2))) + */
/* 	      (1 - eta) * (x - m) / pow(s,2) * exp(-0.5 * pow( (x - m) / s,2)))); */
/* } */

/* double dPVds(double x, double m, double s, double A, double eta) { */
/*   return(A * (eta / (2 * M_PI) * ( pow( pow(x - m, 2) + pow(s / 2.0,2),-1) - */
/* 				   0.5 * pow(s,2) * pow( pow(x - m, 2) + pow(s / 2.0,2),-2) ) + */
/* 	      (1 - eta) * pow(x - m,2) * pow(s,-3) * exp(-0.5 * pow( (x - m) / s,2)))); */
/* } */

/* double dPdeta(double x, double m, double s, double A, double eta) { */
/*   return(A * (s / (2 * M_PI * ( pow(x - m,2) + pow(s / 2,2))) - */
/* 	      exp(-0.5 * pow( (x - m) / s, 2)))); */
/* } */

void pseudovoigt_model(double x, double m, double s, double A, double eta,
			 double *f,
			 double *dfdm, double *dfds, double *dfdA, double *dfdeta) {
  double G,dGdm,dGds,dGdA,dGdeta;
  double L,dLdm,dLds,dLdA,dLdeta;

  gaussian_model(x,m,s,A,eta,&G,&dGdm,&dGds,&dGdA,&dGdeta);
  lorentzian_model(x,m,s,A,eta,&L,&dLdm,&dLds,&dLdA,&dLdeta);

  *f    = eta * L + (1 - eta) * G;
  *dfdm = eta * dLdm + (1 - eta) * dGdm;
  *dfds = eta * dLds + (1 - eta) * dGds;
  *dfdA = eta * dLds + (1 - eta) * dGdA;
  *dfdeta= L - G;
}
  
double pseudovoigt_line(double x, lines *L, int i) {
  return(pseudovoigt(x,
		     L->m[i],
		     L->s[i],
		     L->F[i] * pow(L->s[i] * sqrt(2 * M_PI),-1),
		     L->eta[i]));
}
