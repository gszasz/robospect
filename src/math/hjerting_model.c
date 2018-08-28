/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"
#include <gsl/gsl_sf_dawson.h>



/* http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?bibcode=1948ApJ...108..112H&db_key=AST&page_ind=1&data_type=GIF&type=SCREEN_VIEW&classic=YES */
#define TWOTHIRDS 0.666666666666666667
inline double hjerting_simple(double x, double m, double s, double A, double eta) {
  double u = ((x - m) / (sqrt(2.0) * s));
  double u2 = pow(u,2);
  double F = gsl_sf_dawson(u);  /* dF/du = 1 - 2 * u F */
  double dFdu = 1 - 2 * u * F;
  double k = -2 / sqrt(M_PI);
  double H0 = exp(-1.0 * u2);
  double H1 = k * (dFdu);
  double H2 = H0 * (1 - 2 * u2);
  double H3 = k * ( (TWOTHIRDS) * (1 - u2) - 2 * u * (1 - (TWOTHIRDS) * u2) * F);
  double H4 = H0 * (0.5 - 2 * u2 + (TWOTHIRDS) * pow(u2,2));
  /*  double norm = 1.0 / (1 - (s * eta)/111); */
  double norm = 1.0;
  if (!isfinite(u)) { return(0.0); }
  return(A * (H0 + eta * H1 + pow(eta,2) * H2 + pow(eta,3) * H3 + pow(eta,4) * H4) * norm);
}
#define SMALL 1e-5
void hjerting_model(double x, double m, double s, double A, double eta,
			 double *f,
			 double *dfdm, double *dfds, double *dfdA, double *dfdeta) {
  double u = ((x - m) / (sqrt(2.0) * s));
  double dudm = -1.0           / (sqrt(2.0) * s);
  double duds = -1.0 * (x - m) / (sqrt(2.0) * pow(s,2));
  double u2 = pow(u,2);
  double F = gsl_sf_dawson(u);  /* dF/du = 1 - 2 * u F */
  double dFdu = 1 - 2 * u * F;
  double d2Fdu2 = -2 * ( F + u * dFdu);
  double k = -2 / sqrt(M_PI);
  double H0 = exp(-1.0 * u2);
  double H1 = k * (dFdu);
  double H2 = H0 * (1 - 2 * u2);
  double H3 = k * ( (TWOTHIRDS) * (1 - u2) - 2 * u * (1 - (TWOTHIRDS) * u2) * F);
  double H4 = H0 * (0.5 - 2 * u2 + (TWOTHIRDS) * pow(u2,2));
  
  double dH0du = -2 * u * H0;
  double dH1du = k * d2Fdu2;
  double dH2du = H0 * (-6 * u + 4 * pow(u,3));
  /*    double dH3du = k * ((dFdu + 1) * (-2.0 * TWOTHIRDS * u) - (d2Fdu2) * (1 - TWOTHIRDS * u2)); */
  double dH3du = k * ((-2.0 * TWOTHIRDS * u) - (
						(2 * u - (TWOTHIRDS) * u2 *u) * dFdu +
						(2 - 2.0 * u2) * F));
  double dH4du = H0 * (-5 * u + (20.0/3.0) * pow(u,3) - 2.0 * TWOTHIRDS * pow(u,5));

  if (!isfinite(u)) {
    *f = 0.0;
    *dfdm = 0.0;
    *dfds = 0.0;
    *dfdA = 0.0;
    *dfdeta = 0.0;
    return;
  }
  if (eta > 1.0) { eta = 1.0; }
  
  *dfdA = H0 + eta * H1 + pow(eta,2) * H2 + pow(eta,3) * H3 + pow(eta,4) * H4;
  *f    = A * *dfdA;
  *dfdm = A * (dH0du + eta * dH1du + pow(eta,2) * dH2du + pow(eta,3) * dH3du + pow(eta,4) * dH4du) * dudm;
  /*    *dfdm = A * dH0du * dudm; */
  *dfdeta = A * (H1 + 2 * eta * H2 + 3 * pow(eta,2) * H3 + 4 * pow(eta,3) * H4) ;
  *dfds = *dfdm * duds / dudm - *dfdeta * eta / s;

}

double hjerting_line(double x, lines *L, int i) {
  return(hjerting_simple(x,L->m[i],L->s[i],L->F[i] * pow(L->s[i] * sqrt(2 * M_PI),-1),L->eta[i]));
}
