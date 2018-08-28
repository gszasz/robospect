/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void planck_model(double x, double m, double s, double A, double eta,
		  double *f,
		  double *dfdm, double *dfds, double *dfdA, double *dfdeta) {
  double E = exp(pow(m * x,-1));
  *dfdA = pow(x,-5) * pow(E - 1,-1);
  *f    = A * *dfdA;
  *dfdm = *f * pow(x,-1) * pow(m,-2) * E / (E - 1);
  *dfds = 0.0;
  *dfdeta = 0.0;
}
