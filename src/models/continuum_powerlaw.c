/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void generate_model_continuum_powerlaw(opts *options, data *D, model *M) {
  int i;

  double E;
  double S = 0,Sx = 0,Sy = 0;
  double Sxx = 0,Syy = 0,Sxy = 0;

  double DD,A,B;
  
  for (i = 0; i <= D->N; i++) {
    E = pow(D->e[i],-2);
    S += E;
    Sx += D->x[i] * E;
    Sy += log10(D->y[i]) * E;
    Sxx += pow(D->x[i],2) * E;
    Syy += pow(log10(D->y[i]),2) * E;
    Sxy += D->x[i] * log10(D->y[i]) * E;
  }

  DD = (S * Sxx - Sx * Sx);
  if (DD == 0) {
    A = 0.0;
    B = 0.0;
  }
  else {
    A = (Sy * Sxx - Sx * Sxy) / DD;
    B = (S  * Sxy - Sx * Sy)  / DD;
  }

  /* Do real fit */
  log_comment(options,ROBO_VERBOSE_CONTINUUM,
	      "gmcpl: (%d/%d) (%g %g %g)",
	      i,D->N,
	      A,B,DD
	      );

  for (i = 0; i < D->N; i++) {
    M->continuum[i] = pow(10,A + B * D->x[i]);
  }
  
  log_comment(options,ROBO_VERBOSE_CONTINUUM,
	      "gmcpl: (%d/%d)  (%g %g %g)",
	      i,D->N,
	      A,B,DD
	      );
}
