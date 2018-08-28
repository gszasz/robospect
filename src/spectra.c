/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "robo.h"

double *continuum_normalized_spectra(data *D, model *M) {
  double *CNS = vector_divide(D->yO,M->continuum,D->N);
  return(CNS);
}
double *line_removed_spectra(data *D, model *M) {
  double *LC  = vector_multiply(M->lines,M->continuum,D->N);
  double *LRS = vector_subtract(D->yO,LC,D->N);
  free(LC);
  return(LRS);
}
double *continuum_normalized_line_removed_spectra(data *D, model *M) {
  double *CNS = continuum_normalized_spectra(D,M);
  double *CNLRS = vector_subtract(CNS,M->lines,D->N);
  free(CNS);
  return(CNLRS);
}
double *continuum_normalized_Zvalue(data *D, model *M) {
  double *CNS = vector_divide(D->yO,M->continuum,D->N);
  double *CNZ = vector_subtract_constant(CNS,1.0,D->N);
  double *Z   = vector_divide(CNZ,D->e,D->N);

  free(CNS);
  free(CNZ);
  return(Z);
}

double *continuum_subtracted_spectra(data *D, model *M) {
  double *CNS = vector_subtract(D->yO,M->continuum,D->N);
  return(CNS);
}
double *line_removed_unnormalized_spectra(data *D, model *M) {
  double *LRS = vector_subtract(D->yO,M->lines,D->N);
  return(LRS);
}
double *continuum_subtracted_line_removed_spectra(data *D, model *M) {
  double *CNS = continuum_subtracted_spectra(D,M);
  double *CNLRS = vector_subtract(CNS,M->lines,D->N);
  free(CNS);
  return(CNLRS);
}
double *continuum_subtracted_Zvalue(data *D, model *M) {
  double *CNZ = continuum_subtracted_spectra(D,M);
  double *Z   = vector_divide(CNZ,D->e,D->N);

  free(CNZ);
  return(Z);
}
  

void set_line_fit_spectra(data *D, model *M) {
  double *CNS;

  if (D->y) {
    free(D->y);
  }

  CNS = continuum_normalized_spectra(D,M);
  D->y = vector_subtract_constant(CNS,1.0,D->N);
  free(CNS);
}
void set_continuum_fit_spectra(data *D, model *M) {
  if (D->y) {
    free(D->y);
  }
  if (D->O->flux_calibrated == 0) {
    D->y = line_removed_spectra(D,M);
  }
  else {
    D->y = line_removed_unnormalized_spectra(D,M);
  }
}
void set_noise_fit_spectra(data *D, model *M) {
  double *CNLRS;
  if (D->y) {
    free(D->y);
  }
  if (D->O->flux_calibrated == 0) {
    CNLRS = continuum_normalized_line_removed_spectra(D,M);
    D->y = vector_subtract_constant(CNLRS,1.0,D->N);
  }
  else {
    CNLRS = continuum_subtracted_line_removed_spectra(D,M);
    D->y = vector_subtract_constant(CNLRS,0.0,D->N);
  }
  free(CNLRS);
}
void set_linefinder_Z_spectra(data *D, model *M) {
  if (D->y) {
    free(D->y);
  }
  if (D->O->flux_calibrated == 0) {
    D->y = continuum_normalized_Zvalue(D,M);
  }
  else {
    D->y = continuum_subtracted_Zvalue(D,M);
  }
}
