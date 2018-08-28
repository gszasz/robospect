/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void generate_model_continuum_boxcar(opts *options, data *D, model *M) {
  int i,j;
  stats *S1;
  double *v1 = malloc(sizeof(double) * D->N);
  int n1;

  for (i = 0; i < D->N; i += 1) {
    n1 = 0;
    for (j = 0; j < D->N; j++) {
      if ((D->x[j] > D->x[i] - options->continuum_box)&&(D->x[j] < D->x[i])) {
	v1[n1] = D->y[j];
	n1++;
      }
      else if ((D->x[j] > D->x[i])&&(D->x[j] < D->x[i] + options->continuum_box)) {
	v1[n1] = D->y[j];
	n1++;
      }
      if (D->x[j] > D->x[i] + options->continuum_box) {
	j = D->N + 10;
      }
    }
    S1 = array_stats(v1,n1);

    M->continuum[i] = S1->med;
    if (options->supplied_errors == 0) {
      if (options->flux_calibrated == 0) {
	D->e[i]         = 1.4826 * S1->MAD / S1->med;
      }
      else {
	D->e[i]         = 1.4826 * S1->MAD;
      }
      if (D->e[i] <= 0.0) {
	D->e[i] = D->e[i-1];
      }
    }
    log_comment(options,ROBO_VERBOSE_CONTINUUM,
		"gmcb: (%d/%d) %g (%f %g %g) (%g %g %ld) %g",
		i,D->N,
		options->continuum_box,
		D->x[i],D->y[i],D->e[i],
		S1->med,1.4826 * S1->MAD,S1->N,
		M->continuum[i]
		);
    free(S1);
  }
  
  free(v1);
}
