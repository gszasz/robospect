/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void generate_model_continuum_histogram(opts *options, data *D, model *M) {

  int i,j;
  stats *S1;
  double *v1 = malloc(sizeof(double) * D->N);
  if (v1 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  int n1;
  for (i = 0; i < D->N; i += 1) {
    n1 = 0;
    for (j = 0; j < D->N; j++) {
      if ((D->x[j] > D->x[i] - options->continuum_box)&&(D->x[j] < D->x[i])) {
	v1[n1] = (D->y[j] );
	n1++;
      }
      else if ((D->x[j] > D->x[i])&&(D->x[j] < D->x[i] + options->continuum_box)) {
	v1[n1] = (D->y[j]);
	n1++;
      }
      if (D->x[j] > D->x[i] + options->continuum_box) {
	j = D->N + 10;
      }
    }
    S1 = histogram_and_gaussfit(v1,n1);
    M->continuum[i] = S1->mean;
    if (options->supplied_errors == 0) {
      D->e[i] = S1->sigma / S1->mean;
    }
    log_comment(options,ROBO_VERBOSE_CONTINUUM,
		"gmch: (%d/%d) (%f %f %f) (%f %f %ld) %f",
		i,D->N,
		D->x[i],D->y[i],D->e[i],
		S1->mean,S1->sigma,S1->N,
		M->continuum[i]
		);
    free(S1);
  }
  free(v1);
}
