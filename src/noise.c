/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "robo.h"

void generate_noise_model_null(opts *options, data *D, model *M) {
  int i;
  log_comment(options,ROBO_VERBOSE_NOISE,
	      "gnmn: NOOP");
  for (i = 0; i < D->N; i++) {
    if (D->e[i] <= 0.0 ) {
      D->e[i] = D->e[i-1];
    }
  }
}

void generate_noise_model_boxcar(opts *options, data *D, model *M) {
  int i,j;
  stats *S1;
  double *v1 = malloc(sizeof(double) * D->N);
  int n1;
  
  for (i = 0; i < D->N; i += 1) {
    n1 = 0;
    for (j = 0; j < D->N; j++) {
      if ((D->x[j] > D->x[i] - options->continuum_box)&&(D->x[j] < D->x[i])) {
	v1[n1] = (D->y[j]);
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
    S1 = array_stats(v1,n1);

    D->e[i] = 1.4826 * S1->MAD;
    if (D->e[i] <= 0.0) {
      D->e[i] = D->e[i-1];
    }
    free(S1);
    log_comment(options,ROBO_VERBOSE_NOISE,
		"gnmb: %g %g %g",
		D->x[i],D->y[i],D->e[i]);
  }
  free(v1);
}

void generate_noise_model_Poisson(opts *options, data *D, model *M) {
  int i;
  double Crest;

  for (i = 0; i < D->N; i += 1) {
    Crest = D->y[i] + M->continuum[i];
    if (Crest >= 0.0) {
      D->e[i] = sqrt(Crest);
    }
    else {
      D->e[i] = D->e[i-1];
    }
    log_comment(options,ROBO_VERBOSE_NOISE,
		"gnmP: %g %g %g",
		D->x[i],D->y[i],D->e[i]);
  }
}
    

void generate_noise_model(opts *options, data *D, model *M) {
  switch (M->noise_model) {
  case NOISE_BOXCAR:
  case NOISE_SLOW_BOXCAR:
    generate_noise_model_boxcar(options,D,M);
    break;
  case NOISE_POISSON:
    generate_noise_model_Poisson(options,D,M);
    break;
  case NOISE_NULL:
  case NOISE_SUPPLIED:
  default:
    generate_noise_model_null(options,D,M);
    break;
  }
}
