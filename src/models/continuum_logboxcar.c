/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void generate_model_continuum_logboxcar(opts *options, data *D, model *M) {
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
      if ((D->x[j] > D->x[i] - options->continuum_box)&&(D->x[j] < D->x[i])&&
	  (D->y[j] > 0)) {
	v1[n1] = log10(D->y[j]);
	n1++;
      }
      else if ((D->x[j] > D->x[i])&&(D->x[j] < D->x[i] + options->continuum_box)&&
	       (D->y[j] > 0)) {
	v1[n1] = log10(D->y[j]);
	n1++;
      }
      if (D->x[j] > D->x[i] + options->continuum_box) {
	j = D->N + 10;
      }
    }

    if (n1 > 1) {
      S1 = array_stats(v1,n1);
      
      M->continuum[i] = pow(10,S1->med);
      if (options->supplied_errors == 0) {
	D->e[i]         = (1.4826 * S1->MAD) * log(10);
	D->e[i] *= M->continuum[i];
	if (D->e[i] <= 0.0) {
	  D->e[i] = D->e[i-1];
	}
      }
      free(S1);
    }
    else {
      M->continuum[i] = M->continuum[i-1];
      if (options->supplied_errors == 0) {
	D->e[i] = D->e[i-1];
      }
    }
    log_comment(options,ROBO_VERBOSE_CONTINUUM,
		"gmclb: (%d/%d) (%f %g %g) (%g %g %ld) %f",
		i,D->N,
		D->x[i],D->y[i],D->e[i],
		M->continuum[i],D->e[i],n1,
		M->continuum[i]
		);


  }
  
  free(v1);
}
