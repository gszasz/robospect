/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void generate_model_continuum_robust_linear(opts *options, data *D, model *M) {
  double m = 0.4 / 1500;
  double b = -0.3;
  double *x; 
  double *y; 

  int *count = malloc(sizeof(int) * D->N);
  if (count == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  int i,j,k = 0,z = 0;
  
  for (i = 0; i < D->N; i++) {
    M->continuum[i] = 0.0;
  }
  
  for (i = 0; i < D->N; i++) {
    x = malloc(sizeof(double) * (options->continuum_box + 1));
    if (x == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    y = malloc(sizeof(double) * (options->continuum_box + 1));
    if (y == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    k = 0;
    for (j = -1 * (int) options->continuum_box / 2; j <= (int) options->continuum_box / 2; j++) {
      z = i + j;
      if ((z >= 0)&&(z < D->N)) {
	x[k] = D->x[z];
	y[k] = D->y[z];
	k++;
      }
    }
    robust_linear_fit(x,y,k,&m,&b,1e-6);

    k = 0;

    for (j = -1 * (int) options->continuum_box / 2; j <= (int) options->continuum_box / 2; j++) {
      z = i + j;
      if ((z >= 0)&&(z < D->N)) {
	M->continuum[z] += (m * x[k] + b);
	count[z] ++;
	k++;
      }

    }
    free(x);
    free(y);
  }

  for (i = 0; i < D->N; i++) {
    if (count[i] != 0) {
      M->continuum[i] /= count[i];
      if (M->continuum[i] < 0) {
	M->continuum[i] = 1.0;
      }
    }
  }
  free(count);
    
}
