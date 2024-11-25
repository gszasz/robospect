/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
     
void generate_model_continuum_fft(opts *options, data *D, model *M) {
  int i;
  double *v = malloc(sizeof(double) * D->N);
  if (v == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  double cutoff;

  gsl_fft_real_wavetable * real;
  gsl_fft_halfcomplex_wavetable * hc;
  gsl_fft_real_workspace * work;

  for (i = 0; i < D->N; i++) {
    v[i] = D->y[i];
  }
  
  work = gsl_fft_real_workspace_alloc(D->N);
  real = gsl_fft_real_wavetable_alloc(D->N);

  gsl_fft_real_transform (v, 1, D->N,
			  real, work);
  gsl_fft_real_wavetable_free (real);
  
  cutoff = (options->continuum_box / options->delta_x);
  cutoff = 1.0 / cutoff;
  for (i = 0; i < D->N; i++) {
    log_comment(options,ROBO_VERBOSE_CONTINUUM,
  		"gmcFT: (%d/%d) (%f %f %f) (%f %f %ld) %f",
  		i,D->N,
  		D->x[i],D->y[i],D->e[i],
		cutoff,0.0,0,
		v[i]
  		);
    if (i > cutoff) {
	v[i] = 0;
    }
  }

  
  hc = gsl_fft_halfcomplex_wavetable_alloc(D->N);
  gsl_fft_halfcomplex_inverse (v, 1, D->N,
			       hc, work);
  gsl_fft_halfcomplex_wavetable_free (hc);
  
  gsl_fft_real_workspace_free (work);

  for (i = 0; i < D->N; i++) {
    M->continuum[i] = v[i];
    log_comment(options,ROBO_VERBOSE_CONTINUUM,
  		"gmcF: (%d/%d) (%f %f %f) (%f %f %ld) %f",
  		i,D->N,
  		D->x[i],D->y[i],D->e[i],
		0.0,0.0,0,
  		M->continuum[i]
  		);

  }
  
  free(v);
}
