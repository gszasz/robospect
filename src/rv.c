/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "robo.h"

#define C_KM_S 299792.458
double radial_velocity_transform(double rv, double L) {
  return(L * (1.0 + rv/C_KM_S));
}
double radial_velocity_invtransform(double rv, double L) {
  return(L / (1.0 + rv/C_KM_S));
}

void generate_simulated_altspectra_from_lines(model *rv_model,lines *rv_inlines) {
  int i,j;
  
  for (j = 0; j <  rv_model->N; j++) {
    rv_model->alternate[j] = 0;
  }
  for (i = 0; i < rv_inlines->l; i++) {
    for (j = 0; j <  rv_model->N; j++) {
      if ((rv_model->x[j] < rv_inlines->m[i] - SIGMA_RANGE * rv_inlines->s[i])||
	  (rv_model->x[j] > rv_inlines->m[i] + SIGMA_RANGE * rv_inlines->s[i])) {
	continue;
      }
      rv_model->alternate[j] +=  gaussian_line(rv_model->x[j],rv_inlines,i);
    }
  }
  /* if (rv_model->O->flux_calibrated != 0) { */
  /*   for (j = 0; j <  rv_model->N; j++) { */
  /*     //      rv_model->alternate[j] *= rv_model->continuum[j]; */
  /*   } */
  /* } */
}

void generate_simulated_spectra_from_lines(model *rv_model,lines *rv_inlines) {
  int i,j;
  
  for (j = 0; j <  rv_model->N; j++) {
    rv_model->lines[j] = 0;
  }
  for (i = 0; i < rv_inlines->l; i++) {
    for (j = 0; j <  rv_model->N; j++) {
      if ((rv_model->x[j] < rv_inlines->m[i] - SIGMA_RANGE * rv_inlines->s[i])||
	  (rv_model->x[j] > rv_inlines->m[i] + SIGMA_RANGE * rv_inlines->s[i])) {
	continue;
      }
      rv_model->lines[j] +=  gaussian_line(rv_model->x[j],rv_inlines,i);
    }
  }
  /* if (rv_model->O->flux_calibrated != 0) { */
  /*   for (j = 0; j <  rv_model->N; j++) { */
  /*     //      rv_model->lines[j] *= rv_model->continuum[j]; */
  /*   } */
  /* } */
}

#define SPLINE_TYPE gsl_interp_akima
double rv_xcorr(model *rv_model, double rv_range, double rv_max_error, int rv_steps) {
  int i;
  double corr = 0.0;
  double xt   = 0.0;
  double corr_best = -99e99;
  double rv_best   = 0.0;
  double rv;
  double rv_min = -1.0 * rv_range;
  double rv_max =  1.0 * rv_range;
  gsl_spline *data_spline = gsl_spline_alloc(SPLINE_TYPE ,rv_model->N);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  int rv_err    = rv_max - rv_min;
  int j = 0;
  int count = 0;
  gsl_spline_init(data_spline,rv_model->x,rv_model->alternate,rv_model->N);

  count = 0;
  for (i = 0; i < rv_model->N; i++) {
    if (rv_model->lines[i] != 0) {
      count++;
      log_comment(rv_model->O,ROBO_VERBOSE_XCORR,
		  "XCORR_lines: %d %d %f",
		  i,count,rv_model->x[i]);
    }
  }
  
  log_comment(rv_model->O,ROBO_VERBOSE_XCORR,
	      "XCORR_key: iter Nlines rv corr min best max error");
  
  do {
    for (rv = rv_min; rv <= rv_max; rv += (rv_max - rv_min) / rv_steps) {
      corr = 0.0;
      for (i = 0; i < rv_model->N; i++) {
	if (rv_model->lines[i] != 0) {
	  xt = radial_velocity_invtransform(rv,rv_model->x[i]);
	  if ((xt > rv_model->x[0])&&
	      (xt < rv_model->x[rv_model->N - 1])) {
	    corr += rv_model->lines[i] * gsl_spline_eval(data_spline,xt,acc);
	  }
	}
      }
      
      if (corr > corr_best) {
	corr_best = corr;
	rv_best = rv;
      }
      log_comment(rv_model->O,ROBO_VERBOSE_XCORR,
		  "XCORR: %d %d %f %f  %f %f %f  %f",
		  j,count,rv,log10(corr),
		  rv_min,rv_best,rv_max,rv_err);
    }
    rv_min = rv_best - (rv_max - rv_min) / rv_steps;
    rv_max = rv_best + (rv_max - rv_min) / rv_steps;
    rv_err = rv_max - rv_min;
    j++;

  } while ((rv_err > rv_max_error)&&(j < 10));

  gsl_spline_free(data_spline);
  gsl_interp_accel_free(acc);
  return(rv_best);
}

void measure_radial_velocity(opts *options, data *D, model *M) {
  opts *rv_options = alloc_options();
  data *rv_data;
  model *rv_model;
  lines *rv_inlines;
  double rv;
  int i;
  double av_flux,av_sigma;
  int av_N = 0;
  /* Create a fake RS environment for calculations */
  rv_options->continuum_box = options->continuum_box;
  rv_options->find_lines    = 1;
  rv_options->find_sigma    = options->rv_sigma;
  rv_options->flux_calibrated = options->flux_calibrated;
  rv_options->min_x           = options->min_x;
  rv_options->max_x           = options->max_x;
  /* Things I don't care about */
  rv_options->psf_width       = options->psf_width;
  rv_options->deblend_radius  = options->deblend_radius;
  rv_options->deblend_ratio   = options->deblend_ratio;
  rv_options->deblend_iterations= options->deblend_iterations;
  rv_options->tolerance       = options->tolerance;
  rv_options->relax           = options->relax;
  rv_options->max_iterations  = options->max_iterations;
  rv_options->iteration       = options->iteration;
  rv_options->wavelength_min_error = options->wavelength_min_error;
  rv_options->wavelength_max_error = options->wavelength_max_error;
  rv_options->wavelength_limit = options->wavelength_limit;
  rv_options->measure_radial_velocity = 0;
  rv_options->save_temp = 0;
  rv_options->plot_all = 0;
  rv_options->help = 0;
  rv_options->verbose = options->verbose;
  rv_options->fault = 0;
  rv_options->infilename = 0;
  rv_options->delta_x = options->delta_x;
  rv_options->order = 0;
  rv_options->max_order = 0;
  rv_options->supplied_errors = 0;
  rv_options->path_base = malloc(32 * sizeof(char));
  rv_options->line_list_filename = NULL;
  rv_options->strict = 0;
  rv_options->fits_IO = 0;
  rv_options->fits_row = 0;
  rv_options->log = options->log;
  rv_options->command_name = NULL;
  rv_options->command_line = NULL;
  
  rv_options->continuum_model_name = malloc(32 * sizeof(char));
  snprintf(rv_options->continuum_model_name,32 * sizeof(char),"%s",options->continuum_model_name);
  rv_options->line_model_name = malloc(32 * sizeof(char));
  snprintf(rv_options->line_model_name,32 * sizeof(char),"pre");
  rv_options->noise_model_name = malloc(32 * sizeof(char));
  snprintf(rv_options->noise_model_name,32 * sizeof(char),"%s",options->noise_model_name);
  rv_options->deblend_model_name = malloc(32 * sizeof(char));
  rv_options->function_model_name = malloc(32 * sizeof(char));
    
  rv_data  = malloc(sizeof(data));
  rv_data->N = D->N;
  rv_data->x = vector_copy(D->x , D->N);
  rv_data->y = vector_copy(D->y , D->N);
  rv_data->e = vector_copy(D->e , D->N);
  rv_data->yO = vector_copy(D->yO , D->N);
  rv_data->O = rv_options;
  rv_model = alloc_models(rv_options,rv_data);
  rv_model->L = NULL;
  rv_model->O = rv_options;
  /* Construct a set of observed lines. */
  set_continuum_fit_spectra(rv_data,rv_model);
  generate_model_continuum(rv_options,rv_data,rv_model);
  set_noise_fit_spectra(rv_data,rv_model);
  generate_noise_model(rv_options,rv_data,rv_model);
  set_linefinder_Z_spectra(rv_data,rv_model);
  
  linefinder(rv_options,rv_data,rv_model);
  set_line_fit_spectra(rv_data,rv_model);
  measure_lines(rv_data,rv_model);

  /* Measure mean properties */
  av_flux = 0;
  av_sigma = 0;
  
  for (i = 0; i < rv_model->L->l; i++) {
    if ((isfinite(rv_model->L->F[i]))&&
	(isfinite(rv_model->L->s[i]))) {
      av_flux += rv_model->L->F[i];
      av_sigma += rv_model->L->s[i];
      av_N ++;
    }
  }
  av_flux /= av_N;
  av_sigma /= av_N;

  /* Construct a simulation of the input line list. */
  rv_inlines = alloc_lines(M->L->l);
  for (i = 0; i < M->L->l; i++) {
    rv_inlines->m[i] = M->L->x0[i];
    rv_inlines->s[i] = av_sigma;
    rv_inlines->F[i] = av_flux;
    rv_inlines->eta[i] = 0.0;

  }
  for (i = 0; i < rv_model->L->l; i++) {
    rv_model->L->s[i] = av_sigma;
    rv_model->L->F[i] = av_flux;
  }

  generate_simulated_altspectra_from_lines(rv_model,rv_inlines);
  generate_simulated_spectra_from_lines(rv_model,rv_model->L);

  /* Do cross correlation */
  rv = rv_xcorr(rv_model,
		options->rv_range,
		options->rv_max_error,
		options->rv_steps
		);
  /* Save results */
  options->radial_velocity = rv;
  for (i = 0; i < M->L->l; i++) {
    M->L->x0[i] = radial_velocity_transform(rv,M->L->x0[i]);
  }
  log_comment(options,ROBO_VERBOSE_ALL,"radial_velocity: Measured radial velocity of %g km/s",rv);
  
  /* Free everything we don't need. */
  free_data(rv_data);
  free_models(rv_model);
  rv_options->log = NULL;
  free_options(rv_options);
  
  
}

  
