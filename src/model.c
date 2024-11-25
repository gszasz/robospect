/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "robo.h"
#define NAMESIZE 32

model *alloc_models(opts *options, data *D) {
  int i;
  model *M = malloc(sizeof(model));
  if (M == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  M->N = D->N;
  M->x = malloc(M->N * sizeof(double));
  if (M->x == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  M->continuum = malloc(M->N * sizeof(double));
  if (M->continuum == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  M->lines = malloc(M->N * sizeof(double));
  if (M->lines == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  M->alternate = malloc(M->N * sizeof(double));
  if (M->alternate == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < D->N; i++) {
    M->x[i] = D->x[i];
    M->continuum[i] = 1.0;
    M->lines[i] = 0.0;
  }
  if (!options->continuum_model_name) {
    options->continuum_model_name = malloc(NAMESIZE * sizeof(char));
    if (options->continuum_model_name == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    snprintf(options->continuum_model_name,NAMESIZE,"default");
  }
  if (!options->line_model_name) {
    options->line_model_name = malloc(NAMESIZE * sizeof(char));
    if (options->line_model_name == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    snprintf(options->line_model_name,NAMESIZE,"default");
  }
  if (!options->noise_model_name) {
    options->noise_model_name = malloc(NAMESIZE * sizeof(char));
    if (options->noise_model_name == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    snprintf(options->noise_model_name,NAMESIZE,"default");
  }
  if (!options->deblend_model_name) {
    options->deblend_model_name = malloc(NAMESIZE * sizeof(char));
    if (options->deblend_model_name == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    snprintf(options->deblend_model_name,NAMESIZE,"default");
  }
  if (!options->function_model_name) {
    options->function_model_name = malloc(NAMESIZE * sizeof(char));
    if (options->function_model_name == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    snprintf(options->function_model_name,NAMESIZE,"default");
  }

  /* NOISE MODES */
  if ((strcasecmp(options->noise_model_name,"null") == 0)) {
    M->noise_model = NOISE_NULL;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using noise model %d NOISE_NULL",
		M->noise_model);
  }
  else if ((strcasecmp(options->noise_model_name,"poisson") == 0)) {
    M->noise_model = NOISE_POISSON;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using noise model %d NOISE_POISSON",
		M->noise_model);
  }
  else if (options->supplied_errors == 1) {
    M->noise_model = NOISE_SUPPLIED;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using noise model %d NOISE_SUPPLIED",
		M->noise_model);
  }
  else if ((strcasecmp(options->noise_model_name,"slow") == 0)) {
    M->noise_model = NOISE_SLOW_BOXCAR;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using noise model %d NOISE_SLOW_BOXCAR",
		M->noise_model);
  }
  else {
    M->noise_model = NOISE_BOXCAR;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using noise model %d NOISE_BOXCAR",
		M->noise_model);
  }

  /* CONTINUUM MODES */
  if (strcasecmp(options->continuum_model_name,"boxcar") == 0) {
    M->continuum_model = CONTINUUM_BOXCAR;
    if (M->noise_model == NOISE_BOXCAR) {
      M->noise_model     = NOISE_NULL;
    }
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_BOXCAR",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"logboxcar") == 0) {
    M->continuum_model = CONTINUUM_LOGBOXCAR;
    if (M->noise_model == NOISE_BOXCAR) {
      M->noise_model     = NOISE_NULL;
    }
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_LOGBOXCAR",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"histogram") == 0) {
    M->continuum_model = CONTINUUM_HISTOGRAM;
    if (M->noise_model == NOISE_BOXCAR) {
      M->noise_model     = NOISE_NULL;
    }
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_HISTOGRAM",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"linear") == 0) {
    M->continuum_model = CONTINUUM_ROBUST_LINEAR;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_ROBUST_LINEAR",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"penalized") == 0) {
    M->continuum_model = CONTINUUM_PENALIZED_SPLINE;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_PENALIZED_SPLINE",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"bspline") == 0) {
    M->continuum_model = CONTINUUM_BSPLINE;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_BSPLINE",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"fft") == 0) {
    M->continuum_model = CONTINUUM_FFT;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_FFT",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"peak") == 0) {
    M->continuum_model = CONTINUUM_PEAK;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_PEAK",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"blackbody") == 0) {
    M->continuum_model = CONTINUUM_BLACKBODY;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_BLACKBODY",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"powerlaw") == 0) {
    M->continuum_model = CONTINUUM_POWERLAW;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_POWERLAW",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"devel") == 0) {
    M->continuum_model = CONTINUUM_DEVEL;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_DEVEL",
		M->continuum_model);
  }
  else if (strcasecmp(options->continuum_model_name,"null") == 0) {
    M->continuum_model = CONTINUUM_NULL;
    M->noise_model = NOISE_BOXCAR;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_NULL",
		M->continuum_model);
  }
  else {
    M->continuum_model = CONTINUUM_BOXCAR;
    M->noise_model     = NOISE_NULL;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using continuum model %d CONTINUUM_BOXCAR",
		M->continuum_model);
  }
  /* DEBLENDING MODES */
  if (strcasecmp(options->deblend_model_name,"nlls") == 0) {
    M->deblend_model = DEBLEND_NLLS;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using deblend model %d DEBLEND_NLLS",
		M->deblend_model);
  }
  else {
    M->deblend_model = DEBLEND_NULL;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using deblend model %d DEBLEND_NULL",
		M->deblend_model);
  }
  
  /* LINE FITTING MODES */
  if (strcasecmp(options->line_model_name,"gauss") == 0) {
    M->line_model = LINE_GAUSSIAN;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_GAUSSIAN",
		M->line_model);
  }
  else if (strcasecmp(options->line_model_name,"fixed") == 0) {
    M->line_model = LINE_FIXEDMEAN_GAUSSIAN;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_FIXEDMEAN_GAUSSIAN",
		M->line_model);
  }
  else if (strcasecmp(options->line_model_name,"nofit") == 0) {
    M->line_model = LINE_NULL;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_NULL",
		M->line_model);
  }
  else if (strcasecmp(options->line_model_name,"isolated") == 0) {
    M->line_model = LINE_FIXEDMEAN_ISOLATED_GAUSSIAN;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_ISOLATED_GAUSSIAN",
		M->line_model);
  }
  else if (strcasecmp(options->line_model_name,"pre") == 0) {
    M->line_model = LINE_PRE;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_PRE",
		M->line_model);
  }
  else if (strcasecmp(options->line_model_name,"nonparametric") == 0) {
    M->line_model = LINE_NONPARAMETRIC;
    M->function_model = FUNCTION_NONPARAMETRIC;
    options->function_model_name = options->line_model_name;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_NONPARAMETRIC",
		M->line_model);
  }
  else if (strcasecmp(options->line_model_name,"preoriginal") == 0) {
    M->line_model = LINE_PRE_ORIGINAL;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_PRE",
		M->line_model);
  }
  else if (strcasecmp(options->line_model_name,"twostage") == 0) {
    M->line_model = LINE_TWO_STAGE;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_TWO_STAGE",
		M->line_model);
  }
  else if (strcasecmp(options->line_model_name,"twosticky") == 0) {
    M->line_model = LINE_TWO_STICKY;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_TWO_STICKY",
		M->line_model);
  }
  else if ((strcasecmp(options->line_model_name,"threestage") == 0)||
	   (strcasecmp(options->line_model_name,"best") == 0)) {
    M->line_model = LINE_THREE_STAGE;
    M->deblend_model = DEBLEND_NULL;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_THREE_STAGE",
		M->line_model);
  }
  else {
    M->line_model = LINE_THREE_STAGE;
    M->deblend_model = DEBLEND_NULL;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line model %d LINE_DEFAULT",
		M->line_model);
  }
  
  /* LINE FUNCTION MODES */
  if (strcasecmp(options->function_model_name,"gaussian") == 0) {
    M->function_model = FUNCTION_GAUSSIAN;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line function %d FUNCTION_GAUSSIAN",
		M->function_model);
  }
  else if (strcasecmp(options->function_model_name,"lorentzian") == 0) {
    M->function_model = FUNCTION_LORENTZIAN;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line function %d FUNCTION_LORENTZIAN",
		M->function_model);
  }
  else if (
	   (strcasecmp(options->function_model_name,"pseudovoigt") == 0)) {
    M->function_model = FUNCTION_PSEUDOVOIGT;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line function %d FUNCTION_PSEUDOVOIGT",
		M->function_model);
  }
  else if ((strcasecmp(options->function_model_name,"hjerting") == 0)) {
    M->function_model = FUNCTION_HJERTING;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line function %d FUNCTION_HJERTING",
		M->function_model);
  }
  else if ((strcasecmp(options->function_model_name,"voigt") == 0)||
	   (strcasecmp(options->function_model_name,"humlicek") == 0)) {
    M->function_model = FUNCTION_HUMLICEK;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line function %d FUNCTION_HUMLICEK",
		M->function_model);
  }
  else if ((strcasecmp(options->function_model_name,"skewgauss") == 0)) {
    M->function_model = FUNCTION_SKEWGAUSS;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line function %d FUNCTION_SKEWGAUSS",
		M->function_model);
  }
  else if ((strcasecmp(options->function_model_name,"nonparametric") == 0)) {
    M->function_model = FUNCTION_NONPARAMETRIC;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line function %d FUNCTION_NONPARAMETRIC",
		M->function_model);
  }
  else {
    M->function_model = FUNCTION_GAUSSIAN;
    log_comment(options,ROBO_VERBOSE_ALL,"alloc_models: Using line function %d FUNCTION_DEFAULT",
		M->function_model);
  }  
  
  return(M);
}

void free_models(model *M) {
  free(M->x);
  free(M->continuum);
  free(M->lines);
  free(M->alternate);
  free(M);
}

void generate_model_continuum(opts *options, data *D, model *M) {
    switch (M->continuum_model) {
    case CONTINUUM_ROBUST_LINEAR:
      generate_model_continuum_robust_linear(options,D,M);
      break;
    case CONTINUUM_PENALIZED_SPLINE:
      generate_model_continuum_penalized_spline(options,D,M);
      break;
    case CONTINUUM_BSPLINE:
      generate_model_continuum_bspline(options,D,M);
      break;
    case CONTINUUM_NULL:
      generate_model_continuum_null(options,D,M);
      break;
    case CONTINUUM_HISTOGRAM:
      generate_model_continuum_histogram(options,D,M);
      break;
    case CONTINUUM_FFT:
      generate_model_continuum_fft(options,D,M);
      break;
    case CONTINUUM_PEAK:
      generate_model_continuum_peak_boxcar(options,D,M);
      break;
    case CONTINUUM_BLACKBODY:
      generate_model_continuum_blackbody(options,D,M);
      break;
    case CONTINUUM_POWERLAW:
      generate_model_continuum_powerlaw(options,D,M);
      break;
    case CONTINUUM_DEVEL:
      generate_model_continuum_devel(options,D,M);
      break;
    case CONTINUUM_LOGBOXCAR:
      generate_model_continuum_logboxcar(options,D,M);
      break;
    case CONTINUUM_BOXCAR:
    default:
      generate_model_continuum_boxcar(options,D,M);
      break;
    }
}

void measure_lines(data *D, model *M) {
  int i,j;
  double v;

  /* We should always do the pre-measurement, as we'll want to know what the lines look like before we try to fit. */
  
  measure_lines_PRE(D,M);

  /* Reset flag status. */
  for (i = 0; i < M->L->l; i++) {
    if ((M->L->flags[i] & ROBO_INIT)||
    	(M->L->flags[i] & ROBO_BAD_FIT)||
	(M->line_model == LINE_PRE)) {
      /* This is the only place we should be assigning values to the pre fit. */
      M->L->m[i] = M->L->mp[i];
      M->L->s[i] = M->L->sp[i];
      M->L->F[i] = M->L->Fp[i];
      M->L->eta[i] = 0.0;

      /* Set the errors to dummy values so they're set, but also so it's clear they aren't measured. */
      M->L->dm[i] = -0.01 * M->L->m[i];
      M->L->ds[i] = -0.1 * M->L->s[i];
      M->L->dF[i] = -0.1 * M->L->F[i];
      M->L->deta[i] = -0.1 * M->L->eta[i];
      M->L->chi[i] = -1.0;
    }
    M->L->flags[i] &= ROBO_RESET;
    if ((M->L->flags[i] & ROBO_INIT)&&
	((M->function_model == FUNCTION_PSEUDOVOIGT)||
	 (M->function_model == FUNCTION_HJERTING)||
	 (M->function_model == FUNCTION_HUMLICEK))) {
      M->L->eta[i] = 0.5;
    }
  }
  

  switch (M->line_model) {
  case LINE_GAUSSIAN: 
    measure_lines_gaussian(D,M);
    break;
  case LINE_FIXEDMEAN_GAUSSIAN:
    measure_lines_fixedmean_gaussian(D,M);
    break;
  case LINE_FIXEDMEAN_ISOLATED_GAUSSIAN:
    measure_lines_fixedmean_isolation_gaussian(D,M);
    break;
  case LINE_PRE:
    for (i = 0; i < M->L->l; i++) {
      M->L->m[i] = M->L->mp[i];
      M->L->s[i] = M->L->sp[i];
      M->L->F[i] = M->L->Fp[i];
      M->L->eta[i] = 0.0;
      M->L->chi[i] = -1.0;
    }
    break;
  case LINE_PRE_ORIGINAL:
    measure_lines_PRE_original(D,M);
    break;
  case LINE_TWO_STAGE:
    measure_lines_twostage(D,M);
    break;
  case LINE_TWO_STICKY:
    measure_lines_twostage_sticky(D,M);
    break;
  case LINE_THREE_STAGE:
    measure_lines_threestage(D,M);
    break;
  case LINE_NONPARAMETRIC:
    measure_lines_nonparametric(D,M);
    break;
  case LINE_NULL:
  default:
    measure_lines_NULL(D,M);
    break;
  }

  switch (M->deblend_model) {
  case DEBLEND_NLLS:
    deblend_lines_nlls(D,M);
    break;
  case DEBLEND_NULL:
  default:
    break;
  }

  /* Validate lines here. Truncate excessive errors. */
  set_flags(D,M);

  if (M->function_model != FUNCTION_NONPARAMETRIC) {
    /* Accumulate here. */
    for (j = 0; j < M->N; j++) {
      M->lines[j] = 0.0;
      M->alternate[j] = 0.0;
    }
    for (i = 0; i < M->L->l; i++) {
      if (!(M->L->flags[i] & (ROBO_FIT_FAIL|ROBO_FIT_REFUSED|ROBO_FIT_IGNORED))) {
	for (j = 0; j < M->N; j++) {
	  if (M->function_model != FUNCTION_HUMLICEK) {
	    if ((M->x[j] < M->L->m[i] - SIGMA_RANGE * M->L->s[i])||
		(M->x[j] > M->L->m[i] + SIGMA_RANGE * M->L->s[i])) {
	      continue;
	    }
	  }
	  v = 0;
	  if (M->function_model == FUNCTION_GAUSSIAN) {
	    v = gaussian_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_LORENTZIAN) {
	    v = lorentzian_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_PSEUDOVOIGT) {
	    v = pseudovoigt_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_HJERTING) {
	    v = hjerting_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_HUMLICEK) {
	    v = humlicek_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_SKEWGAUSS) {
	    v = skewgauss_line(M->x[j],M->L,i);
	  }
	  log_comment(D->O,ROBO_VERBOSE_INTERMEDIATE,"%d %g : %g %g %d",i,M->L->x0[i],v,v,j );
	  if (v == v) {  /* This should always be true, as we have previous validated the lines. */
	    M->lines[j] += v;
	  }
	}
      }
      if (!(M->L->flags[i] & (ROBO_ALT_REFUSED))) {
	for (j = 0; j < M->N; j++) {
	  if ((M->x[j] < M->L->mp[i] - SIGMA_RANGE * M->L->sp[i])||
	      (M->x[j] > M->L->mp[i] + SIGMA_RANGE * M->L->sp[i])) {
	    continue;
	  }
	  v = gaussian_line_alt(M->x[j],M->L,i);
	  if (v == v) {
	    M->alternate[j] += v;
	  }
	}
      }
    }
    /* Check that all lines contribute to an improved chi^2 value. */
    chi_check_lines(D,M);
  }
  /* Correct accumulated model spectra for flux calibrated data back to flux calibrated (temporary). */
  if (D->O->flux_calibrated) {
    double *t = vector_multiply(M->lines,M->continuum,M->N);
    free(M->lines);
    M->lines = t;
  }
}  

void set_flags(data *D, model *M) {
  int i;
  for (i = 0; i < M->L->l; i++) {
    /* Values must be finite. */
    if ((!isfinite(M->L->m[i]))||
	(!isfinite(M->L->s[i]))||
	(!isfinite(M->L->F[i]))) {
      M->L->flags[i] |= ROBO_FIT_FAIL;
    }
    if (M->function_model == FUNCTION_PSEUDOVOIGT) {
      if (!isfinite(M->L->eta[i])) {
	M->L->flags[i] |= ROBO_FIT_FAIL;
      }
    }
    /* Lines must be within the spectra */
    if ((M->L->m[i] < D->O->min_x)||
	(M->L->m[i] > D->O->max_x)) {
      M->L->m[i] = M->L->x0[i];
      M->L->flags[i] |= ROBO_FIT_FAIL;
    }
    /* Errors may not exceed a certain level. */
    if (fabs(M->L->dm[i]) > ERROR_MAX_SN * fabs(M->L->m[i])) {
      M->L->dm[i] = -1.0;
      M->L->flags[i] |= ROBO_BAD_ERROR_VAL;
    }
    if (fabs(M->L->ds[i]) > ERROR_MAX_SN * fabs(M->L->s[i])) {
      M->L->ds[i] = -1.0;
      M->L->flags[i] |= ROBO_BAD_ERROR_VAL;
    }
    if (fabs(M->L->dF[i]) > ERROR_MAX_SN * fabs(M->L->F[i])) {
      M->L->dF[i] = -1.0;
      M->L->flags[i] |= ROBO_BAD_ERROR_VAL;
    }
    if (fabs(M->L->deta[i]) > ERROR_MAX_SN * fabs(M->L->eta[i])) {
      M->L->deta[i] = -1.0;
      if (M->function_model == FUNCTION_PSEUDOVOIGT) {
	M->L->flags[i] |= ROBO_BAD_ERROR_VAL;
      }
    }
    /* This is a hack that I should rework to better determine if a fit is good or not. */
    if ((D->O->strict == 1)&&(D->O->flux_calibrated == 0)) {
      if ((fabs(M->L->m[i] - M->L->x0[i]) > D->O->strict_center)||
	  (fabs(M->L->s[i]) > D->O->strict_width)||
	  ((fabs(M->L->F[i]) > D->O->strict_flux)&&(fabs(M->L->F[i] / M->L->s[i]) > fabs(M->L->s[i])))) {
	M->L->flags[i] |= ROBO_FIT_IGNORED;
      }
    }
    /* See if the alternate model is reasonable. */
    if ((!isfinite(M->L->mp[i]))||
	(!isfinite(M->L->sp[i]))||
	(!isfinite(M->L->Fp[i]))) {
      M->L->flags[i] |= ROBO_ALT_REFUSED;
    }
    /* If alternate model d isagrees with the fit model, reject the alternate. */
    if ((M->L->sp[i] == 0.0)||
	(M->L->sp[i] > 100 * M->L->s[i])||
	(M->L->sp[i] > 0.01 * (D->O->max_x - D->O->min_x))) {
      M->L->flags[i] |= ROBO_ALT_REFUSED;
    }
    /* Set the model function used, if used. */
#if MULTIFUNCTION
    if (M->function_model == FUNCTION_GAUSSIAN) {
      M->L->flags[i] |= ROBO_FUNC_GAUSSIAN;
    }
    if (M->function_model == FUNCTION_LORENTZIAN) {
      M->L->flags[i] |= ROBO_FUNC_LORENTZIAN;
    }
    if (M->function_model == FUNCTION_HJERTING) {
      M->L->flags[i] |= ROBO_FUNC_HJERTING;
    }
    if (M->function_model == FUNCTION_HUMLICEK) {
      M->L->flags[i] |= ROBO_FUNC_HUMLICEK;
    }
#endif
  }  
}

void chi_check_lines(data *D, model *M) {
  int i,j;
  double chi,chiP;
  int N;
  double v = 0;

  for (i = 0; i < M->L->l; i++) {
    chi  = 0;
    chiP = 0;
    N = 0;
    if (!(M->L->flags[i] & (ROBO_FIT_FAIL|ROBO_FIT_IGNORED))) {
      for (j = 0; j < M->N; j++) {
	if (M->function_model != FUNCTION_HUMLICEK) {
	  if ((M->x[j] < M->L->m[i] - SIGMA_RANGE * M->L->s[i])||
	      (M->x[j] > M->L->m[i] + SIGMA_RANGE * M->L->s[i])) {
	    continue;
	  }
	}
	if (M->function_model == FUNCTION_GAUSSIAN) {
	  v = gaussian_line(M->x[j],M->L,i);
	}
	else if (M->function_model == FUNCTION_LORENTZIAN) {
	  v = lorentzian_line(M->x[j],M->L,i);
	}
	else if (M->function_model == FUNCTION_PSEUDOVOIGT) {
	  v = pseudovoigt_line(M->x[j],M->L,i);
	}
	else if (M->function_model == FUNCTION_HJERTING) {
	  v = hjerting_line(M->x[j],M->L,i);
	}
	else if (M->function_model == FUNCTION_HUMLICEK) {
	  v = humlicek_line(M->x[j],M->L,i);
	}
	else if (M->function_model == FUNCTION_SKEWGAUSS) {
	  v = skewgauss_line(M->x[j],M->L,i);
	}
	log_comment(D->O,ROBO_VERBOSE_INTERMEDIATE,"%d %g : %g %g %d",i,M->L->x0[i],v,v,j );

	if ((v == v)&&(isfinite(D->e[j]))&&(D->e[j] != 0)) {
	  /*	if (v == v) { */
	  chi  += pow((D->y[j] - M->lines[j]) / D->e[j],2);
	  chiP += pow((D->y[j] - (M->lines[j] - v)) / D->e[j],2);
	  N += 1;
	}
      }
      log_comment(D->O,ROBO_VERBOSE_INTERMEDIATE,"CHIX: %d %g : %g %g %d",i,M->L->x0[i],chi,chiP,N);
      /*      if (chiP / (D->N - (3.0 * (M->L->l - 1))) < chi / (D->N - (3.0 * M->L->l))) {  */
      if (chiP == chi) { /* We had zero contribution in this line anyway.  */
      }
      else if (((chiP > 3)&&(chiP < chi + 3))||
      	       ((chiP < 3)&&(chi > 0.5 * chiP))) { /* Using chi^2 for of AIC.  chiP  */
                           	 /* represents teh "simple" model excluding this line.  chi is teh complex model  */
                                 /* Therefore, we only should accept chiP if it has a lower AIC than chi. */
      	                         /* The +3 adds the complexity of an extra set of parameters involved in this */
      	                         /* complex model. */
	for (j = 0; j < M->N; j++) {
	  if (M->function_model != FUNCTION_HUMLICEK) {
	    if ((M->x[j] < M->L->m[i] - SIGMA_RANGE * M->L->s[i])||
		(M->x[j] > M->L->m[i] + SIGMA_RANGE * M->L->s[i])) {
	      continue;
	    }
	  }
	  if (M->function_model == FUNCTION_GAUSSIAN) {
	    v = gaussian_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_LORENTZIAN) {
	    v = lorentzian_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_PSEUDOVOIGT) {
	    v = pseudovoigt_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_HJERTING) {
	    v = hjerting_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_HUMLICEK) {
	    v = humlicek_line(M->x[j],M->L,i);
	  }
	  else if (M->function_model == FUNCTION_SKEWGAUSS) {
	    v = skewgauss_line(M->x[j],M->L,i);
	  }
	  log_comment(D->O,ROBO_VERBOSE_INTERMEDIATE,"%d %g : %g %g %d",i,M->L->x0[i],v,v,j );

	  if (v == v) {
	    M->lines[j] -= v;
	  }
	}
      }
    }
  }
}


