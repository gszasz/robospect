/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

/* This function needs to calculate the best nonparametric model under */
/* the following conditions: */
/*  1. The line shape is symmetric around the mean, m. */
/*  2. The line is constrained to be monotonically decreasing as (x - m) increases. */
/*  3. A Gaussianity measure needs to be calculated to determine how much this */
/*     line deviates from a true Gaussian (likely to be easier in log-space). */
/*  4. Defined as the minimum envelope that matches the observed result. */
#define SPLINE_TYPE gsl_interp_akima

void measure_lines_nonparametric(data *D, model *M) {
  int i,j;
  double v;

  gsl_spline *data_spline = gsl_spline_alloc(SPLINE_TYPE ,D->N);
  gsl_spline *error_spline = gsl_spline_alloc(SPLINE_TYPE ,D->N);
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  
  gsl_spline *np_model;
  gsl_interp_accel *npacc = gsl_interp_accel_alloc();

  double chi = 0;
  double ochi = 99e99;

  int     failed_update = 0;
  int     envelope_N = 5;
  int     envelope_i = 0;
  double *envelope_X = malloc(sizeof(double) * envelope_N);
  if (envelope_X == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  double *envelope_Y = malloc(sizeof(double) * envelope_N);
  if (envelope_Y == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  double x1,x2;
  double y1,y2;
  double dy1,dy2;
  double ddy1,ddy2;
  double gradient;
  
  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }

  gsl_spline_init(data_spline,D->x,D->y,D->N);
  gsl_spline_init(error_spline,D->x,D->e,D->N);

  for (i = 0; i < M->L->l; i++) {
      /* Do math here. */
      /* Identify center sample */
    v = gsl_spline_eval(data_spline,M->L->x0[i],acc);

    /* Identify best mean using symmetry */
    /* chi^2 = sum( (Y(x - m) - Y(m - x)) / E) */
    /* dchi/dm = dY/dm(x-m) - dY/dm(m - x) */
    chi = 0;
    ochi = 99e99;
    /* Determine envelope */

    /* This needs to do three things: */
    /*  1. Average left and right flux values and choose the smaller. */
    /*  2. Average left and right flux derivative and choose the larger. */
    /*  3. Force flux to monotonically decrease. */
    /*  4. Force flux derivative to rise to a maximum, and then monotonically decrease */

    if ((M->L->m[i] < D->O->min_x)||
	(M->L->m[i] > D->O->max_x)) {
      M->L->m[i] = M->L->x0[i];
    }
    y1 = gsl_spline_eval(data_spline,M->L->m[i],acc);
    dy1= gsl_spline_eval_deriv(data_spline,M->L->m[i],acc);
    ddy1=gsl_spline_eval_deriv2(data_spline,M->L->m[i],acc);
    y2 = y1; dy2 = dy1; ddy2 = ddy1;
    j = 0;
    envelope_Y[envelope_i] = y1;

    /*    while (envelope_Y[envelope_i] / gsl_spline_eval(data_spline,M->L->m[i],acc) > 0.0) { */
    while (j < 1000) {
      /* Check that we're within bounds */
      x1 = M->L->m[i] - j * D->O->delta_x;
      x2 = M->L->m[i] + j * D->O->delta_x;
      if ((x1 < D->O->min_x)||(x2 < D->O->min_x)||
	  (x1 > D->O->max_x)||(x2 > D->O->max_x)) {
	break;
      }

      /* Calculate splines */
      y1 = gsl_spline_eval(data_spline,x1,acc);
      y2 = gsl_spline_eval(data_spline,x2,acc);

      /* Select an initial guess. */
      if (fabs(y1) < fabs(y2)) {
	envelope_Y[envelope_i] = y1;
	envelope_X[envelope_i] = 1.0 * j * D->O->delta_x;
      }
      else {
	envelope_Y[envelope_i] = y2;
	envelope_X[envelope_i] = 1.0 * j * D->O->delta_x;
      }

      /* Check that we've decreased from the last sample */
      if ((envelope_i > 0) && (fabs(envelope_Y[envelope_i]) > fabs(envelope_Y[envelope_i - 1]))) {
	/* We failed this check, and should use gradient data */
	failed_update++;
	gradient = fmax(fabs(dy1),fabs(dy2));
	envelope_Y[envelope_i] = envelope_Y[envelope_i - 1] +
	  /*	  failed_update */
	  D->O->delta_x * gradient;
	gradient -= D->O->delta_x * fmax(fabs(ddy1),fabs(ddy2));
      }
      else {
	failed_update = 0;
	dy1 = gsl_spline_eval_deriv(data_spline,x1,acc);
	dy2 = gsl_spline_eval_deriv(data_spline,x2,acc);
	ddy1 = gsl_spline_eval_deriv2(data_spline,x1,acc);
	ddy2 = gsl_spline_eval_deriv2(data_spline,x2,acc);
      }
      log_comment(D->O,ROBO_VERBOSE_LINE,
		  "mlNP_iter: %d %d %d %g %g (%g %g) (%g %g %g %g) %g %g",
		  i,j,envelope_i,
		  envelope_X[envelope_i],envelope_Y[envelope_i],y1,y2,
		  dy1,dy2,ddy1,ddy2,
		  failed_update * D->O->delta_x * fmin(fabs(dy1),fabs(dy2)) -
		  pow(failed_update * D->O->delta_x,2) * 0.5 * fmin(fabs(ddy1),fabs(ddy2)),
		  envelope_Y[envelope_i] / gsl_spline_eval(data_spline,M->L->m[i],acc));
      if (envelope_Y[envelope_i] / gsl_spline_eval(data_spline,M->L->m[i],acc) < 0.0) {
	envelope_Y[envelope_i] = 0.0;
	break;
      }
      envelope_i++;
      if (envelope_i + 1 >= envelope_N) {
	envelope_N += 5;
	envelope_X = realloc(envelope_X,sizeof(double) * envelope_N);
	envelope_Y = realloc(envelope_Y,sizeof(double) * envelope_N);
      }
      j++;
    }

    envelope_N = envelope_i;
    /* Allocate the envelope spline */
    if (envelope_N <= 1) {
      M->L->F[i] = envelope_Y[0];
      M->L->s[i] = 0.0;
      M->L->eta[i] = 0.0;
      M->L->dm[i] = 0.0;
      M->L->ds[i] = 0.0;
      M->L->dF[i] = 0.0;
      continue;
    }
    else if (envelope_N <= 5) {
      np_model = gsl_spline_alloc(gsl_interp_linear,envelope_N);
    }
    else {
      np_model = gsl_spline_alloc(SPLINE_TYPE,envelope_N);
    }
    for (j = 1; j < envelope_N; j++) {
      if (envelope_X[j] <= envelope_X[j-1]) {
	envelope_X[j] = envelope_X[j-1] + 1e-4;
      }
    }
    gsl_spline_init(np_model,envelope_X,envelope_Y,envelope_N);
    /* Determine continuum level/zero gradient point */

    /* If the flux level hits the continuum level, we should terminate there. */
    /* The flux derivative constraint should help keep noise from bothering too much. */
    /* We should also measure the flux here, and save that. */

    /* Calculate non-Gaussianity */

    /* This should be stored in the eta value, and should represent how discrepant this envelope is */
    /* I'm not totally sure on the math for that yet. */
    M->L->F[i] = 0.0;
    M->L->s[i] = 0.0;
    M->L->eta[i] = 0.0;
    M->L->dm[i] = 0.0;
    M->L->ds[i] = 0.0;
    M->L->dF[i] = 0.0;

    /* This can be done easier using the spline code */
    M->L->F[i] = 2 * gsl_spline_eval_integ(np_model,0,envelope_X[envelope_N - 1],npacc);
    for (j = 0; j < envelope_N - 1; j++) {
      /* Sigma is defined from the FWHM, for lack of a better method. */
      if (M->L->s[i] == 0.0) {
	if (fabs(envelope_Y[j]) < 0.5 * fabs(envelope_Y[0])) {
	  M->L->s[i] = envelope_X[j] / (sqrt(2 * log(2)));
	  M->L->ds[i] = (envelope_X[j] - envelope_X[j-1]) / (sqrt(2 * log(2)));
	}
      }
      /* Calculate an eta proxy */
      if ((M->L->s[i] != 0.0)&&(M->L->eta[i] == 0.0)) {
	if (fabs(envelope_Y[j]) < 0.1 * fabs(envelope_Y[0])) {
	  M->L->eta[i] = envelope_X[j] / (envelope_X[j] + (sqrt(-0.5 * log(0.1)) * M->L->s[i]));
	}
      }
    }

    /* Add this line to the line accumulator */
    chi = 0;
    ochi = 0;

    for (j = 0; j < M->N; j++) {
      if (fabs(D->x[j] - M->L->m[i]) < envelope_X[envelope_N - 1]) {
	v = gsl_spline_eval(np_model,fabs(D->x[j] - M->L->m[i]),npacc);
	chi+= pow((D->y[j] - (M->lines[j] + v))/D->e[j],2); 
	ochi += pow((D->y[j] - M->lines[j])/D->e[j],2);
      }
    }
    if (chi < ochi) {
      for (j = 0; j < M->N; j++) {
	if (fabs(D->x[j] - M->L->m[i]) < envelope_X[envelope_N - 1]) {
	  v = gsl_spline_eval(np_model,fabs(D->x[j] - M->L->m[i]),npacc);
	  M->lines[j] += v;
	}
      }
    }
    log_comment(D->O,ROBO_VERBOSE_LINE,
		"mlNP: (%d/%d) %f (%f %f %f %f) +- (%f %f %f %f)",
		i,M->L->l,
		M->L->x0[i],
		M->L->m[i],M->L->s[i],M->L->F[i],M->L->eta[i],
		M->L->dm[i],M->L->ds[i],M->L->dF[i],M->L->deta[i]);

    free(envelope_X);
    free(envelope_Y);
    gsl_spline_free(np_model);
    gsl_interp_accel_reset(npacc);
    gsl_interp_accel_reset(acc);
  }

/* Free the splines used here. */
  gsl_spline_free(data_spline);
  gsl_spline_free(error_spline);
  gsl_interp_accel_free(acc);
}
  
  
