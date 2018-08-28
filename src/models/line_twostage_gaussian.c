/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void measure_lines_twostage(data *D, model *M) {
  int i,j;
  double tolerance = D->O->tolerance;
  double relax     = D->O->relax;
  double v = 0.0;
  int iterations;
  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }

  /* Pre-measure the lines from the data and an estimate of the fwhm */
  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }

  measure_lines_PRE(D,M);

  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }
  
  /* Loop over each known line. */
  for (i = 0; i < M->L->l; i++) {
    M->L->dm[i] = 0; M->L->ds[i] = 0; M->L->dF[i] = 0;
    M->L->chi[i] = 0;
    iterations = gaussfit(D->x,D->y,D->e,D->N,
			  &(M->L->m[i]),&(M->L->s[i]),&(M->L->F[i]),
			  &(M->L->dm[i]),&(M->L->ds[i]),&(M->L->dF[i]),
			  1,1,1,relax,
			  &(M->L->chi[i]),tolerance,MAX_ITERATIONS);
    log_comment(D->O,ROBO_VERBOSE_LINE,
		"ml2S: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f)",
		i,M->L->l,iterations,M->L->chi[i],M->L->x0[i],
		M->L->m[i],M->L->s[i],M->L->F[i],
		M->L->dm[i],M->L->ds[i],M->L->dF[i]);
	    
    /* Set flags */
    if ((!isfinite(M->L->m[i]))||
	(!isfinite(M->L->s[i]))||
	(!isfinite(M->L->F[i]))) {
      M->L->flags[i] |= ROBO_FIT_FAIL;
    }
    if (iterations >= MAX_ITERATIONS) {
      M->L->flags[i] |= ROBO_MAX_ITER;
    }
    if (fabs(M->L->m[i] - M->L->mp[i]) > LARGE_SHIFT * M->L->sp[i]) {
      M->L->flags[i] |= ROBO_FIT_LARGE_SHIFT;
    }
    /* End flags */
    /* Add to line model */
    if (((fabs(M->L->m[i] - M->L->x0[i]) < 100))&&
	(((fabs(M->L->s[i]) < 0.1)||
	  ((fabs(M->L->F[i]) < 10)&&(fabs(M->L->F[i] / M->L->s[i]) > fabs(M->L->s[i])))))) {
      for (j = 0; j < M->N; j++) {
	v = gaussian_line(M->x[j],M->L,i);
	if (v == v) {
	  M->lines[j] += v;
	}
      }      
    }
    else {
      M->L->flags[i] |= ROBO_FIT_IGNORED;
    }
    /* End line model addition */
  }
}
