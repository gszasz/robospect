/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void measure_lines_fixedmean_gaussian(data *D, model *M) {
  int i,j;
  double dm,ds,dF;
  double m1,s1;
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
    dm = 0; ds = 0; dF = 0; M->L->chi[i] = 0;
    m1 = M->L->m[i]; s1 = M->L->s[i]; 
    iterations = gaussfit(D->x,D->y,D->e,D->N,
			  &(M->L->m[i]),&(M->L->s[i]),&(M->L->F[i]),
			  &dm,&ds,&dF,
			  0,1,1,relax,
			  &(M->L->chi[i]),tolerance,MAX_ITERATIONS);
    fprintf(stderr,"mlfmG: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f)\n",
	    i,M->L->l,iterations,M->L->chi[i],M->L->x0[i],M->L->m[i],M->L->s[i],M->L->F[i],
	    dm,ds,dF);
    
    /* Set flags */
    if ((!isfinite(M->L->m[i]))||
	(!isfinite(M->L->s[i]))||
	(!isfinite(M->L->F[i]))) {
      M->L->flags[i] |= ROBO_FIT_FAIL;
    }
    if (iterations >= MAX_ITERATIONS) {
      M->L->flags[i] |= ROBO_MAX_ITER;
    }
    if (fabs(M->L->m[i] - m1) > LARGE_SHIFT * s1) {
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
