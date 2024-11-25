/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"


void measure_lines_fixedmean_isolation_gaussian(data *D, model *M) {
  int i,j;
  double tolerance = D->O->tolerance;
  double relax     = D->O->relax;
  double v = 0.0;
  int iterations;

  double *y = malloc(sizeof(double) * D->N);
  if (y == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  
  /* Subtract the current model from the data. */
  for (j = 0; j < M->N; j++) {
    y[j] = D->y[j] - M->alternate[j];
  }
  
  /* Loop over each known line. */
  for (i = 0; i < M->L->l; i++) {
    /* Add this line's pre-model back into the data. */

    for (j = 0; j < M->N; j++) {
      if ((M->x[j] < M->L->mp[i] - SIGMA_RANGE * M->L->sp[i])|| 
	  (M->x[j] > M->L->mp[i] + SIGMA_RANGE * M->L->sp[i])) {
	continue;
      }
      v = gaussian_line_alt(M->x[j],M->L,i);
      if (v == v) {
	y[j] -= v;
      }
    }
    
    M->L->dm[i] = 0; M->L->ds[i] = 0; M->L->dF[i] = 0;
    M->L->chi[i] = 0;
    iterations = gaussfit(D->x,y,D->e,D->N,
			  &(M->L->m[i]),&(M->L->s[i]),&(M->L->F[i]),
			  &(M->L->dm[i]),&(M->L->ds[i]),&(M->L->dF[i]),
			  1,1,1,relax,
			  &(M->L->chi[i]),tolerance,100);

    log_comment(D->O,ROBO_VERBOSE_LINE,"mlI2Ss: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f)",
		i,M->L->l,iterations,M->L->chi[i],M->L->x0[i],
		M->L->m[i],M->L->s[i],M->L->F[i],
		M->L->dm[i],M->L->ds[i],M->L->dF[i]);
    /* Set flags */
    if (iterations >= MAX_ITERATIONS) {
      M->L->flags[i] |= ROBO_MAX_ITER;
    }
    if (fabs(M->L->m[i] - M->L->mp[i]) > LARGE_SHIFT * M->L->sp[i]) {
      M->L->flags[i] |= ROBO_FIT_LARGE_SHIFT;
    }
    /* End flags */
    if ((fabs((M->L->m[i] - M->L->mp[i]) / M->L->sp[i]) > 2.0)) {
      M->L->m[i] = M->L->mp[i];
      M->L->s[i] = M->L->sp[i];
      M->L->F[i] = M->L->Fp[i];
      M->L->flags[i] |= ROBO_FIT_RECENTER;
      iterations = gaussfit(D->x,y,D->e,D->N,
			    &(M->L->m[i]),&(M->L->s[i]),&(M->L->F[i]),
			    &(M->L->dm[i]),&(M->L->ds[i]),&(M->L->dF[i]),
			    0,1,1,relax,
			    &(M->L->chi[i]),tolerance,100);
      /* Set flags */
      if (iterations >= MAX_ITERATIONS) {
	M->L->flags[i] |= ROBO_MAX_ITER;
      }
      if (fabs(M->L->m[i] - M->L->mp[i]) > LARGE_SHIFT * M->L->sp[i]) {
	M->L->flags[i] |= ROBO_FIT_LARGE_SHIFT;
      }
      /* End flags */

      log_comment(D->O,ROBO_VERBOSE_LINE,"mlIF: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f)",
		  i,M->L->l,iterations,M->L->chi[i],M->L->x0[i],
		  M->L->m[i],M->L->s[i],M->L->F[i],
		  M->L->dm[i],M->L->ds[i],M->L->dF[i]);

    }
    for (j = 0; j < M->N; j++) {
      if ((M->x[j] < M->L->mp[i] - SIGMA_RANGE * M->L->sp[i])|| 
	  (M->x[j] > M->L->mp[i] + SIGMA_RANGE * M->L->sp[i])) {
	continue;
      }
      v = gaussian_line_alt(M->x[j],M->L,i);
      if (v == v) {
	y[j] += v;
      }
    }
  }
}
