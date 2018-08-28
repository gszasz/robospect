/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void measure_lines_twostage_sticky(data *D, model *M) {
  int i;
  double tolerance = D->O->tolerance;
  double relax     = D->O->relax;
  int iterations;

  int vm = 1, vs = 1, vF = 1;
  /* Loop over each known line. */
  for (i = 0; i < M->L->l; i++) {
    vm = 1;
    M->L->dm[i] = 0; M->L->ds[i] = 0; M->L->dF[i] = 0;
    M->L->chi[i] = 0;
    iterations = multigauss(D->x,D->y,D->e,D->N,
			    &(M->L->m[i]),&(M->L->s[i]),&(M->L->F[i]),1,
			  &(M->L->dm[i]),&(M->L->ds[i]),&(M->L->dF[i]),
			  &vm,&vs,&vF,relax,
			  &(M->L->chi[i]),tolerance,100);
    log_comment(D->O,ROBO_VERBOSE_LINE,"ml2Ss: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f)",
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
      vm = 0;
      iterations = multigauss(D->x,D->y,D->e,D->N,
			      &(M->L->m[i]),&(M->L->s[i]),&(M->L->F[i]),1,
			      &(M->L->dm[i]),&(M->L->ds[i]),&(M->L->dF[i]),
			      &vm,&vs,&vF,relax,
			      &(M->L->chi[i]),tolerance,100);
      /* Set flags */
      if (iterations >= MAX_ITERATIONS) {
	M->L->flags[i] |= ROBO_MAX_ITER;
      }
      if (fabs(M->L->m[i] - M->L->mp[i]) > LARGE_SHIFT * M->L->sp[i]) {
	M->L->flags[i] |= ROBO_FIT_LARGE_SHIFT;
      }
      /* End flags */

      log_comment(D->O,ROBO_VERBOSE_LINE,"mlF: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f)",
		  i,M->L->l,iterations,M->L->chi[i],M->L->x0[i],
		  M->L->m[i],M->L->s[i],M->L->F[i],
		  M->L->dm[i],M->L->ds[i],M->L->dF[i]);

    }
  }
}
