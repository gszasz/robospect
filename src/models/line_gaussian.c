/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void measure_lines_gaussian(data *D, model *M) {
  int i,j;
  double dm,ds,dF;
  double tolerance = D->O->tolerance;
  double relax     = D->O->relax;
  double v = 0.0;
  int iterations;
  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }
  
  /* Loop over each known line. */
  for (i = 0; i < M->L->l; i++) {
    dm = 0; ds = 0; dF = 0; M->L->chi[i] = 0;
    M->L->s[i] = 0.08;
    M->L->F[i] = -0.05;
    iterations = gaussfit(D->x,D->y,D->e,D->N,
			  &(M->L->m[i]),&(M->L->s[i]),&(M->L->F[i]),
			  &dm,&ds,&dF,
			  1,1,1,relax,
			  &(M->L->chi[i]),tolerance,100);
    fprintf(stderr,"mlG: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f)\n",
	    i,M->L->l,iterations,M->L->chi[i],M->L->x0[i],M->L->m[i],M->L->s[i],M->L->F[i],
	    dm,ds,dF);

    if ((fabs(M->L->m[i] - M->L->x0[i]) < 100)&&
	(fabs(M->L->s[i]) < 50)) 
      {
      for (j = 0; j < M->N; j++) {
	v = gaussian_line(M->x[j],M->L,i);
	if (v == v) {
	  M->lines[j] += v;
	}
      }
    }
  }
}
