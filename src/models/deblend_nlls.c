/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <ctype.h>

#include "../robo.h"

void deblend_lines_nlls(data *D, model *M) {
  int i,j,k;
  int B;
  double *m, *s, *F;
  double *dm,*ds,*dF;
  double chi;
  int *vm,*vs,*vF;
  double tolerance = 1e-3;
  int iterations;


  for (i = 1; i < M->L->b; i++) {
    B = 1;
    m = malloc(B * sizeof(double));
    if (m == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    s = malloc(B * sizeof(double));
    if (s == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    F = malloc(B * sizeof(double));
    if (F == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    dm= malloc(B * sizeof(double));
    if (dm == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    ds= malloc(B * sizeof(double));
    if (ds == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    dF= malloc(B * sizeof(double));
    if (dF == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    vm= malloc(B * sizeof(int));
    if (vm == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    vs= malloc(B * sizeof(int));
    if (vs == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    vF= malloc(B * sizeof(int));
    if (vF == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    	       
    for (j = 0; j < M->L->l; j++) {
      if (M->L->blend_group[j] == i) {
	m[B-1] = M->L->mp[j];
	s[B-1] = M->L->sp[j];
	F[B-1] = M->L->Fp[j];

	dm[B-1] = 0.0;
	ds[B-1] = 0.0;
	dF[B-1] = 0.0;

	vm[B-1] = 1;
	vs[B-1] = 1;
	vF[B-1] = 1;
	B++;
	
	m = realloc(m,B * sizeof(double));
	s = realloc(s,B * sizeof(double));
	F = realloc(F,B * sizeof(double));
	dm= realloc(dm,B * sizeof(double));
	ds= realloc(ds,B * sizeof(double));
	dF= realloc(dF,B * sizeof(double));
	vm= realloc(vm,B * sizeof(int));
	vs= realloc(vs,B * sizeof(int));
	vF= realloc(vF,B * sizeof(int));
      }
    }
    if (B - 1 < 1) {
      continue;
    }
    iterations = multigauss(D->x,D->y,D->e,D->N,
			    m,s,F,B - 1,
			    dm,ds,dF,
			    vm,vs,vF,
			    1.0,
			    &chi,
			    tolerance,100);
    for (j = 0; j < B; j++) {
      log_comment(D->O,ROBO_VERBOSE_LINE,
		  "dbN: (%d/%ld) (%d/%d) %d %f (%f,%f,%f) += (%f,%f,%f)",
		  i,M->L->b,
		  j,B,
		  iterations,
		  chi,
		  m[j],s[j],F[j],dm[j],ds[j],dF[j]);
    }
    k = 0;
    for (j = 0; j < M->L->l; j++) {
      if (M->L->blend_group[j] == i) {
	M->L->m[j] = m[k];
	M->L->s[j] = s[k];
	M->L->F[j] = F[k];
	M->L->dm[j] = dm[k];
	M->L->ds[j] = ds[k];
	M->L->dF[j] = dF[k];
	M->L->flags[j] |= ROBO_FIT_DEBLEND;
	k++;
      }
    }
  }
}
		      
