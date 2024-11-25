/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void measure_lines_threestage(data *D, model *M) {
  int i,j;
  double tolerance = D->O->tolerance;
  double relax     = D->O->relax;
  int iterations;

  /* Loop over each known line. */
  int B;
  double *m, *s, *F,*eta;
  double *dm,*ds,*dF,*deta;
  double chi;
  int *vm,*vs,*vF,*veta;

  for (j = 0; j < M->L->l; j++) {
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
    eta = malloc(B * sizeof(double));
    if (eta == NULL) {
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
    deta = malloc(B * sizeof(double));
    if (deta == NULL) {
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
    veta = malloc(B * sizeof(int));
    if (veta == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }

    m[B-1] = M->L->mp[j];
    s[B-1] = M->L->sp[j];
    F[B-1] = M->L->Fp[j];
    eta[B-1] = M->L->eta[j];
    
    dm[B-1] = 0.0;
    ds[B-1] = 0.0;
    dF[B-1] = 0.0;
    deta[B-1] = 0.0;
    
    vm[B-1] = 1;
    vs[B-1] = 1;
    vF[B-1] = 1;
    veta[B-1] = 1;
    if (M->L->blend_group[j] != 0) {
      for (i = j+1; i < M->L->l - 1; i++) {
      	log_comment(D->O,ROBO_VERBOSE_DEBUG,
		    "334: %d %d %d %d %d",
		    j,i,M->L->blend_group[i],M->L->blend_group[j],B);
	
	if (M->L->blend_group[i] == M->L->blend_group[j]) {
	  B++;
	  m = realloc(m,B * sizeof(double));
	  s = realloc(s,B * sizeof(double));
	  F = realloc(F,B * sizeof(double));
	  eta = realloc(eta,B * sizeof(double));
	  dm= realloc(dm,B * sizeof(double));
	  ds= realloc(ds,B * sizeof(double));
	  dF= realloc(dF,B * sizeof(double));
	  deta=realloc(deta,B * sizeof(double));
	  vm= realloc(vm,B * sizeof(int));
	  vs= realloc(vs,B * sizeof(int));
	  vF= realloc(vF,B * sizeof(int));
	  veta = realloc(veta,B * sizeof(int));
	  
	  m[B-1] = M->L->mp[i];
	  s[B-1] = M->L->sp[i];
	  F[B-1] = M->L->Fp[i];
	  eta[B-1] = M->L->eta[j];
	  
	  dm[B-1] = 0.0;
	  ds[B-1] = 0.0;
	  dF[B-1] = 0.0;
	  deta[B-1] = 0.0;
	  
	  vm[B-1] = 1;
	  vs[B-1] = 1;
	  vF[B-1] = 1;
	  veta[B-1] = 1;
	  log_comment(D->O,ROBO_VERBOSE_DEBUG,
		      "333: %d %d %d %d %d",
		      j,i,M->L->blend_group[i],M->L->blend_group[j],B);
	}
	else {
	  i = M->L->l + 1;
	}
      }
    }
#if MULTIFUNCTION
    iterations = multifunction(D->x,D->y,D->e,D->N,
    			       M->function_model,
    			       m,s,F,eta,B,
    			       dm,ds,dF,deta,
    			       vm,vs,vF,veta,
    			       relax,
    			       &chi,
    			       tolerance,100);
#else
    iterations = multigauss(D->x,D->y,D->e,D->N,
    			    m,s,F,B ,
    			    dm,ds,dF,
    			    vm,vs,vF,
    			    relax,
    			    &chi,
    			    tolerance,100);
#endif
    for (i = 0; i < B; i++) {
      log_comment(D->O,ROBO_VERBOSE_LINE,
		  "db3: (%d/%ld) (%d/%d) %d %f (%f,%f,%f,%f) += (%f,%f,%f,%f)",
		  j+i,M->L->l,
		  i+1,B,
		  iterations,
		  chi,
		  m[i],s[i],F[i],eta[i],dm[i],ds[i],dF[i],deta[i]);
      M->L->m[j+i] = m[i];
      M->L->s[j+i] = s[i];
      M->L->F[j+i] = F[i];
      M->L->eta[j+i] = eta[i];
      M->L->dm[j+i] = dm[i];
      M->L->ds[j+i] = ds[i];
      M->L->dF[j+i] = dF[i];
      M->L->deta[j+i] = deta[i];
      if (B > 1) {
	M->L->flags[j+i] |= ROBO_FIT_DEBLEND;
      }
      M->L->chi[j+i] = chi;
      if (iterations >= MAX_ITERATIONS) {
	M->L->flags[j+i] |= ROBO_MAX_ITER;
      }
    }
    free(m);
    free(s);
    free(F);
    free(eta);
    free(dm);
    free(ds);
    free(dF);
    free(deta);
    free(vm);
    free(vs);
    free(vF);
    free(veta);
    
    j = j + i - 1;
  }
}
