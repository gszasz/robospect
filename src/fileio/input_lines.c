/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

lines *alloc_lines(int size) {
  lines *L = malloc(sizeof(lines));
  if (L == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->l = size;
  L->b = 1;
  L->m = malloc(size * sizeof(double));
  if (L->m == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->s = malloc(size * sizeof(double));
  if (L->s == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->F = malloc(size * sizeof(double));
  if (L->F == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->eta = malloc(size * sizeof(double));
  if (L->eta == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  
  L->dm = malloc(size * sizeof(double));
  if (L->dm == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->ds = malloc(size * sizeof(double));
  if (L->ds == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->dF = malloc(size * sizeof(double));
  if (L->dF == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->deta = malloc(size * sizeof(double));
  if (L->deta == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  
  L->mp = malloc(size * sizeof(double));
  if (L->mp == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->sp = malloc(size * sizeof(double));
  if (L->sp == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->Fp = malloc(size * sizeof(double));
  if (L->Fp == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->etap = malloc(size * sizeof(double));
  if (L->etap == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  
  L->x0 = malloc(size * sizeof(double));
  if (L->x0 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->chi = malloc(size * sizeof(double));
  if (L->chi == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  L->comment = malloc(size * sizeof(char **));
  if (L->comment == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->manual  = malloc(size * sizeof(int));
  if (L->manual == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->flags   = malloc(size * sizeof(int));
  if (L->flags == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->blend_group = malloc(size * sizeof(int));
  if (L->blend_group == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->x0[0] = 0.0;
  return(L);  
}

void realloc_lines(lines *L, int size) {
  L->l = size;
  L->x0 = realloc(L->x0,size * sizeof(double));
  if (L->x0 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->chi = realloc(L->chi,size * sizeof(double));
  if (L->chi == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->m = realloc(L->m,size * sizeof(double));
  if (L->m == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->s = realloc(L->s,size * sizeof(double));
  if (L->s == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->F = realloc(L->F,size * sizeof(double));
  if (L->F == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->eta = realloc(L->eta,size * sizeof(double));
  if (L->eta == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->dm = realloc(L->dm,size * sizeof(double));
  if (L->dm == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->ds = realloc(L->ds,size * sizeof(double));
  if (L->ds == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->dF = realloc(L->dF,size * sizeof(double));
  if (L->dF == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->deta = realloc(L->deta,size * sizeof(double));
  if (L->deta == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->mp = realloc(L->mp,size * sizeof(double));
  if (L->mp == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->sp = realloc(L->sp,size * sizeof(double));
  if (L->sp == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->Fp = realloc(L->Fp,size * sizeof(double));
  if (L->Fp == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->etap = realloc(L->etap,size * sizeof(double));
  if (L->etap == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  L->comment = realloc(L->comment,size * sizeof(char **));
  if (L->comment == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->manual  = realloc(L->manual, size * sizeof(int));
  if (L->manual == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->flags   = realloc(L->flags,  size * sizeof(int));
  if (L->flags == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  L->blend_group = realloc(L->blend_group, size * sizeof(int));
  if (L->blend_group == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

void free_lines(lines *L) {
  free(L->x0);
  free(L->chi);
  free(L->m);
  free(L->s);
  free(L->F);
  free(L->eta);
  free(L->dm);
  free(L->ds);
  free(L->dF);
  free(L->deta);
  free(L->mp);
  free(L->sp);
  free(L->Fp);
  free(L->etap);
  free(L->comment);
  free(L->manual);
  free(L->flags);
  free(L->blend_group);
  free(L);
}

lines *merge_lines(lines *A, lines *B) {
  lines *L = alloc_lines(A->l + B->l);
  int a,b,i;

  a = 0;
  b = 0;
  for (i = 0; i < L->l; i++) {
    
    if ((a < A->l)&&(b < B->l)) {
      if (A->x0[a] < B->x0[b]) {
	L->x0[i]  = A->x0[a];
	L->chi[i] = A->chi[a];
	L->m[i]   = A->m[a];
	L->s[i]   = A->s[a];
	L->F[i]   = A->F[a];
	L->eta[i] = A->eta[a];
	L->dm[i]  = A->dm[a];
	L->ds[i]  = A->ds[a];
	L->dF[i]  = A->dF[a];
	L->deta[i] = A->deta[a];
	L->mp[i]  = A->mp[a];
	L->sp[i]  = A->sp[a];
	L->Fp[i]  = A->Fp[a];
	L->etap[i]  = A->etap[a];
	
	L->manual[i]      = A->manual[a];
	if (A->comment[a]) {
	  L->comment[i]     = strdup(A->comment[a]);
	}
	else {
	  L->comment[i]     = 0;
	}
	L->flags[i]       = A->flags[a];
	L->blend_group[i] = A->blend_group[a];
	a++;
      }
      else {
	L->x0[i]  = B->x0[b];
	L->chi[i] = B->chi[b];
	L->m[i]   = B->m[b];
	L->s[i]   = B->s[b];
	L->F[i]   = B->F[b];
	L->eta[i] = B->eta[b];
	L->dm[i]  = B->dm[b];
	L->ds[i]  = B->ds[b];
	L->dF[i]  = B->dF[b];
	L->deta[i] = B->deta[b];
	L->mp[i]  = B->mp[b];
	L->sp[i]  = B->sp[b];
	L->Fp[i]  = B->Fp[b];
	L->etap[i]  = B->etap[b];
	
	L->manual[i]      = B->manual[b];
	if (B->comment[b]) {
	  L->comment[i]     = strdup(B->comment[b]);
	}
	else {
	  L->comment[i]     = 0;
	}

	L->flags[i]       = B->flags[b];
	L->blend_group[i] = B->blend_group[b];
	b++;
      }
    }
    else if (a < A->l) {
      L->x0[i]  = A->x0[a];
      L->chi[i] = A->chi[a];
      L->m[i]   = A->m[a];
      L->s[i]   = A->s[a];
      L->F[i]   = A->F[a];
      L->eta[i] = A->eta[a];
      L->dm[i]  = A->dm[a];
      L->ds[i]  = A->ds[a];
      L->dF[i]  = A->dF[a];
      L->deta[i] = A->deta[a];
      L->mp[i]  = A->mp[a];
      L->sp[i]  = A->sp[a];
      L->Fp[i]  = A->Fp[a];
      L->etap[i]  = A->etap[a];
      
      L->manual[i]      = A->manual[a];
      if (A->comment[a]) {
	L->comment[i]     = strdup(A->comment[a]);
      }
      	else {
	  L->comment[i]     = 0;
	}

      L->flags[i]       = A->flags[a];
      L->blend_group[i] = A->blend_group[a];
      a++;
    }
    else if (b < B->l) {
      L->x0[i]  = B->x0[b];
      L->chi[i] = B->chi[b];
      L->m[i]   = B->m[b];
      L->s[i]   = B->s[b];
      L->F[i]   = B->F[b];
      L->eta[i] = B->eta[b];
      L->dm[i]  = B->dm[b];
      L->ds[i]  = B->ds[b];
      L->dF[i]  = B->dF[b];
      L->deta[i] = B->deta[b];
      L->mp[i]  = B->mp[b];
      L->sp[i]  = B->sp[b];
      L->Fp[i]  = B->Fp[b];
      L->etap[i]  = B->etap[b];
      
      L->manual[i]      = B->manual[b];
      if (B->comment[b]) {
	L->comment[i]     = strdup(B->comment[b]);
      }
      	else {
	  L->comment[i]     = 0;
	}

      L->flags[i]       = B->flags[b];
      L->blend_group[i] = B->blend_group[b];
      b++;
    }
  }
  free_lines(A);
  free_lines(B);
  return(L);
}

lines *append_lines(lines *L, double x0) {
  int N;
  if (!L) {
    L = alloc_lines(1);
  }
  else {
    realloc_lines(L,L->l + 1);
  }
  N = L->l - 1;
  L->x0[N] = x0;
  L->m[N]  = x0;
  L->s[N]  = 1.0;
  L->F[N]  = -0.05;
  L->dm[N] = 0.0;
  L->ds[N] = 0.0;
  L->dF[N] = 0.0;
  L->mp[N] = x0;
  L->sp[N] = 0.0;
  L->Fp[N] = 0.0;
  L->etap[N] = 0.0;
  L->eta[N] = 0.0;
  L->deta[N] = 0.0;
  L->chi[N] = 0.0;
  L->manual[N] = 0;
  L->flags[N] = ROBO_INIT;
  L->comment[N] = 0;
  L->blend_group[N] = 0;
  return(L);
}

lines *read_line_list(opts *options, int *N) {
  FILE *linesfile;
  char line[240];
  int i = 0;
  int v;
  lines *L;
  double X;
  *N = 10;
  L = alloc_lines(*N);
  linesfile = fopen(options->line_list_filename,"r");

  while(fgets(line,240,linesfile)) {
    for (v = 0; v < strlen(line); v++) {
      if ( line[v] == '\n' || line[v] == '\r' )
	line[v] = '\0';
    }
    v = sscanf(line,"%lf\n",&X);
    log_comment(options,ROBO_VERBOSE_DEBUG,"read_line_list: %d %d %f",i,*N,X);
    if ((X < options->min_x)||(X > options->max_x)) {
      continue;
    }
    if (v != -1) {
      L->x0[i] = X;
      L->m[i] = X;
      L->s[i] = 1.0;
      L->F[i] = -0.05;
      L->eta[i] = 0.0;
      L->dm[i] = 0.0;
      L->ds[i] = 0.0;
      L->dF[i] = 0.0;
      L->deta[i] = 0.0;
      L->mp[i] = X;
      L->sp[i] = 0.0;
      L->Fp[i] = 0.0;
      L->etap[i] = 0.0;

      L->manual[i] = 1;
      L->comment[i] = strdup(line);
      L->flags[i] = ROBO_INIT;
      L->blend_group[i] = 0;
      i++;
      if (i >= *N) {
	*N += 10;
	realloc_lines(L,*N);
      }
    }
  }

  *N = i;
  realloc_lines(L,*N);
  fclose(linesfile);
  log_comment(options,ROBO_VERBOSE_IO,"read_line_list: read %d lines",*N);
  
  return(L);
}

void validate_line_peaks(opts *options, data *D, model *M) {
  int i,j;
  int j_peak;
  double *x;
  double *dw;
  double R = 99e99;
  int    rx = 0;
  int    Nmax = 0;
  j_peak = 0;
  log_comment(options,ROBO_VERBOSE_DEBUG,"validating peaks");
  /* Scan for lines that may be offset */
  x = malloc(M->L->l * sizeof(double));
  if (x == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  dw = malloc(M->L->l * sizeof(double));
  if (dw == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  for (i = 0; i < M->L->l; i++) {
    /* Find the current peak bin. */
    for (j = j_peak; D->x[j] < M->L->x0[i]; j++) {
      j_peak = j;
      rx = j_peak;
      R = 99e99;
    }

    /* Look to the left */
    for (j = j_peak; j > 0; j--) {
      if (((D->y[j] - D->y[j - 1]) < 0)&&
	  ((D->y[j] - D->y[j + 1]) < 0)) {
	/* This is a peak */
	log_comment(options,ROBO_VERBOSE_DEBUG,"found left peak %g %g",M->L->x0[i],D->x[j]);
	rx = j;
	R = (D->y[rx] * D->x[rx] + D->y[rx-1] * D->x[rx-1] + D->y[rx+1] * D->x[rx+1]) /
	  (D->y[rx] + D->y[rx-1] + D->y[rx+1]);
	R -= M->L->x0[i];
	j = -1;
	
      }
    }
    /* Look to the right */
    for (j = j_peak; j < D->N; j++) {
      if (((D->y[j] - D->y[j - 1]) < 0)&&
	  ((D->y[j] - D->y[j + 1]) < 0)) {
	/* This is a peak */
	if (fabs(D->x[j] - M->L->x0[i]) < fabs(R)) {
	  log_comment(options,ROBO_VERBOSE_DEBUG,"found right peak %g %g",M->L->x0[i],D->x[j]);
	  rx = j;
	  R = (D->y[rx] * D->x[rx] + D->y[rx-1] * D->x[rx-1] + D->y[rx+1] * D->x[rx+1]) /
	    (D->y[rx] + D->y[rx-1] + D->y[rx+1]);
	  R -= M->L->x0[i];
	  j = D->N + 10;
	}
      }
    }
    /* Assign this set to the things to fit and fix */
    if ((rx != j_peak)&&
	(fabs(D->x[rx] - M->L->x0[i]) > options->wavelength_min_error)) {
      log_comment(options,ROBO_VERBOSE_DEBUG,"found offset peak %g %g %d",M->L->x0[i],D->x[rx],Nmax);
      M->L->flags[i] |= ROBO_BAD_WAVELENGTH;
      x[i] = M->L->x0[i];
      dw[i] = (D->y[rx] * D->x[rx] + D->y[rx-1] * D->x[rx-1] + D->y[rx+1] * D->x[rx+1]) /
	(D->y[rx] + D->y[rx-1] + D->y[rx+1]) - M->L->x0[i];
      Nmax++;
    }
    else {
      if (Nmax < options->wavelength_limit) {
	for (j = 0; j <= i; j++) {
	  M->L->flags[j] &= ~ROBO_BAD_WAVELENGTH;
	}
	Nmax = 0;
      }
      else {
	Nmax++;
      }
    }
  }
  fit_peak_offsets(x,dw,M,options->wavelength_max_error);
  free(x);
  free(dw);
}

void fit_peak_offsets(double *x, double *dw, model *M, double max) {
  double S = 0,Sx = 0,Sy = 0;
  double Sxx = 0,Syy = 0,Sxy = 0;
  double D = 0;
  double m = 0,b = 0;
  int i;

  for (i = 0; i < M->L->l; i++) {
    if (M->L->flags[i] & ROBO_BAD_WAVELENGTH) {
      S   += 1.0;
      Sx  += x[i];
      Sy  += dw[i];
      Sxx += x[i] * x[i];
      Syy += dw[i]*dw[i];
      Sxy += x[i] * dw[i];
    }
  }
  D = S * Sxx - Sx * Sx;
  
  if (D != 0) {
    m = (S * Sxy - Sx * Sy) / D;
    b = (Sy * Sxx - Sx * Sxy) / D;
    for (i = 0; i < M->L->l; i++) {
      if (M->L->flags[i] & ROBO_BAD_WAVELENGTH) {
	if (fabs(m * M->L->x0[i] + b) < max) {
	  M->L->x0[i] += m * M->L->x0[i] + b;
	  M->L->flags[i] |= ROBO_FIX_WAVELENGTH;
	}
      }
    }
  }
}
	
	
  
    
