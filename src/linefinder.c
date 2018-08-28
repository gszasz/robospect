/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "robo.h"

#define AUTO_COMMENT_LENGTH 80
#define AUTO_COMMENT_FORMAT "Auto-found on iteration %ld at %.2f with %.2f sigma significance."
#define AUTO_COMMENT_VALUES options->iteration,L->x0[j],fabs(y_peak)

void linefinder(opts *options, data *D, model *M) {
  lines *L;
  if (options->find_lines & 1) {
    if ((options->find_lines & 2)&&(M->L)) {
      log_comment(options,ROBO_VERBOSE_ID,"lf: Using Prior version");
      L = linefinder_with_prior(options,D,M);
    }
    else {
      log_comment(options,ROBO_VERBOSE_ID,"lf: Using Naive version");
      L = linefinder_naive(options,D,M);
    }

    if (M->L) {
      free(M->L);
    }

    M->L = L;
  }
}

lines *linefinder_naive(opts *options, data *D, model *M) {
  lines *L;
  int i;
  int in_a_line = 0;
  int Nlines = 10;
  int j = 0;

  double x_peak = 99e99;
  double y_peak = 99e99;

  L = alloc_lines(Nlines);
  
  for (i = 0; i < D->N; i++) {
    log_comment(options,ROBO_VERBOSE_DEBUG,"lfn: %g %g",D->x[i],D->y[i]);
    if (in_a_line == 0) {
      if (fabs(D->y[i]) > options->find_sigma) {
	in_a_line = 1;
	x_peak = D->x[i];
	y_peak = D->y[i];
      }
    }
    if (in_a_line == 1) {
      if (fabs(D->y[i]) <= options->find_sigma) {
	L->x0[j] = x_peak;
	L->m[j]  = x_peak;
	L->s[j]  = 1.0;
	L->F[j]  = -0.05;
	L->eta[j] = 0.0;
	L->manual[j] = 0;
	L->comment[j] = malloc(sizeof(char) * AUTO_COMMENT_LENGTH);
	snprintf(L->comment[j],AUTO_COMMENT_LENGTH,AUTO_COMMENT_FORMAT,AUTO_COMMENT_VALUES);
	L->flags[j] = ROBO_INIT;
	L->blend_group[j] = 0;
	j++;
	if (j >= Nlines) {
	  Nlines += 10;
	  realloc_lines(L,Nlines);
	}
	in_a_line = 0;
      }
      else {
	if (fabs(D->y[i]) > fabs(y_peak)) {
	  x_peak = D->x[i];
	  y_peak = D->y[i];
	}
      }
    }
  }

  Nlines = j;
  realloc_lines(L,Nlines);
  log_comment(options,ROBO_VERBOSE_ID,"lfn: Found %d lines",Nlines);
  return(L);
}

lines *linefinder_with_prior(opts *options, data *D, model *M) {
  lines *L;
  int i;
  int in_a_line = 0;
  int Nlines = 10;
  int j = 0;
  int l = 0;
  
  double x_peak = 99e99;
  double y_peak = 99e99;

  L = alloc_lines(Nlines);
  for (i = 0; i < D->N; i++) {
    if (in_a_line == 0) {
      if (fabs(D->y[i]) >= options->find_sigma) {
	in_a_line = 1;
	x_peak = D->x[i];
	y_peak = D->y[i];
      }
      else if (((D->x[i] >= M->L->x0[l]))&&(l < M->L->l)) {
	L->x0[j] = M->L->x0[l];
	L->m[j] = M->L->m[l];
	L->s[j] = M->L->s[l];
	L->F[j] = M->L->F[l];
	L->eta[j] = M->L->eta[l];
	L->comment[j] = M->L->comment[l];
	L->manual[j]  = M->L->manual[l];
	L->flags[j] = M->L->flags[l];
	L->blend_group[j] = 0;
	j++;
	l++;
	if (j >= Nlines) {
	  Nlines += 10;
	  realloc_lines(L,Nlines);
	}
      }	
    }
    if (in_a_line == 1) {
      if (fabs(D->y[i]) <= options->find_sigma) {
	if ((D->x[i] > M->L->x0[l])&&(l < M->L->l)) {
	  L->x0[j] = M->L->x0[l];
	  L->m[j] = M->L->m[l];
	  L->s[j] = M->L->s[l];
	  L->F[j] = M->L->F[l];
	  L->eta[j] = M->L->eta[l];
	  L->comment[j] = M->L->comment[l];
	  L->manual[j]  = M->L->manual[l];
	  L->flags[j] = M->L->flags[l];
	  L->blend_group[j] = 0;
	  j++;
	  l++;
	}
	else {
	  L->x0[j] = x_peak;
	  L->m[j]  = x_peak;
	  L->s[j]  = 1.0;
	  L->F[j]  = -0.05;
	  L->eta[j] = 0.0;
	  L->manual[j] = 0;
	  L->comment[j] = malloc(sizeof(char) * AUTO_COMMENT_LENGTH);
	  snprintf(L->comment[j],AUTO_COMMENT_LENGTH,AUTO_COMMENT_FORMAT,AUTO_COMMENT_VALUES);
	  L->flags[j] = ROBO_INIT;
	  L->blend_group[j] = 0;
	  j++;
	}
	if (j >= Nlines) {
	  Nlines += 10;
	  realloc_lines(L,Nlines);
	}
	in_a_line = 0;
      }
      else {
	if (fabs(D->y[i]) > fabs(y_peak)) {
	  x_peak = D->x[i];
	  y_peak = D->y[i];
	}
      }
    }
  }

  Nlines = j;
  realloc_lines(L,Nlines);
  log_comment(options,ROBO_VERBOSE_ID,"lfp: Found %d lines (%d in prior)",Nlines,M->L->l);
  return(L);
}
