/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void generate_model_continuum_null(opts *options, data *D, model *M) {
  int i;
  log_comment(options,ROBO_VERBOSE_CONTINUUM,
	      "gmcn: Setting continuum level to 1.0");
  for (i = 0; i < D->N; i += 1) {
    M->continuum[i] = 1.0;
  }
}
