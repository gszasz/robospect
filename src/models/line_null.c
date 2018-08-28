/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void measure_lines_NULL(data *D, model *M) {
  int j;
  log_comment(D->O,ROBO_VERBOSE_LINE,
	      "mlN: Setting line flux to 0.0");
  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }
}
