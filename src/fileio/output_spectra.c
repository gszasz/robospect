/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

void write_spectra_data(opts *options, data *D, model *M) {
  FILE *spectfile = fopen(generate_filename(options,"SPECTRA"),"w");
  int i;
  fprintf(spectfile,"## %s\n",options->command_line);
  fprintf(spectfile,"#wavelength flux err_est continuum lineflux model_spect fitflux \n");
  for (i = 0; i < D->N; i++) {
    fprintf(spectfile,"%5.8g %g %g %g %g %g %g\n",
	    D->x[i],D->yO[i],D->e[i],
	    M->continuum[i],M->lines[i],
	    options->flux_calibrated ?
	    M->continuum[i] + M->lines[i] :
	    M->continuum[i] * (1.0 + M->lines[i]),
	    D->y[i]);
  }
  fclose(spectfile);
}
