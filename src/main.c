/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
/* This space for versioning comments. */
/* v2.16 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>

#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_blas.h>

#include "robo.h"

int main(int argc, char *argv[]) {
  data *D;
  model *M;
  opts *options = set_options_from_command_line(argc,argv);

  int i;
  int Nsamples;

  if (options->fault != 0) {
    usage_block(options);
    free_options(options);
    return(options->fault);
  }

  /* ********************* */
  /* SETUP DATA AND MODELS */
  /* ********************* */
  for (options->order = options->fits_row; options->order <= options->max_order; options->order++) {
    /* Load the data, and make a copy of the original data so we can work in place. */
    D = read_spectrum_data(options,&Nsamples);
    log_comment(options,ROBO_VERBOSE_ALL,"main: Datafile %s read. %d samples.",options->infilename,D->N);
    /* Allocate the data holding the models. */
    M = alloc_models(options,D);
    log_comment(options,ROBO_VERBOSE_ALL,"main: Models allocated.");
    /* Load line list file. */
    if (options->line_list_filename) {
      M->L = read_line_list(options,&i);
      validate_line_peaks(options,D,M);
      if (options->measure_radial_velocity) {
	measure_radial_velocity(options,D,M);
      }
    }
    else {
      M->L = NULL;
    }
    
    D->O = options;
    M->O = options;
    /* ********* */
    /* MAIN LOOP */
    /* ********* */
    for (options->iteration = 0; options->iteration < options->max_iterations; options->iteration++) {
      log_comment(options,ROBO_VERBOSE_ALL,"main: Begin iteration: %d/%d %s",options->iteration + 1,options->max_iterations,get_time());
      /* Smooth the data to generate a continuum model. */
      log_comment(options,ROBO_VERBOSE_ALL,"main: i%d Model continuum",options->iteration + 1);
      set_continuum_fit_spectra(D,M);
      generate_model_continuum(options,D,M);
      
      /* Make an estimate of the noise by scanning the data and saving the binned sigma */
      log_comment(options,ROBO_VERBOSE_ALL,"main: i%d Model noise",options->iteration + 1);
      set_noise_fit_spectra(D,M);
      generate_noise_model(options,D,M);
      
      /* Optionally attempt to find the lines. */
      if (options->find_lines & 1) {
	log_comment(options,ROBO_VERBOSE_ALL,"main: i%d Finding lines",options->iteration + 1);
	set_linefinder_Z_spectra(D,M);
	linefinder(options,D,M);
      }
      
      /* Fit lines */
      set_line_fit_spectra(D,M);
      
      log_comment(options,ROBO_VERBOSE_ALL,"main: i%d Measure lines",options->iteration + 1);
      measure_lines(D,M);
      
      if (options->save_temp == 1) {
	/* Save intermediate data. */
	write_line_data(options,D,M);
	write_spectra_data(options,D,M);
#if HAVE_LIBPLPLOTD
	render_line_plots(options,D,M);
#endif
      }

      log_comment(options,ROBO_VERBOSE_ALL,"main: i%d End iteration",options->iteration + 1);
    }
    
    /* *********** */
    /* SAVE OUTPUT */
    /* *********** */
    
    /* Write out to files named something like our input file. */
    if (strcmp(options->infilename,"/dev/stdin") != 0) {
      log_comment(options,ROBO_VERBOSE_ALL,"main: Writing output files");
      write_line_data(options,D,M);
      write_spectra_data(options,D,M);
#if HAVE_LIBPLPLOTD
      render_line_plots(options,D,M);
#endif
    }

    /* **** */
    /* FREE */
    /* **** */
    free_data(D);
    free_models(M);
  }
  log_comment(options,ROBO_VERBOSE_ALL,"main: Program complete");
  free_options(options);
  return(0);
}

