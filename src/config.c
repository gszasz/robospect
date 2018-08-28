/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <time.h>
#include "robo.h"

/* Set up a default option structure. */
opts *alloc_options() {
  opts *O = malloc(sizeof(opts));

  O->psf_width = 5.0;                /* to remove */
  O->continuum_box = 40.0;
  O->deblend_radius = 4.0;           /* to remove */
  O->deblend_ratio  = 0.1;           /* to remove */
  O->deblend_iterations = 3;           /* to remove */
  O->tolerance      = 1e-3;
  O->relax          = 1.0;           /* replace with better handling */
  
  O->max_iterations = 1;
  O->iteration      = 1;

  O->order          = 0;
  O->max_order      = 0;

  O->line_list_filename = NULL;
  O->find_lines     = 0;
  O->find_sigma     = 3.0;
  O->wavelength_min_error = 0.1;
  O->wavelength_max_error = 0.5;
  O->wavelength_limit = 2000;
  O->radial_velocity = 0.0;
  O->measure_radial_velocity = 0;
  O->rv_range       = 300;
  O->rv_max_error   = 1e-2;
  O->rv_steps       = 100;
  O->rv_sigma       = 10.0;
    
  O->save_temp      = 0;
  O->plot_all       = 0;
  O->help           = 0;
  O->verbose        = ROBO_VERBOSE_DEFAULT;
  O->fault          = 0;           /* to remove */

  O->strict         = 1;           /* replace with better handling */
  O->strict_center  = 100.0;
  O->strict_width   = 0.5;
  O->strict_flux    = 10.0;
  O->flux_calibrated = 0;
  O->fits_IO        = 0;
  O->fits_row       = 0;

  O->min_x          = 0;
  O->max_x          = 0;
  O->log            = stderr;
  O->command_name   = NULL;
  O->command_line   = NULL;
  O->continuum_model_name = NULL;
  O->line_model_name = NULL;
  O->noise_model_name = NULL;
  O->deblend_model_name = NULL;
  O->function_model_name = NULL;

  return(O);
}

void free_options (opts *O) {
  if (O->log != stderr) {
    log_comment(O,ROBO_VERBOSE_ALL,"# %s END",get_time());
    if (O->log) {
      fclose(O->log);
    }
    O->log = stderr;
  }
  if (O->command_name) {
    free(O->command_name);
  }
  if (O->command_line) {
    free(O->command_line);
  }
  free(O);
}

opts *set_options_from_command_line(int argc, char *argv[]) {
  opts *O = alloc_options();
  set_options(O,argc,argv);
  return(O);
}

void set_verbosity(opts *O, char *arg) {
  if (strcasecmp(arg,"io") == 0) {
    O->verbose |= ROBO_VERBOSE_IO;
  }
  else if (strcasecmp(arg,"id") == 0) {
    O->verbose |= ROBO_VERBOSE_ID;
  }
  else if (strcasecmp(arg,"intermediate") == 0) {
    O->verbose |= ROBO_VERBOSE_INTERMEDIATE;
  }
  else if (strcasecmp(arg,"line") == 0) {
    O->verbose |= ROBO_VERBOSE_LINE;
  }
  else if (strcasecmp(arg,"continuum") == 0) {
    O->verbose |= ROBO_VERBOSE_CONTINUUM;
  }
  else if (strcasecmp(arg,"noise") == 0) {
    O->verbose |= ROBO_VERBOSE_NOISE;
  }
  else if (strcasecmp(arg,"math") == 0) {
    O->verbose |= ROBO_VERBOSE_MATH;
  }
  else if (strcasecmp(arg,"screen") == 0) {
    O->verbose |= ROBO_VERBOSE_SCREEN;
  }
  else if (strcasecmp(arg,"debug") == 0) {
    O->verbose |= ROBO_VERBOSE_DEBUG;
  }
  else if (strcasecmp(arg,"xcorr") == 0) {
    O->verbose |= ROBO_VERBOSE_XCORR;
  }
  else if (strcasecmp(arg,"all") == 0) {
    O->verbose |= ROBO_VERBOSE_ALL;
  }
  else if (strcasecmp(arg,"silent") == 0) {
    O->verbose = 0;
  }
  else {
    O->verbose |= ROBO_VERBOSE_DEFAULT;
    fprintf(stderr,"Unknown verbosity option: %s. Stop.\n",arg);
  }
}  

void usage_block(opts *O) {
  fprintf(stderr,
	  "Usage: %s <options> <input.spectrum>\n",O->command_name);
  /*       1234567891123456789212345678931234567894123456789512345678961234567898*/
  fprintf(stderr,
	  "  -h, --help                     This help                            (%s)\n",
	  O->help ? "found" : "");

  /*       1234567891123456789212345678931234567894123456789512345678961234567898*/
  fprintf(stderr,
	  "\n                        FITTING PARAMETERS\n");
  fprintf(stderr,
	  "  -V, --continuum_box <N_AA>     Width of box used for continuum\n"
	  "                                  and noise estimates.                (%f)\n",
	  O->continuum_box);
  fprintf(stderr,
	  "  -r, --deblend_radius <N_sig>   Number of a line's stdev to look\n"
	  "                                  for possibly blended lines.         (%f)\n",
	  O->deblend_radius);
  fprintf(stderr,
	  "  -R, --deblend_ratio <F>        Fraction of a line's peak to bother\n"
	  "                                  with possibly blended lines.        (%f)\n",
	  O->deblend_ratio);
  fprintf(stderr,
	  "  -d, --deblend_iterations <N>   Number of times to attempt to find\n"
	  "                                  neighboring lines.                  (%d)\n",
	  O->deblend_iterations);
  fprintf(stderr,
	  "  -T, --tolerance <eps>          Tolerance value for fitting.         (%g)\n",
	  O->tolerance);

  /*       1234567891123456789212345678931234567894123456789512345678961234567898*/
  fprintf(stderr,
	  "\n                        MODEL SELECTION\n");
  fprintf(stderr,
	  "  -C, --continuum_model <NAME>   Select continuum model to use:       (%s)\n",
	  O->continuum_model_name ? O->continuum_model_name : "boxcar");
  fprintf(stderr,
	  "        boxcar (default)            Simple boxcar model\n"
	  "        logboxcar                   log10 boxcar model\n"
	  "        peak                        Fit boxcar across maxima\n"
	  "        blackbody                   Fit blackbody curve to continuum\n"
	  "        powerlaw                    Fit powerlaw to continuum\n"
	  "        histogram                   Fit Gaussian to data histogram\n"); /* to remove */
  fprintf(stderr,
	  "        fft                         Use FFT frequency-space filtering\n" /* to remove */
	  "        null                        Do no fit, assume continuum = 1.0\n"
	  "        devel                       Development continuum model\n"       /* to remove */
	  "        bspline                     Devel: better with discontinuities\n"); /* to remove */
  fprintf(stderr,
	  "  -N, --noise_model <NAME>       Select noise model to use:           (%s)\n",
	  O->noise_model_name ? O->noise_model_name : "boxcar");
  fprintf(stderr,
	  "        boxcar (default)            Simple boxcar model\n"
	  "        slow                        Slower two-pass boxcar model\n"
	  "        Poisson                     Poissonian noise model\n"
	  "        null                        Do no fit, inherit from continuum fit.\n"); 
  fprintf(stderr,
	  "  -M, --line_model <NAME>        Select line model to use:            (%s)\n",
	  O->line_model_name ? O->line_model_name : "best");
  fprintf(stderr,
	  "        null                        Do no fit, assume lines = 0.0\n"
	  "        gauss                       Simple NLLS Gaussian fitter\n"        
	  "        fixed                       Simple NLLS Gaussian, fixed mean\n"        
	  "        nofit                       Do not fit lines\n"                        
	  "        pre                         Estimate lines from FWHM\n");
  fprintf(stderr,
	  "        nonparametric               Estimate lines from non-parametric model\n"
	  "        twostage                    Pre + Gauss\n"                             /* to remove */
	  "        twosticky                   Pre + Gauss + fixed\n"                     /* to remove */
	  "        best (default)              Pre + automatic deblending (best)\n");     
  fprintf(stderr,
	  "  -Q, --function_model <NAME>    Select line function to use:         (%s)\n",
	  O->function_model_name ? O->function_model_name : "Gaussian");
  fprintf(stderr,
	  "        Gaussian                    Gaussian model.\n"
	  "        Lorentzian                  Lorentzian model (experimental).\n"
	  "        Hjerting                    Voigt approximation (experimental).\n"
	  "        Humlicek                    Voigt approximation (experimental).\n"
	  "        Nonparametric               Non-parametric model (experimental).\n"
	  "        SkewGauss                   Skew-Gaussian model  (experimental).\n"
	  );
  fprintf(stderr,
	  "  -D, --deblend_model <NAME>     Select deblending model to use:      (%s)\n", /* to remove */
	  O->deblend_model_name ? O->deblend_model_name : "null");
  fprintf(stderr,
	  "        null (default)              Do not deblend lines.\n"
	  "        nlls                        NLLS Gaussian method.\n");

  /*       1234567891123456789212345678931234567894123456789512345678961234567898*/
  fprintf(stderr,
	  "\n                        LINE IDENTIFICATION\n");
  fprintf(stderr,
	  "  -L, --line_list <file>         Fit the lines specified in this file (%s)\n",
	  O->line_list_filename ? O->line_list_filename : "");
  fprintf(stderr,
	  "  -F, --naive_find_lines         Attempt to find additional lines to\n"
	  "                                  fit. Uses sigma threshold.          (%s)\n",
	  (O->find_lines & 1) ? "found" : "");
  fprintf(stderr,
	  "  -f, --find_sigma <s>           Set the threshold for line finding.  (%f)\n",
	  O->find_sigma);
  fprintf(stderr,
	  "      --wavelength_min_error <dw>Set the limits for acceptable error \n"
	  "                                  in the peak finding code.           (%f)\n",
	  O->wavelength_min_error);
  fprintf(stderr,
	  "      --wavelength_max_error <dw>Set the limits for acceptable error \n"
	  "                                  in the peak finding code.           (%f)\n",
	  O->wavelength_max_error);
  fprintf(stderr,
	  "      --wavelength_limit <N>     Set the limit for the number of lines\n"
	  "                                  with bad wavelengths before the\n"
	  "                                  code attempts to correct.           (%d)\n",
	  O->wavelength_limit);
  fprintf(stderr,
	  "      --radial_velocity <v>      Apply a known radial velocity\n"
	  "                                  correction to wavelengths in units\n"
	  "                                  of (wavelength_unit * 1e13)/sec.    (%f)\n",
	  O->radial_velocity);
  fprintf(stderr,
	  "      --measure_radial_velocity  Attempt to match lines to linelist to\n"
	  "                                  measure the radial velocity.        (%d)\n",
	  O->measure_radial_velocity);
  fprintf(stderr,
	  "      --radial_velocity_range    Search for radial velocities +/- this\n"
	  "                                  value.                              (%f)\n",
	  O->rv_range);
  fprintf(stderr,
	  "      --radial_velocity_error    Set the maximum rv error allowed.    (%f)\n",
	  O->rv_max_error);
  fprintf(stderr,
	  "      --radial_velocity_step     Set the number of correlation steps. (%d)\n",
	  O->rv_steps);
  fprintf(stderr,
	  "      --radial_velocity_sigma    Set the detection threshold for RV.  (%f)\n",
	  O->rv_sigma);
  /*       1234567891123456789212345678931234567894123456789512345678961234567898*/
  fprintf(stderr,
	  "\n                        PROGRAM CONTROL\n");
  fprintf(stderr,
	  "  -i, --iterations <i>           Number of continuum/line iterations  (%ld)\n",
	  O->max_iterations);
  fprintf(stderr,
	  "      --loosen                   Loosen strict line fit checks.       (%d)\n"
	  "                                  Use --strict_ options to set check\n"
	  "                                  limits for individual checks.\n",
	  (O->strict == 1) ? 0 : 1);
  fprintf(stderr,
	  "      --strict_center            Set strict center constraint on      (%f)\n"
	  "                                  |m_fit - m_in| check.\n",
	  O->strict_center);
  fprintf(stderr,
	  "      --strict_width             Set strict width limit on |s_fit|.   (%f)\n",
	  O->strict_width);
  fprintf(stderr,
	  "      --strict_flux              Set strict flux limit on |F_fit|.    (%f)\n",
	  O->strict_flux);
  fprintf(stderr,
	  "      --flux_calibrated          Input is in flux units(erg/s/cm^2/A) (%d)\n",
	  O->flux_calibrated);
  fprintf(stderr,
	  "      --fits_row <N>             Fit row N of 2-d fits spectrum.      (%d)\n",
	  O->fits_row);
  fprintf(stderr,
	  "  -I, --save_temp                Save intermediate results            (%s)\n",
	  O->save_temp ? "found" : "");
  fprintf(stderr,
	  "  -A, --plot_all                 Default plot only contains lines\n"
	  "                                  from the line list. This option\n"
	  "                                  allows found lines in the plot.     (%s)\n",
	  O->plot_all ? "found" : "");
  fprintf(stderr,
	  "  -P, --path_base <output_root>  Output file rootname.                (%s)\n",
	  O->path_base ? O->path_base : "");
  fprintf(stderr,
  	  "  -v, --verbose <verbose_level>  Change verbosity level:              (0x%04x)\n"
	  "                                  A logfile is used for most levels.\n",
  	  O->verbose);
  fprintf(stderr,
	  "        silent                      No messages              (%s)\n"
  	  "        default                     Normal messages          (%s)\n"
  	  "        io                          File IO                  (%s)\n"
	  "        xcorr                       Radial Velocity          (%s)\n"
  	  "        id                          Line identification      (%s)\n",
	  (O->verbose == 0) ? "found" : "",
  	  (O->verbose & ROBO_VERBOSE_DEFAULT) ? "found" : "",
  	  (O->verbose & ROBO_VERBOSE_IO) ? "found" : "",
	  (O->verbose & ROBO_VERBOSE_XCORR) ? "found" : "",
  	  (O->verbose & ROBO_VERBOSE_ID) ? "found" : "");
  fprintf(stderr,
  	  "        line                        Line fitting             (%s)\n"
  	  "        continuum                   Continuum fitting        (%s)\n"
  	  "        noise                       Noise fitting            (%s)\n",
	  (O->verbose & ROBO_VERBOSE_LINE) ? "found" : "",
  	  (O->verbose & ROBO_VERBOSE_CONTINUUM) ? "found" : "",
  	  (O->verbose & ROBO_VERBOSE_NOISE) ? "found" : "");
  fprintf(stderr,
  	  "        math                        Math operations          (%s)\n"
  	  "        screen                      Print messages to screen (%s)\n"
  	  "        debug                       Code debug messages      (%s)\n"
  	  "        all                         Enable all messages      (%s)\n",
  	  (O->verbose & ROBO_VERBOSE_MATH) ? "found" : "",
  	  (O->verbose & ROBO_VERBOSE_SCREEN) ? "found" : "",
  	  (O->verbose & ROBO_VERBOSE_DEBUG) ? "found" : "",
	  (O->verbose == ROBO_VERBOSE_ALL) ? "found" : "");
  fprintf(stderr,"\n");
#if HAVE_LIBPLPLOTD
  fprintf(stderr,
	  "Compiled with libplplotd support for plotting.\n");
#endif
#if HAVE_LIBCFITSIO
  fprintf(stderr,
	  "Compiled with libcfitsio support for fits io.\n");
#endif
  fprintf(stderr,
	  "\n"
	  "Report bugs to: watersc1@ifa.hawaii.edu\n"
	  "(C) 2010-2013 Chris Waters and Julie Hollek\n"
	  "This is free software; see the source for copying conditions.\n"
	  "There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
}


char * get_time() {
  time_t curtime = time(NULL);
  char *output = malloc(32 * sizeof(char));
  strftime(output,32,DATE_FORMAT,localtime(&curtime));
  return(output);
}

void set_options(opts *O, int argc, char *argv[]) {
  char c = 0;
  int i = 0;
  int j = 0;
  int opt_index;

  /* Rework this to use a better argument list for long options */
  static struct option long_options[] = {
    {"help", no_argument, 0, 'h'},
    {"psf_width", required_argument, 0, 'p'},
    {"continuum_box", required_argument, 0, 'V'},
    {"deblend_radius", required_argument, 0, 'r'},
    {"deblend_ratio", required_argument, 0, 'R'},
    {"deblend_iterations",required_argument, 0, 'd'},
    {"tolerance", required_argument, 0, 'T'},
    {"iterations", required_argument, 0, 'i'},
    {"save_temp", no_argument, 0, 'I'},
    {"continuum_model", required_argument, 0, 'C'},
    {"noise_model", required_argument, 0, 'N'},
    {"line_model", required_argument, 0, 'M'},
    {"deblend_model", required_argument, 0, 'D'},
    {"function_model", required_argument, 0, 'Q'},
    
    {"line_list", required_argument, 0, 'L'},
    {"naive_find_lines", no_argument, 0, 'F'},
    {"find_sigma", required_argument, 0, 'f'},
    {"wavelength_min_error",required_argument, 0, 'e'},
    {"wavelength_max_error",required_argument, 0, 'E'},
    {"wavelength_limit",required_argument, 0, 'w'},
    {"radial_velocity",required_argument,  0, '3'},
    {"measure_radial_velocity",no_argument, 0, '4'},
    {"radial_velocity_sigma", required_argument, 0, '8'},
    {"radial_velocity_range", required_argument, 0, '5'},
    {"radial_velocity_error", required_argument, 0, '6'},
    {"radial_velocity_step", required_argument, 0, '7'},
    {"loosen", no_argument, 0, 'l'},
    {"strict_center", required_argument, 0, 'z'},
    {"strict_width", required_argument, 0, 'x'},
    {"strict_flux", required_argument, 0, 'y'},
    {"flux_calibrated", no_argument, 0, '1'},
    {"fits_row",required_argument, 0, '2'},
    {"plot_all",   no_argument, 0, 'A'},
    {"path_base", required_argument, 0, 'P'},

    {"verbose", optional_argument, 0, 'v'},
    {0,0,0,0}
  };
  
  /* Store command line into the option struct. */
  i = 0;
  O->command_name = realloc(O->command_name,sizeof(char) * (i + strlen(argv[j]) + 2));
  O->command_name = strcat(O->command_name,argv[j]);  
  for (j = 0; j < argc; j++) {
    O->command_line = realloc(O->command_line,sizeof(char) * (i + strlen(argv[j]) + 2));
    i = i + strlen(argv[j]) + 2;
    if (j > 0) {
      O->command_line = strcat(O->command_line," ");
    }
    O->command_line = strcat(O->command_line,argv[j]);
  }
  /* Parse options. */
  while ((c = getopt_long(argc,argv,"hFf:e:E:w:p:V:r:R:T:i:d:L:C:N:M:D:v:P:lz:x:y:1AI2:Q:3:45:6:7:8:",
			  long_options,&opt_index)) != -1) {
    switch (c) {
    case 'h':
      O->help = 1;
      O->fault = 1;
      break;
    case 'p':
      O->psf_width = atof(optarg);
      break;
    case 'V':
      O->continuum_box = atof(optarg);
      break;
    case 'r':
      O->deblend_radius = atof(optarg);
      break;
    case 'R':
      O->deblend_ratio = atof(optarg);
      break;
    case 'd':
      O->deblend_iterations = atoi(optarg);
      break;
    case 'T':
      O->tolerance = atof(optarg);
      break;
    case 'v':
      if (optarg) {
	set_verbosity(O,optarg);
      }
      else {
      	O->verbose++;
      }
      break;
    case 'i':
      O->max_iterations = atoi(optarg);
      break;
    case 'I':
      O->save_temp = 1;
      break;
    case 'F':
      O->find_lines += 1;
      break;
    case 'f':
      O->find_sigma = atof(optarg);
      break;
    case 'e':
      O->wavelength_min_error = atof(optarg);
      break;
    case 'E':
      O->wavelength_max_error = atof(optarg);
      break;
    case 'w':
      O->wavelength_limit = atoi(optarg);
      break;
    case '3':
      O->radial_velocity = atof(optarg);
      break;
    case '4':
      O->measure_radial_velocity = 1;
      break;
    case '5':
      O->rv_range = atof(optarg);
      break;
    case '6':
      O->rv_max_error = atof(optarg);
      break;
    case '7':
      O->rv_steps = atoi(optarg);
      break;
    case '8':
      O->rv_sigma = atof(optarg);
      break;
    case 'l':
      O->strict = 0;
      break;
    case 'z':
      O->strict_center = atof(optarg);
      break;
    case 'x':
      O->strict_width = atof(optarg);
      break;
    case 'y':
      O->strict_flux = atof(optarg);
      break;
    case '1':
      O->flux_calibrated = 1;
      break;
    case '2':
      O->fits_row = atoi(optarg);
      O->max_order = atoi(optarg);
      break;
    case 'A':
      O->plot_all = 1;
      break;
    case 'L':
      O->find_lines += 2;
      O->line_list_filename = optarg;
      break;
    case 'C':
      O->continuum_model_name = optarg;
      break;
    case 'N':
      O->noise_model_name = optarg;
      break;
    case 'M':
      O->line_model_name = optarg;
      break;
    case 'D':
      O->deblend_model_name = optarg;
      break;
    case 'Q':
      O->function_model_name = optarg;
      break;
    case 'P':
      O->path_base = optarg;
      break;
    case '?':
      if (isprint(optopt)) {
	fprintf(stderr, "Unknown option `-%c'.\n",optopt);
      }
      else {
	fprintf(stderr, "Unknown option character `\\x%x'.\n",optopt);
      }
      O->fault = EX_USAGE;
      break;
    default:
      fprintf(stderr, "Unknown option `-%c'.\n",optopt);
      O->fault = EX_USAGE;
      break;
    }
  }
  /* Handle odd cases of the user not specifying what they want.  Not that I'm looking at you as I write this, JULIE! */
  if (O->find_lines == 0) {
    O->help = 1;
    fprintf(stderr,
	    "!! You must specify either an input linelist (-L) or ask for lines to be found (-F).\n"
	    );
    O->fault = EX_NOINPUT;
  }
  
  /* Load input file or read from stdin */
  if (argv[optind] == NULL) {
    O->infilename = malloc(32 * sizeof(char));
    snprintf(O->infilename,32,"/dev/stdin");
  }
  else {
    O->infilename = argv[optind];
  }

  /* If verbosity is more than zero, open a logfile. */
  if (O->verbose > ROBO_VERBOSE_DEFAULT) {
    O->log = fopen(generate_filename(O,"LOG"),"w");
    log_comment(O,ROBO_VERBOSE_ALL,"# %s %s",get_time(),O->command_line);
  }
}

