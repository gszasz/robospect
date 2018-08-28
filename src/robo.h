/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#define __USE_GNU
#define _GNU_SOURCE

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <time.h>
#include <sysexits.h>
#include "config.h"


#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics.h>

#if HAVE_LIBPLPLOTD
#include <plplot/plplot.h>
#endif

#if HAVE_LIBCFITSIO
#ifdef ODD_CFITSIO
#include <cfitsio/fitsio.h>
#else
#include <fitsio.h>
#endif
#endif

#if __APPLE__

#endif

#ifndef __ROBO_H__
#define __ROBO_H__


#ifndef isfinite
#  define isfinite(x) \
  (sizeof (x) == sizeof (float)                                            \
    ? __finitef (x)                                                        \
    : sizeof (x) == sizeof (double)                                        \
     ? __finite (x) : __finitel (x))
#endif

/* BEGIN LINE FLAG BITMASKS */
#define ROBO_MAX_ITER         0x000001     /* Solver hit maximum number of iterations without returning a fit within tolerance. */
#define ROBO_FIT_FAIL         0x000002     /* Solver returned impossible value and aborted further computation. */
#define ROBO_FIT_IGNORED      0x000004     /* Line rejected and ignored from consideration due to concerns with fit parameters. */
#define ROBO_FIT_REFUSED      0x000008     /* Line rejected and refused from consideration due to lack of chi^2 improvement with inclusion. */
#define ROBO_FIT_LARGE_SHIFT  0x000010     /* Fit mean shifted significantly from expected initial value. */
#define ROBO_FIT_RECENTER     0x000020     /* Fit mean reset to initial value and refit with fixed mean. */
#define ROBO_BAD_WAVELENGTH   0x000040     /* Suspected bad wavelength solution around this line. */
#define ROBO_FIX_WAVELENGTH   0x000080     /* Automatically corrected supplied line peak to fit bad wavelength solution. */

#define ROBO_FIT_BLEND        0x000100     /* Line believed to be part of a blend. */
#define ROBO_FIT_DEBLEND      0x000200     /* Line solution based on deblending modeling. */

#define ROBO_ALT_REFUSED      0x001000     /* Alternate line model rejected and refused due to lack of chi^2 improvement with inclusion. */
#define ROBO_BAD_ERROR_VAL    0x002000     /* Error value out of range and ignored. */   
#define ROBO_RESET            0x0000a0     /* Flags to reset on each solve iteration. */
#define ROBO_BAD_FIT          0x00000f     /* Generic concerns with line.  Not used in final model. */
#define ROBO_INIT             0x010000     /* This signifies we haven't done anything to the line yet. */
#define ROBO_FUNC_GAUSSIAN    0x100000     /* Line fit using a Gaussian model. */
#define ROBO_FUNC_LORENTZIAN  0x200000     /* Line fit using a Lorentzian model. */
#define ROBO_FUNC_HUMLICEK    0x400000
#define ROBO_FUNC_HJERTING    0x800000
#define ROBO_FUNC_NONPARAM    0x020000     /* Not happy that this is here.  Need to extend bitmask length? */
/* END LINE FLAG BITMASKS */

/* Verbosity information bitmasks. */
#define ROBO_VERBOSE_DEFAULT   0x0001
#define ROBO_VERBOSE_IO        0x0002
#define ROBO_VERBOSE_ID        0x0004
#define ROBO_VERBOSE_INTERMEDIATE 0x0008
#define ROBO_VERBOSE_LINE      0x0010
#define ROBO_VERBOSE_CONTINUUM 0x0020
#define ROBO_VERBOSE_NOISE     0x0040
#define ROBO_VERBOSE_MATH      0x0080
#define ROBO_VERBOSE_SCREEN    0x0100
#define ROBO_VERBOSE_DEBUG     0x0200
#define ROBO_VERBOSE_XCORR     0x0400
#define ROBO_VERBOSE_ALL       0xffff

#define DATE_FORMAT "%a %b %d %H:%M:%S %Y"
#define SIGMA_RANGE 8.0
#define ERROR_MAX_SN 10
#define MULTIFUNCTION 1

#define MIN_CONTINUUM_COUNT 25
#define MAX_CONTINUUM_COUNT 500


typedef struct {
  double psf_width;
  double continuum_box;
  double deblend_radius;
  double deblend_ratio;
  int    deblend_iterations;
  double tolerance;
  double relax;
  
  long max_iterations;
  long iteration;

  int    find_lines;
  double find_sigma;
  double wavelength_min_error;
  double wavelength_max_error;
  int    wavelength_limit;
  double radial_velocity;
  int    measure_radial_velocity;
  double rv_range;
  double rv_max_error;
  int    rv_steps;
  double rv_sigma;
  
  int save_temp;
  int plot_all;
  int help;
  int verbose;
  int fault;
  
  char *infilename;
  double min_x;
  double max_x;
  double delta_x;
  int order;
  int max_order;
  int supplied_errors;
  
  char *path_base;
  char *line_list_filename;

  char *line_model_name;
  char *continuum_model_name;
  char *noise_model_name;
  char *deblend_model_name;
  char *function_model_name;
  
  int strict;
  double strict_center;
  double strict_width;
  double strict_flux;
  int flux_calibrated;
  int fits_IO;
  int fits_row;
  
  FILE *log;
  char *command_name;
  char *command_line;
} opts;

typedef struct {
  double *m;
  double *s;
  double *F;
  double *eta;
  double *chi;
  int    Niter;
  
  double *dm;
  double *ds;
  double *dF;
  double *deta;
  
  double *mp;
  double *sp;
  double *Fp;
  double *etap;
  
  double *x0;

  char   **comment;
  int    *manual;
  int    *flags;
  int    *blend_group;
  long b;
  long l;
} lines;

typedef struct {
  double *x;
  double *y;
  double *e;
  int    *order;

  double *yO;
  long N;
  int max_order;
  opts *O;
} data;

typedef struct {
  double *x;
  double *continuum;
  double *lines;
  double *alternate;
  long N;

  int continuum_model;
  int line_model;
  int noise_model;
  int deblend_model;
  int function_model;
  
  lines *L;
  opts *O;
} model;

typedef struct {
  long N;
  double mean;
  double sigma;
  double min;
  double max;
  double med;
  double MAD;
} stats;

typedef enum {
  CONTINUUM_BOXCAR,
  CONTINUUM_LOGBOXCAR,
  CONTINUUM_ROBUST_LINEAR,
  CONTINUUM_PENALIZED_SPLINE,
  CONTINUUM_BSPLINE,
  CONTINUUM_HISTOGRAM,
  CONTINUUM_FFT,
  CONTINUUM_PEAK,
  CONTINUUM_BLACKBODY,
  CONTINUUM_POWERLAW,
  CONTINUUM_DEVEL,
  CONTINUUM_NULL
} continuum_model;

typedef enum {
  LINE_GAUSSIAN,
  LINE_FIXEDMEAN_GAUSSIAN,
  LINE_FIXEDMEAN_ISOLATED_GAUSSIAN,
  LINE_ADAPTIVE_MM,
  LINE_PRE,
  LINE_PRE_ORIGINAL,
  LINE_TWO_STAGE,
  LINE_TWO_STICKY,
  LINE_THREE_STAGE,
  LINE_NONPARAMETRIC,
  LINE_NULL
} line_model;

typedef enum {
  NOISE_NULL,
  NOISE_SUPPLIED,
  NOISE_POISSON,
  NOISE_SLOW_BOXCAR,
  NOISE_BOXCAR
} noise_model;

typedef enum {
  DEBLEND_NULL,
  DEBLEND_NLLS
} deblend_model;

typedef enum {
  FUNCTION_GAUSSIAN,
  FUNCTION_LORENTZIAN,
  FUNCTION_PSEUDOVOIGT,
  FUNCTION_HJERTING,
  FUNCTION_HUMLICEK,
  FUNCTION_PLANCK,
  FUNCTION_NONPARAMETRIC,
  FUNCTION_SKEWGAUSS
} function_model;

/* config.c */
opts *alloc_options();
void free_options (opts *O);
opts *set_options_from_command_line(int argc, char *argv[]);
void set_options(opts *O, int argc, char *argv[]);
void usage_block(opts *O);
char * get_time();
/* fileio.c */
data *alloc_data(int size); 
void realloc_data(data *D, int size);
void free_data(data *D);
data *read_spectrum_data(opts *options, int *N);
data *read_data_ascii(opts *options, int *N); /* O */
#if HAVE_LIBCFITSIO
/* Quick macro to do error reporting. */
#define FRE fits_report_error(stderr,status); if (status != 0) { fprintf(stderr, "  AT %s:%d\n",__FILE__,__LINE__); } ; status = 0;

data *read_data_fits(opts *options, int *N);

typedef struct {
  double *x0;
  double *xd;
  int dispersion_axis;
  int N;
} fits_WCS;

fits_WCS *read_fits_WCS(fitsfile *ff);
double   calculate_fits_WCS(fits_WCS *WCS,int x, int y);
void     free_fits_WCS(fits_WCS *WCS);

#endif

lines *alloc_lines(int size);
void realloc_lines(lines *L, int size);
void free_lines(lines *L);
lines *append_lines(lines *L, double x0);
lines *merge_lines(lines *A, lines *B);
lines *read_line_list(opts *options, int *N);

char *generate_filename(opts *options, char *rule);
void write_spectra_data(opts *options, data *D, model *M);
void write_line_data(opts *options,  data *D, model *M);
#if HAVE_LIBPLPLOTD
void render_line_plots(opts *options,  data *D, model *M);
#endif
void log_comment(opts *options, int level, char *format, ...);

void validate_line_peaks(opts *options, data *D, model *M);
void fit_peak_offsets(double *x, double *dw, model *M, double max);
/* linefinder.c */
void linefinder(opts *options, data *D, model *M);
lines *linefinder_naive(opts *options, data *D, model *M);
lines *linefinder_with_prior(opts *options, data *D, model *M);

/* rv.c */
void measure_radial_velocity(opts *options, data *D, model *M);

/* math.c */
double gaussian(double x, double m, double s, double F);
stats *array_stats(double *x, int N);
stats *array_stats_safe(double *x, int N);
stats *histogram_and_gaussfit(double *x, int N);
/*double gaussian_line(double x, lines *L, int i); */
double gaussian_line_alt(double x, lines *L, int i);
double equivalent_width(lines *L, int i);
double equivalent_width_alt(lines *L, int i);
double equivalent_width_error(lines *L, int i);
void vector_convolve(double *x1, double *x2, double *conv, int N);
double *make_psf(double s, int N);
void robust_linear_fit(double *x, double *y, int N, double *m, double *b, double tolerance);
double interp(double x0, double *X, double *Y, int N);
/*double lorentzian_line(double x, lines *L, int i); */
/*double pseudovoigt_line(double x, lines *L, int i); */
/*double hjerting_line(double x, lines *L, int i); */
/* nlls-gaussfit.c */
/* inline double f(double x, double m, double s, double A); */
/* inline double dfdA(double x, double m, double s, double A); */
/* inline double dfdm(double x, double m, double s, double A); */
/* inline double dfds(double x, double m, double s, double A); */
/* Function to be fit, in this case, a Gaussian.  Also, the derivatives */
/* of the function with respect to all fit parameters.                  */
/* inline double gaussian(double x, double m, double s, double F) { */
/*   return(f_gauss(x,m,s,F)); */
/* } */

int invert3x3(double **in, double **out);
  
int gaussfit(double *X, double *Y, double *E, int N,
	     double *m, double *s, double *A,
	     double *dm, double *ds, double *dA,
	     int vm, int vs, int vA, double relax,
	     double *chi, double tolerance, int max_iter);

/* nlls-multigauss.c */
int multigauss(double *X, double *Y, double *E, int N,
	       double *m, double *s, double *A, int M,
	       double *dm, double *ds, double *dA,
	       int *vm, int *vs, int *vA,
	       double relax,
	       double *chi, double tolerance, int max_iter);
/* math */
int multifunction(double *X, double *Y, double *E, int N,
		  int fit_function,
		  double *m, double *s, double *A, double *eta, int M,
		  double *dm, double *ds, double *dA, double *deta,
		  int *vm, int *vs, int *vA, int *veta,
		  double relax,
		  double *chi, double tolerance, int max_iter);
double gaussian_line(double x, lines *L, int i);
void gaussian_model(double x, double m, double s, double A, double eta,
		    double *f,
		    double *dfdm, double *dfds, double *dfdA, double *dfdeta);
double lorentzian_line(double x, lines *L, int i);
void lorentzian_model(double x, double m, double s, double A, double eta,
		      double *f,
		      double *dfdm, double *dfds, double *dfdA, double *dfdeta);
double pseudovoigt_line(double x, lines *L, int i);
void pseudovoigt_model(double x, double m, double s, double A, double eta,
			 double *f,
			 double *dfdm, double *dfds, double *dfdA, double *dfdeta);
double hjerting_line(double x, lines *L, int i);
void hjerting_model(double x, double m, double s, double A, double eta,
			 double *f,
			 double *dfdm, double *dfds, double *dfdA, double *dfdeta);

double humlicek_line(double x, lines *L, int i);
void humlicek_model(double x, double m, double s, double A, double eta,
		    double *f,
		    double *dfdm, double *dfds, double *dfdA, double *dfdeta);
		   

double skewgauss_line(double x, lines *L, int i);
void skewgauss_model(double x, double m, double s, double A, double eta,
		     double *f,
		     double *dfdm, double *dfds, double *dfdA, double *dfdeta);
		   
void planck_model(double x, double m, double s, double A, double eta,
		  double *f,
		  double *dfdm, double *dfds, double *dfdA, double *dfdeta);
  
  
  
  
/* gsllm-multigauss.c */
int lm_multigauss(double *X, double *Y, double *E, int N,
		  double *m, double *s, double *A, int M,
		  double *dm, double *ds, double *dA,
		  int *vm, int *vs, int *vA,
		  double relax,
		  double *chi, double tolerance, int max_iter);

/* noise.c */
void generate_noise_model(opts *options, data *D, model *M);

/* models.c */
#define MAX_ITERATIONS 100
#define LARGE_SHIFT    2.0

model *alloc_models(opts *options, data *D);
void free_models(model *M);
void generate_model_continuum(opts *options, data *D, model *M);
void generate_model_continuum_boxcar(opts *options, data *D, model *M);
void generate_model_continuum_logboxcar(opts *options, data *D, model *M);
void generate_model_continuum_peak_boxcar(opts *options, data *D, model *M);
void generate_model_continuum_histogram(opts *options, data *D, model *M);
void generate_model_continuum_robust_linear(opts *options, data *D, model *M);
void generate_model_continuum_penalized_spline(opts *options, data *D, model *M);
void generate_model_continuum_bspline(opts *options, data *D, model *M);
void generate_model_continuum_fft(opts *options, data *D, model *M);
void generate_model_continuum_blackbody(opts *options, data *D, model *M);
void generate_model_continuum_powerlaw(opts *options, data *D, model *M);
void generate_model_continuum_devel(opts *options, data *D, model *M);
void generate_model_continuum_null(opts *options, data *D, model *M);

void measure_lines(data *D, model *M);
void chi_check_lines(data *D, model *M);
void set_flags(data *D, model *M);
void measure_lines_gaussian(data *D, model *M);
void measure_lines_fixedmean_gaussian(data *D, model *M);
void measure_lines_fixedmean_isolation_gaussian(data *D, model *M);
void measure_lines_adaptive_mm(data *D, model *M);
void measure_lines_NULL(data *D, model *M);
void measure_lines_PRE(data *D, model *M);
void measure_lines_PRE_original(data *D, model *M);
void measure_lines_twostage(data *D, model *M);
void measure_lines_twostage_sticky(data *D, model *M);
void measure_lines_threestage(data *D, model *M);

void measure_lines_nonparametric(data *D, model *M);

void deblend_lines_nlls(data *D, model *M);

/* spectra.c */
double *continuum_normalized_spectra(data *D, model *M);
double *line_removed_spectra(data *D, model *M);
double *continuum_normalized_line_removed_spectra(data *D, model *M);
double *continuum_normalized_Zvalue(data *D, model *M);
void set_line_fit_spectra(data *D, model *M);
void set_continuum_fit_spectra(data *D, model *M);
void set_noise_fit_spectra(data *D, model *M);
void set_linefinder_Z_spectra(data *D, model *M);		

/* vectors.c */
double *vector_add(double *A, double *B, int N);
double *vector_subtract(double *A, double *B, int N);
double *vector_subtract_constant(double *A, double B, int N);
double *vector_multiply(double *A, double *B, int N);
double *vector_divide(double *A, double *B, int N);
double *vector_constant(double V, int N);
double *vector_copy(double *A, int N);


#endif
