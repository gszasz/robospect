AM_CFLAGS = -Wall -pedantic -Wunreachable-code -march=native -msse2 -mmmx -mfpmath=sse -falign-loops -falign-functions -maccumulate-outgoing-args -time -ggdb3 -I/opt/local/include -I/usr/include/plplot  -O3
VPATH = math:fileio:models

I_SRC = input_data.c input_lines.c input_fitsWCS.c
O_SRC = output_lines.c output_plots.c output_spectra.c output_logs.c
IO_SRC = ${I_SRC} ${O_SRC}
C_SRC = continuum_boxcar.c continuum_logboxcar.c continuum_bspline.c continuum_null.c continuum_penalized_spline.c continuum_robust_linear.c continuum_histogram.c continuum_fft.c continuum_peak_boxcar.c continuum_blackbody.c continuum_devel.c continuum_powerlaw.c
L_SRC = line_fixedmean_gaussian.c line_gaussian.c line_isolated_gaussian.c line_null.c line_pre.c line_threestage_gaussian.c line_twostage_gaussian.c line_twostage_sticky_gaussian.c line_nonparametric.c
D_SRC = deblend_nlls.c
M_SRC = gaussian_model.c lorentzian_model.c pseudovoigt_model.c hjerting_model.c planck_model.c nlls-multifunction.c humlicek_model.c skewgauss_model.c
# There are test functions.
# zaghloul_model.c abrarov_zaghloul.c
MODEL_SRC = ${C_SRC} ${L_SRC} ${D_SRC} ${M_SRC}
SRC   = main.c config.c fileio.c math.c noise.c model.c vectors.c spectra.c nlls-gaussfit.c linefinder.c  nlls-multigauss.c gsllm-multigauss.c rv.c ${IO_SRC} ${MODEL_SRC}


bin_PROGRAMS = robospect
robospect_SOURCES = robo.h config.h $(SRC)