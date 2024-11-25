/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

#if HAVE_LIBPLPLOTD
void render_line_plots(opts *options, data *D, model *M) {
  int i,j;
  double *x;
  double *data;
  double *model;
  double *alt_model;
  double *continuum;
  int Nsamples;
  int Noffset;

  double x_min;
  double x_max;
  double y_min;
  double y_max;
  double y_continuum = 1.0;
  double y_data = 0.0;

  char outfile[1024];
  char label[256];
  double xline[2];
  double yline[2];

  /* Load and setup the output terminal. */
  plsdev("psc");
  /* plspage(72,72,700,1000,39,30); */
  plspage(72,72,612,792,0,0);
  plsetopt("portrait","3");

  /* Fix the stupid way they set the default colors */
  plscol0(0,255,255,255);
  plscol0(1,0,0,0);
  plscol0(2,255,0,0);
  plscol0(3,0,255,0);
  plscol0(4,0,0,255);

  snprintf(outfile,sizeof(outfile),"%s",generate_filename(options,"PLOT"));
  plsfnam(outfile);

  plinit();

  plssub(2,5);
  
  for (i = 0; i < M->L->l; i++) {
    if ((M->L->x0[i] < D->x[0])) {
      continue;
    }
    if (M->L->x0[i] > D->x[D->N - 1]) {
      continue;
    }
    if (!(M->L->manual[i] || options->plot_all)) {
      continue;
    }

    /* x-limits are easy. */
    x_min = M->L->x0[i] - 5 * M->L->s[i];
    x_max = M->L->x0[i] + 5 * M->L->s[i];

    /* Set y-limits. */
    /* Find the continuum level and the data level at the line center. */
    for (j = 0; j < D->N - 1; j++) {
      if ((D->x[j] <= M->L->m[i])&&(D->x[j+1] >= M->L->m[i])) {
	y_continuum = M->continuum[j];
	y_data      = D->yO[j];
	break;
      }
    }

    if (y_data < y_continuum) {
      y_min = 0.8 * y_data;
      y_max = 1.2 * y_continuum;
    }
    else {
      y_min = 0.8 * y_continuum;
      y_max = 1.2 * y_data;
    }
    
    xline[0] = M->L->x0[i];
    xline[1] = M->L->x0[i];
    yline[0] = y_min;
    yline[1] = y_max;    
    
    if ((!isfinite(x_min))||(!isfinite(x_max))||(!isfinite(y_min))||(!isfinite(y_max))) {
      continue;
    }

    /* Draw box */
    plcol0(1);
    plenv(x_min,x_max,y_min,y_max,0,0);
    pllab("wavelength","flux","");

    /* Label of best fit parameters. */
    snprintf(label,sizeof(label),"m = %.4f",
	     M->L->m[i]);
    plptex(0.05 * (x_max - x_min) + x_min,
	   0.25 * (y_max - y_min) + y_min,
	   9999,
	   0,
	   0.0,
	   label);
    snprintf(label,sizeof(label),"s = %.3f",
	     fabs(M->L->s[i]));
    plptex(0.05 * (x_max - x_min) + x_min,
	   0.20 * (y_max - y_min) + y_min,
	   9999,
	   0,
	   0.0,
	   label);
    if (options->flux_calibrated == 0) {
      snprintf(label,sizeof(label),"EQ = %.2f",
	       (equivalent_width(M->L,i)));
    }
    else {
      snprintf(label,sizeof(label),"Flux = %.6g",
	       M->L->F[i] * interp(M->L->m[i],D->x,M->continuum,D->N));
    }      
    plptex(0.05 * (x_max - x_min) + x_min,
	   0.15 * (y_max - y_min) + y_min,
	   9999,
	   0,
	   0.0,
	   label);
    snprintf(label,sizeof(label),"chi = %.3f",
	     (M->L->chi[i] > 1e6) ?
	     NAN :
	     M->L->chi[i]);
    plptex(0.05 * (x_max - x_min) + x_min,
	   0.10 * (y_max - y_min) + y_min,
	   9999,
	   0,
	   0.0,
	   label);
    snprintf(label,sizeof(label),"mP = %.4f",
	     M->L->mp[i]);
    plptex(0.80 * (x_max - x_min) + x_min,
	   0.25 * (y_max - y_min) + y_min,
	   9999,
	   0,
	   0.0,
	   label);
    snprintf(label,sizeof(label),"sP = %.3f",
	     fabs(M->L->sp[i]));
    plptex(0.80 * (x_max - x_min) + x_min,
	   0.20 * (y_max - y_min) + y_min,
	   9999,
	   0,
	   0.0,
	   label);
    snprintf(label,sizeof(label),"EQ = %.2f",
	     fabs(equivalent_width_alt(M->L,i)));
    plptex(0.80 * (x_max - x_min) + x_min,
	   0.15 * (y_max - y_min) + y_min,
	   9999,
	   0,
	   0.0,
	   label);
    snprintf(label,sizeof(label),"flags = 0x%06x",M->L->flags[i]);
    plptex(0.80 * (x_max - x_min) + x_min,
	   0.10 * (y_max - y_min) + y_min,
	   9999,
	   0,
	   0.0,
	   label);
    if (M->L->comment[i]) {
      snprintf(label,sizeof(label),"%s",
	       M->L->comment[i]);
      plptex(0.05 * (x_max - x_min) + x_min,
	     0.05 * (y_max - y_min) + y_min,
	     9999,
	     0,
	     0.0,
	     label);
    }      
    /* Dashed vertical line at line center */
    plcol0(4);
    pllsty(2);
    plline(2,xline,yline);

    if (options->verbose & ROBO_VERBOSE_SCREEN) {
      plcol0(5);
      xline[0] = M->L->mp[i] - options->deblend_radius * M->L->sp[i];
      xline[1] = M->L->mp[i] - options->deblend_radius * M->L->sp[i];
      yline[0] = y_min;
      yline[1] = y_max;    
      plline(2,xline,yline);
      xline[0] = M->L->mp[i] + options->deblend_radius * M->L->sp[i];
      xline[1] = M->L->mp[i] + options->deblend_radius * M->L->sp[i];
      yline[0] = y_min;
      yline[1] = y_max;    
      plline(2,xline,yline);
      plwidth(1.0);
    }
    plcol0(1);
    pllsty(1);
    
    /* Determine what data to plot in this box */
    Noffset = 0;
    Nsamples = 0;
    for (j = 0; j < D->N; j++) {
      if (Noffset > 0) {
	Nsamples++;
      }
      else {
	if (D->x[j] >= x_min) {
	  Noffset = j;
	}
      }
      if (D->x[j] >= x_max) {
	continue;
      }
    }
    x = malloc(sizeof(double) * Nsamples);
    if (x == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    data = malloc(sizeof(double) * Nsamples);
    if (data == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    model = malloc(sizeof(double) * Nsamples);
    if (model == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    alt_model = malloc(sizeof(double) * Nsamples);
    if (alt_model == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    continuum = malloc(sizeof(double) * Nsamples);
    if (continuum == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    for (j = 0; j < Nsamples; j++) {
      x[j] = D->x[j + Noffset];
      data[j] = D->yO[j + Noffset];
      if (options->flux_calibrated == 0) {
	model[j] = (1 + M->lines[j + Noffset]) * M->continuum[j + Noffset];
	alt_model[j] = (1 + M->alternate[j + Noffset]) * M->continuum[j + Noffset];
      }
      else {
	model[j] = (M->lines[j + Noffset]) + M->continuum[j + Noffset];
	alt_model[j] = (M->alternate[j + Noffset]) + M->continuum[j + Noffset];
      }
      continuum[j] = M->continuum[j + Noffset];

    }	       

    /* Plot data and model */
    plcol0(3);    
    plline(Nsamples,x,data);
    plcol0(6);
    plline(Nsamples,x,alt_model);
    plcol0(2);
    plline(Nsamples,x,model);
    plcol0(4);
    plline(Nsamples,x,continuum);
    
    free(x);
    free(data);
    free(model);
    free(alt_model);
    free(continuum);
  }
  
  plend();
}
#endif
