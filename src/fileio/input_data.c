/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

data *alloc_data(int size) {
  data *D = malloc(sizeof(data));
  if (D == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  D->N = size;
  D->x = malloc(size * sizeof(double));
  if (D->x == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  D->y = malloc(size * sizeof(double));
  if (D->y == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  D->e = malloc(size * sizeof(double));
  if (D->e == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  D->yO= malloc(size * sizeof(double));
  if (D->yO == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  return(D);
}

void realloc_data(data *D, int size) {
  D->N = size;
  D->x = realloc(D->x,size * sizeof(double));
  if (D->x == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  D->y = realloc(D->y,size * sizeof(double));
  if (D->y == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  D->e = realloc(D->e,size * sizeof(double));
  if (D->e == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  D->yO= realloc(D->yO,size * sizeof(double));
  if (D->yO == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
}

void free_data(data *D) {
  free(D->x);
  free(D->y);
  free(D->e);
  free(D->yO);
  free(D);
}

#define C_KM_S 299792.458
double radial_velocity(opts *options, double L) {
  return(L * (1.0 + options->radial_velocity/C_KM_S));
}

data *read_spectrum_data(opts *options, int *N) {
#if HAVE_LIBCFITSIO
  fitsfile *ff;
  int status;
  fits_open_file(&ff,options->infilename,READONLY,&status);
  if (status == 0) {
    fits_close_file(ff,&status);
    options->fits_IO = 1;
    return(read_data_fits(options,N));
  }
#endif
  return(read_data_ascii(options,N));
}

data *read_data_ascii(opts *options, int *N) {
  FILE *datafile;
  char line[240];
  int i = 0;
  int v;
  data *F;
  double X,Y;
  int order = 0;
  double E = 0;
  
  int last_was = 0;
  options->min_x = 99e99;
  options->max_x = -99e99;
  options->delta_x = 0;
  *N = 10;
  F = alloc_data(*N);
  options->supplied_errors = 0;
  datafile = fopen(options->infilename,"r");
  while(fgets(line,240,datafile)) {
    v = sscanf(line,"%lf %lf %lf %d",&X,&Y,&E,&order);
    if (v > 0) {
      X = radial_velocity(options,X);
      if (Y == 0.0) {
	continue;
      }
      if (v == 4) {
	if (order >= options->max_order) {
	  options->max_order = order;
	}
      }
      else {
	order = 0;
      }
      if (order == options->order) { /* Working order */
	F->x[i] = X;
	if (X < options->min_x) {
	  options->min_x = X;
	  last_was = 1;
	}
	if ((last_was == 1)&&(Y == 0.0)) {
	  options->min_x = X;
	  last_was = 1;
	  continue;
	}
	else {
	  last_was = 0;
	}
	if (X > options->max_x) {
	  options->max_x = X;
	  if (i > 0) {
	    options->delta_x = X - F->x[i-1];
	  }
	}
	if (v == 1) {
	  F->y[i]  = 0.0;
	  F->yO[i] = 0.0;
	  options->supplied_errors = 0;
	}
	if (v >= 2) {
	  F->y[i]  = Y;
	  F->yO[i] = Y;
	  options->supplied_errors = 0;
	}
	if (v >= 3) {
	  F->e[i] = E;
	  if ((E)&&(E != 0.0)) {
	    options->supplied_errors = 1;
	  }
	  else {
	    options->supplied_errors = 0;
	  }
	}
	i++;
	if (i >= *N) {
	  *N += 10;
	  realloc_data(F,*N);
	}
      } /* End working order */
    }
  }
  *N = i;
  realloc_data(F,*N);
  fclose(datafile);
  if (*N == 0) {
    return(F);
  }
  log_comment(options,ROBO_VERBOSE_IO,
	      "read_data_ascii: Read %d samples between %f and %f (spacing %f) for order %d/%d with errors: %d",
	      *N,options->min_x,options->max_x,options->delta_x,options->order,options->max_order,options->supplied_errors);
  if (options->continuum_box < MIN_CONTINUUM_COUNT * options->delta_x) {
    log_comment(options,ROBO_VERBOSE_DEFAULT,
		"read_data_ascii: Updating continuum box from %g to %g because pixel spacing (%g) is too large (need more samples)",
		options->continuum_box,options->delta_x * MIN_CONTINUUM_COUNT,options->delta_x);
    options->continuum_box = options->delta_x * MIN_CONTINUUM_COUNT;
  }
  if ((*N / (1.0 * MAX_CONTINUUM_COUNT) < 100)&&
      (options->continuum_box > MAX_CONTINUUM_COUNT * options->delta_x)) {
    log_comment(options,ROBO_VERBOSE_DEFAULT,
		"read_data_ascii: Updating continuum box from %g to %g because pixel spacing (%g) is too small",
		options->continuum_box,options->delta_x * MAX_CONTINUUM_COUNT,options->delta_x);
    options->continuum_box = options->delta_x * MAX_CONTINUUM_COUNT;
  }
  return(F);
}

#if HAVE_LIBCFITSIO
data *read_data_fits(opts *options, int *N) {
  fitsfile *ff = 0;
  int status = 0;
  int bitpix;
  int naxis;
  long lpixel[2];
  long fpixel[2] = {1,1};
  int x,y;
  double xrefval,yrefval,xrefpix,yrefpix,xinc,yinc,rot;
  char coordtype[80];
  int i,j,last_was = 0;
  double X,Y;
  double *image = 0;

  data *F;
  fits_WCS *WCS;
  /* Set robospect info */
  options->min_x = 99e99;
  options->max_x = -99e99;
  options->delta_x = 0;
  *N = 10;
  F = alloc_data(*N);

  /* Open fits file. */
  fits_open_file(&ff,options->infilename,READONLY,&status); FRE;
  fits_get_img_type(ff,&bitpix,&status); FRE;
  fits_get_img_dim(ff,&naxis,&status); FRE;
  log_comment(options,ROBO_VERBOSE_IO,"read_data_fits: File open: naxis: %d bitpix %d",naxis,bitpix);

  fits_read_img_coord(ff,
		      &xrefval,&yrefval,
		      &xrefpix,&yrefpix,
		      &xinc,&yinc,&rot,
		      coordtype,&status);
  fits_get_img_size(ff,naxis,lpixel,&status); FRE;
  if (naxis == 1) {
    image = malloc(lpixel[0] * sizeof(double));
    if (image == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    fits_read_pix(ff,TDOUBLE,fpixel,lpixel[0],NULL,image,&naxis,&status); FRE;
    lpixel[1] = 1;
  }
  else if (naxis == 2) {
    fpixel[1] = options->order + 1;  /* This probably should be correctly handled by a DISP-AXIS header. */
    image = malloc(lpixel[0] * sizeof(double));
    if (image == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    fits_read_pix(ff,TDOUBLE,fpixel,lpixel[0] ,NULL,image,&naxis,&status); FRE;
    if (options->max_order == 0) {
      options->max_order = lpixel[1] - 1;
    }
  }
  log_comment(options,ROBO_VERBOSE_IO,"read_data_fits: span: %d %d %d %d",
	      fpixel[0],fpixel[1],lpixel[0],lpixel[1]);
  j = 0;

  WCS = read_fits_WCS(ff);

  y = options->order;
  for (x = 0; x < lpixel[0]; x++) {
    i = x;
    X = calculate_fits_WCS(WCS,x,y);
    X = radial_velocity(options,X);
    Y = image[i];
    
    F->x[j] = X;
    if (X < options->min_x) {
      options->min_x = X;
      last_was = 1;
    }
    if ((last_was == 1)&&(Y == 0.0)) {
      options->min_x = X;
      last_was = 1;
      continue;
    }
    else {
      last_was = 0;
    }
    if (X > options->max_x) {
      options->max_x = X;
      if (i > 0) {
	options->delta_x = X - F->x[i-1];
      }
    }
    F->y[j]  = Y;
    F->yO[j] = Y;
    j++;
    if (j >= *N) {
      *N += 10;
      realloc_data(F,*N);
    }
  }

  *N = j;
  realloc_data(F,*N);

  free(image);
  fits_close_file(ff,&status); FRE;

  log_comment(options,ROBO_VERBOSE_IO,"read_data_fits: Read %d samples between %f and %f for order %d/%d\n\n",*N,options->min_x,options->max_x,options->order,options->max_order);

  return(F);
}
#endif

  
