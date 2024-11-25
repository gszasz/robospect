/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"
#if HAVE_LIBCFITSIO
typedef enum {
  SPECT_LINEAR,
  SPECT_EQUISPEC,
  SPECT_MULTISPEC,
  SPECT_NULL
} spect_type;

fits_WCS *read_fits_WCS(fitsfile *ff) {
  fits_WCS *WCS = malloc(sizeof(fits_WCS));
  if (WCS == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
  }

  int i;
  char value[80];
  char comment[80];
  int status = 0;
  spect_type S = SPECT_NULL;
  int naxis = 0;
  double scale_factor = 1.0;
  double crval,crpix,cdelt;

  /* MULTISPEC params */
  char keyword[9];
  char *fullspect;
  char *token;
  int  k, pos;
  int  ret;
  
  fits_read_keyword(ff,"WAT0_001",value,comment,&status);  FRE;
  if ((strcasecmp(value,"'system=world'") == 0)||(status != 0)) {
    S = SPECT_LINEAR;
    WCS->dispersion_axis = 0;
    status = 0;
  }
  else if (strcasecmp(value,"'system=equispec'") == 0) {
    S = SPECT_EQUISPEC;
    WCS->dispersion_axis = 0;
  }
  else if (strcasecmp(value,"'system=multispec'") == 0) {
    WCS->dispersion_axis = 0;
    S = SPECT_MULTISPEC;
    
  }
  if (S == SPECT_NULL) {
    S = SPECT_LINEAR; /* Try to do something */
  }
  fits_read_keyword(ff,"NAXIS", value,NULL,&status); FRE;
  naxis = atoi(value);
  if (naxis == 1) {
    fits_read_keyword(ff,"NAXIS1",value,NULL,&status); FRE;
  }
  else {
    fits_read_keyword(ff,"NAXIS2",value,NULL,&status); FRE;
  }
  WCS->N = atoi(value);
  WCS->x0 = malloc(WCS->N * sizeof(double));
  if (WCS->x0 == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  WCS->xd = malloc(WCS->N * sizeof(double));
  if (WCS->xd == NULL) {
    fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  fits_read_keyword(ff,"WAT1_001",value,NULL,&status); FRE;
  
  if (strcasecmp(value,"'wtype=linear label=Wavelength units=angstroms'") == 0) {
    scale_factor = 1.0;
  }
  if (strcasecmp(value,"'wtype=linear label=Wavelength units=nanometers'") == 0) {
    scale_factor = 10.0;
  }
  if (strcasecmp(value,"'wtype=linear label=Wavelength units=meters'") == 0) {
    scale_factor = 1e10;
  }

  if ((S == SPECT_LINEAR)||(S == SPECT_EQUISPEC)) {
    fits_read_keyword(ff,"CRVAL1",value,NULL,&status); FRE;
    crval = atof(value);
    fits_read_keyword(ff,"CRPIX1",value,NULL,&status); FRE;
    crpix = atof(value);
    fits_read_keyword(ff,"CD1_1", value,NULL,&status); /* no FRE here, as we need to catch status. */
    cdelt = atof(value);
    if (status != 0) {
      fits_read_keyword(ff,"CDELT1", value,NULL,&status); FRE;
      cdelt = atof(value);
    }
    
    for (i = 0; i < WCS->N; i++) {
      WCS->x0[i] = (crval - cdelt * (crpix - 1.0)) * scale_factor;
      WCS->xd[i] = cdelt * scale_factor;
    }
  }
  else if (S == SPECT_MULTISPEC) {
    status = 0;
    fullspect = malloc(WCS->N * 2 * 80);
    if (fullspect == NULL) {
      fprintf(stderr, "robospect: %s: %d: Cannot allocate memory\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    
    for (k = 1, pos = 0; k < WCS->N * 2; k++) {
      ret = snprintf(keyword, sizeof(keyword), "WAT2_%03d", k);
      if (ret < 0) {
        fprintf(stderr, "robospect: %s: %d: Buffer overflow\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
      }
      fits_read_key(ff,TSTRING,keyword,value,NULL,&status);
      if (status == 0) {
        pos += sprintf(&fullspect[pos], "%-68s", value);
      }
      else {
	k = WCS->N * 3;
	status = 0;
      }
    }

    i = 0;
    token = strtok(fullspect," ");
    do {
      if (strncmp(token,"spec",4) == 0) {
	token = strtok(NULL," "); /* equals sign */
	token = strtok(NULL," \""); /* ap */
	token = strtok(NULL," "); /* beam */
	token = strtok(NULL," "); /* dtype */
	token = strtok(NULL," "); /* lambda */
	crval = atof(token);
	token = strtok(NULL," "); /* dlambda */
	cdelt = atof(token);
	token = strtok(NULL," "); /* nlambda */
	token = strtok(NULL," "); /* z */
	crpix = atof(token); /* lazy */
	token = strtok(NULL," "); /* aplow */
	token = strtok(NULL," \""); /* aphigh */

	WCS->x0[i] = crval / (1 + crpix);
	WCS->xd[i] = cdelt / (1 + crpix);

	i++;
      }
      token = strtok(NULL," \"");	
    } while (token != NULL);
  }
  return(WCS);
      
}

double calculate_fits_WCS(fits_WCS *WCS, int x, int y) {
  if (WCS->dispersion_axis == 0) {
    return(WCS->x0[y] + x * WCS->xd[y]);
  }
  else {
    return(WCS->x0[x] + y * WCS->xd[x]);
  }
}

void free_fits_WCS(fits_WCS *WCS) {
  free(WCS->x0);
  free(WCS->xd);
  free(WCS);
}
#endif

