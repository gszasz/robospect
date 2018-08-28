/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

int compD (const void *a, const void *b) {
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}


void generate_model_continuum_devel(opts *options, data *D, model *M) {
  generate_model_continuum_boxcar(options,D,M);

#ifdef DEVEL
  int i,j;
  stats *S1;
  double *v1 = malloc(sizeof(double) * D->N);
  int n1;

  for (i = 0; i < D->N; i += 1) {
    n1 = 0;
    n2 = 0;
    for (j = 0; j < D->N; j++) {
      if ((D->x[j] > D->x[i] - options->continuum_box)&&(D->x[j] < D->x[i])) {
	v1[n1] = (D->y[j] - M->lines[j]);
	n1++;
      }
      else if ((D->x[j] > D->x[i])&&(D->x[j] < D->x[i] + options->continuum_box)) {
	v1[n1] = (D->y[j] - M->lines[j]);
	n1++;
      }
      if (D->x[j] > D->x[i] + options->continuum_box) {
	j = D->N + 10;
      }
    }
    S1 = array_stats(v1,n1);
    /* Create scratch array we can sort; */
    double *v2 = malloc(sizeof(double) * n1);
    for (j = 0; j < n1; j++) {
      v2[j] = v1[j] - S1->med;
    }
    qsort(v2,n1,sizeof(double),compD);

    double sigma = 0.0;
    double ns = 0.0;
    for (j = 0; j < n1; j++) {
      if (v2[j] > 0.0) {
    	sigma += pow(v2[j],2);
    	ns += 1.0;
      }
    }
    sigma = sqrt(sigma / ns);
    
    /* Construct cdf */
    double *cdf = malloc(sizeof(double) * n1);

    for (j = 0; j < n1; j++) {
      cdf[j] = 1.0 * (j) / (1.0 * n1);
    }

    /* Construct a weight */
    double *w = malloc(sizeof(double) * n1);

    M->continuum[i] = 0.0;
    ns = 0.0;
    for (j = 0; j < n1; j++) {
      double kexp = 0.5 + 0.5 * erf((v2[j])/(sqrt(2) * sigma));
      double kcdf = cdf[j];
      if (kcdf == 0.0) { kcdf = 1.0; }
      w[j] = kexp / kcdf;

      if (v2[j] > 0.0) { w[j] = 1.0;}
      M->continuum[i] += (v2[j] + S1->med) * w[j];
      ns += w[j];
    }
    M->continuum[i] /= ns;

    int count,k;
    double *pdf = malloc(sizeof(double) * 20);
    for (k = 0; k < 20; k++) {
      pdf[k] = 0;
    }
    double pdf_min = 99e99;
    double pdf_max = -99e99;
    double box_min,box_max;
    
    for (j = 0; j < n1; j++) {
      if (v1[j] < pdf_min) { pdf_min = v1[j]; }
      if (v1[j] > pdf_max) { pdf_max = v1[j]; }
    }

    for (j = 0; j < n1; j++) {
      k = (int) (20 * ((v1[j] - pdf_min)/(pdf_max - pdf_min)));
      pdf[k] ++;
      if (i == 21145) {
	fprintf(stderr,"%d %d %g %g %g %g\n",j,k,pdf_min,pdf_max,v1[j],pdf[k]);
      }
    }
    count = 0;
    for (k = 0; k < 20; k++) {
      if (pdf[k] > count) {
	box_min = k/20.0 * (pdf_max - pdf_min) + pdf_min;
	box_max = (k+1)/20.0 * (pdf_max - pdf_min) + pdf_min;
	count = pdf[k];
	if (i == 21145) {
	  fprintf(stderr,"%d %d %g %g\n",k,count,box_min,box_max);
	}

      }
    }
    
    pdf_min = box_min;
    pdf_max = box_max;
    for (k = 0; k < 20; k++) {
      pdf[k] = 0;
    }
    for (j = 0; j < n1; j++) {
      k = (int) (20 * ((v1[j] - pdf_min)/(pdf_max - pdf_min)));
      if ((k >= 0)&&(k < 20)) {
	pdf[k] ++;
      }
      if (i == 21145) {
	fprintf(stderr,"%d %d %g %g %g %g\n",j,k,pdf_min,pdf_max,v1[j],pdf[k]);
      }
    }
    count = 0;
    for (k = 0; k < 20; k++) {
      if (pdf[k] > count) {
	box_min = k/20.0 * (pdf_max - pdf_min) + pdf_min;
	box_max = (k+1)/20.0 * (pdf_max - pdf_min) + pdf_min;
	count = pdf[k];
	if (i == 21145) {
	  fprintf(stderr,"%d %d %g %g\n",k,count,box_min,box_max);
	}
      }
    }
    M->continuum[i] = (box_max + box_min) / 2.0;

    /* Free values; */
    free(w);
    free(cdf);
    free(v2);

    if (options->supplied_errors == 0) {
      if (options->flux_calibrated == 0) {
	D->e[i]         = sigma;
      }
      else {
	D->e[i]         = sigma;
      }
      if (D->e[i] <= 0.0) {
	D->e[i] = D->e[i-1];
      }
    }
    log_comment(options,ROBO_VERBOSE_CONTINUUM,
		"gmcb: (%d/%d) %g (%f %g %g) (%g %g %ld) %g",
		i,D->N,
		options->continuum_box,
		D->x[i],D->y[i],D->e[i],
		S1->med,1.4826 * S1->MAD,S1->N,
		M->continuum[i]
		);
    free(S1);
  }
  
  free(v1);
#endif
}

