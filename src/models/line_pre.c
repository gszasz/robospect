/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

/* Discussion on FWHM for Gaussian lines:
   Given a gaussian distribution defined by PDF = 1 / sqrt(2 * pi) * exp(-0.5 * u^2),
   then you can define a fractional max location by:
   PDF / PDF_MAX = exp(-0.5 * u^2) / exp(-0.5 * 0^2);
   PDF / PDF_MAX = exp(-0.5 * u^2)
   sqrt(-2.0 * log(R)) = u
   So, for a given ratio of the maximum R, the half-width R-max is located at u.
   This leads to the full-width R-max location being at
   FWRM = 2.0 * sqrt(-2.0 * log(R))

   Since we've used u, we can add the sigma dependence back in by noting that u = x / s:
   FWRM = 2.0 * sqrt(-2.0 * log(R)) * sigma
   
   The standard value is the FWHM, or R = 0.5:
   FWHM = 2.0 * sqrt(-2.0 * log(0.5)) = 2.0 * sqrt(2.0 * log(2.0))

   Do this for quarter max, R = 0.25:
   FWQM = 2.0 * sqrt(-2.0 * log(0.25)) = 2.0 * sqrt(2.0 * log(4))
        = 2.0 * sqrt(2.0 * 2.0 * log(2.0)) = sqrt(2) * FWHM
   And for 3/4 max, R = 0.75:
   FW3QM = 2.0 * sqrt(-2.0 * log(0.75)) = 2.0 * sqrt(2.0 * log(4/3))
         = 2.0 * sqrt(2.0 * (log(4) - log(3)))
	 = 2.0 * sqrt(2.0 * (2.0 * log(2.0) - log(3)))
	 =~ 0.64423 * FWHM

   If a distribution is "wingier" than a Gaussian, then FWQM will be inflated, such that
   FWQM_wingy > FWQM_Gaussian
	 
*/   

int measure_lines_PRE_individual(data *D, lines *L, int i,int j0) {
  double R = D->O->delta_x;
  int x = 0;
  int j;
  double hwhm1 = 0.0,hwhm2 = 0.0;
  double hwqm1 = 0.0,hwqm2 = 0.0;
  double hw3qm1 = 0.0,hw3qm2 = 0.0;
  double hwhm1_x = 0.0,hwhm2_x = 0.0;
  double hwqm1_x = 0.0,hwqm2_x = 0.0;
  double hw3qm1_x = 0.0,hw3qm2_x = 0.0;
  /* Measure center by finding the sample closest to the expected mean, and taking the centroid around the nearest pixels */
  for (j = j0 ; j < D->N; j++) {
    if (fabs(D->x[j] - L->x0[i]) < R) {
      R = fabs(D->x[j] - L->x0[i]);
      x = j;
      if ((x > 0)&&(x < D->N - 1)) {
	L->mp[i] = (D->y[x] * D->x[x] + D->y[x-1] * D->x[x-1] + D->y[x+1] * D->x[x+1]) / (D->y[x] + D->y[x-1] + D->y[x+1]);
      }
      else {
	L->mp[i] = D->x[x];
      }
    }
  }
  
  /* Determine the lower bound for the FWHM via linear interpolation */
  for (j = x; j > j0; j--) {
    if (fabs(D->y[j - 1]) < 0.5 * fabs(D->y[x])) {
      hwhm1_x = ((D->x[j] - D->x[j - 1]) / (D->y[j] - D->y[j - 1])) * (0.5 * D->y[x] - D->y[j-1]) + D->x[j-1];
      log_comment(D->O,ROBO_VERBOSE_DEBUG,"X1 %d %g %g %g %g %g %g %g",
		  j,hwhm1_x,D->x[j],D->x[j-1],fabs(D->y[j]),fabs(D->y[j-1]),D->x[x],fabs(D->y[x]));
      j = 0;
    }
  }
  for (j = x; j > j0; j--) {
    if (fabs(D->y[j - 1]) < 0.75 * fabs(D->y[x])) {
      hw3qm1_x = ((D->x[j] - D->x[j - 1]) / (D->y[j] - D->y[j - 1])) * (0.75 * D->y[x] - D->y[j-1]) + D->x[j-1];
      log_comment(D->O,ROBO_VERBOSE_DEBUG,"X1a %d %g %g %g %g %g %g %g",
		  j,hw3qm1_x,D->x[j],D->x[j-1],fabs(D->y[j]),fabs(D->y[j-1]),D->x[x],fabs(D->y[x]));
      j = 0;
    }
  }
  for (j = x; j > j0; j--) {
    if (fabs(D->y[j - 1]) < 0.25 * fabs(D->y[x])) {
      hwqm1_x = ((D->x[j] - D->x[j - 1]) / (D->y[j] - D->y[j - 1])) * (0.25 * D->y[x] - D->y[j-1]) + D->x[j-1];
      log_comment(D->O,ROBO_VERBOSE_DEBUG,"X1a %d %g %g %g %g %g %g %g",
		  j,hwqm1_x,D->x[j],D->x[j-1],fabs(D->y[j]),fabs(D->y[j-1]),D->x[x],fabs(D->y[x]));
      j = 0;
    }
  }
  /* Find the upper bound for the FWHM via linear interpolation */
  for (j = x; j < D->N; j++) {
    if (fabs(D->y[j]) < 0.5 * fabs(D->y[x])) {
      hwhm2_x = ((D->x[j] - D->x[j - 1]) / (D->y[j] - D->y[j - 1])) * (0.5 * D->y[x] - D->y[j-1]) + D->x[j-1];
      log_comment(D->O,ROBO_VERBOSE_DEBUG,"X2 %d %g %g %g %g %g %g %g",
		  j,hwhm2_x,D->x[j],D->x[j-1],fabs(D->y[j]),fabs(D->y[j-1]),D->x[x],fabs(D->y[x]));
      j = D->N + 100;
    }
  }
  for (j = x; j < D->N; j++) {
    if (fabs(D->y[j]) < 0.75 * fabs(D->y[x])) {
      hw3qm2_x = ((D->x[j] - D->x[j - 1]) / (D->y[j] - D->y[j - 1])) * (0.75 * D->y[x] - D->y[j-1]) + D->x[j-1];
      log_comment(D->O,ROBO_VERBOSE_DEBUG,"X2a %d %g %g %g %g %g %g %g",
		  j,hw3qm2_x,D->x[j],D->x[j-1],fabs(D->y[j]),fabs(D->y[j-1]),D->x[x],fabs(D->y[x]));
      j = D->N + 100;
    }
  }
  for (j = x; j < D->N; j++) {
    if (fabs(D->y[j]) < 0.25 * fabs(D->y[x])) {
      hwqm2_x = ((D->x[j] - D->x[j - 1]) / (D->y[j] - D->y[j - 1])) * (0.25 * D->y[x] - D->y[j-1]) + D->x[j-1];
      log_comment(D->O,ROBO_VERBOSE_DEBUG,"X2a %d %g %g %g %g %g %g %g",
		  j,hwqm2_x,D->x[j],D->x[j-1],fabs(D->y[j]),fabs(D->y[j-1]),D->x[x],fabs(D->y[x]));
      j = D->N + 100;
    }
  }
  /* Convert the upper and lower bounds into true widths   */
  if (hwhm1_x == 0.0) { /* No solution found in this direction */
    hwhm1 = 0.0;
  }
  else { /* Calculate the offset from the line center */
    hwhm1 = fabs(hwhm1_x - L->mp[i]);
  }
  if (hwhm2_x == 0.0) {
    hwhm2 = 0.0;
  }
  else {
    hwhm2 = fabs(hwhm2_x - L->mp[i]);
  }
  /* 3-quarter max */
  if (hw3qm1_x == 0.0) {
    hw3qm1 = 0.0;
  }
  else {
    hw3qm1 = fabs(hw3qm1_x - L->mp[i]);
  }
  if (hw3qm2_x == 0.0) {
    hw3qm2 = 0.0;
  }
  else {
    hw3qm2 = fabs(hw3qm2_x - L->mp[i]);
  }
  /* quarter max */
  if (hwqm1_x == 0.0) {
    hwqm1 = 0.0;
  }
  else {
    hwqm1 = fabs(hwqm1_x - L->mp[i]);
  }
  if (hwqm2_x == 0.0) {
    hwqm2 = 0.0;
  }
  else {
    hwqm2 = fabs(hwqm2_x - L->mp[i]);
  }
  log_comment(D->O,ROBO_VERBOSE_DEBUG,"FWHM: %g: %g %g %g %g %g %g",
	      L->x0[i],hwhm1,hwhm2,hw3qm1,hw3qm2,hwqm1,hwqm2);
  
  /* Average hwhm values to get a sigma value, and calculate the gaussian flux */
  if ((hwhm1 == 0.0)&&(hwhm2 == 0.0)) {
    L->sp[i] = (hw3qm2 + hw3qm1) / 1.55223;/* inverse of the value from above. */
  }
  else if ((hwhm1 == 0.0)||(hwhm1 > 2.0 * hwhm2)) {
    L->sp[i] = (hwhm2) / (sqrt(2.0 * log(2.0)));
  }
  else if ((hwhm2 == 0.0)||(hwhm2 > 2.0 * hwhm1)) {
    L->sp[i] = (hwhm1) / (sqrt(2.0 * log(2.0)));
  }
  else {
    L->sp[i] = (hwhm2 + hwhm1) / (2.0 * sqrt(2.0 * log(2.0)));
  }
  if (L->sp[i] == 0.0) {
    L->sp[i] = D->x[x+1] - D->x[x];
  }
  L->Fp[i] = D->y[x] * pow(L->sp[i] * sqrt(2 * M_PI),1);

  /* Make a guess at the eta value */
  if (((hwhm2 + hwhm1) / (hw3qm2 + hw3qm1)) > 1.68) {
    L->etap[i] = -132.4711 + 79.3913 * ((hwhm2 + hwhm1) / (hw3qm2 + hw3qm1));
  }
  else {
    L->etap[i] = -18.7118 + 11.9942 * ((hwhm2 + hwhm1) / (hw3qm2 + hw3qm1));
  }
  if (L->etap[i] < 0.0) { L->etap[i] = 0.0; }
  
  j0 = x;

  log_comment(D->O,ROBO_VERBOSE_LINE,
	      "mlP: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f) eta (%f)",
	      i,L->l,0,-1.0,L->x0[i],L->mp[i],L->sp[i],L->Fp[i],
	      hwhm2 + 0.0 *hwhm1,hw3qm2 + 0.0 * hw3qm1,hwqm2 + 0.0 * hwqm1,L->etap[i]);
  return(j0);
}

void measure_lines_PRE(data *D, model *M) {
  int i,j,j0;
  j0 = 0;

  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }
  for (i = 0; i < M->L->l; i++) {
    j0 = measure_lines_PRE_individual(D,M->L,i,j0);
    if ((fabs(M->L->mp[i] - M->L->x0[i]) > M->L->sp[i])||
	(fabs(M->L->mp[i] - M->L->x0[i]) > 100)) {
      M->L->flags[i] |= ROBO_ALT_REFUSED;
    }
  }
  M->L->b = 1;
  for (i = 0; i < M->L->l; i++) {
    if (M->L->flags[i] & ROBO_ALT_REFUSED) {
      continue;
    }
    for (j = i+1; j< M->L->l; j++) {
      if (M->L->flags[j] & ROBO_ALT_REFUSED) {
	continue;
      }
      if (((M->L->mp[j] - M->L->mp[i]) < (D->O->deblend_radius * (fabs(M->L->sp[j]) + fabs(M->L->sp[i]))))) { 
	if (M->L->blend_group[i] == 0) {
	  M->L->flags[i] |= ROBO_FIT_BLEND;
	  M->L->blend_group[i] = M->L->b;
	  M->L->b++;
	}
	M->L->blend_group[j] = M->L->blend_group[i];
	M->L->flags[j] |= ROBO_FIT_BLEND;
	log_comment(D->O,ROBO_VERBOSE_DEBUG,
		    "mlPG: (%d %d) %d %f %f %f %f",
		    i,j,M->L->blend_group[i],
		    M->L->mp[i],M->L->mp[j],
		    (M->L->mp[j] - M->L->mp[i]),
		    (D->O->deblend_radius * (fabs(M->L->sp[j]) + fabs(M->L->sp[i]))));
      }
      else {
	log_comment(D->O,ROBO_VERBOSE_DEBUG,
		    "mlPG: (%d %d) %d %f %f %f %f",
		    i,j,0,
		    M->L->mp[i],M->L->mp[j],
		    (M->L->mp[j] - M->L->mp[i]),
		    (D->O->deblend_radius * (fabs(M->L->sp[j]) + fabs(M->L->sp[i]))));
      }

    }
  }
}
	


void measure_lines_PRE_second_version(data *D, model *M) {
  int i,j,j0,k;
  int x;
  lines *dwarf_lines = NULL;
  j0 = 0;

  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }
  
  for (k = 0; k <= (D->O->deblend_iterations - 1) * (M->deblend_model != DEBLEND_NULL); k++) {
    for (i = 0; i < M->L->l; i++) {
      j0 = measure_lines_PRE_individual(D,M->L,i,j0);
      /* End basic line parameter measurement. */
      /* Now, attempt to see if any of these lines are likely to be blended. */
      if ((M->deblend_model != DEBLEND_NULL)&&(M->L->manual[i] == 1)) {
	/*      dwarf_lines->l = 0; */
	double min, max;
	int min_x,max_x;
	int found_match;
	x = j0;
	/* Look to the left of the line. */
    	min = fabs(D->y[j0]);
    	max = -99e99;
    	min_x = -1;
    	max_x = -1;
    	for (j = x - 1; (D->x[j] > M->L->mp[i] - D->O->deblend_radius * M->L->sp[i])&&(j > 0); j--) {
    	  if ((fabs(D->y[j]) < min)&&(max_x == -1)) {
    	    min = fabs(D->y[j]);
    	    min_x = j;
    	  }
    	  else {
    	    if (fabs(D->y[j]) - min > D->e[j]) {
    	      max_x = j;
    	      max = fabs(D->y[j]);
    	    }
    	  }
    	  log_comment(D->O,ROBO_VERBOSE_DEBUG,"dc1: (%g %g %g) (%g %g %g) (%g %g %d) (%g %g %d)",
    		      D->x[j],D->y[j],D->e[j],
    		      D->x[j0],D->y[j0],D->e[j0],
    		      min,min_x != -1 ?  D->x[min_x] : -1.0 ,min_x,
    		      max,max_x != -1 ?  D->x[max_x] : -1.0 ,max_x);
    	}
    	if ((max > fabs(D->O->deblend_ratio * D->y[j0]))&&(max > min)&&(max > D->e[j0])) {
    	  /* We've found something. */
    	  M->L->flags[i] |= ROBO_FIT_BLEND;
    	  found_match = 0;
    	  /* Find the actual peak (which may be outside the Nsigma box) */
    	  for (j = max_x; j >= 0; j--) {
    	    if (fabs(D->y[j]) > fabs(D->y[max_x])) {
    	      max_x = j;
    	    }
    	    else {
    	      break;
    	    }
    	  }
    	  /* Check that this isn't another known line by making sure it's not within two pixels of the known line */
    	  for (j = i ; (j >= i - 50)&&(j >= 0); j--) {
    	    log_comment(D->O,ROBO_VERBOSE_DEBUG,"ASDF: %d %d %d %g %g %g %g",
    			j,max_x,j0,D->x[max_x],M->L->x0[j],fabs(M->L->x0[j] - D->x[max_x]),fabs(D->x[max_x] - D->x[max_x - 2]));
    	    if (fabs(M->L->x0[j] - D->x[max_x]) < fabs(D->x[max_x] - D->x[max_x - 2])) {
    	      found_match = 1;
    	      M->L->flags[j] |= ROBO_FIT_BLEND;
    	    }
    	  }
	  if (dwarf_lines) {
	    for (j = dwarf_lines->l; (j > dwarf_lines->l - 50)&&(j >= 0); j--) {
	      log_comment(D->O,ROBO_VERBOSE_DEBUG,"ASDF3a: %d %d %d %g %g %g %g",
			  j,max_x,j0,D->x[max_x],dwarf_lines->x0[j],fabs(dwarf_lines->x0[j] - D->x[max_x]),fabs(D->x[max_x] - D->x[max_x - 2]));
	      
	      if (fabs(dwarf_lines->x0[j] - D->x[max_x]) < fabs(D->x[max_x] - D->x[max_x - 2])) {
		found_match = 1;
		M->L->flags[j] |= ROBO_FIT_BLEND;
		break;
	      }
	    }
	  }
	
    	  if (found_match == 0) {
    	    /* Add this line to the list of new lines */
	    if (dwarf_lines) {
	      if (D->x[max_x] < dwarf_lines->x0[dwarf_lines->l - 1]) {
		log_comment(D->O,ROBO_VERBOSE_LINE,"A: I just inserted a line less than the last line's value: line: %g new: %g last: %g",
			    M->L->x0[i],D->x[max_x],dwarf_lines->x0[dwarf_lines->l - 1]);
	      }
	    }
    	    dwarf_lines = append_lines(dwarf_lines,D->x[max_x]);
    	    dwarf_lines->flags[dwarf_lines->l - 1] |= ROBO_FIT_BLEND;
	    dwarf_lines->comment[dwarf_lines->l - 1] = malloc(sizeof(char) * 100);
	    snprintf(dwarf_lines->comment[dwarf_lines->l - 1],100,
		     "Dwarf line at %g adjacent to line at %g(%d)",
		     D->x[max_x],M->L->mp[i],i);		    
	   
	    measure_lines_PRE_individual(D,dwarf_lines,dwarf_lines->l - 1,0);
    	    log_comment(D->O,ROBO_VERBOSE_LINE,"DC1: Probable blend line peak around %g (%g %g) %d %g",
    			D->x[max_x],D->y[max_x],D->O->deblend_ratio * D->y[j0],found_match,M->L->x0[j]);
    	  }
    	}
	/* Look to the right of the line */
    	min = fabs(D->y[j0]);
    	max = -99e99;
    	min_x = -1;
    	max_x = -1;
    	for (j = x + 1; (D->x[j] < M->L->mp[i] + D->O->deblend_radius * M->L->sp[i])&&(j < D->N); j++) {
    	  if ((fabs(D->y[j]) < min)&&(max_x == -1)) {
    	    min = fabs(D->y[j]);
    	    min_x = j;
    	  }
    	  else {
    	    if (fabs(D->y[j]) - min > D->e[j]) {
    	      max_x = j;
    	      max = fabs(D->y[j]);
    	    }
    	  }
    	  log_comment(D->O,ROBO_VERBOSE_DEBUG,"dc2: (%g %g %g) (%g %g %g) (%g %g %d) (%g %g %d)",
    		      D->x[j],D->y[j],D->e[j],
    		      D->x[j0],D->y[j0],D->e[j0],
    		      min,min_x != -1 ?  D->x[min_x] : -1.0 ,min_x,
    		      max,max_x != -1 ?  D->x[max_x] : -1.0 ,max_x);
    	}
    	if ((max > fabs(D->O->deblend_ratio * D->y[j0]))&&(max > min)&&(max > D->e[j0])) {
    	  /* We've found something */
    	  M->L->flags[i] |= ROBO_FIT_BLEND;
    	  found_match = 0;
    	  /* Find the actual peak (which may be outside the Nsigma box) */
    	  for (j = max_x; j < D->N; j++) {
    	    if (fabs(D->y[j]) > fabs(D->y[max_x])) {
    	      max_x = j;
    	    }
    	    else {
    	      break;
    	    }
    	  }
    	  /* Check that this isn't another known line */
    	  for (j = i + 1; (j < i + 50)&&(j < M->L->l); j++) {
    	    log_comment(D->O,ROBO_VERBOSE_LINE,"ASDF2: %d %d %d %g %g %g %g",
    			j,max_x,j0,D->x[max_x],M->L->x0[j],fabs(M->L->x0[j] - D->x[max_x]),fabs(D->x[max_x] - D->x[max_x - 2]));
	    
    	    if (fabs(M->L->x0[j] - D->x[max_x]) < fabs(D->x[max_x] - D->x[max_x - 2])) {
    	      found_match = 1;
    	      M->L->flags[j] |= ROBO_FIT_BLEND;
    	      break;
    	    }
    	  }
    	  /* Check that we haven't found this line already */
	  if (dwarf_lines) {
	    for (j = dwarf_lines->l; (j > dwarf_lines->l - 50)&&(j >= 0); j--) {
	      log_comment(D->O,ROBO_VERBOSE_LINE,"ASDF3: %d %d %d %g %g %g %g",
			  j,max_x,j0,D->x[max_x],dwarf_lines->x0[j],fabs(dwarf_lines->x0[j] - D->x[max_x]),fabs(D->x[max_x] - D->x[max_x - 2]));
	      
	      if (fabs(dwarf_lines->x0[j] - D->x[max_x]) < fabs(D->x[max_x] - D->x[max_x - 2])) {
		found_match = 1;
		M->L->flags[j] |= ROBO_FIT_BLEND;
		break;
	      }
	    }
	  }
    	  log_comment(D->O,ROBO_VERBOSE_LINE,"DC2: Probable blend line peak around %g (%g %g) %d %g",
    		      D->x[max_x],D->y[max_x],D->O->deblend_ratio * D->y[j0],found_match,M->L->x0[j]);
	  
    	  if (found_match == 0) {
    	    /* Add this line to the list of new lines. */
	    
    	    dwarf_lines = append_lines(dwarf_lines,D->x[max_x]);
	    
    	    dwarf_lines->flags[dwarf_lines->l - 1] |= ROBO_FIT_BLEND;
	    dwarf_lines->comment[dwarf_lines->l - 1] = malloc(sizeof(char) * 100);
	    snprintf(dwarf_lines->comment[dwarf_lines->l - 1],100,
		     "Dwarf line at %g adjacent to line at %g(%d)",
		     D->x[max_x],M->L->mp[i],i);		    

	    measure_lines_PRE_individual(D,dwarf_lines,dwarf_lines->l - 1,0);
    	    log_comment(D->O,ROBO_VERBOSE_LINE,"DC2: Probable blend line peak around %g (%g %g) %d %g",
    			D->x[max_x],D->y[max_x],D->O->deblend_ratio * D->y[j0],found_match,M->L->x0[j]);
    	  }
    	}
      }
    }
    if (dwarf_lines) {
      M->L = merge_lines(M->L,dwarf_lines);
    }
    dwarf_lines = 0;
  }
  M->L->b = 1;
  x = 0;
  for (j = 0; j < M->L->l; j++) {
    if (M->L->flags[j] & ROBO_FIT_BLEND) {
      /* This claims to be blended. */
      if ((j > 0)&&(M->L->flags[j-1] & ROBO_FIT_BLEND)) {
	/* My previous neighbor claims to be blended as well */
	if (M->L->m[j] > M->L->m[j-1] + D->O->deblend_radius * M->L->s[j-1]) {
	  /* My neighbor is too far away, so I'm not in that group. */
	  M->L->b++;
	  x = 1;
	}
      }
      if ((j < M->L->l - 1)&&(M->L->flags[j+1] & ROBO_FIT_BLEND)) {
	/* My subsequent neighbor claims to be blended as well */
	if (M->L->m[j] < M->L->m[j+1] - D->O->deblend_radius * M->L->s[j+1]) {
	  /* That neighbor is too far away, so that's not in my group. */
	}
	else {
	}
      }
      M->L->blend_group[j] = M->L->b;
      M->L->flags[j] |= ROBO_FIT_BLEND;
      x++;
      if (x >= 30) {
	M->L->b++;
	x = 0;
      }
    }
    else if ((j > 0)&&(M->L->flags[j-1] & ROBO_FIT_BLEND)) {
      M->L->b++;
      x = 1;
    }
  }
      
    
}

void measure_lines_PRE_original(data *D, model *M) {
  int i,j,j0;
  double R;
  int x;
  double hwhm1 = 0,hwhm2 = 0;
  j0 = 0;

  for (j = 0; j < M->N; j++) {
    M->lines[j] = 0.0;
  }
  
  for (i = 0; i < M->L->l; i++) {
    R = 1.0;
    x = 0;

    /* Measure center by finding the sample closest to the expected mean */
    for (j = j0 ; j < D->N; j++) {
      if (fabs(D->x[j] - M->L->x0[i]) < R) {
	R = fabs(D->x[j] - M->L->x0[i]);
	x = j;
	M->L->m[i] = D->x[x];
      }
    }
    if (x == 0) {
      continue; /* We failed to find a valid match */
    }

    /* Determine the lower bound for the FWHM. */
    for (j = j0; j < x; j++) {
      if (fabs(D->y[j]) < 0.5 * fabs(D->y[x])) {
	hwhm1 = fabs(D->x[j] - D->x[x]);
      }
    }

    /* Find the upper bound for the FWHM */
    for (j = x + 1; j < D->N; j++) {
      if (fabs(D->y[j]) < 0.5 * fabs(D->y[x])) {
	hwhm2 = fabs(D->x[j] - D->x[x]);
	j = D->N + 100;
      }
    }
    /* Average hwhm values to get a sigma value, and calculate the gaussian flux */
    M->L->s[i] = (hwhm1 + hwhm2) / (2.0 * sqrt(2.0 * log(2.0)));
    M->L->F[i] = D->y[x] * pow(M->L->s[i] * sqrt(2 * M_PI),1);
    j0 = x;

    M->L->mp[i] = M->L->m[i];
    M->L->sp[i] = M->L->s[i];
    M->L->Fp[i] = M->L->F[i];
    
    log_comment(D->O,ROBO_VERBOSE_LINE,
		"mlP: (%d/%ld) %d %f %f (%f,%f,%f) +- (%f,%f,%f)",
		i,M->L->l,0,-1.0,M->L->x0[i],M->L->m[i],M->L->s[i],M->L->F[i],
		-1.0,-1.0,-1.0);

  }
}
