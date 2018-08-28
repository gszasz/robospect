/* Copyright 2010-2013 Christopher Waters and Julie Hollek */
#include "../robo.h"

  
void write_line_data(opts *options, data *D, model *M) {
  FILE *linefile  = fopen(generate_filename(options,"LINES"),"w");
  int i;
  fprintf(linefile,"## %s\n",options->command_line);
  fprintf(linefile,"## Flags\n");
  fprintf(linefile,"## ROBO_MAX_ITER:        0x%06x\n",ROBO_MAX_ITER);
  fprintf(linefile,"## ROBO_FIT_FAIL:        0x%06x\n",ROBO_FIT_FAIL);
  fprintf(linefile,"## ROBO_FIT_IGNORED:     0x%06x\n",ROBO_FIT_IGNORED);
  fprintf(linefile,"## ROBO_FIT_REFUSED:     0x%06x\n",ROBO_FIT_REFUSED);
  fprintf(linefile,"## ROBO_ALT_REFUSED:     0x%06x\n",ROBO_ALT_REFUSED);
  
  fprintf(linefile,"## ROBO_FIT_LARGE_SHIFT: 0x%06x\n",ROBO_FIT_LARGE_SHIFT);
  fprintf(linefile,"## ROBO_FIT_RECENTER:    0x%06x\n",ROBO_FIT_RECENTER);

  fprintf(linefile,"## ROBO_FIT_BLEND:       0x%06x\n",ROBO_FIT_BLEND);
  fprintf(linefile,"## ROBO_FIT_DEBLEND:     0x%06x\n",ROBO_FIT_DEBLEND);
  fprintf(linefile,"## PARAMETER_UNITS:\n");
  if (options->flux_calibrated == 0) {
    fprintf(linefile,"## Ang     Ang        Ang    Ang              Ang     Ang    Ang            Ang        Ang     Ang            mAng      mAng\n");
    fprintf(linefile,"#x0        mean       sigma  flux     eta     dmean   dsigma dflux  deta    FWHM_mean  FWHM_s  FWHM_F         EQW       dEQW  chi    flags    group comment\n");
  }
  else {
    fprintf(linefile,"## Ang     Ang        Ang    Ang              Ang     Ang    Ang            Ang        Ang     Ang            erg/s/cm^2 erg/s/cm^2\n");
    fprintf(linefile,"#x0        mean       sigma  flux     eta     dmean   dsigma dflux  deta    FWHM_mean  FWHM_s  FWHM_F         FLUX       dFLUX  chi    flags    group comment\n");
  }    
  for (i = 0; i < M->L->l; i++) {
    if (M->L->comment[i]) {
      fprintf(linefile,"%9.4f % 9.4f % 7.4f % 7.4f % 7.4f % 9.4f % 7.4f % 7.4f % 7.4f % 9.4f % 7.4f % 7.4f  %10.4g %10.4g % 7.4f 0x%06x %5d # %s\n",
	      M->L->x0[i],
	      M->L->m[i],M->L->s[i],M->L->F[i],M->L->eta[i],
	      M->L->dm[i],M->L->ds[i],M->L->dF[i],M->L->deta[i],
	      M->L->mp[i],M->L->sp[i],M->L->Fp[i],
	      options->flux_calibrated ? M->L->F[i] * interp(M->L->m[i],D->x,M->continuum,D->N) : equivalent_width(M->L,i),
	      options->flux_calibrated ? M->L->dF[i] * interp(M->L->m[i],D->x,M->continuum,D->N) : equivalent_width_error(M->L,i),
	      M->L->chi[i] > 1e6 ? NAN : M->L->chi[i],
	      M->L->flags[i],
	      M->L->blend_group[i],
	      M->L->comment[i]);
    }
    else {
      fprintf(linefile,"%9.4f % 9.4f % 7.4f % 7.4f % 7.4f % 9.4f % 7.4f % 7.4f % 7.4f % 9.4f % 7.4f % 7.4f  %10.4g %10.4g % 7.4f 0x%06x %5d # %s\n",
	      M->L->x0[i],
	      M->L->m[i],M->L->s[i],M->L->F[i],M->L->eta[i],
	      M->L->dm[i],M->L->ds[i],M->L->dF[i],M->L->deta[i],
	      M->L->mp[i],M->L->sp[i],M->L->Fp[i],
	      options->flux_calibrated ? M->L->F[i] * interp(M->L->m[i],D->x,M->continuum,D->N) : equivalent_width(M->L,i),
	      options->flux_calibrated ? M->L->dF[i] * interp(M->L->m[i],D->x,M->continuum,D->N) : equivalent_width_error(M->L,i),
	      M->L->chi[i] > 1e6 ? NAN : M->L->chi[i],
	      M->L->flags[i],
	      M->L->blend_group[i],
	      "");
    }      
  }
  fclose(linefile);  
}
