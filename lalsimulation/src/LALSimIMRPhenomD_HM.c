/*
 * Copyright (C) 2016 Sebastian Khan
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif


#include <math.h>
#include <stdio.h>

#include <lal/LALSimIMR.h>

/* This allows us to reuse internal IMRPhenomD functions without making those functions XLAL */
#include "LALSimIMRPhenomD_internals.c"
#include "LALSimRingdownCW.c"

#include "LALSimIMRPhenomD_HM.h"

#ifndef _OPENMP
#define omp ignore
#endif



// THIS IS MY TESTING FUNCTION
void XLALSimIMRPhenomDHMTEST(void) {
   printf("Hi! I come from XLALSimInspiralPhenomDHMTEST!\n");
}

/**
 * Lionel's QNM higher mode fits
 */
double XLALSimIMRPhenomDHMfring(const REAL8 eta,     /**< symmetric mass-ratio */
                                const REAL8 chi1z,   /**< aligned-spin on larger BH */
                                const REAL8 chi2z,   /**< aligned-spin on smaller BH */
                                const REAL8 finspin, /**< dimensionless final spin */
                                const UINT4 ell,     /**< ell mode */
                                const UINT4 mm       /**< m mode */
                            )
{

    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    const REAL8 Mw = CW07102016( KAPPA(finspin, ell, mm), ell, mm, 0 ) / 2. / LAL_PI; /* convert from angular GW frequency */
    const REAL8 return_val = Mw / (1.0 - EradRational0815(eta, chi1z, chi2z)); /* scale by predicted final mass */

  return return_val;
}


double XLALSimIMRPhenomDHMInspiralFreqScale( const REAL8 f, const UINT4 m )
{
    return 2.0 * f / m;
}

/**
 * XLALSimIMRPhenomDHMFreqDomainMapHM
 * Input waveform frequency in Geometric units (Mf_wf)
 * and computes what frequency this corresponds
 * to scaled to the 22 mode.
 */
double XLALSimIMRPhenomDHMFreqDomainMapHM( const REAL8 Mf_wf,
                                           const UINT4 ell,
                                           const UINT4 mm,
                                           const REAL8 eta,
                                           const REAL8 chi1z,
                                           const REAL8 chi2z,
                                           const UINT4 AmpFlag
                              )
{

    /* compute predicted final spin */
    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);
    REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z);
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    // initialise
    REAL8 Mf_22   = 0.; /* the geometric frequency scaled to the 22 mode, NOTE: called f_map in notes */
    REAL8 Mf_1_22  = AMP_fJoin_INS; /* inspiral joining frequency from PhenomD */
    REAL8 Mf_RD_22 = XLALSimIMRPhenomDHMfring(eta, chi1z, chi2z, finspin, ell, mm);

    REAL8 Mf_1_lm  = XLALSimIMRPhenomDHMInspiralFreqScale( Mf_1_22, mm );
    REAL8 Mf_RD_lm = XLALSimIMRPhenomDHMfring(eta, chi1z, chi2z, finspin, ell, mm);


    /* from technical notes: if AmpFlag is True */
   if ( Mf_wf < Mf_1_lm ) {
       /* inspiral */
       Mf_22 = XLALSimIMRPhenomDHMInspiralFreqScale( Mf_wf, mm );
   } else if ( Mf_1_lm <= Mf_wf && Mf_wf < Mf_RD_lm  ) {
       /* intermediate */
       const REAL8 S = (Mf_RD_22 - Mf_1_22) / (Mf_RD_lm - Mf_1_lm) ;
       Mf_22 = S * ( Mf_wf - Mf_1_lm ) + Mf_1_22 ;
   } else if ( Mf_RD_lm <= Mf_wf ) {
       /* ringdown */
       if ( AmpFlag==1 ) {
           Mf_22 = Mf_wf - Mf_RD_lm + Mf_RD_22 ;
       } else if ( AmpFlag==0 ) {
           const REAL8 rho_lm = Mf_RD_22 / Mf_RD_lm ;
           Mf_22 = rho_lm * Mf_wf;
       }
   }

    return Mf_22;
}
