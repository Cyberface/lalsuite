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
 * returns real part of the ringdown frequency
 */
double XLALSimIMRPhenomDHMfring(const REAL8 eta,     /**< symmetric mass-ratio */
                                const REAL8 chi1z,   /**< aligned-spin on larger BH */
                                const REAL8 chi2z,   /**< aligned-spin on smaller BH */
                                const REAL8 finspin, /**< dimensionless final spin */
                                const INT4 ell,     /**< ell mode */
                                const INT4 mm       /**< m mode */
                            )
{

    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");
    const complex double ZZ = CW07102016( KAPPA(finspin, ell, mm), ell, mm, 0 );
    const REAL8 Mf_RD = creal(ZZ) / 2. / LAL_PI; /* GW ringdown frequency, converted from angular frequency */
    const REAL8 return_val = Mf_RD / (1.0 - EradRational0815(eta, chi1z, chi2z)); /* scale by predicted final mass */

    return return_val;
}

/**
 * Lionel's QNM higher mode fits
 * returns imaginary part of the ringdown frequency
 */
double XLALSimIMRPhenomDHMfdamp(const REAL8 eta,     /**< symmetric mass-ratio */
                                const REAL8 chi1z,   /**< aligned-spin on larger BH */
                                const REAL8 chi2z,   /**< aligned-spin on smaller BH */
                                const REAL8 finspin, /**< dimensionless final spin */
                                const INT4 ell,     /**< ell mode */
                                const INT4 mm       /**< m mode */
                            )
{

    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");
    const complex double ZZ = CW07102016( KAPPA(finspin, ell, mm), ell, mm, 0 );
    const REAL8 fdamp = cimag(ZZ)  / 2. / LAL_PI ; /* this is the 1./tau in the complex QNM */
    const REAL8 return_val = fdamp / (1.0 - EradRational0815(eta, chi1z, chi2z)); /* scale by predicted final mass */

    return return_val;
}


/**
 * Equation 20 arXiv:1508.07253 (called f_peak in paper)
 * analytic location of maximum of AmpMRDAnsatz
 * generalised for PhenomD_HM model.
 * Uses different values for fRD and fDM depending
 * on ell and m spherical harmonic mode.
 */

/*NOTE: I THINK THIS FUNCTION IS NOT NEEDED ACTUALLY*/
double XLALSimIMRPhenomDHMfmaxCalc(
                                   const REAL8 eta,
                                   const REAL8 chi1z,
                                   const REAL8 chi2z,
                                   const INT4 ell,
                                   const INT4 mm
                               )
{
    const REAL8 chi = chiPN(eta, chi1z, chi2z);
    const REAL8 gamma2 = gamma2_fun(eta, chi);
    const REAL8 gamma3 = gamma3_fun(eta, chi);

    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);
    const REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z);
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    const REAL8 fRD = XLALSimIMRPhenomDHMfring(eta, chi1z, chi2z, finspin, ell, mm);
    const REAL8 fDM = XLALSimIMRPhenomDHMfdamp(eta, chi1z, chi2z, finspin, ell, mm);


    // NOTE: There's a problem with this expression from the paper becoming imaginary if gamma2>=1
    // Fix: if gamma2 >= 1 then set the square root term to zero.
    if (gamma2 <= 1)
        return fabs(fRD + (fDM*(-1 + sqrt(1 - pow_2_of(gamma2)))*gamma3)/gamma2);
    else
        return fabs(fRD + (fDM*(-1)*gamma3)/gamma2);
}

/*
 * XLALSimIMRPhenomDHMInspiralFreqScale
 * This function takes as input the frequency 'f' at a given 'm'
 * spherical harmonic mode and returns a scaled frequency,
 * scaled to the m=2 mode assuming a PN, inspiral only like scaling.
 */
double XLALSimIMRPhenomDHMInspiralFreqScale( const REAL8 f, const INT4 mm )
{
    return 2.0 * f / mm;
}

/**
 * XLALSimIMRPhenomDHMFreqDomainMapHM
 * Input waveform frequency in Geometric units (Mf_wf)
 * and computes what frequency this corresponds
 * to scaled to the 22 mode.
 */
double XLALSimIMRPhenomDHMFreqDomainMapHM( const REAL8 Mf_wf,
                                           const INT4 ell,
                                           const INT4 mm,
                                           const REAL8 eta,
                                           const REAL8 chi1z,
                                           const REAL8 chi2z,
                                           const INT4 AmpFlag
                              )
{

    /* compute predicted final spin */
    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);
    REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z); /* dimensionless final spin */
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    // initialise
    REAL8 Mf_22   = 0.; /* the geometric frequency scaled to the 22 mode, NOTE: called f_map in notes */
    REAL8 Mf_1_22  = 0.; /* initalise variable */
    if ( AmpFlag==1 ) {
        /* For amplitude */
        Mf_1_22  = AMP_fJoin_INS; /* inspiral joining frequency from PhenomD [amplitude model], for the 22 mode */
    } else if ( AmpFlag==0 ) {
        /* For phase */
        Mf_1_22  = PHI_fJoin_INS; /* inspiral joining frequency from PhenomD [phase model], for the 22 mode */
    }

    /* FIXME: added this fudge factor so that the discontinuity in the phase derivative occurs after the peak of the phase derivative. */
    /* FIXME PLEASE :) */
    REAL8 FUDGE_FACTOR = 1.3;

    REAL8 Mf_RD_22 = FUDGE_FACTOR * XLALSimIMRPhenomDHMfring(eta, chi1z, chi2z, finspin, 2, 2); /* 22 mode ringdown frequency, geometric units */

    REAL8 Mf_1_lm  = Mf_1_22 * mm / 2.0; /* Convert from 22 to lm, opposite to what XLALSimIMRPhenomDHMInspiralFreqScale does */
    REAL8 Mf_RD_lm = FUDGE_FACTOR * XLALSimIMRPhenomDHMfring(eta, chi1z, chi2z, finspin, ell, mm); /* (ell, mm) ringdown frequency, geometric units */

    /* The following if ladder determines which frequencies get scaled in which way */
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
           /* For amplitude */
           Mf_22 = Mf_wf - Mf_RD_lm + Mf_RD_22 ;
       } else if ( AmpFlag==0 ) {
           /* For phase */
           const REAL8 rho_lm = Mf_RD_22 / Mf_RD_lm ;
           Mf_22 = rho_lm * Mf_wf;
       }
   }

    /* This alternaitve to the above block of code just assumes inspiral PN frequency scaling throughout */
    // Mf_22 = XLALSimIMRPhenomDHMInspiralFreqScale( Mf_wf, mm );

    return Mf_22;
}

/*
 * For a given frequency and ell and m spherical harmonic mode
 * return the frequency scaled to give the leading order PN
 * amplitude for the given ell and m modes.
 * Calcuated from mathematica function: FrequencyPower[f, {ell, m}] / FrequencyPower[f, {2, 2}]
 * FrequencyPower function just returns the leading order PN term in the amplitude.
 */
double XLALSimIMRPhenomDHMPNFrequencyScale( REAL8 Mf,
                                        INT4 ell,
                                        INT4 mm ) {

    /*FIXME: Precompute these powers*/

    /*initialise answer*/
    REAL8 ans = 0.0;

    if ( ell==2 && mm==2 ) {
        ans = 1.0;
    } else if ( ell==2 && mm==1 ) {
        ans = pow(Mf, 1.0/3.0);
    } else if ( ell==3 && mm==3 ) {
        ans = pow(Mf, 1.0/3.0);
    } else if ( ell==3 && mm==2 ) {
        ans = pow(Mf, 2.0/3.0);
    } else if ( ell==4 && mm==4 ) {
        ans = pow(Mf, 2.0/3.0);
    } else if ( ell==4 && mm==3 ) {
        ans = Mf;
    } else if ( ell==5 && mm==5 ) {
        ans = Mf;
    } else if ( ell==5 && mm==4 ) {
        ans = pow(Mf, 4.0/3.0);
    } else if ( ell==6 && mm==6 ) {
        ans = pow(Mf, 4.0/3.0);
    } else if ( ell==6 && mm==5 ) {
        ans = pow(Mf, 5.0/3.0);
    } else {
        XLALPrintError("XLAL Error - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }

    return ans;

}

/* FIXME: returns leading order PN amplitude for given ell and m mode.
 * This is from mma notebook 'leadingPNamp.nb' in /work/projects/PhenomHM
 */
double XLALSimIMRPhenomDHMPNAmplitudeLeadingOrder( REAL8 Mf_wf,
                                                REAL8 eta,
                                                INT4 ell,
                                                INT4 mm ) {

    //taking the absolute value of complex terms

    /*initialise answer*/
    REAL8 ans = 0.0;

    REAL8 dummy = 0.5;

    if ( ell==2 && mm==2 ) {
        ans = 0.674677 * sqrt(eta) * pow(Mf_wf, -7.0/6.0);
    } else if ( ell==2 && mm==1 ) {
        ans = 0.329376 * sqrt( 1.0 - 4.0 * eta ) * sqrt( eta ) * pow(Mf_wf, -5.0/6.0);
    } else if ( ell==3 && mm==3 ) {
        ans = 0.767106 * sqrt( 1.0 - 4.0 * eta ) * sqrt( eta ) * pow(Mf_wf, -5.0/6.0);
    } else if ( ell==3 && mm==2 ) {
        ans = 0.407703 * sqrt( eta ) * pow( Mf_wf, -1.0/2.0);
    } else if ( ell==4 && mm==4 ) {
        ans = ( 1.08721 - 3.26162*eta ) * sqrt( eta ) * pow( Mf_wf, -1.0/2.0);
    } else if ( ell==4 && mm==3 ) {
        ans =  ( 0.570006 - 1.14001*eta ) * sqrt( 1.0 - 4.0 * eta ) * sqrt( eta ) * pow( Mf_wf, -1.0/6.0 );
    } else if ( ell==5 && mm==5 ) {
        ans = 3.39713 * (0.5 - eta) * sqrt( 1.0 - 4.0 * eta ) * sqrt( eta ) * pow(Mf_wf, -1.0/6.0);
    } else if ( ell==5 && mm==4 ) {
        ans = 0.859267 * sqrt(eta) * pow(Mf_wf, 1.0/6.0);
    } else if ( ell==6 && mm==6 ) {
        ans = 2.80361 * sqrt(eta) * pow(Mf_wf, 1.0/6.0);
    } else if ( ell==6 && mm==5 ) {
        ans = 1.36104 * sqrt( 1.0 - 4.0 * eta ) * sqrt(eta) * pow(Mf_wf, 1.0/2.0);
    } else {
        XLALPrintError("XLAL Error - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }

    return ans;
}

double XLALSimIMRPhenomDHMAmplitude( double Mf_wf,
                                    double eta,
                                    double chi1z,
                                    double chi2z,
                                    int ell,
                                    int mm
                                )
{

    int errcode = XLAL_SUCCESS;
    errcode = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_useful_powers() failed.");
    // Calculate phenomenological parameters
    const REAL8 finspin = FinalSpin0815(eta, chi1z, chi2z); // dimensionless final spin used in PhenomD

    if (finspin < MIN_FINAL_SPIN)
            XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                            the model might misbehave here.", finspin);

    IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1z, chi2z, finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    AmpInsPrefactors amp_prefactors;
    errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");

    const INT4 AmpFlagTrue = 1; /* FIXME: Could make this a global variable too */
    double Mf_22 = XLALSimIMRPhenomDHMFreqDomainMapHM( Mf_wf, ell, mm, eta, chi1z, chi2z, AmpFlagTrue );

    UsefulPowers powers_of_Mf_22;
    errcode = init_useful_powers(&powers_of_Mf_22, Mf_22);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_powers failed for Mf_22");

    double PhenDamp = IMRPhenDAmplitude(Mf_22, pAmp, &powers_of_Mf_22, &amp_prefactors);

    /* compute amplitude ratio correction to take 22 mode in to (ell, mm) mode amplitude */
    double MfAtScale_wf = 0.0001; /* FIXME: This should be made a global variable in header. */
    double MfAtScale_22 = XLALSimIMRPhenomDHMFreqDomainMapHM( MfAtScale_wf, ell, mm, eta, chi1z, chi2z, AmpFlagTrue );

    UsefulPowers powers_of_MfAtScale_22;
    errcode = init_useful_powers(&powers_of_MfAtScale_22, MfAtScale_22);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_powers failed for MfAtScale_22");

    /* see technical document for description of below lines with A_R and R */
    double A_R_num = XLALSimIMRPhenomDHMPNAmplitudeLeadingOrder( MfAtScale_wf, eta, ell, mm );
    double A_R_den = XLALSimIMRPhenomDHMPNFrequencyScale(MfAtScale_22, ell, mm) * IMRPhenDAmplitude(MfAtScale_22, pAmp, &powers_of_MfAtScale_22, &amp_prefactors);
    double ampRatio = A_R_num/A_R_den;
    double R = ampRatio * XLALSimIMRPhenomDHMPNFrequencyScale(Mf_22, ell, mm);

    double HMamp = PhenDamp * R;

    LALFree(pAmp);

    return HMamp;
}

double XLALSimIMRPhenomDHMPhase( double Mf_wf,
                                    double eta,
                                    double chi1z,
                                    double chi2z,
                                    int ell,
                                    int mm
                                )
{

    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);
    REAL8 M = m1 + m2;

    LALDict *extraParams = NULL;

    int status = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initiate useful powers of pi.");


    // Calculate phenomenological parameters
    const REAL8 finspin = FinalSpin0815(eta, chi1z, chi2z); //FinalSpin0815 - 0815 is like a version number

    if (finspin < MIN_FINAL_SPIN)
          XLAL_PRINT_WARNING("Final spin (Mf=%g) and ISCO frequency of this system are small, \
                          the model might misbehave here.", finspin);

    if (extraParams==NULL)
    extraParams=XLALCreateDict();
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    IMRPhenomDPhaseCoefficients *pPhi = ComputeIMRPhenomDPhaseCoefficients(eta, chi1z, chi2z, finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1, m2, chi1z, chi2z, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSS(m1, m2, M, chi1z, chi2z) * pn->v[0]);

    PhiInsPrefactors phi_prefactors;
    status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

    // Compute coefficients to make phase C^1 continuous (phase and first derivative)
    ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors);

    const INT4 AmpFlagFalse = 0; /* FIXME: Could make this a global variable too */
    double Mf_22 = XLALSimIMRPhenomDHMFreqDomainMapHM( Mf_wf, ell, mm, eta, chi1z, chi2z, AmpFlagFalse );

    UsefulPowers powers_of_Mf_22;
    status = init_useful_powers(&powers_of_Mf_22, Mf_22);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf_22 failed");

    /* phi_lm(f) = m * phi_22(f_22) / 2.0 */
    double PhenDphase = mm * IMRPhenDPhase(Mf_22, pPhi, pn, &powers_of_Mf_22, &phi_prefactors) / 2.0;

    LALFree(pPhi);
    LALFree(pn);

    return PhenDphase;
}

int XLALSimIMRPhenomDHMCoreOneMode(
    COMPLEX16FrequencySeries **hptilde,
    COMPLEX16FrequencySeries **hctilde,
    const REAL8 deltaF,                /**< frequency resolution */
    const REAL8 f_min,                 /**< start frequency */
    const REAL8 f_max,                 /**< end frequency */
    const REAL8 m1_in,                 /**< Mass of companion 1 (kg) */
    const REAL8 m2_in,                 /**< Mass of companion 2 (kg) */
    double chi1z_in,
    double chi2z_in,
    int ell,
    int mm
)   {

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* sanity checks on input parameters */
    XLAL_CHECK(0 != hptilde, XLAL_EFAULT, "hptilde is null");
    XLAL_CHECK(0 != hctilde, XLAL_EFAULT, "hctilde is null");
    if (*hptilde) XLAL_ERROR(XLAL_EFAULT);
    if (*hctilde) XLAL_ERROR(XLAL_EFAULT);

    if (m1_in <= 0) XLAL_ERROR(XLAL_EDOM, "m1_in must be positive\n");
    if (m2_in <= 0) XLAL_ERROR(XLAL_EDOM, "m2_in must be positive\n");

    if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
    if (f_min <= 0) XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
    if (f_max < 0) XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");

    /* ensure that m1 is the larger mass and chi1z is the spin on m1 */
    REAL8 chi1z, chi2z, m1_tmp, m2_tmp;
    if (m1_in>m2_in) {
       chi1z = chi1z_in;
       chi2z = chi2z_in;
       m1_tmp   = m1_in;
       m2_tmp   = m2_in;
   } else { /* swap spins and masses */
       chi1z = chi2z_in;
       chi2z = chi1z_in;
       m1_tmp   = m2_in;
       m2_tmp   = m1_in;
    }

    /* external: SI; internal: solar masses */
    const REAL8 m1 = m1_tmp / LAL_MSUN_SI;
    const REAL8 m2 = m2_tmp / LAL_MSUN_SI;

    const REAL8 q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

    if (q > MAX_ALLOWED_MASS_RATIO)
      XLAL_PRINT_WARNING("Warning: The model is not supported for high mass ratio, see MAX_ALLOWED_MASS_RATIO\n");

    const REAL8 M = m1 + m2;
    const REAL8 M_sec = M * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
    const REAL8 eta = m1 * m2 / (M * M);

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    if (chi1z > 1.0 || chi1z < -1.0 || chi2z > 1.0 || chi2z < -1.0)
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    /* Coalesce at t=0 */
    // shift by overall length in time
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, -1. / deltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    /* Allocate htilde */
    size_t n = NextPow2(f_max / deltaF) + 1;

    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
    XLAL_CHECK ( *hptilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu for f_max=%f, deltaF=%g.", n, f_max, deltaF);

    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
    XLAL_CHECK ( *hctilde, XLAL_ENOMEM, "Failed to allocated waveform COMPLEX16FrequencySeries of length %zu for f_max=%f, deltaF=%g.", n, f_max, deltaF);


    memset((*hptilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitMultiply(&((*hptilde)->sampleUnits), &((*hptilde)->sampleUnits), &lalSecondUnit);

    memset((*hctilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitMultiply(&((*hctilde)->sampleUnits), &((*hctilde)->sampleUnits), &lalSecondUnit);

    /* range that will have actual non-zero waveform values generated */
    size_t ind_min = (size_t) (f_min / deltaF);
    size_t ind_max = (size_t) (f_max / deltaF);
    XLAL_CHECK ( (ind_max<=n) && (ind_min<=ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=htilde->data>length=%zu.", ind_min, ind_max, n);


    /* Now generate the waveform */
    #pragma omp parallel for
    for (size_t i = ind_min; i < ind_max; i++)
    {
      REAL8 Mf_wf = M_sec * i * deltaF; // geometric frequency

    //   REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
    //   REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors);
      double HMamp = XLALSimIMRPhenomDHMAmplitude( Mf_wf, eta, chi1z, chi2z, ell, mm );
      double HMphase = XLALSimIMRPhenomDHMPhase( Mf_wf, eta, chi1z, chi2z, ell, mm );

      // phi -= t0*(Mf-MfRef) + phi_precalc;
      // ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
      ((*hptilde)->data->data)[i] = HMamp * cexp(-I * HMphase);
      ((*hctilde)->data->data)[i] = ((*hptilde)->data->data)[i];

    }

    return XLAL_SUCCESS;
}
