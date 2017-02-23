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
#include <lal/SphericalHarmonics.h>

/* This allows us to reuse internal IMRPhenomD functions without making those functions XLAL */
#include "LALSimIMRPhenomD_internals.c"
#include "LALSimRingdownCW.c"


#include "LALSimIMRPhenomHM.h"

#ifndef _OPENMP
#define omp ignore
#endif

// #define NMODES_MAX 8
#define NMODES_MAX 5
// #define NMODES_MAX 4

static int NMODES = NMODES_MAX;
static const int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {4,4}, {5,5} };


/* dimensionless frequency of last data point in waveform */
#define Mf_CUT_HM 0.5




/**
 * Lionel's QNM higher mode fits
 * returns real part of the ringdown frequency
 */
double XLALSimIMRPhenomHMfring(const REAL8 eta,     /**< symmetric mass-ratio */
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
double XLALSimIMRPhenomHMfdamp(const REAL8 eta,     /**< symmetric mass-ratio */
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



/* START: newer functions  */

//mathematica function postfRDflm
double XLALIMRPhenomHMpostfRDflm(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag){
    double ans = 0.0;
    if ( AmpFlag==1 ) {
        /* For amplitude */
        ans = Mf-Mf_RD_lm+Mf_RD_22; /*Used for the Amplitude as an approx fix for post merger powerlaw slope */
    } else if ( AmpFlag==0 ) {
        /* For phase */
        REAL8 Rholm = Mf_RD_22/Mf_RD_lm;
        ans = Rholm * Mf;          /* Used for the Phase */
    }

    return ans;
}

// mathematica function Ti
// domain mapping function - inspiral
double XLALIMRPhenomHMTi(REAL8 Mf, const INT4 mm){
    return 2.0 * Mf / mm;
}

// mathematica function Trd
// domain mapping function - ringdown
double XLALIMRPhenomHMTrd(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag){
    return XLALIMRPhenomHMpostfRDflm(Mf, Mf_RD_22, Mf_RD_lm, AmpFlag);
}

// mathematica function Tm
// domain mapping function - intermediate
double XLALIMRPhenomHMTm(REAL8 Mf, const INT4 mm, REAL8 fi, REAL8 fr, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag){
    REAL8 Trd = XLALIMRPhenomHMTrd(fr, Mf_RD_22, Mf_RD_lm, AmpFlag);
    REAL8 Ti = XLALIMRPhenomHMTi(fi, mm);
    return (Mf - fi) * ( Trd - Ti )/( fr - fi ) + Ti;
}

double XLALIMRPhenomHMSlopeAm(const INT4 mm, REAL8 fi, REAL8 fr, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag){
    //Am = ( Trd[fr]-Ti[fi] )/( fr - fi );
    REAL8 Trd = XLALIMRPhenomHMTrd(fr, Mf_RD_22, Mf_RD_lm, AmpFlag);
    REAL8 Ti = XLALIMRPhenomHMTi(fi, mm);

    return ( Trd-Ti )/( fr - fi );
}

double XLALIMRPhenomHMSlopeBm(const INT4 mm, REAL8 fi, REAL8 fr, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag){
    //Bm = Ti[fi] - fi*Am;
    REAL8 Ti = XLALIMRPhenomHMTi(fi, mm);
    REAL8 Am = XLALIMRPhenomHMSlopeAm(mm, fi, fr, Mf_RD_22, Mf_RD_lm, AmpFlag);
    return Ti - fi*Am;
}

int XLALIMRPhenomHMMapParams(REAL8 *a, REAL8 *b, REAL8 flm, REAL8 fi, REAL8 fr, REAL8 Ai, REAL8 Bi, REAL8 Am, REAL8 Bm, REAL8 Ar, REAL8 Br){
    // Defne function to output map params used depending on
    if ( flm <= fi ){
        *a = Ai;
        *b = Bi;
    } else if ( fi < flm && flm <= fr ){
        *a = Am;
        *b = Bm;
    } else if ( fr < flm ){
        *a = Ar;
        *b = Br;
    };
    return XLAL_SUCCESS;
}

/**
 * XLALSimIMRPhenomHMRholm = MfRD22/MfRDlm.
 * Ratio of the (l,m)=(2,2) ringdown frequency to input (l.m) mode.
 */
double XLALSimIMRPhenomHMRholm(REAL8 eta, REAL8 chi1z, REAL8 chi2z, const INT4 ell, const INT4 mm){
    /* compute predicted final spin */
    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);
    REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z); /* dimensionless final spin */
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");
    REAL8 Mf_RD_22 = XLALSimIMRPhenomHMfring(eta, chi1z, chi2z, finspin, 2, 2); /* 22 mode ringdown frequency (real part of ringdown), geometric units */
    REAL8 Mf_RD_lm = XLALSimIMRPhenomHMfring(eta, chi1z, chi2z, finspin, ell, mm);
    REAL8 Rholm = Mf_RD_22/Mf_RD_lm;
    return Rholm;
}

/**
 * XLALSimIMRPhenomHMTaulm = MfDMlm/MfDM22.
 * Ratio of the input (l.m) mode ringdown damping time to (l,m)=(2,2).
 */
double XLALSimIMRPhenomHMTaulm(REAL8 eta, REAL8 chi1z, REAL8 chi2z, const INT4 ell, const INT4 mm){
    /* compute predicted final spin */
    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);
    REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z); /* dimensionless final spin */
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");
    REAL8 Mf_DM_22 = XLALSimIMRPhenomHMfdamp(eta, chi1z, chi2z, finspin, 2, 2); /* 22 mode ringdown frequency (real part of ringdown), geometric units */
    REAL8 Mf_DM_lm = XLALSimIMRPhenomHMfdamp(eta, chi1z, chi2z, finspin, ell, mm);
    REAL8 Taulm = Mf_DM_lm/Mf_DM_22;
    return Taulm;
}

int XLALIMRPhenomHMFreqDomainMapParams( REAL8 *a,/**< [Out]  */
                                         REAL8 *b,/**< [Out]  */
                                         REAL8 *fi,/**< [Out]  */
                                         REAL8 *fr,/**< [Out]  */
                                         REAL8 *f1,/**< [Out]  */
                                         REAL8 *f2lm,/**< [Out]  */
                                         const REAL8 flm, /**< input waveform frequency */
                                         const INT4 ell, /**< spherical harmonics ell mode */
                                         const INT4 mm, /**< spherical harmonics m mode */
                                         const REAL8 eta, /**< symmetric mass ratio */
                                         const REAL8 chi1z, /**< aligned spin on larger body */
                                         const REAL8 chi2z, /**< aligned spin on smaller body */
                                         const int AmpFlag /**< is ==1 then computes for amplitude, if ==0 then computes for phase */)
{

    /*check output points are NULL*/
    XLAL_CHECK(a != NULL, XLAL_EFAULT);
    XLAL_CHECK(b != NULL, XLAL_EFAULT);
    XLAL_CHECK(fi != NULL, XLAL_EFAULT);
    XLAL_CHECK(fr != NULL, XLAL_EFAULT);
    XLAL_CHECK(f1 != NULL, XLAL_EFAULT);
    XLAL_CHECK(f2lm != NULL, XLAL_EFAULT);

    /* compute predicted final spin */
    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);
    REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z); /* dimensionless final spin */
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    // Account for different f1 definition between PhenomD Amplitude and Phase derivative models
    // initialise
    REAL8 Mf_1_22  = 0.; /* initalise variable */
    if ( AmpFlag==1 ) {
        /* For amplitude */
        Mf_1_22  = AMP_fJoin_INS; /* inspiral joining frequency from PhenomD [amplitude model], for the 22 mode */
    } else if ( AmpFlag==0 ) {
        /* For phase */
        Mf_1_22  = PHI_fJoin_INS; /* inspiral joining frequency from PhenomD [phase model], for the 22 mode */
    }

    *f1 = Mf_1_22;

    REAL8 Mf_RD_22 = XLALSimIMRPhenomHMfring(eta, chi1z, chi2z, finspin, 2, 2); /* 22 mode ringdown frequency (real part of ringdown), geometric units */
    // REAL8 Mf_DM_22 = XLALSimIMRPhenomHMfdamp(eta, chi1z, chi2z, finspin, 2, 2); /* (2, 2) damping time (complex part of ringdown), geometric units */

    REAL8 Mf_RD_lm = XLALSimIMRPhenomHMfring(eta, chi1z, chi2z, finspin, ell, mm);
    // REAL8 Mf_DM_lm = XLALSimIMRPhenomHMfdamp(eta, chi1z, chi2z, finspin, ell, mm);

    // Define a ratio of QNM frequencies to be used for scaling various quantities
	REAL8 Rholm = Mf_RD_22/Mf_RD_lm;

    // Given experiments with the l!=m modes, it appears that the QNM scaling rather than the PN scaling may be optimal for mapping f1
    REAL8 Mf_1_lm = Mf_1_22 / Rholm;

    // (* Handle cases for post-f3 mapping *)
    // REAL8 postfRDflm = XLALIMRPhenomHMpostfRDflm(flm, Mf_RD_22, Mf_RD_lm, AmpFlag);

    /* Define transition frequencies */
	*fi = Mf_1_lm;
	*fr = Mf_RD_lm;

    // (* Define functions to be applied in each domain *)
    // REAL8 Ti = XLALIMRPhenomHMTi(flm, mm);
    // REAL8 Trd = XLALIMRPhenomHMTrd(flm, Mf_RD_22, Mf_RD_lm, AmpFlag);
    // REAL8 Tm = XLALIMRPhenomHMTm(flm, mm, *fi, *fr, Mf_RD_22, Mf_RD_lm, AmpFlag);

    /*Define the slope and intercepts of the linear transformation used*/
	REAL8 Ai = 2.0/mm;
    REAL8 Bi = 0.0;
	REAL8 Am = XLALIMRPhenomHMSlopeAm(mm, *fi, *fr, Mf_RD_22, Mf_RD_lm, AmpFlag);
    REAL8 Bm = XLALIMRPhenomHMSlopeBm(mm, *fi, *fr, Mf_RD_22, Mf_RD_lm, AmpFlag);


    REAL8 Ar = 0.0;
    REAL8 Br = 0.0;
    if ( AmpFlag==1 ) {
        /* For amplitude */
        Ar = 1.0;
        Br = -Mf_RD_lm+Mf_RD_22;
    } else if ( AmpFlag==0 ) {
        /* For phase */
        Ar = Rholm;
        Br = 0.0;
    }

    // Defne function to output map params used depending on
    // *a = 0.0;
    // *b = 0.0;
    int ret = XLALIMRPhenomHMMapParams(a, b, flm, *fi, *fr, Ai, Bi, Am, Bm, Ar, Br);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMMapParams failed in XLALIMRPhenomHMFreqDomainMapParams (1)\n");
        XLAL_ERROR(XLAL_EDOM);
    }


    REAL8 a2 = 0.0;
    REAL8 b2 = 0.0;
    REAL8 frfi = ( *fr + *fi ) / 2.0;
    ret = XLALIMRPhenomHMMapParams(&a2, &b2, frfi, *fi, *fr, Ai, Bi, Am, Bm, Ar, Br);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMMapParams failed in XLALIMRPhenomHMFreqDomainMapParams (2)\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    *f2lm = ( (Mf_RD_22 / 2.0) - b2 ) / a2;

    return XLAL_SUCCESS;
}

/**
 * XLALSimIMRPhenomHMFreqDomainMap
 * Input waveform frequency in Geometric units (Mflm)
 * and computes what frequency this corresponds
 * to scaled to the 22 mode.
 */
double XLALSimIMRPhenomHMFreqDomainMap(REAL8 Mflm,
                                        const INT4 ell,
                                        const INT4 mm,
                                        const REAL8 eta,
                                        const REAL8 chi1z,
                                        const REAL8 chi2z,
                                        const int AmpFlag)
{

    //Mflm here has the same meaning as Mf_wf in XLALSimIMRPhenomHMFreqDomainMapHM (old deleted function).

    REAL8 a = 0.;
    REAL8 b = 0.;
    /*following variables not used in this funciton but are returned in XLALIMRPhenomHMFreqDomainMapParams*/
    REAL8 fi = 0.;
    REAL8 fr = 0.;
    REAL8 f1 = 0.;
    REAL8 f2lm = 0.;
    int ret = XLALIMRPhenomHMFreqDomainMapParams(&a, &b, &fi, &fr, &f1, &f2lm, Mflm, ell, mm, eta, chi1z, chi2z, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALSimIMRPhenomHMFreqDomainMap\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    REAL8 Mf22 = a * Mflm + b;
    return Mf22;
}

/* END: newer functions  */




/*
 * For a given frequency and ell and m spherical harmonic mode
 * return the frequency scaled to give the leading order PN
 * amplitude for the given ell and m modes.
 * Calcuated from mathematica function: FrequencyPower[f, {ell, m}] / FrequencyPower[f, {2, 2}]
 * FrequencyPower function just returns the leading order PN term in the amplitude.
 */
double XLALSimIMRPhenomHMPNFrequencyScale( REAL8 Mf,
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
double XLALSimIMRPhenomHMPNAmplitudeLeadingOrder( REAL8 Mf_wf,
                                                REAL8 eta,
                                                INT4 ell,
                                                INT4 mm ) {

    //taking the absolute value of complex terms

    /*initialise answer*/
    REAL8 ans = 0.0;

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


static double ComputeAmpRatio(REAL8 eta, REAL8 chi1z, REAL8 chi2z, INT4 ell, INT4 mm, AmpInsPrefactors amp_prefactors, IMRPhenomDAmplitudeCoefficients *pAmp);
static double ComputeAmpRatio(REAL8 eta, REAL8 chi1z, REAL8 chi2z, INT4 ell, INT4 mm, AmpInsPrefactors amp_prefactors, IMRPhenomDAmplitudeCoefficients *pAmp){
    /* compute amplitude ratio correction to take 22 mode in to (ell, mm) mode amplitude */
    const INT4 AmpFlagTrue = 1; /* FIXME: Could make this a global variable too */
    double MfAtScale_wf_amp = 0.0001; /* FIXME: This should be made a global variable in header. */
    double MfAtScale_22_amp = XLALSimIMRPhenomHMFreqDomainMap( MfAtScale_wf_amp, ell, mm, eta, chi1z, chi2z, AmpFlagTrue );

    UsefulPowers powers_of_MfAtScale_22_amp;
    int errcode = init_useful_powers(&powers_of_MfAtScale_22_amp, MfAtScale_22_amp);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_powers failed for MfAtScale_22_amp");

    /* see technical document for description of below lines with A_R and R */
    double A_R_num = XLALSimIMRPhenomHMPNAmplitudeLeadingOrder( MfAtScale_wf_amp, eta, ell, mm );
    double A_R_den = XLALSimIMRPhenomHMPNFrequencyScale(MfAtScale_22_amp, ell, mm) * IMRPhenDAmplitude(MfAtScale_22_amp, pAmp, &powers_of_MfAtScale_22_amp, &amp_prefactors);
    double ampRatio = A_R_num/A_R_den;

    return ampRatio;
}

double XLALSimIMRPhenomHMAmplitude( double Mf_wf,
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
    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);

    const REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z);
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1z, chi2z, finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    AmpInsPrefactors amp_prefactors;
    errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");

    const INT4 AmpFlagTrue = 1; /* FIXME: Could make this a global variable too */
    double Mf_22 =  XLALSimIMRPhenomHMFreqDomainMap(Mf_wf, ell, mm, eta, chi1z, chi2z, AmpFlagTrue);

    UsefulPowers powers_of_Mf_22;
    errcode = init_useful_powers(&powers_of_Mf_22, Mf_22);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_powers failed for Mf_22");

    double PhenDamp = IMRPhenDAmplitude(Mf_22, pAmp, &powers_of_Mf_22, &amp_prefactors);

    double ampRatio = ComputeAmpRatio(eta, chi1z, chi2z, ell, mm, amp_prefactors, pAmp);

    double R = ampRatio * XLALSimIMRPhenomHMPNFrequencyScale(Mf_22, ell, mm);

    double HMamp = PhenDamp * R;

    LALFree(pAmp);

    return HMamp;
}

double XLALSimIMRPhenomHMAmplitudeOpt( double Mf_wf, double eta, double chi1z, double chi2z, int ell, int mm, IMRPhenomDAmplitudeCoefficients *pAmp, AmpInsPrefactors * amp_prefactors );

double XLALSimIMRPhenomHMAmplitudeOpt( double Mf_wf,
                                    double eta,
                                    double chi1z,
                                    double chi2z,
                                    int ell,
                                    int mm,
                                    IMRPhenomDAmplitudeCoefficients *pAmp,
                                    AmpInsPrefactors * amp_prefactors
                                )
{

    const INT4 AmpFlagTrue = 1; /* FIXME: Could make this a global variable too */
    double Mf_22 =  XLALSimIMRPhenomHMFreqDomainMap(Mf_wf, ell, mm, eta, chi1z, chi2z, AmpFlagTrue);

    UsefulPowers powers_of_Mf_22;
    int errcode = init_useful_powers(&powers_of_Mf_22, Mf_22);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_powers failed for Mf_22");

    double PhenDamp = IMRPhenDAmplitude(Mf_22, pAmp, &powers_of_Mf_22, amp_prefactors);

    double ampRatio = ComputeAmpRatio(eta, chi1z, chi2z, ell, mm, *amp_prefactors, pAmp);

    double R = ampRatio * XLALSimIMRPhenomHMPNFrequencyScale(Mf_22, ell, mm);

    double HMamp = PhenDamp * R;

    // LALFree(pAmp);

    return HMamp;
}


/*TODO Should also probably add LALDict as a argument here.*/
int XLALSimIMRPhenomHMPhasePreComp(HMPhasePreComp *q, const INT4 ell, const INT4 mm, const REAL8 eta, const REAL8 chi1z, const REAL8 chi2z, const double finspin)
{
    /*
    { ai,bi,fi,fr,f1,f2lm,fRD22,fRDlm,Ti,Trd,Tm } = FreqDomainMapParams[ 0.0001, mode, eta,chi1,chi2,False ];
	{ a2lm,b2lm,fi,fr,f1,f2lm,fRD22,fRDlm,Ti,Trd,Tm } = FreqDomainMapParams[ f2lm+0.0001, mode, eta,chi1,chi2,False ];
	{ ar,br,fi,fr,f1,f2lm,fRD22,fRDlm,Ti,Trd,Tm } = FreqDomainMapParams[ fr+0.0001, mode, eta,chi1,chi2,False ];

	PhDBconst = IMRPhenDPhasev2[ a2lm*fi   + b2lm, eta,chi1,chi2, mode ]/a2lm;
	PhDCconst = IMRPhenDPhasev2[ a2lm*f2lm + b2lm, eta,chi1,chi2, mode ]/a2lm;
	PhDDconst = IMRPhenDPhasev2[ ar*fr + br, eta,chi1,chi2, mode]/ar;

	PhDBAterm = IMRPhenDPhasev2[ ai*fi + bi, eta,chi1,chi2, mode ]/ai;

    Return[{ai,bi,a2lm,b2lm,ar,br,fi,f2lm,fr,PhDBconst,PhDCconst,PhDDconst,PhDBAterm}];
    */

    REAL8 ai = 0.0;
    REAL8 bi = 0.0;
    REAL8 a2lm = 0.0;
    REAL8 b2lm = 0.0;
    REAL8 ar = 0.0;
    REAL8 br = 0.0;
    REAL8 fi = 0.0;
    REAL8 f1 = 0.0;
    REAL8 f2lm = 0.0;
    REAL8 fr = 0.0;

    const INT4 AmpFlag = 0;

    /* NOTE: As long as Mfshit + f2lm isn't >= fr then the value of the shift is arbitrary. */
    const REAL8 Mfshift = 0.0001;

    int ret = XLALIMRPhenomHMFreqDomainMapParams(&ai, &bi, &fi, &fr, &f1, &f2lm, Mfshift, ell, mm, eta, chi1z, chi2z, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALSimIMRPhenomHMPhasePreComp - inspiral\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->ai = ai;
    q->bi = bi;


    ret = XLALIMRPhenomHMFreqDomainMapParams(&a2lm, &b2lm, &fi, &fr, &f1, &f2lm, f2lm+Mfshift, ell, mm, eta, chi1z, chi2z, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALSimIMRPhenomHMPhasePreComp - intermediate\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->a2lm = a2lm;
    q->b2lm = b2lm;

    ret = XLALIMRPhenomHMFreqDomainMapParams(&ar, &br, &fi, &fr, &f1, &f2lm, fr+Mfshift, ell, mm, eta, chi1z, chi2z, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALSimIMRPhenomHMPhasePreComp - merger-ringdown\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->ar = ar;
    q->br = br;

    q->fi = fi;
    q->fr = fr;
    q->f2lm = f2lm;

    LALDict *extraParams = NULL;

    int status = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "Failed to initiate useful powers of pi.");

    if (extraParams==NULL)
    extraParams=XLALCreateDict();
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    IMRPhenomDPhaseCoefficients *pPhi = ComputeIMRPhenomDPhaseCoefficients(eta, chi1z, chi2z, finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);

    // Convention m1 >= m2
    // FIXME: change input function args to always be m1, m2, chi1z, chi2z and never eta!
    REAL8 Seta = sqrt(1.0 - 4.0*eta);
    REAL8 m1 = 0.5 * (1.0 + Seta);
    REAL8 m2 = 0.5 * (1.0 - Seta);
    REAL8 M = m1 + m2;

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

    double Rholm = XLALSimIMRPhenomHMRholm(eta, chi1z, chi2z, ell, mm);
    double Taulm = XLALSimIMRPhenomHMTaulm(eta, chi1z, chi2z, ell, mm);

    // Compute coefficients to make phase C^1 continuous (phase and first derivative)
    ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);

    REAL8 PhDBMf = a2lm*fi + b2lm;
    UsefulPowers powers_of_PhDBMf;
    status = init_useful_powers(&powers_of_PhDBMf, PhDBMf);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_PhDBMf failed");
    q->PhDBconst = IMRPhenDPhase(PhDBMf, pPhi, pn, &powers_of_PhDBMf, &phi_prefactors, Rholm, Taulm)/a2lm;

    REAL8 PhDCMf = a2lm*f2lm + b2lm;
    UsefulPowers powers_of_PhDCMf;
    status = init_useful_powers(&powers_of_PhDCMf, PhDCMf);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_PhDCMf failed");
    q->PhDCconst = IMRPhenDPhase(PhDCMf, pPhi, pn, &powers_of_PhDCMf, &phi_prefactors, Rholm, Taulm)/a2lm;

    REAL8 PhDDMf = ar*fr + br;
    UsefulPowers powers_of_PhDDMf;
    status = init_useful_powers(&powers_of_PhDDMf, PhDDMf);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_PhDDMf failed");
    q->PhDDconst = IMRPhenDPhase(PhDDMf, pPhi, pn, &powers_of_PhDDMf, &phi_prefactors, Rholm, Taulm)/ar;

    REAL8 PhDBAMf = ai*fi + bi;
    UsefulPowers powers_of_PhDBAMf;
    status = init_useful_powers(&powers_of_PhDBAMf, PhDBAMf);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_PhDBAMf failed");
    q->PhDBAterm = IMRPhenDPhase(PhDBAMf, pPhi, pn, &powers_of_PhDBAMf, &phi_prefactors, Rholm, Taulm)/ai;

    return XLAL_SUCCESS;

}
double XLALSimIMRPhenomHMPhase( double Mf_wf, double eta, double chi1z, double chi2z, int ell, int mm, HMPhasePreComp *q );

double XLALSimIMRPhenomHMPhase( double Mf_wf, /**< input frequency in geometric units*/
                                    double eta,
                                    double chi1z,
                                    double chi2z,
                                    int ell,
                                    int mm,
                                    HMPhasePreComp *q
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
    const REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z);
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

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

    double Rholm = XLALSimIMRPhenomHMRholm(eta, chi1z, chi2z, ell, mm);
    double Taulm = XLALSimIMRPhenomHMTaulm(eta, chi1z, chi2z, ell, mm);

    // Compute coefficients to make phase C^1 continuous (phase and first derivative)
    ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);

    // const INT4 AmpFlagFalse = 0; /* FIXME: Could make this a global variable too */
    // double Mf_22 = XLALSimIMRPhenomHMFreqDomainMap( Mf_wf, ell, mm, eta, chi1z, chi2z, AmpFlagFalse );


    // UsefulPowers powers_of_Mf_22;
    // status = init_useful_powers(&powers_of_Mf_22, Mf_22);
    // XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf_22 failed");

    /* phi_lm(f) = m * phi_22(f_22) / 2.0 */
    // double PhenDphase = mm * IMRPhenDPhase(Mf_22, pPhi, pn, &powers_of_Mf_22, &phi_prefactors, Rholm, Taulm) / 2.0;


    REAL8 Mf = 0.0;
    REAL8 Mf2lm = 0.0;
    REAL8 Mfr = 0.0;
    REAL8 retphase = 0.0;
    REAL8 tmpphaseB = 0.0;
    REAL8 tmpphaseC = 0.0;

    // This if ladder is in the mathematica function HMPhase. PhenomHMDev.nb

    if ( Mf_wf <= q->fi ){
        Mf = q->ai * Mf_wf + q->bi;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, &phi_prefactors, Rholm, Taulm) / q->ai;
    } else if ( q->fi < Mf_wf && Mf_wf <= q->f2lm ){
        Mf = q->a2lm*Mf_wf + q->b2lm;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, &phi_prefactors, Rholm, Taulm) / q->a2lm - q->PhDBconst + q->PhDBAterm;
    } else if ( q->f2lm < Mf_wf && Mf_wf <= q->fr ){
        Mf2lm = q->a2lm*q->f2lm + q->b2lm;
        UsefulPowers powers_of_Mf2lm;
        status = init_useful_powers(&powers_of_Mf2lm, Mf2lm);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf2lm failed");
        tmpphaseB = IMRPhenDPhase(Mf2lm, pPhi, pn, &powers_of_Mf2lm, &phi_prefactors, Rholm, Taulm) / q->a2lm;
        Mf = q->a2lm*Mf_wf + q->b2lm;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, &phi_prefactors, Rholm, Taulm) / q->a2lm - q->PhDCconst + tmpphaseB - q->PhDBconst + q->PhDBAterm;
    } else if ( Mf_wf > q->fr ) {
        Mf2lm = q->a2lm*q->f2lm + q->b2lm;
        UsefulPowers powers_of_Mf2lm;
        status = init_useful_powers(&powers_of_Mf2lm, Mf2lm);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf2lm failed");
        tmpphaseB = IMRPhenDPhase(Mf2lm, pPhi, pn, &powers_of_Mf2lm, &phi_prefactors, Rholm, Taulm) / q->a2lm;

        Mfr = q->a2lm*q->fr + q->b2lm;
        UsefulPowers powers_of_Mfr;
        status = init_useful_powers(&powers_of_Mfr, Mfr);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mfr failed");
        tmpphaseC = IMRPhenDPhase(Mfr, pPhi, pn, &powers_of_Mfr, &phi_prefactors, Rholm, Taulm) / q->a2lm - q->PhDCconst + tmpphaseB - q->PhDBconst + q->PhDBAterm;;

        Mf = q->ar*Mf_wf + q->br;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, &phi_prefactors, Rholm, Taulm) / q->ar - q->PhDDconst + tmpphaseC;
    } else {
        XLALPrintError("XLAL_ERROR - should not get here - in function XLALSimIMRPhenomHMPhase");
        XLAL_ERROR(XLAL_EDOM);
    }

    /* phase shift due to leading order complex amplitude
    Luc Blancet 1310.1528. Section 9.5
    "Spherical hrmonic modes for numerical relativity" */
    /*initialise answer*/
    REAL8 cShift = 0.0;

    if ( ell==2 && mm==1 ) {
        cShift = LAL_PI/2.0; /* i shift */
    } else if ( ell==3 && mm==3 ) {
        cShift = -LAL_PI/2.0; /* -i shift */
    } else if ( ell==4 && mm==4 ) {
        cShift = LAL_PI; /* -1 shift */
    } else if ( ell==4 && mm==3 ) {
        cShift = -LAL_PI/2.0; /* -i shift */
    } else if ( ell==5 && mm==5 ) {
        cShift = LAL_PI/2.0; /* i shift */
    } else if ( ell==5 && mm==4 ) {
        cShift = LAL_PI; /* -1 shift */
    } else if ( ell==6 && mm==5 ) {
        cShift = LAL_PI/2.0; /* i shift */
    } /*else {
        XLALPrintError("XLAL Error (in function XLALSimIMRPhenomHMPhase) - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }*/
    /*TODO: Need to have a list at the begining and a function to check the input
    lm mode to see if it is one that is included in the model.*/

    retphase += cShift;

    LALFree(pPhi);
    LALFree(pn);

    return retphase;
}

double XLALSimIMRPhenomHMPhaseOpt( double Mf_wf, int ell, int mm, HMPhasePreComp *q, PNPhasingSeries *pn, IMRPhenomDPhaseCoefficients *pPhi, PhiInsPrefactors * phi_prefactors, double Rholm, double Taulm );

double XLALSimIMRPhenomHMPhaseOpt( double Mf_wf, /**< input frequency in geometric units*/
                                    int ell,
                                    int mm,
                                    HMPhasePreComp *q,
                                    PNPhasingSeries *pn,
                                    IMRPhenomDPhaseCoefficients *pPhi,
                                    PhiInsPrefactors * phi_prefactors,
                                    double Rholm,
                                    double Taulm
                                )
{



    // const INT4 AmpFlagFalse = 0; /* FIXME: Could make this a global variable too */
    // double Mf_22 = XLALSimIMRPhenomHMFreqDomainMap( Mf_wf, ell, mm, eta, chi1z, chi2z, AmpFlagFalse );


    // UsefulPowers powers_of_Mf_22;
    // status = init_useful_powers(&powers_of_Mf_22, Mf_22);
    // XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf_22 failed");

    /* phi_lm(f) = m * phi_22(f_22) / 2.0 */
    // double PhenDphase = mm * IMRPhenDPhase(Mf_22, pPhi, pn, &powers_of_Mf_22, &phi_prefactors, Rholm, Taulm) / 2.0;


    REAL8 Mf = 0.0;
    REAL8 Mf2lm = 0.0;
    REAL8 Mfr = 0.0;
    REAL8 retphase = 0.0;
    REAL8 tmpphaseB = 0.0;
    REAL8 tmpphaseC = 0.0;
    int status = 0;

    // This if ladder is in the mathematica function HMPhase. PhenomHMDev.nb

    if ( Mf_wf <= q->fi ){
        Mf = q->ai * Mf_wf + q->bi;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->ai;
    } else if ( q->fi < Mf_wf && Mf_wf <= q->f2lm ){
        Mf = q->a2lm*Mf_wf + q->b2lm;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->a2lm - q->PhDBconst + q->PhDBAterm;
    } else if ( q->f2lm < Mf_wf && Mf_wf <= q->fr ){
        Mf2lm = q->a2lm*q->f2lm + q->b2lm;
        UsefulPowers powers_of_Mf2lm;
        status = init_useful_powers(&powers_of_Mf2lm, Mf2lm);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf2lm failed");
        tmpphaseB = IMRPhenDPhase(Mf2lm, pPhi, pn, &powers_of_Mf2lm, phi_prefactors, Rholm, Taulm) / q->a2lm;
        Mf = q->a2lm*Mf_wf + q->b2lm;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->a2lm - q->PhDCconst + tmpphaseB - q->PhDBconst + q->PhDBAterm;
    } else if ( Mf_wf > q->fr ) {
        Mf2lm = q->a2lm*q->f2lm + q->b2lm;
        UsefulPowers powers_of_Mf2lm;
        status = init_useful_powers(&powers_of_Mf2lm, Mf2lm);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf2lm failed");
        tmpphaseB = IMRPhenDPhase(Mf2lm, pPhi, pn, &powers_of_Mf2lm, phi_prefactors, Rholm, Taulm) / q->a2lm;

        Mfr = q->a2lm*q->fr + q->b2lm;
        UsefulPowers powers_of_Mfr;
        status = init_useful_powers(&powers_of_Mfr, Mfr);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mfr failed");
        tmpphaseC = IMRPhenDPhase(Mfr, pPhi, pn, &powers_of_Mfr, phi_prefactors, Rholm, Taulm) / q->a2lm - q->PhDCconst + tmpphaseB - q->PhDBconst + q->PhDBAterm;;

        Mf = q->ar*Mf_wf + q->br;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->ar - q->PhDDconst + tmpphaseC;
    } else {
        XLALPrintError("XLAL_ERROR - should not get here - in function XLALSimIMRPhenomHMPhase");
        XLAL_ERROR(XLAL_EDOM);
    }

    /* phase shift due to leading order complex amplitude
    Luc Blancet 1310.1528. Section 9.5
    "Spherical hrmonic modes for numerical relativity" */
    /*initialise answer*/
    REAL8 cShift = 0.0;

    if ( ell==2 && mm==1 ) {
        cShift = LAL_PI/2.0; /* i shift */
    } else if ( ell==3 && mm==3 ) {
        cShift = -LAL_PI/2.0; /* -i shift */
    } else if ( ell==4 && mm==4 ) {
        cShift = LAL_PI; /* -1 shift */
    } else if ( ell==4 && mm==3 ) {
        cShift = -LAL_PI/2.0; /* -i shift */
    } else if ( ell==5 && mm==5 ) {
        cShift = LAL_PI/2.0; /* i shift */
    } else if ( ell==5 && mm==4 ) {
        cShift = LAL_PI; /* -1 shift */
    } else if ( ell==6 && mm==5 ) {
        cShift = LAL_PI/2.0; /* i shift */
    } /*else {
        XLALPrintError("XLAL Error (in function XLALSimIMRPhenomHMPhase) - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }*/
    /*TODO: Need to have a list at the begining and a function to check the input
    lm mode to see if it is one that is included in the model.*/

    retphase += cShift;

    // LALFree(pPhi);
    // LALFree(pn);

    return retphase;
}

/* Function to add modes for frequency-domain structures */
/* This function was lifted from the EOBNRv2HM_ROM code */
static INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym);

/********************* Function to add modes for frequency-domain structures ********************/

/* Helper function to add a mode to hplus, hcross in Fourier domain
 * - copies the function XLALSimAddMode, which was done only for TD structures */
 /* This function was lifted from the EOBNRv2HM_ROM code */
static INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym) {
  /* Deleted the definition of the string 'func': usage ? */
  COMPLEX16 Y;
  UINT4 j;
  COMPLEX16 hlmtildevalue;

  /* Checks LAL_CHECK_VALID_SERIES and LAL_CHECK_CONSISTENT_TIME_SERIES removed
   * - they do not seem available for frequency series ? */

  INT4 minus1l; /* (-1)^l */
  if ( l%2 ) minus1l = -1;
  else minus1l = 1;
  if ( sym ) { /* equatorial symmetry: add in -m mode */
    Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
    COMPLEX16 Ymstar = conj(XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m));
    COMPLEX16 factorp = 1./2*(Y + minus1l*Ymstar);
    COMPLEX16 factorc = I/2*(Y - minus1l*Ymstar);
    // COMPLEX16* datap = hptilde->data->data;
    // COMPLEX16* datac = hctilde->data->data;
    for ( j = 0; j < hlmtilde->data->length; ++j ) {
      hlmtildevalue = (hlmtilde->data->data[j]);
    //   datap[j] += factorp*hlmtildevalue;
    //   datac[j] += factorc*hlmtildevalue;
      hptilde->data->data[j] += factorp*hlmtildevalue;
      hctilde->data->data[j] += factorc*hlmtildevalue;
    }
  }
  else { /* not adding in the -m mode */
    Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
    COMPLEX16 factorp = 1./2*Y;
    COMPLEX16 factorc = I/2*Y;
    // COMPLEX16* datap = hptilde->data->data;
    // COMPLEX16* datac = hctilde->data->data;
    for ( j = 0; j < hlmtilde->data->length; ++j ) {
      hlmtildevalue = (hlmtilde->data->data[j]);
    //   datap[j] += factorp*hlmtildevalue;
    //   datac[j] += factorc*hlmtildevalue;
      hptilde->data->data[j] += factorp*hlmtildevalue;
      hctilde->data->data[j] += factorc*hlmtildevalue;
    }
  }

  return 0;
}

/*
 *
 * New functions in for the 'clean' code
 *
 *
 *
 */

static COMPLEX16 IMRPhenomHMSingleModehlmOpt2(REAL8 eta, REAL8 chi1z, REAL8 chi2z, int ell, int mm, double Mf, HMPhasePreComp *z, IMRPhenomDAmplitudeCoefficients *pAmp, AmpInsPrefactors * amp_prefactors, PNPhasingSeries *pn, IMRPhenomDPhaseCoefficients *pPhi, PhiInsPrefactors * phi_prefactors, double Rholm, double Taulm, double phi_precalc );
/*
 * This returns the hlm of a single (l,m) mode at a single frequency Mf
 */
static COMPLEX16 IMRPhenomHMSingleModehlmOpt2(
        REAL8 eta, REAL8 chi1z, REAL8 chi2z,
        int ell,
        int mm,
        double Mf,
        HMPhasePreComp *z,
        IMRPhenomDAmplitudeCoefficients *pAmp,
        AmpInsPrefactors * amp_prefactors,
        PNPhasingSeries *pn,
        IMRPhenomDPhaseCoefficients *pPhi,
        PhiInsPrefactors * phi_prefactors,
        double Rholm,
        double Taulm,
        double phi_precalc /**< m*phi0 + HMphaseRef*/
) {

    /*
    * In this function we should pass the
    * phenom model parameters and the (l,m)
    * mode to return a given hlm, not summed
    * with spherical harmonics.
    * Can be evaluated at a single geometric frequency (Mf).
    */

    double HMamp = XLALSimIMRPhenomHMAmplitudeOpt( Mf, eta, chi1z, chi2z, ell, mm, pAmp, amp_prefactors );
    // printf("Mf = %.10f ,    HMamp Opt = %.16f\n", Mf, HMamp);
    double HMphase = XLALSimIMRPhenomHMPhaseOpt( Mf, ell, mm, z, pn, pPhi, phi_prefactors, Rholm, Taulm );

    /* compute reference phase at reference frequency */

    // factor of m spherical harmonic mode b/c phi0 is orbital phase
    /* NOTE: Does HMphaseRef already have the mm scaling? as it's m*(phi0 + phiref) */
    HMphase -= phi_precalc;

    COMPLEX16 hlm = HMamp * cexp(-I * HMphase);
    // printf("Mf = %f     hlm  =  %f + i %f\nx",Mf, creal(hlm), cimag(hlm));

    return hlm;
}

/* given the final frequency in Mf and the total mass, calculate the final frequency in Hz */
static REAL8 ComputeIMRPhenomHMfmax(REAL8 Mf, REAL8 f_min, REAL8 f_max, REAL8 M);
static REAL8 ComputeIMRPhenomHMfmax(REAL8 Mf /**< geometric frequency*/,
                                     REAL8 f_min /**< low frequency in Hz*/,
                                     REAL8 f_max /**< end frequency in Hz*/,
                                     REAL8 M /**< total mass (Msun)*/){

      const REAL8 M_sec = M * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency

      const REAL8 fCut = Mf/M_sec; // convert Mf -> Hz
      // Somewhat arbitrary end point for the waveform.
      // Chosen so that the end of the waveform is well after the ringdown.
      if (fCut <= f_min)
          XLAL_ERROR(XLAL_EDOM, "(fCut = %g Hz) <= f_min = %g\n", fCut, f_min);

      /* default f_max to Cut */
      REAL8 f_max_prime = f_max;
      f_max_prime = f_max ? f_max : fCut;
      f_max_prime = (f_max_prime > fCut) ? fCut : f_max_prime;
      if (f_max_prime <= f_min)
          XLAL_ERROR(XLAL_EDOM, "f_max <= f_min\n");

      return f_max_prime;
}

// static REAL8 Computet0(REAL8 eta, REAL8 chi1z, REAL8 chi2z, REAL8 finspin, INT4 ell, INT4 mm);
static REAL8 Computet0(REAL8 eta, REAL8 chi1z, REAL8 chi2z, REAL8 finspin);
static REAL8 Computet0(REAL8 eta, REAL8 chi1z, REAL8 chi2z, REAL8 finspin){

    /*This should be OK to pass NULL for LALDict here.*/
    LALDict *extraParams = NULL;
    IMRPhenomDPhaseCoefficients *pPhi = ComputeIMRPhenomDPhaseCoefficients(eta, chi1z, chi2z, finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);

    IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1z, chi2z, finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);

    // double Rholm = XLALSimIMRPhenomHMRholm(eta, chi1z, chi2z, ell, mm);
    // double Taulm = XLALSimIMRPhenomHMTaulm(eta, chi1z, chi2z, ell, mm);

    //time shift so that peak amplitude is approximately at t=0
    //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
    //NOTE: All modes will have the same time offset. So we use the 22 mode.
    //If we just use the 22 mode then we pass 1.0, 1.0 into DPhiMRD.
    const REAL8 t0 = DPhiMRD(pAmp->fmaxCalc, pPhi, 1.0, 1.0);
    return t0;
}



int EnforcePrimaryIsm1(REAL8 *m1, REAL8 *m2, REAL8 *chi1z, REAL8 *chi2z);
int EnforcePrimaryIsm1(REAL8 *m1, REAL8 *m2, REAL8 *chi1z, REAL8 *chi2z){
    REAL8 chi1z_tmp, chi2z_tmp, m1_tmp, m2_tmp;
    if (*m1>*m2) {
       chi1z_tmp = *chi1z;
       chi2z_tmp = *chi2z;
       m1_tmp   = *m1;
       m2_tmp   = *m2;
   } else { /* swap spins and masses */
       chi1z_tmp = *chi2z;
       chi2z_tmp = *chi1z;
       m1_tmp   = *m2;
       m2_tmp   = *m1;
    }
    *m1 = m1_tmp;
    *m2 = m2_tmp;
    *chi1z = chi1z_tmp;
    *chi2z = chi2z_tmp;

    if (*m1 < *m2)
        XLAL_ERROR(XLAL_EDOM, "XLAL_ERROR in EnforcePrimaryIsm1. When trying\
 to enfore that m1 should be the larger mass.\
 After trying to enforce this m1 = %f and m2 = %f\n", *m1, *m2);

    return XLAL_SUCCESS;
}


/* the goal of this function is to compute hlm for a list of modes
 * and output them in the data type
 * SphHarmFrequencySeries
 */

/* NOTE: I've change the name from hlmsphharmfreqseries to hlms */

/* NOTE: Some code duplication here because I'd like to keep this function XLAL */
int XLALIMRPhenomHMMultiModehlmOpt(
    SphHarmFrequencySeries **hlms, /**< [out] can access multiple modes with units of Mf */
    REAL8 m1Msun,                        /**< primary mass in Msun*/
    REAL8 m2Msun,                        /**< secondary mass in Msun*/
    REAL8 chi1z,
    REAL8 chi2z,
    REAL8 deltaF,/**< Hz */
    REAL8 f_min,/**< Hz */
    REAL8 f_max, /**< Hz */
    REAL8 fRef_in,  /**< reference frequency in Hz */
    REAL8 phi0, /**< reference orbital phase */
    REAL8 distance /**< distance to source in SI */
) {

    LALDict *extraParams = NULL;

    //here masses are in Msun
    int ret = EnforcePrimaryIsm1(&m1Msun, &m2Msun, &chi1z, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "EnforcePrimaryIsm1 failed");

    /*Unfortunately duplication of M, eta and M_sec here and in XLALIMRPhenomHMMultiModeStrain*/
    const REAL8 M = m1Msun + m2Msun;
    const REAL8 eta = m1Msun * m2Msun / (M * M);
    const REAL8 M_sec = M * LAL_MTSUN_SI;

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    if (chi1z > 1.0 || chi1z < -1.0 || chi2z > 1.0 || chi2z < -1.0)
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    // if no reference frequency given, set it to the starting GW frequency
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    REAL8 f_max_prime = ComputeIMRPhenomHMfmax(Mf_CUT_HM, f_min, f_max, M);

    /* This does not work :( */
    // #define NMODES_MAX LMAX
    // int NMODES = NMODES_MAX;
    // int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1} };
    //
    //
    // if ( LMAX==2 ) {
    //     int NMODES = NMODES_MAX;
    //     int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1} };
    // } else if ( LMAX==3 ) {
    //     int NMODES = NMODES_MAX;
    //     int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {3,2} };
    // } else {
    //     XLAL_ERROR(XLAL_EDOM, "LMAX = %d not supported.\n", LMAX);
    // }

    // Convention m1 >= m2
    const REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1Msun, m2Msun, chi1z, chi2z);
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    /*compute phenomD amp and phase coefficients*/

    IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1z, chi2z, finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    AmpInsPrefactors amp_prefactors;
    int errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");
    if (extraParams==NULL)
    extraParams=XLALCreateDict();
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    IMRPhenomDPhaseCoefficients *pPhi = ComputeIMRPhenomDPhaseCoefficients(eta, chi1z, chi2z, finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1Msun, m2Msun, chi1z, chi2z, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    REAL8 testGRcor=1.0;
    testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);

    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSS(m1Msun, m2Msun, M, chi1z, chi2z) * pn->v[0])* testGRcor;

    PhiInsPrefactors phi_prefactors;
    int status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

    printf("HI IM FROM WORKING FUNC phi_prefactors.two_thirds = %f\n", phi_prefactors.two_thirds);
    /* Compute the amplitude pre-factor */
    const REAL8 amp0 = M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    // shift by overall length in time
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, -1. / deltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    size_t n = NextPow2(f_max_prime / deltaF) + 1;
    /* range that will have actual non-zero waveform values generated */
    size_t ind_min = (size_t) (f_min / deltaF);
    size_t ind_max = (size_t) (f_max_prime / deltaF);
    XLAL_CHECK ( (ind_max<=n) && (ind_min<=ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=hptilde->data>length=%zu.", ind_min, ind_max, n);

    const REAL8 MfRef = M_sec * fRef;





    /*
    * In this function we should setup
    * all the variables that we can precompute
    * to generate the PhenomD amplitude and phase functions
    * We then loop over the static function 'IMRPhenomHMSingleModehlm'
    * to generate the modes.
    * We then sum them up at the end.
    */

    /* Now we have all the PhenomD model parameters, which actually correspond
     * to the (l,m)=(2,2) parameters, we can feed these into the function
     * IMRPhenomHMSingleModehlm to generate any mode we want. */

     //TODO: Turn t0 computation into a function
     /*NOTE: We could compute the t0 here as the time of the peak of the
     22 mode and make that the assumtion.*/


         /*
         * NOTE: I'm not sure what Mf should be used for the reference time... I think it should be the scaled one. And it should come from the amplitude
         * NOTE: I'm not sure what Mf should be used for the reference time... I think it should be the scaled one. And it should come from the amplitude
         */

     const REAL8 t0 = Computet0(eta, chi1z, chi2z, finspin);




     for( int j=0; j<NMODES; j++ ){

         int ell = ModeArray[j][0];
         int mm = ModeArray[j][1];
         printf("ell = %i\n",ell);
         printf("mm = %i\n",mm);

         double Rholm = XLALSimIMRPhenomHMRholm(eta, chi1z, chi2z, ell, mm);
         double Taulm = XLALSimIMRPhenomHMTaulm(eta, chi1z, chi2z, ell, mm);

         // Compute coefficients to make phase C^1 continuous (phase and first derivative)
         ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);


         /* compute phenomHM pre computations */
         /* NOTE: Need to make this an input and NOT part of the frequency loop! */
         HMPhasePreComp z;
         ret = XLALSimIMRPhenomHMPhasePreComp(&z, ell, mm, eta, chi1z, chi2z, finspin);
         if (ret != XLAL_SUCCESS){
             XLALPrintError("XLAL Error - XLALSimIMRPhenomHMPhasePreComp failed\n");
             XLAL_ERROR(XLAL_EDOM);
         }

         double HMphaseRef = XLALSimIMRPhenomHMPhaseOpt( MfRef, ell, mm, &z, pn, pPhi, &phi_prefactors, Rholm, Taulm );

         /* compute reference phase at reference frequency */

         // factor of m spherical harmonic mode b/c phi0 is orbital phase
         /* NOTE: Does HMphaseRef already have the mm scaling? as it's m*(phi0 + phiref) */
         REAL8 phi_precalc = mm*phi0 + HMphaseRef;


         /* we loop over (l,m) and use a temporary hlm frequency series to store the results of a single mode */
         COMPLEX16FrequencySeries *hlm = NULL;
         hlm = XLALCreateCOMPLEX16FrequencySeries("hlm: single mode", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
         memset(hlm->data->data, 0, n * sizeof(COMPLEX16));

        /* NOTE: Do I need this bit? */
        /* XLALUnitMultiply(hlm->sampleUnits, hlm->sampleUnits, &lalSecondUnit); */

         /* loop over frequency */
         /* Now generate the waveform for a single (l,m) mode */
         #pragma omp parallel for
         for (size_t i = ind_min; i < ind_max; i++)
         {
            REAL8 Mf = M_sec * i * deltaF; /* geometric frequency */
            /* now we can compute the hlm */
            /* TODO: fix phase and time offsets */
            /* TODO: inclusion of reference frequency */
            // phi -= t0*(Mf-MfRef) + phi_precalc;
            // ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
            /* construct hlm at single frequency point and return */
            // (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(eta, chi1z, chi2z, ell, mm, Mf, MfRef, phi0, &z);

            //optimised version
            (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlmOpt2(eta, chi1z, chi2z, ell, mm, Mf, &z, pAmp, &amp_prefactors, pn , pPhi, &phi_prefactors, Rholm, Taulm, phi_precalc );
            // (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlmOpt(eta, chi1z, chi2z, ell, mm, Mf, MfRef, phi0, &z, pAmp);

            // printf("f = %f    Mf = %f      (hlm->data->data)[i]  =  %f + i %f\nx",i * deltaF, Mf, creal((hlm->data->data)[i]), cimag((hlm->data->data)[i]));
            (hlm->data->data)[i] *= cexp(-I * t0*(Mf-MfRef) );
            /*from phenomD for referene*/
            //  REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
            //  REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors);
            //  phi -= t0*(Mf-MfRef) + phi_precalc;
            //  ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
         }
        *hlms = XLALSphHarmFrequencySeriesAddMode(*hlms, hlm, ell, mm);
         /* Destroy hlm in (l,m) loop */
         XLALDestroyCOMPLEX16FrequencySeries(hlm);
     }

     LALFree(pAmp);
     LALFree(pPhi);
     LALFree(pn);

    return XLAL_SUCCESS;
}

/* This function will be a wrapper of
 * XLALIMRPhenomHMMultiModehlm
 * It will take the
 * SphHarmFrequencySeries **hlms as inputs
 * and compbine them with the Spherical harmonics
 * at a given inclination angle and azimuthal angle
 * and return the hptilde and hctilde
 */
 /*TODO: Add LALDict as an argument to XLALIMRPhenomHMMultiModeStrain */
 /*TODO: Change name from XLALIMRPhenomHMMultiModeStrain to XLALIMRPhenomHM */
int XLALIMRPhenomHMMultiModeStrain(
    COMPLEX16FrequencySeries **hptilde, /**< [out] */
    COMPLEX16FrequencySeries **hctilde, /**< [out] */
    REAL8 m1,                        /**< primary mass in SI*/
    REAL8 m2,                        /**< secondary mass in SI*/
    REAL8 chi1z,
    REAL8 chi2z,
    REAL8 deltaF, /**< Hz*/
    REAL8 f_min, /**< Hz*/
    REAL8 f_max, /**< Hz*/
    REAL8 fRef_in, /**< reference frequency in Hz */
    REAL8 phi0,/**< reference orbital phase */
    REAL8 inclination,
    REAL8 distance /**< distance to source in SI */
) {

    /* sanity checks on input parameters */

    if (m1 <= 0) XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
    if (m2 <= 0) XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");

    if (deltaF <= 0) XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
    if (f_min <= 0) XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
    if (f_max < 0) XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");

   /* external: SI; internal: solar masses */
   m1 /= LAL_MSUN_SI;
   m2 /= LAL_MSUN_SI;

   int ret = EnforcePrimaryIsm1(&m1, &m2, &chi1z, &chi2z);
   XLAL_CHECK(XLAL_SUCCESS == ret, ret, "EnforcePrimaryIsm1 failed");

    const REAL8 q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

    if (q > MAX_ALLOWED_MASS_RATIO)
      XLAL_PRINT_WARNING("Warning: The model is not supported for high mass ratio, see MAX_ALLOWED_MASS_RATIO\n");

    const REAL8 M = m1 + m2; /* total mass (Msun) */
    const REAL8 eta = m1 * m2 / (M * M);

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    if (chi1z > 1.0 || chi1z < -1.0 || chi2z > 1.0 || chi2z < -1.0)
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    // if no reference frequency given, set it to the starting GW frequency
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    REAL8 f_max_prime = ComputeIMRPhenomHMfmax(Mf_CUT_HM, f_min, f_max, M);

    /* evaluate XLALIMRPhenomHMMultiModehlm */

    // SphHarmFrequencySeries *hlms=NULL;
    SphHarmFrequencySeries **hlms=XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlms=NULL;

    /*TODO: Add LALDict as an argument to XLALIMRPhenomHMMultiModehlm */
    ret = XLALIMRPhenomHMMultiModehlmOpt(hlms, m1, m2, chi1z, chi2z, deltaF, f_min, f_max_prime, fRef, phi0, distance);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALIMRPhenomHMMultiModehlm(&hlms) failed");

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    // shift by overall length in time
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, -1. / deltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    /* compute array sizes */
    size_t n = NextPow2(f_max_prime / deltaF) + 1;
    /* range that will have actual non-zero waveform values generated */
    size_t ind_min = (size_t) (f_min / deltaF);
    size_t ind_max = (size_t) (f_max_prime / deltaF);
    XLAL_CHECK ( (ind_max<=n) && (ind_min<=ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=hptilde->data>length=%zu.", ind_min, ind_max, n);

    /* Allocate hptilde and hctilde */
    *hptilde = XLALCreateCOMPLEX16FrequencySeries("hptilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
    if (!(hptilde) ) XLAL_ERROR(XLAL_EFUNC);
    memset((*hptilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&(*hptilde)->sampleUnits, &(*hptilde)->sampleUnits, &lalSecondUnit);

    *hctilde = XLALCreateCOMPLEX16FrequencySeries("hctilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
    if (!(hctilde) ) XLAL_ERROR(XLAL_EFUNC);
    memset((*hctilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&(*hctilde)->sampleUnits, &(*hctilde)->sampleUnits, &lalSecondUnit);


    /* Adding the modes to form hplus, hcross
     * - use of a function that copies XLALSimAddMode but for Fourier domain structures */
    INT4 sym; /* sym will decide whether to add the -m mode (when equatorial symmetry is present) */
    for( int i=0; i<NMODES; i++){
      INT4 ell = ModeArray[i][0];
      INT4 mm = ModeArray[i][1];
    //   printf("computing hlm for mode l=%d, m=%d\n", ell, mm); /*NOTE: Remove this*/
      COMPLEX16FrequencySeries* hlm = XLALSphHarmFrequencySeriesGetMode(*hlms, ell, mm);
      if (!(hlm)) XLAL_ERROR(XLAL_EFUNC);

      if ( mm==0 ) {
          sym = 0; /* We test for hypothetical m=0 modes */
      } else {
          sym = 1;
      }
      FDAddMode( *hptilde, *hctilde, hlm, inclination, 0., ell, mm, sym); /* The phase \Phi is set to 0 - assumes phiRef is defined as half the phase of the 22 mode h22 (or the first mode in the list), not for h = hplus-I hcross */
    }

    XLALDestroySphHarmFrequencySeries(*hlms);
    XLALFree(hlms);

    return XLAL_SUCCESS;
}


/* convenience function to compute hlm for a single mode - from
which the Alm and philm (amplitude and phase) of a particular
mode can be obtained. Useful to compare to NR */
int XLALSimIMRPhenomHMSingleModehlm(COMPLEX16FrequencySeries **hlmtilde, /**< [out] */
                                    REAL8 m1Msun,
                                    REAL8 m2Msun,
                                    REAL8 chi1z,
                                    REAL8 chi2z,
                                    REAL8 deltaF,
                                    REAL8 f_min,
                                    REAL8 f_max,
                                    REAL8 fRef_in,
                                    REAL8 phi0,
                                    REAL8 distance,
                                    INT4 ell,
                                    INT4 mm){
    LALDict *extraParams = NULL;

    //here masses are in Msun
    int ret = EnforcePrimaryIsm1(&m1Msun, &m2Msun, &chi1z, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "EnforcePrimaryIsm1 failed");

printf(" m1Msun= %f\n",m1Msun);
printf(" m2Msun= %f\n",m2Msun);
printf(" chi1z= %f\n",chi1z);
printf(" chi2z= %f\n",chi2z);

    /*Unfortunately duplication of M, eta and M_sec here and in XLALIMRPhenomHMMultiModeStrain*/
    const REAL8 M = m1Msun + m2Msun;
    const REAL8 eta = m1Msun * m2Msun / (M * M);
    const REAL8 M_sec = M * LAL_MTSUN_SI;

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    if (chi1z > 1.0 || chi1z < -1.0 || chi2z > 1.0 || chi2z < -1.0)
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    // if no reference frequency given, set it to the starting GW frequency
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    REAL8 f_max_prime = ComputeIMRPhenomHMfmax(Mf_CUT_HM, f_min, f_max, M);


    /* This does not work :( */
    // #define NMODES_MAX LMAX
    // int NMODES = NMODES_MAX;
    // int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1} };
    //
    //
    // if ( LMAX==2 ) {
    //     int NMODES = NMODES_MAX;
    //     int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1} };
    // } else if ( LMAX==3 ) {
    //     int NMODES = NMODES_MAX;
    //     int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {3,2} };
    // } else {
    //     XLAL_ERROR(XLAL_EDOM, "LMAX = %d not supported.\n", LMAX);
    // }

    // Convention m1 >= m2
    const REAL8 finspin = XLALSimIMRPhenomDFinalSpin(m1Msun, m2Msun, chi1z, chi2z);
    if (finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

    /*compute phenomD amp and phase coefficients*/

    IMRPhenomDAmplitudeCoefficients *pAmp = ComputeIMRPhenomDAmplitudeCoefficients(eta, chi1z, chi2z, finspin);

    printf("pAmp->chi1 = %f\n", pAmp->chi1);

    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    AmpInsPrefactors amp_prefactors;
    int errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");
    if (extraParams==NULL)
    extraParams=XLALCreateDict();
    XLALSimInspiralWaveformParamsInsertPNSpinOrder(extraParams,LAL_SIM_INSPIRAL_SPIN_ORDER_35PN);
    IMRPhenomDPhaseCoefficients *pPhi = ComputeIMRPhenomDPhaseCoefficients(eta, chi1z, chi2z, finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1Msun, m2Msun, chi1z, chi2z, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);


    printf("pPhi->beta1 = %f\n", pPhi->beta1);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    REAL8 testGRcor=1.0;
    testGRcor += XLALSimInspiralWaveformParamsLookupNonGRDChi6(extraParams);

    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSS(m1Msun, m2Msun, M, chi1z, chi2z) * pn->v[0])* testGRcor;

    PhiInsPrefactors phi_prefactors;
    int status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");
    /* NOTE: There seems to be a problem here, with the phi_prefactors as in IMRPhenomHMMultiModehlmOpt they work fine here but in this
    function they do not. */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */

    printf("phi_prefactors.two_thirds = %f\n", phi_prefactors.two_thirds);

    /* Compute the amplitude pre-factor */
    const REAL8 amp0 = M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    // shift by overall length in time
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, -1. / deltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    /* compute array sizes */
    size_t n = NextPow2(f_max_prime / deltaF) + 1;
    /* range that will have actual non-zero waveform values generated */
    size_t ind_min = (size_t) (f_min / deltaF);
    size_t ind_max = (size_t) (f_max_prime / deltaF);
    XLAL_CHECK ( (ind_max<=n) && (ind_min<=ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=hptilde->data>length=%zu.", ind_min, ind_max, n);


    /* Allocate hptilde and hctilde */
    *hlmtilde = XLALCreateCOMPLEX16FrequencySeries("hlmtilde: FD waveform", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
    if (!(hlmtilde) ) XLAL_ERROR(XLAL_EFUNC);
    memset((*hlmtilde)->data->data, 0, n * sizeof(COMPLEX16));
    XLALUnitDivide(&(*hlmtilde)->sampleUnits, &(*hlmtilde)->sampleUnits, &lalSecondUnit);

    const REAL8 MfRef = M_sec * fRef;

    /*
    * In this function we should setup
    * all the variables that we can precompute
    * to generate the PhenomD amplitude and phase functions
    * We then loop over the static function 'IMRPhenomHMSingleModehlm'
    * to generate the modes.
    * We then sum them up at the end.
    */

    /* Now we have all the PhenomD model parameters, which actually correspond
     * to the (l,m)=(2,2) parameters, we can feed these into the function
     * IMRPhenomHMSingleModehlm to generate any mode we want. */

     //TODO: Turn t0 computation into a function
     /*NOTE: We could compute the t0 here as the time of the peak of the
     22 mode and make that the assumtion.*/


         /*
         * NOTE: I'm not sure what Mf should be used for the reference time... I think it should be the scaled one. And it should come from the amplitude
         * NOTE: I'm not sure what Mf should be used for the reference time... I think it should be the scaled one. And it should come from the amplitude
         */

     const REAL8 t0 = Computet0(eta, chi1z, chi2z, finspin);

     printf("ell = %i\n",ell);
     printf("mm = %i\n",mm);

     double Rholm = XLALSimIMRPhenomHMRholm(eta, chi1z, chi2z, ell, mm);
     double Taulm = XLALSimIMRPhenomHMTaulm(eta, chi1z, chi2z, ell, mm);

     printf("Rholm = %f\n",Rholm);
     printf("Taulm = %f\n",Taulm);

     // Compute coefficients to make phase C^1 continuous (phase and first derivative)
     ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);
     printf("pPhi->C1Int = %f\n", pPhi->C1Int);
     /* compute phenomHM pre computations */
     /* NOTE: Need to make this an input and NOT part of the frequency loop! */
     HMPhasePreComp z;
     ret = XLALSimIMRPhenomHMPhasePreComp(&z, ell, mm, eta, chi1z, chi2z, finspin);
     if (ret != XLAL_SUCCESS){
         XLALPrintError("XLAL Error - XLALSimIMRPhenomHMPhasePreComp failed\n");
         XLAL_ERROR(XLAL_EDOM);
     }

     double HMphaseRef = XLALSimIMRPhenomHMPhaseOpt( MfRef, ell, mm, &z, pn, pPhi, &phi_prefactors, Rholm, Taulm );
     printf("HMphaseRef = %f\n",HMphaseRef);
     /* compute reference phase at reference frequency */

     // factor of m spherical harmonic mode b/c phi0 is orbital phase
     /* NOTE: Does HMphaseRef already have the mm scaling? as it's m*(phi0 + phiref) */
     REAL8 phi_precalc = mm*phi0 + HMphaseRef;

    /* NOTE: Do I need this bit? */
    /* XLALUnitMultiply(hlm->sampleUnits, hlm->sampleUnits, &lalSecondUnit); */

     /* loop over frequency */
     /* Now generate the waveform for a single (l,m) mode */
     #pragma omp parallel for
     for (size_t i = ind_min; i < ind_max; i++)
     {
        REAL8 Mf = M_sec * i * deltaF; /* geometric frequency */
        /* now we can compute the hlm */
        /* TODO: fix phase and time offsets */
        /* TODO: inclusion of reference frequency */
        // phi -= t0*(Mf-MfRef) + phi_precalc;
        // ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
        /* construct hlm at single frequency point and return */
        // (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(eta, chi1z, chi2z, ell, mm, Mf, MfRef, phi0, &z);

        //optimised version
        ((*hlmtilde)->data->data)[i] = amp0 * IMRPhenomHMSingleModehlmOpt2(eta, chi1z, chi2z, ell, mm, Mf, &z, pAmp, &amp_prefactors, pn , pPhi, &phi_prefactors, Rholm, Taulm, phi_precalc );
        // (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlmOpt(eta, chi1z, chi2z, ell, mm, Mf, MfRef, phi0, &z, pAmp);

        // printf("f = %f    Mf = %f      (hlm->data->data)[i]  =  %f + i %f\nx",i * deltaF, Mf, creal((hlm->data->data)[i]), cimag((hlm->data->data)[i]));
        ((*hlmtilde)->data->data)[i] *= cexp(-I * t0*(Mf-MfRef) );
        /*from phenomD for referene*/
        //  REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
        //  REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors);
        //  phi -= t0*(Mf-MfRef) + phi_precalc;
        //  ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
     }

     LALFree(pAmp);
     LALFree(pPhi);
     LALFree(pn);

    return XLAL_SUCCESS;
}
