/*
 * Copyright (C) 2017 Sebastian Khan, Francesco Pannarale
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

//FP: make PhenomDQuantities, powers_of_MfAtScale_22_amp, etc. EXTERNAL?

/* List of modelled modes */
static int NMODES = NMODES_MAX;
static const int ModeArray[NMODES_MAX][2] = { {2,2}, {2,1}, {3,3}, {4,4}, {5,5} };

/* Dimensionless frequency of last data point in waveform */
#define Mf_CUT_HM 0.5
/* Activates amplitude part of the model */
#define AmpFlagTrue 1
#define MfAtScale_wf_amp 0.0001

/* START: newer functions  */

int init_useful_mf_powers(UsefulMfPowers *p, REAL8 number)
{
	XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
	XLAL_CHECK(!(number < 0) , XLAL_EDOM, "number must be non-negative");

	// consider changing pow(x,1/6.0) to cbrt(x) and sqrt(x) - might be faster
	p->itself = number;
        p->sqrt = sqrt(number);
	//p->sixth = pow(number, 1/6.0); //FP
	p->sixth = cbrt(p->sqrt);
        p->m_sixth = 1.0/p->sixth;
	p->third = p->sixth * p->sixth;
	p->two_thirds = p->third * p->third;
	p->four_thirds = number * p->third;
	p->five_thirds = number * p->two_thirds;
	p->two = number * number;
	p->seven_thirds = p->third * p->two;
	p->eight_thirds = p->two_thirds * p->two;
        p->m_seven_sixths = p->m_sixth/number;
        p->m_five_sixths = p->m_seven_sixths*p->third;
        p->m_sqrt = 1./p->sqrt;

	return XLAL_SUCCESS;
}

/* TODO: at this point it is probably unnecessary to have UsefulPowers and
 * UsefulMfPowers: one struct should suffice */

/* Given a full UsefulMfPowers variable, extract a UsefulPowers one from it */
int downsize_useful_mf_powers(UsefulPowers *out, UsefulMfPowers *in);
int downsize_useful_mf_powers(UsefulPowers *out, UsefulMfPowers *in)
{
	XLAL_CHECK(0 != in, XLAL_EFAULT, "in is NULL");

	out->third          = in->third;
	out->two_thirds     = in->two_thirds;
	out->four_thirds    = in->four_thirds;
	out->five_thirds    = in->five_thirds;
	out->two            = in->two;
	out->seven_thirds   = in->seven_thirds;
	out->eight_thirds   = in->eight_thirds;
        out->inv            = 1./in->itself; 
        out->m_seven_sixths = in->m_seven_sixths;
        out->m_third        = 1./in->third;
        out->m_two_thirds   = out->m_third*out->m_third;
        out->m_five_thirds  = out->inv*out->m_two_thirds;

	return XLAL_SUCCESS;
}

/* Precompute a bunch of PhenomD related quantities and store them filling in a PhenomDStorage variable */
int init_PhenomD_Storage(PhenomDStorage* p, const REAL8 m1, const REAL8 m2, const REAL8 chi1z, const REAL8 chi2z)
{
  XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");

  p->m1 = m1; /* Inherits units from m1 in function arguments */
  p->m2 = m2; /* Inherits units from m2 in function arguments */
  p->Mtot = m1+m2; /* Inherits units from m1 and m2 in function arguments */
  p->eta = m1*m2/(p->Mtot*p->Mtot);
  p->Inv1MinusEradRational0815 = 1.0/(1.0-EradRational0815(p->eta, chi1z, chi2z));
  p->finspin = XLALSimIMRPhenomDFinalSpin(m1, m2, chi1z, chi2z); /* dimensionless final spin */
  if (p->finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

  /* Precompute Lionel's QNM higher mode fits (both real and imaginary parts of
     the complex ringdown frequency).
     Annoying, but faster than looking up ell's and mm's.
     WARNING: bug prone, as one may access PhenomHMfring and PhenomHMfdamp
     for ell's and mm's not in ModeArray. */
  complex double ZZ;
  const REAL8 inv2Pi = 0.5/LAL_PI;
  for( int j=0; j<NMODES; j++ ){
    int ell = ModeArray[j][0];
    int mm = ModeArray[j][1];
    ZZ = CW07102016( KAPPA(p->finspin, ell, mm), ell, mm, 0 );
    /* lm mode ringdown frequency (real part of ringdown), geometric units */
    const REAL8 Mf_RD = inv2Pi * creal(ZZ); /* GW ringdown frequency, converted from angular frequency */
    p->PhenomHMfring[ell][mm] = Mf_RD * p->Inv1MinusEradRational0815; /* scale by predicted final mass */
    /* lm mode ringdown damping time (imaginary part of ringdown), geometric units */
    const REAL8 fdamp = inv2Pi * cimag(ZZ); /* this is the 1./tau in the complex QNM */
    p->PhenomHMfdamp[ell][mm] = fdamp * p->Inv1MinusEradRational0815; /* scale by predicted final mass */
  }
  p->Mf_RD_22 = p->PhenomHMfring[2][2];
  p->Mf_DM_22 = p->PhenomHMfdamp[2][2];

  for( int j=0; j<NMODES; j++ ){
    int ell = ModeArray[j][0];
    int mm = ModeArray[j][1];
    p->Rholm[ell][mm] = p->Mf_RD_22/p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm]/p->Mf_DM_22;

  }

  /* A bunch of useful powers used in XLALSimIMRPhenomHMPNAmplitudeLeadingOrder */
  REAL8 sqrteta = sqrt(p->eta);
  REAL8 Seta = sqrt( 1.0 - 4.0 * p->eta );
  REAL8 ans = sqrteta;
  for( int j=0; j<NMODES; j++ ){
    int ell = ModeArray[j][0];
    int mm = ModeArray[j][1];

    if ( ell==2 ) {
        if (mm==2 ) {
          ans *= 0.674677;
        } else { // mm==1
          ans *= 0.329376 * Seta;
        }
    } else if ( ell==3 ) {
        if ( mm==3 ) {
          ans *= 0.767106 * Seta;
        }
        else { // mm==2
          ans *= 0.407703;
        }
    } else if ( ell==4 ) {
        if ( mm==4 ) {
          ans *= ( 1.08721 - 3.26162*p->eta );
        } else { // mm==3
          ans *=  ( 0.570006 - 1.14001*p->eta ) * Seta;
        }
    } else if ( ell==5 ) {
        if ( mm==5 ) {
          ans *= 3.39713 * (0.5 - p->eta) * Seta;
        }
        else { // mm==4
          ans *= 0.859267;
        }
    } else if ( ell==6 ) {
        if ( mm==6 ) {
          ans *= 2.80361;
        }
        else { // mm==5
          ans *= 1.36104 * Seta;
        }
    } else {
        ans = 0.0;
    }
    p->pow_Mf_wf_prefactor[ell][mm] = ans;
  }

  return XLAL_SUCCESS;
}

/**
 * Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
 * (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but was not
 * available when PhenomD was tuned.  This is similar to Subtract3PNSS in
 * LALSimIMRPhenomD_internals.c, but avoids recomputing eta and squaring the chi's.
 */
static double Subtract3PNSpinSpin(double m1, double m2, double M, double eta, double chi1, double chi2);
static double Subtract3PNSpinSpin(double m1, double m2, double M, double eta, double chi1, double chi2){
  REAL8 m1M = m1 / M;
  REAL8 m2M = m2 / M;
  REAL8 pn_ss3 =  (326.75L/1.12L + 557.5L/1.8L*eta)*eta*chi1*chi2;
  pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m1M-120.L*m1M*m1M) + (-4108.25L/6.72L-108.5L/1.2L*m1M+125.5L/3.6L*m1M*m1M)) *m1M*m1M * chi1*chi1;
  pn_ss3 += ((4703.5L/8.4L+2935.L/6.L*m2M-120.L*m2M*m2M) + (-4108.25L/6.72L-108.5L/1.2L*m2M+125.5L/3.6L*m2M*m2M)) *m2M*m2M * chi2*chi2;
  return pn_ss3;
}

//mathematica function postfRDflm
//double XLALIMRPhenomHMpostfRDflm(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities){
double XLALIMRPhenomHMTrd(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities);
double XLALIMRPhenomHMTrd(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities){
    double ans = 0.0;
    if ( AmpFlag==1 ) {
        /* For amplitude */
        ans = Mf-Mf_RD_lm+Mf_RD_22; /*Used for the Amplitude as an approx fix for post merger powerlaw slope */
    } else {
        /* For phase */
        REAL8 Rholm = PhenomDQuantities->Rholm[ell][mm];
        ans = Rholm * Mf;          /* Used for the Phase */
    }

    return ans;
}
// mathematica function Trd
// domain mapping function - ringdown
//UNUSED double XLALIMRPhenomHMTrd(REAL8 Mf, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities){
//    return XLALIMRPhenomHMpostfRDflm(Mf, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, mm, PhenomDQuantities);
//}

// mathematica function Ti
// domain mapping function - inspiral
double XLALIMRPhenomHMTi(REAL8 Mf, const INT4 mm){
    return 2.0 * Mf / mm;
}

void XLALIMRPhenomHMSlopeAmAndBm(double *Am, double *Bm, const INT4 mm, REAL8 fi, REAL8 fr, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, PhenomDStorage* PhenomDQuantities);
void XLALIMRPhenomHMSlopeAmAndBm(double *Am, double *Bm, const INT4 mm, REAL8 fi, REAL8 fr, REAL8 Mf_RD_22, REAL8 Mf_RD_lm, const INT4 AmpFlag, const INT4 ell, PhenomDStorage* PhenomDQuantities){
    REAL8 Trd = XLALIMRPhenomHMTrd(fr, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, mm, PhenomDQuantities);
    REAL8 Ti = XLALIMRPhenomHMTi(fi, mm);

    //Am = ( Trd[fr]-Ti[fi] )/( fr - fi );
    *Am = ( Trd-Ti )/( fr - fi );

    //Bm = Ti[fi] - fi*Am;
    *Bm = Ti - fi*(*Am);
}

int XLALIMRPhenomHMMapParams(REAL8 *a, REAL8 *b, REAL8 flm, REAL8 fi, REAL8 fr, REAL8 Ai, REAL8 Bi, REAL8 Am, REAL8 Bm, REAL8 Ar, REAL8 Br){
    // Define function to output map params used depending on
    if ( flm > fi ){
        if ( flm > fr ){
            *a = Ar;
            *b = Br;
        } else {
            *a = Am;
            *b = Bm;
        }
    }
    else {
        *a = Ai;
        *b = Bi;
    };
    return XLAL_SUCCESS;
}

int XLALIMRPhenomHMFreqDomainMapParams(REAL8 *a, REAL8 *b, REAL8 *fi, REAL8 *fr, REAL8 *f1, REAL8 *f2lm, const REAL8 flm, const INT4 ell, const INT4 mm, PhenomDStorage *PhenomDQuantities, const int AmpFlag); 
int XLALIMRPhenomHMFreqDomainMapParams( REAL8 *a,/**< [Out]  */
                                        REAL8 *b,/**< [Out]  */
                                        REAL8 *fi,/**< [Out]  */
                                        REAL8 *fr,/**< [Out]  */
                                        REAL8 *f1,/**< [Out]  */
                                        REAL8 *f2lm,/**< [Out]  */
                                        const REAL8 flm, /**< input waveform frequency */
                                        const INT4 ell, /**< spherical harmonics ell mode */
                                        const INT4 mm, /**< spherical harmonics m mode */
                                        PhenomDStorage *PhenomDQuantities, /**< Stores quantities in order to calculate them only once */
                                        const int AmpFlag /**< is ==1 then computes for amplitude, if ==0 then computes for phase */)
{

    /*check output points are NULL*/
    XLAL_CHECK(a != NULL, XLAL_EFAULT);
    XLAL_CHECK(b != NULL, XLAL_EFAULT);
    XLAL_CHECK(fi != NULL, XLAL_EFAULT);
    XLAL_CHECK(fr != NULL, XLAL_EFAULT);
    XLAL_CHECK(f1 != NULL, XLAL_EFAULT);
    XLAL_CHECK(f2lm != NULL, XLAL_EFAULT);

    /* Account for different f1 definition between PhenomD Amplitude and Phase derivative models */
    REAL8 Mf_1_22  = 0.; /* initalise variable */
    if ( AmpFlag==1 ) {
        /* For amplitude */
        Mf_1_22  = AMP_fJoin_INS; /* inspiral joining frequency from PhenomD [amplitude model], for the 22 mode */
    } else {
        /* For phase */
        Mf_1_22  = PHI_fJoin_INS; /* inspiral joining frequency from PhenomD [phase model], for the 22 mode */
    }

    *f1 = Mf_1_22;

    REAL8 Mf_RD_22 = PhenomDQuantities->Mf_RD_22;
    REAL8 Mf_RD_lm = PhenomDQuantities->PhenomHMfring[ell][mm];

    // Define a ratio of QNM frequencies to be used for scaling various quantities
    REAL8 Rholm = PhenomDQuantities->Rholm[ell][mm];

    // Given experiments with the l!=m modes, it appears that the QNM scaling rather than the PN scaling may be optimal for mapping f1
    REAL8 Mf_1_lm = Mf_1_22 / Rholm;

    /* Define transition frequencies */
    *fi = Mf_1_lm;
    *fr = Mf_RD_lm;

    /*Define the slope and intercepts of the linear transformation used*/
    REAL8 Ai = 2.0/mm;
    REAL8 Bi = 0.0;
    REAL8 Am;
    REAL8 Bm;
    XLALIMRPhenomHMSlopeAmAndBm(&Am, &Bm, mm, *fi, *fr, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, PhenomDQuantities);

    REAL8 Ar = 1.0;
    REAL8 Br = 0.0;
    if ( AmpFlag==1 ) {
        /* For amplitude */
        Br = -Mf_RD_lm+Mf_RD_22;
    } else {
        /* For phase */
        Ar = Rholm;
    }

    /* Define function to output map params used depending on */
    int ret = XLALIMRPhenomHMMapParams(a, b, flm, *fi, *fr, Ai, Bi, Am, Bm, Ar, Br);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMMapParams failed in XLALIMRPhenomHMMapParams (1)\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    REAL8 frfi = 0.5 * ( *fr + *fi );
    *f2lm = frfi;

    return XLAL_SUCCESS;
}

/**
 * XLALSimIMRPhenomHMFreqDomainMap
 * Input waveform frequency in Geometric units (Mflm)
 * and computes what frequency this corresponds
 * to scaled to the 22 mode.
 */
double XLALSimIMRPhenomHMFreqDomainMap(REAL8 Mflm, const INT4 ell, const INT4 mm, PhenomDStorage* PhenomDQuantities, const int AmpFlag);
double XLALSimIMRPhenomHMFreqDomainMap(REAL8 Mflm,
                                        const INT4 ell,
                                        const INT4 mm,
                                        PhenomDStorage* PhenomDQuantities,
                                        const int AmpFlag)
{

    /* Mflm here has the same meaning as Mf_wf in XLALSimIMRPhenomHMFreqDomainMapHM (old deleted function). */
    REAL8 a = 0.;
    REAL8 b = 0.;
    /* Following variables not used in this funciton but are returned in XLALIMRPhenomHMFreqDomainMapParams */
    REAL8 fi = 0.;
    REAL8 fr = 0.;
    REAL8 f1 = 0.;
    REAL8 f2lm = 0.;
    int ret = XLALIMRPhenomHMFreqDomainMapParams(&a, &b, &fi, &fr, &f1, &f2lm, Mflm, ell, mm, PhenomDQuantities, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALSimIMRPhenomHMFreqDomainMapParams\n");
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
double XLALSimIMRPhenomHMPNFrequencyScale( UsefulPowers *p, REAL8 Mf, INT4 ell, INT4 mm);
double XLALSimIMRPhenomHMPNFrequencyScale( UsefulPowers *p,
                                           REAL8 Mf,
                                           INT4 ell,
                                           INT4 mm ) {

    /* Initialise answer */
    REAL8 ans = 0.0;

    //FP: just compute the power in here 
    if ( ell==2 ) {
        if ( mm==2 ) {
            ans = 1.0;
        }
        else { //mm==1
            ans = p->third;
        }
    } else if ( ell==3 ) {
        if ( mm==3 ) {
          ans = p->third;
        }
        else { //mm==2
          ans = p->two_thirds;
        }
    } else if ( ell==4 ) {
        if ( mm==4 ) {
          ans = p->two_thirds;
        }
        else { //mm==3
          ans = Mf;
        }
    } else if ( ell==5 ) {
        if ( mm==5 ) {
          ans = Mf;
        }
        else { //mm==4
          ans = p->four_thirds;
        }
    } else if ( ell==6 ) {
        if ( mm==6 ) {
          ans = p->four_thirds;
        }
        else { //mm==5
          ans = p->five_thirds;
        }
    } else {
        XLALPrintError("XLAL Error - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }

    return ans;

}

/* FIXME: returns leading order PN amplitude for given ell and m mode.
 * This is from mma notebook 'leadingPNamp.nb' in /work/projects/PhenomHM
 */
double XLALSimIMRPhenomHMPNAmplitudeLeadingOrder(INT4 ell, INT4 mm, PhenomDStorage *PhenomDQuantities, UsefulMfPowers *powers_of_Mf_wf);
double XLALSimIMRPhenomHMPNAmplitudeLeadingOrder(INT4 ell, INT4 mm, PhenomDStorage *PhenomDQuantities, UsefulMfPowers *powers_of_Mf_wf) {
    /* Initialise answer */
    REAL8 pow_Mf_wf_prefactor = PhenomDQuantities->pow_Mf_wf_prefactor[ell][mm];
    REAL8 ans = 0.0;

    //FP: do a Rholm style thing here for speed up. 
    if ( ell==2 ) {
        if ( mm==2 ) {
          ans = powers_of_Mf_wf->m_seven_sixths;
        } else { //mm==1
          ans = powers_of_Mf_wf->m_five_sixths;
        }
    } else if ( ell==3 ) {
        if ( mm==3 ) {
          ans = powers_of_Mf_wf->m_five_sixths;
        }
        else { //mm==2
          ans = powers_of_Mf_wf->m_sqrt;
        }
    } else if ( ell==4 ) {
        if ( mm==4 ) {
          ans = powers_of_Mf_wf->m_sqrt;
        }
        else { //mm==3
          ans = powers_of_Mf_wf->m_sixth;
        }
    } else if ( ell==5 ) {
        if ( mm==5 ) {
          ans = powers_of_Mf_wf->m_sixth;
        }
        else { //mm==4
          ans = powers_of_Mf_wf->sixth;
        }
    } else if ( ell==6 ) {
        if ( mm==6 ) {
          ans = powers_of_Mf_wf->sixth;
        }
        else { //mm==5
          ans = powers_of_Mf_wf->sqrt;
        }
    } else {
        XLALPrintError("XLAL Error - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }
    ans *= pow_Mf_wf_prefactor;

    return ans;
}

static double ComputeAmpRatio(INT4 ell, INT4 mm, AmpInsPrefactors amp_prefactors, IMRPhenomDAmplitudeCoefficients *pAmp, PhenomDStorage *PhenomDQuantities, UsefulMfPowers *powers_of_MfAtScale_22_amp, UsefulPowers *downsized_powers_of_MfAtScale_22_amp);
static double ComputeAmpRatio(INT4 ell, INT4 mm, AmpInsPrefactors amp_prefactors, IMRPhenomDAmplitudeCoefficients *pAmp, PhenomDStorage *PhenomDQuantities, UsefulMfPowers *powers_of_MfAtScale_22_amp, UsefulPowers *downsized_powers_of_MfAtScale_22_amp){

    /* See technical document for description of below lines with A_R and R */
    double A_R_num = XLALSimIMRPhenomHMPNAmplitudeLeadingOrder(ell, mm, PhenomDQuantities, powers_of_MfAtScale_22_amp);
    double A_R_den = XLALSimIMRPhenomHMPNFrequencyScale(downsized_powers_of_MfAtScale_22_amp, powers_of_MfAtScale_22_amp->itself, ell, mm) * IMRPhenDAmplitude(powers_of_MfAtScale_22_amp->itself, pAmp, downsized_powers_of_MfAtScale_22_amp, &amp_prefactors);
    double ampRatio = A_R_num/A_R_den;

    return ampRatio;
}

double XLALSimIMRPhenomHMAmplitude(double Mf_wf, int ell, int mm, IMRPhenomDAmplitudeCoefficients *pAmp, AmpInsPrefactors * amp_prefactors, PhenomDStorage * PhenomDQuantities, UsefulMfPowers *powers_of_MfAtScale_22_amp, UsefulPowers *downsized_powers_of_MfAtScale_22_amp);
double XLALSimIMRPhenomHMAmplitude( double Mf_wf,
                                    int ell,
                                    int mm,
                                    IMRPhenomDAmplitudeCoefficients *pAmp,
                                    AmpInsPrefactors * amp_prefactors,
                                    PhenomDStorage * PhenomDQuantities,
                                    UsefulMfPowers *powers_of_MfAtScale_22_amp,
                                    UsefulPowers *downsized_powers_of_MfAtScale_22_amp
                                  )
{
    double Mf_22 =  XLALSimIMRPhenomHMFreqDomainMap(Mf_wf, ell, mm, PhenomDQuantities, AmpFlagTrue);

    UsefulPowers powers_of_Mf_22;
    int errcode = XLAL_SUCCESS;
    errcode = init_useful_powers(&powers_of_Mf_22, Mf_22);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_powers failed for Mf_22");

    double PhenDamp = IMRPhenDAmplitude(Mf_22, pAmp, &powers_of_Mf_22, amp_prefactors);

    double ampRatio = ComputeAmpRatio(ell, mm, *amp_prefactors, pAmp, PhenomDQuantities, powers_of_MfAtScale_22_amp, downsized_powers_of_MfAtScale_22_amp);

    double R = ampRatio * XLALSimIMRPhenomHMPNFrequencyScale(&powers_of_Mf_22, Mf_22, ell, mm);

    double HMamp = PhenDamp * R;

    return HMamp;
}


/*TODO Should also probably add LALDict as a argument here.*/
int XLALSimIMRPhenomHMPhasePreComp(HMPhasePreComp *q, const INT4 ell, const INT4 mm, const REAL8 eta, const REAL8 chi1z, const REAL8 chi2z, PhenomDStorage *PhenomDQuantities);
int XLALSimIMRPhenomHMPhasePreComp(HMPhasePreComp *q, const INT4 ell, const INT4 mm, const REAL8 eta, const REAL8 chi1z, const REAL8 chi2z, PhenomDStorage *PhenomDQuantities)
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

    int ret = XLALIMRPhenomHMFreqDomainMapParams(&ai, &bi, &fi, &fr, &f1, &f2lm, Mfshift, ell, mm, PhenomDQuantities, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALIMRPhenomHMFreqDomainMapParams - inspiral\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->ai = ai;
    q->bi = bi;

    ret = XLALIMRPhenomHMFreqDomainMapParams(&a2lm, &b2lm, &fi, &fr, &f1, &f2lm, f2lm+Mfshift, ell, mm, PhenomDQuantities, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALIMRPhenomHMFreqDomainMapParams - intermediate\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->a2lm = a2lm;
    q->b2lm = b2lm;

    ret = XLALIMRPhenomHMFreqDomainMapParams(&ar, &br, &fi, &fr, &f1, &f2lm, fr+Mfshift, ell, mm, PhenomDQuantities, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - XLALIMRPhenomHMFreqDomainMapParams failed in XLALIMRPhenomHMFreqDomainMapParams - merger-ringdown\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->ar = ar;
    q->br = br;

    q->fi = fi;
    q->fr = fr;
    q->f2lm = f2lm;

    const LALSimInspiralTestGRParam *extraParams = NULL;
    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, PhenomDQuantities->finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);

    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, PhenomDQuantities->m1, PhenomDQuantities->m2, chi1z, chi2z, 1.0, 1.0, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSpinSpin(PhenomDQuantities->m1, PhenomDQuantities->m2, PhenomDQuantities->Mtot, eta, chi1z, chi2z) * pn->v[0]);

    PhiInsPrefactors phi_prefactors;
    int status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

    double Rholm = PhenomDQuantities->Rholm[ell][mm];
    double Taulm = PhenomDQuantities->Taulm[ell][mm];

    /* Compute coefficients to make phase C^1 continuous (phase and first derivative) */
    ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);

    REAL8 PhDBMf = a2lm*fi + b2lm;
    UsefulPowers powers_of_PhDBMf;
    status = init_useful_powers(&powers_of_PhDBMf, PhDBMf);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_PhDBMf failed");
    q->PhDBconst = IMRPhenDPhase(PhDBMf, pPhi, pn, &powers_of_PhDBMf, &phi_prefactors, Rholm, Taulm)/a2lm;

    REAL8 PhDCMf = a2lm*f2lm + b2lm;

    if (PhDCMf < 0.){

        printf("eta = %f\n", eta);
        printf("chi1z = %f\n", chi1z);
        printf("chi2z = %f\n", chi2z);
        printf("ell = %i\n", ell);
        printf("mm = %i\n", mm);

        printf("a2lm = %f\n", a2lm);
        printf("f2lm = %f\n", f2lm);
        printf("b2lm = %f\n", b2lm);
        printf("PhDCMf = %f\n", PhDCMf);

    }

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

    LALFree(pPhi);
    LALFree(pn);

    return XLAL_SUCCESS;

}

double XLALSimIMRPhenomHMPhase( double Mf_wf, int ell, int mm, HMPhasePreComp *q, PNPhasingSeries *pn, IMRPhenomDPhaseCoefficients *pPhi, PhiInsPrefactors *phi_prefactors, double Rholm, double Taulm );
double XLALSimIMRPhenomHMPhase( double Mf_wf, /**< input frequency in geometric units*/
                                int ell,
                                int mm,
                                HMPhasePreComp *q,
                                PNPhasingSeries *pn,
                                IMRPhenomDPhaseCoefficients *pPhi,
                                PhiInsPrefactors *phi_prefactors,
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

    if ( !(Mf_wf > q->fi) ){
        Mf = q->ai * Mf_wf + q->bi;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->ai;
    } else if ( !(Mf_wf > q->f2lm) ){
        Mf = q->a2lm*Mf_wf + q->b2lm;
        UsefulPowers powers_of_Mf;
        status = init_useful_powers(&powers_of_Mf, Mf);
        XLAL_CHECK(XLAL_SUCCESS == status, status, "init_useful_powers for powers_of_Mf failed");
        retphase = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_Mf, phi_prefactors, Rholm, Taulm) / q->a2lm - q->PhDBconst + q->PhDBAterm;
    } else if ( !(Mf_wf > q->fr) ){
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

    /* 
     * Phase shift due to leading order complex amplitude
     * [L.Blancet, arXiv:1310.1528 (Sec. 9.5)]
     * "Spherical hrmonic modes for numerical relativity"
     */
    /*TODO: Need to have a list at the begining and a function to check the input
    lm mode to see if it is one that is included in the model.*/
    /* Initialise answer */
    REAL8 cShift = 0.0;
    if ( mm==1 ) {
        cShift = LAL_PI_2; /* i shift */
    } else if ( mm==3 ) {
        cShift = -LAL_PI_2; /* -i shift */
    } else if ( mm==4 ) {
        cShift = LAL_PI; /* -1 shift */
    } else if ( mm==5 ) {
        cShift = LAL_PI_2; /* i shift */
    }
    retphase += cShift;

    // LALFree(pPhi);
    // LALFree(pn);

    return retphase;
}

/********************* Function to add modes for frequency-domain structures ********************/

/*
 * Helper function to add a mode to hplus, hcross in Fourier domain
 * - copies the function XLALSimAddMode, which was done only for TD structure
 * This function was lifted from the EOBNRv2HM_ROM code 
 */
static INT4 FDAddMode(COMPLEX16FrequencySeries *hptilde, COMPLEX16FrequencySeries *hctilde, COMPLEX16FrequencySeries *hlmtilde, REAL8 theta, REAL8 phi, INT4 l, INT4 m, INT4 sym);
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
  if ( sym ) { /* Equatorial symmetry: add in -m mode */
    Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
    COMPLEX16 Ymstar = conj(XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, -m));
    COMPLEX16 factorp = 0.5*(Y + minus1l*Ymstar);
    COMPLEX16 factorc = I*0.5*(Y - minus1l*Ymstar);
    //COMPLEX16* datap = hptilde->data->data;
    //COMPLEX16* datac = hctilde->data->data;
    for ( j = 0; j < hlmtilde->data->length; ++j ) {
      hlmtildevalue = (hlmtilde->data->data[j]);
      //datap[j] += factorp*hlmtildevalue;
      //datac[j] += factorc*hlmtildevalue;
      hptilde->data->data[j] += factorp*hlmtildevalue;
      hctilde->data->data[j] += factorc*hlmtildevalue;
    }
  }
  else { /* not adding in the -m mode */
    Y = XLALSpinWeightedSphericalHarmonic(theta, phi, -2, l, m);
    COMPLEX16 factorp = 0.5*Y;
    COMPLEX16 factorc = I*factorp;
    //COMPLEX16* datap = hptilde->data->data;
    //COMPLEX16* datac = hctilde->data->data;
    for ( j = 0; j < hlmtilde->data->length; ++j ) {
      hlmtildevalue = (hlmtilde->data->data[j]);
      //datap[j] += factorp*hlmtildevalue;
      //datac[j] += factorc*hlmtildevalue;
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

/*
 * This returns the hlm of a single (l,m) mode at a single frequency Mf
 */
static COMPLEX16 IMRPhenomHMSingleModehlm(
        int ell,
        int mm,
        double Mf,
        HMPhasePreComp *z,
        IMRPhenomDAmplitudeCoefficients *pAmp,
        AmpInsPrefactors *amp_prefactors,
        PNPhasingSeries *pn,
        IMRPhenomDPhaseCoefficients *pPhi,
        PhiInsPrefactors *phi_prefactors,
        double Rholm,
        double Taulm,
        double phi_precalc, /**< m*phi0 + HMphaseRef*/
        PhenomDStorage *PhenomDQuantities,
        UsefulMfPowers *powers_of_MfAtScale_22_amp,
        UsefulPowers *downsized_powers_of_MfAtScale_22_amp
) {

    /*
     * In this function we should pass the phenom model parameters
     * and the (l,m) mode to return a given hlm, not summed with
     * spherical harmonics.
     * Can be evaluated at a single geometric frequency (Mf).
     */

    double HMamp = XLALSimIMRPhenomHMAmplitude( Mf, ell, mm, pAmp, amp_prefactors, PhenomDQuantities, powers_of_MfAtScale_22_amp, downsized_powers_of_MfAtScale_22_amp );
    double HMphase = XLALSimIMRPhenomHMPhase( Mf, ell, mm, z, pn, pPhi, phi_prefactors, Rholm, Taulm );

    /* Compute reference phase at reference frequency */

    /* Factor of m spherical harmonic mode b/c phi0 is orbital phase */
    /* NOTE: Does HMphaseRef already have the mm scaling? as it's m*(phi0 + phiref) */
    HMphase -= phi_precalc;

    // Debug lines
    // double Seta = sqrt(1.0 - 4.0*pPhi->eta);
    // double m1 = 0.5 * (1.0 + Seta);
    // double m2 = 0.5 * (1.0 - Seta);
    // const REAL8 M_sec = (m1+m2) * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
    // if ( Mf > 0.05 && Mf < 0.051 ){
    //     printf("Mf = %.9f Hz = %f mode l = %i, m = %i       HMphase -= phi_precalc = %.9f\n", Mf, Mf/M_sec, ell, mm, HMphase);
    // }

    COMPLEX16 hlm = HMamp * cexp(-I * HMphase);
    // printf("Mf = %f     hlm  =  %f + i %f\nx",Mf, creal(hlm), cimag(hlm));

    return hlm;
}

/* Given the final frequency in Mf and the total mass, calculate the final frequency in Hz */
static REAL8 ComputeIMRPhenomHMfmax(REAL8 Mf, REAL8 f_min, REAL8 f_max, REAL8 M);
static REAL8 ComputeIMRPhenomHMfmax(REAL8 Mf    /**< geometric frequency */,
                                    REAL8 f_min /**< low frequency in Hz */,
                                    REAL8 f_max /**< end frequency in Hz */,
                                    REAL8 M     /**< total mass (Msun) */
                                   ){

      const REAL8 M_sec = M * LAL_MTSUN_SI; // Conversion factor Hz -> dimensionless frequency
      const REAL8 fCut = Mf/M_sec; // convert Mf -> Hz

      /*
       * Somewhat arbitrary end point for the waveform.
       * Chosen so that the end of the waveform is well after the ringdown.
       */
      if (!(fCut > f_min))
          XLAL_ERROR(XLAL_EDOM, "(fCut = %g Hz) <= f_min = %g\n", fCut, f_min);

      /* Default f_max to Cut */
      REAL8 f_max_prime = f_max;
      f_max_prime = f_max ? f_max : fCut;
      f_max_prime = (f_max_prime > fCut) ? fCut : f_max_prime;
      if (!(f_max_prime > f_min))
          XLAL_ERROR(XLAL_EDOM, "f_max <= f_min\n");

      return f_max_prime;
}

/* Compute t0 as the time of the peak of the 22 mode */
static REAL8 Computet0(REAL8 eta, REAL8 chi1z, REAL8 chi2z, REAL8 finspin);
static REAL8 Computet0(REAL8 eta, REAL8 chi1z, REAL8 chi2z, REAL8 finspin){

    const LALSimInspiralTestGRParam *extraParams = NULL;
    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);

    IMRPhenomDAmplitudeCoefficients *pAmp;
    pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
    ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);

    // double Rholm = XLALSimIMRPhenomHMRholm(eta, chi1z, chi2z, ell, mm);
    // double Taulm = XLALSimIMRPhenomHMTaulm(eta, chi1z, chi2z, ell, mm);

    //time shift so that peak amplitude is approximately at t=0
    //For details see https://www.lsc-group.phys.uwm.edu/ligovirgo/cbcnote/WaveformsReview/IMRPhenomDCodeReview/timedomain
    //NOTE: All modes will have the same time offset. So we use the 22 mode.
    //If we just use the 22 mode then we pass 1.0, 1.0 into DPhiMRD.
    const REAL8 t0 = DPhiMRD(pAmp->fmaxCalc, pPhi, 1.0, 1.0);

    LALFree(pPhi);
    LALFree(pAmp);

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


/*
 * The goal of this function is to compute hlm for a list of modes
 * and output them in the data type SphHarmFrequencySeries.
 */

/* NOTE: name changed from hlmsphharmfreqseries to hlms */
/* NOTE: Some code duplication in order to keep this function XLAL */
int XLALIMRPhenomHMMultiModehlm(SphHarmFrequencySeries **hlms, REAL8 m1Msun, REAL8 m2Msun, REAL8 chi1z, REAL8 chi2z, REAL8 deltaF, REAL8 f_min, REAL8 f_max, REAL8 fRef_in, REAL8 phi0, REAL8 distance);
int XLALIMRPhenomHMMultiModehlm(
    SphHarmFrequencySeries **hlms, /**< [out] can access multiple modes with units of Mf */
    REAL8 m1Msun,                  /**< primary mass in Msun */
    REAL8 m2Msun,                  /**< secondary mass in Msun */
    REAL8 chi1z,                   /**< primary spin parameter */
    REAL8 chi2z,                   /**< secondary spin parameter */
    REAL8 deltaF,                  /**< Hz */
    REAL8 f_min,                   /**< Hz */
    REAL8 f_max,                   /**< Hz */
    REAL8 fRef_in,                 /**< reference frequency in Hz */
    REAL8 phi0,                    /**< reference orbital phase */
    REAL8 distance                 /**< distance to source in SI */
) {

    /* Powers of pi */
    int errcode = XLAL_SUCCESS;
    //FP: these are known, so turn them in to #define's
    errcode = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_useful_powers() failed: failed to initiate useful powers of pi.");

    /* Here masses are in Msun */
    errcode = EnforcePrimaryIsm1(&m1Msun, &m2Msun, &chi1z, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryIsm1 failed");

    /* Compute quantities/parameters related to PhenomD only once and store them */
    PhenomDStorage PhenomDQuantities;
    errcode = init_PhenomD_Storage(&PhenomDQuantities, m1Msun/LAL_MSUN_SI, m2Msun/LAL_MSUN_SI, chi1z, chi2z);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_PhenomD_Storage failed");
    
    const REAL8 M = PhenomDQuantities.Mtot*LAL_MSUN_SI;
    const REAL8 eta = PhenomDQuantities.eta;
    const REAL8 M_sec = M * LAL_MTSUN_SI;

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    if (abs(chi1z) > 1.0 || abs(chi2z) > 1.0)
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    /* If no reference frequency given, set it to the starting GW frequency */
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    REAL8 f_max_prime = ComputeIMRPhenomHMfmax(Mf_CUT_HM, f_min, f_max, M);

    /* Compute phenomD amp coefficients */
    IMRPhenomDAmplitudeCoefficients *pAmp;
    pAmp = (IMRPhenomDAmplitudeCoefficients *) XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
    ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, PhenomDQuantities.finspin);
    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    AmpInsPrefactors amp_prefactors;
    errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");

    /* Compute phenomD phase coefficients */
    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = (IMRPhenomDPhaseCoefficients *) XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    const LALSimInspiralTestGRParam *extraParams = NULL;
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, PhenomDQuantities.finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);

    PNPhasingSeries *pn = NULL;
    // FP: Malloc?
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1Msun, m2Msun, chi1z, chi2z, 1.0, 1.0, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    REAL8 testGRcor=1.0;
    if (extraParams!=NULL)
    {
        if (XLALSimInspiralTestGRParamExists(extraParams,"dchi6"))  testGRcor += XLALSimInspiralGetTestGRParam(extraParams,"dchi6");
    }

    // Was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSpinSpin(m1Msun, m2Msun, M, eta, chi1z, chi2z) * pn->v[0])* testGRcor;

    //FP: Malloc and free this too?
    PhiInsPrefactors phi_prefactors;
    //phi_prefactors = (PhiInsPrefactors *) XLALMalloc(sizeof(PhiInsPrefactors));
    int status = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == status, status, "init_phi_ins_prefactors failed");

    /* Compute the amplitude pre-factor */
    const REAL8 amp0 = M * LAL_MRSUN_SI * M_sec / distance;

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    REAL8 InvDeltaF = 1./deltaF; 
    /* Shift by overall length in time */
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, - InvDeltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    size_t n = NextPow2(f_max_prime * InvDeltaF) + 1;
    /* Range that will have actual non-zero waveform values generated */
    size_t ind_min = (size_t) (f_min * InvDeltaF);
    size_t ind_max = (size_t) (f_max_prime * InvDeltaF);
    XLAL_CHECK ( !(ind_max>n) && !(ind_min>ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=hptilde->data>length=%zu.", ind_min, ind_max, n);

    const REAL8 MfRef = M_sec * fRef;

    /*
     * In this function we should setup all the variables that we can
     * precompute to generate the PhenomD amplitude and phase functions
     * We then loop over the static function 'IMRPhenomHMSingleModehlm'
     * to generate the modes. We then sum them up at the end.
     */

    /* Now we have all the PhenomD model parameters, which actually correspond
     * to the (l,m)=(2,2) parameters, we can feed these into the function
     * IMRPhenomHMSingleModehlm to generate any mode we want. */

    /*
     * NOTE: I'm not sure what Mf should be used for the reference time...
     * I think it should be the scaled one. And it should come from the amplitude
     */

     const REAL8 t0 = Computet0(eta, chi1z, chi2z, PhenomDQuantities.finspin);

     for( int j=0; j<NMODES; j++ ){

         int ell = ModeArray[j][0];
         int mm = ModeArray[j][1];

         double Rholm = PhenomDQuantities.Rholm[ell][mm];
         double Taulm = PhenomDQuantities.Taulm[ell][mm];

         /* Compute coefficients to make phase C^1 continuous (phase and first derivative) */
         ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);

         /* PhenomHM pre-computations */
         /* NOTE: Need to make this an input and NOT part of the frequency loop! */
         HMPhasePreComp z;
         int ret = XLALSimIMRPhenomHMPhasePreComp(&z, ell, mm, eta, chi1z, chi2z, &PhenomDQuantities);
         if (ret != XLAL_SUCCESS){
             XLALPrintError("XLAL Error - XLALSimIMRPhenomHMPhasePreComp failed\n");
             XLAL_ERROR(XLAL_EDOM);
         }

         double HMphaseRef = XLALSimIMRPhenomHMPhase( MfRef, ell, mm, &z, pn, pPhi, &phi_prefactors, Rholm, Taulm );
         HMphaseRef = HMphaseRef * 0.; //FIXME: This seems unnecessary but maybe it's a placeholder?

         /* Compute reference phase at reference frequency */

         // Factor of m spherical harmonic mode b/c phi0 is orbital phase
         /* NOTE: Does HMphaseRef already have the mm scaling? as it's m*(phi0 + phiref) */
         /* NOTE: Because phenomD and phenomHM use the function XLALSimInspiralTaylorF2AlignedPhasing
         to generate the inspiral SPA TF2 phase it does NOT contain the -LAL_PI / 4. phase shift.
         if it did there would be an extra -(2.0-m) * LAL_PI/8.0 term here.*/
         // REAL8 phi_precalc = mm*phi0 + HMphaseRef;
         /* FIXME: */
         /* NOTE: Quick and dirty fix to get a reference phase working.*/
         /* NOTE: Moving the reference phase to the spherical harmonic.*/
         REAL8 phi_precalc = phi0 * 0.;

         //FP: malloc and then free?
         /* compute amplitude ratio correction to take 22 mode in to (ell, mm) mode amplitude */
         double MfAtScale_22_amp = XLALSimIMRPhenomHMFreqDomainMap( MfAtScale_wf_amp, ell, mm, &PhenomDQuantities, AmpFlagTrue );
         UsefulMfPowers powers_of_MfAtScale_22_amp;
         errcode = init_useful_mf_powers(&powers_of_MfAtScale_22_amp, MfAtScale_22_amp);
         XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_mf_powers failed for MfAtScale_22_amp");
         UsefulPowers downsized_powers_of_MfAtScale_22_amp;
         errcode = downsize_useful_mf_powers(&downsized_powers_of_MfAtScale_22_amp, &powers_of_MfAtScale_22_amp);
         XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "downsized_init_useful_powers failed");

         /* We loop over (l,m) and use a temporary hlm frequency series to store the results of a single mode */
         COMPLEX16FrequencySeries *hlm = NULL;
         hlm = XLALCreateCOMPLEX16FrequencySeries("hlm: single mode", &ligotimegps_zero, 0.0, deltaF, &lalStrainUnit, n);
         memset(hlm->data->data, 0, n * sizeof(COMPLEX16));

         /* NOTE: Do I need this bit? */
         /* XLALUnitMultiply(hlm->sampleUnits, hlm->sampleUnits, &lalSecondUnit); */

         /* Now generate the waveform for a single (l,m) mode, i.e. compute hlm*/
         /* Loop over frequency */
         REAL8 M_sec_dF = M_sec * deltaF;
         COMPLEX16 It0 = I*t0; 
         #pragma omp parallel for
         for (size_t i = ind_min; i < ind_max; i++)
         {
            REAL8 Mf = i * M_sec_dF; /* geometric frequency */
            /* TODO: fix phase and time offsets */
            /* TODO: inclusion of reference frequency */
            // phi -= t0*(Mf-MfRef) + phi_precalc;
            // ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
            /* construct hlm at single frequency point and return */
            // (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(eta, chi1z, chi2z, ell, mm, Mf, MfRef, phi0, &z);

            (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(ell, mm, Mf, &z, pAmp, &amp_prefactors, pn, pPhi, &phi_prefactors, Rholm, Taulm, phi_precalc, &PhenomDQuantities, &powers_of_MfAtScale_22_amp, &downsized_powers_of_MfAtScale_22_amp);
            /* NOTE: The frequency used in the time shift term is the fourier variable of the gravitational wave frequency. i.e., Not rescaled. */
            /* NOTE: normally the t0 term is multiplied by 2pi but the 2pi has been absorbed into the t0. */
            // (hlm->data->data)[i] *= cexp(-I * LAL_PI*t0*(Mf-MfRef)*(2.0-mm) );
            // (hlm->data->data)[i] *= cexp(-I * t0*(Mf-MfRef)*(2.0-mm) );
            (hlm->data->data)[i] *= cexp(-It0*(Mf-MfRef)); //FIXME: 2.0-mm is gone, is this intended? If yes, delete this comment.
            /* From phenomD for referene*/
            // REAL8 amp = IMRPhenDAmplitude(Mf, pAmp, &powers_of_f, &amp_prefactors);
            // REAL8 phi = IMRPhenDPhase(Mf, pPhi, pn, &powers_of_f, &phi_prefactors);
            // phi -= t0*(Mf-MfRef) + phi_precalc;
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

/* This function is a wrapper of XLALIMRPhenomHMMultiModehlm.
 * It takes the SphHarmFrequencySeries **hlms as inputs and
 * compbines them with the Spherical harmonics at a given
 * inclination angle and azimuthal angle. It returns the
 * hptilde and hctilde.
 */
 /*TODO: Add LALDict as an argument to XLALIMRPhenomHMMultiModeStrain */
 /*TODO: Change name from XLALIMRPhenomHMMultiModeStrain to XLALIMRPhenomHM */
int XLALIMRPhenomHMMultiModeStrain(
    COMPLEX16FrequencySeries **hptilde, /**< [out] */
    COMPLEX16FrequencySeries **hctilde, /**< [out] */
    REAL8 m1,                           /**< primary mass in SI */
    REAL8 m2,                           /**< secondary mass in SI */
    REAL8 chi1z,                        /**< primary spin parameter */
    REAL8 chi2z,                        /**< secondary spin parameter */
    REAL8 deltaF,                       /**< Hz */
    REAL8 f_min,                        /**< Hz */
    REAL8 f_max,                        /**< Hz */
    REAL8 fRef_in,                      /**< reference frequency in Hz */
    REAL8 phi0,                         /**< reference orbital phase */
    REAL8 inclination,                  /**< inclination... */
    REAL8 distance                      /**< distance to source in SI */
) {

    /* Sanity checks on input parameters */
    if (m1 < 0) XLAL_ERROR(XLAL_EDOM, "m1 must be positive\n");
    if (m2 < 0) XLAL_ERROR(XLAL_EDOM, "m2 must be positive\n");
    if (abs(chi1z) > 1.0 || abs(chi2z) > 1.0 )
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");
    if (deltaF < 0) XLAL_ERROR(XLAL_EDOM, "deltaF must be positive\n");
    if (f_min < 0) XLAL_ERROR(XLAL_EDOM, "f_min must be positive\n");
    if (f_max < 0) XLAL_ERROR(XLAL_EDOM, "f_max must be greater than 0\n");

    /* External: SI; internal: solar masses */
    m1 /= LAL_MSUN_SI;
    m2 /= LAL_MSUN_SI;

    int ret = EnforcePrimaryIsm1(&m1, &m2, &chi1z, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "EnforcePrimaryIsm1 failed");

    /* Mass ratio >= 1 convention */
    const REAL8 q = (m1 > m2) ? (m1 / m2) : (m2 / m1);

    if (q > MAX_ALLOWED_MASS_RATIO)
      XLAL_PRINT_WARNING("Warning: The model is not supported for high mass ratio, see MAX_ALLOWED_MASS_RATIO\n");

    const REAL8 M = m1 + m2; /* total mass (Msun) */
    const REAL8 eta = m1 * m2 / (M * M);

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    /* If no reference frequency given, set it to the starting GW frequency */
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    REAL8 f_max_prime = ComputeIMRPhenomHMfmax(Mf_CUT_HM, f_min, f_max, M);

    /* Evaluate XLALIMRPhenomHMMultiModehlm */

    // SphHarmFrequencySeries *hlms=NULL;
    SphHarmFrequencySeries **hlms=XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlms=NULL;

    /*TODO: Add LALDict as an argument to XLALIMRPhenomHMMultiModehlm */
    ret = XLALIMRPhenomHMMultiModehlm(hlms, m1, m2, chi1z, chi2z, deltaF, f_min, f_max_prime, fRef, phi0, distance);
    XLAL_CHECK(XLAL_SUCCESS == ret, ret, "XLALIMRPhenomHMMultiModehlm(&hlms) failed");

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    /* Shift by overall length in time */
    REAL8 InvDeltaF = 1./deltaF;
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, - InvDeltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    /* Compute array sizes */
    size_t n = NextPow2(f_max_prime * InvDeltaF) + 1;
    /* Range that will have actual non-zero waveform values generated */
    size_t ind_min = (size_t) (f_min * InvDeltaF);
    size_t ind_max = (size_t) (f_max_prime * InvDeltaF);
    XLAL_CHECK ( !(ind_max>n) && !(ind_min>ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=hptilde->data>length=%zu.", ind_min, ind_max, n);

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
      COMPLEX16FrequencySeries* hlm = XLALSphHarmFrequencySeriesGetMode(*hlms, ell, mm);
      if (!(hlm)) XLAL_ERROR(XLAL_EFUNC);

      /* We test for hypothetical m=0 modes */
      if ( mm==0 ) {
          sym = 0;
      } else {
          sym = 1;
      }
    //   FDAddMode( *hptilde, *hctilde, hlm, inclination, 0., ell, mm, sym); /* The phase \Phi is set to 0 - assumes phiRef is defined as half the phase of the 22 mode h22 (or the first mode in the list), not for h = hplus-I hcross */
    FDAddMode( *hptilde, *hctilde, hlm, inclination, phi0, ell, mm, sym); /* Added phi0 here as a quick fix for the reference phase. not sure if it should be m * phi0 or m/2*phi0 . */
    }

    XLALDestroySphHarmFrequencySeries(*hlms);
    XLALFree(hlms);

    return XLAL_SUCCESS;
}


/* Convenience function to compute hlm for a single mode - from which
 * the Alm and philm (amplitude and phase) of a particular mode can be
 * obtained. Useful to compare to NR */
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

    int errcode = XLAL_SUCCESS;
    errcode = init_useful_powers(&powers_of_pi, LAL_PI);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_useful_powers() failed: failed to initiate useful powers of pi.");

    /* Here masses are in Msun */
    errcode = EnforcePrimaryIsm1(&m1Msun, &m2Msun, &chi1z, &chi2z);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "EnforcePrimaryIsm1 failed");

    // printf(" m1Msun= %f\n",m1Msun);
    // printf(" m2Msun= %f\n",m2Msun);
    // printf(" chi1z= %f\n",chi1z);
    // printf(" chi2z= %f\n",chi2z);

    PhenomDStorage PhenomDQuantities;
    errcode = init_PhenomD_Storage(&PhenomDQuantities, m1Msun/LAL_MSUN_SI, m2Msun/LAL_MSUN_SI, chi1z, chi2z);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_PhenomD_Storage failed");

    const REAL8 M = PhenomDQuantities.Mtot*LAL_MSUN_SI;//m1Msun + m2Msun;
    const REAL8 eta = PhenomDQuantities.eta;
    const REAL8 M_sec = M * LAL_MTSUN_SI;

    if (eta > 0.25 || eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");

    if (abs(chi1z) > 1.0 || abs(chi2z) > 1.0)
        XLAL_ERROR(XLAL_EDOM, "Spins outside the range [-1,1] are not supported\n");

    // if no reference frequency given, set it to the starting GW frequency
    REAL8 fRef = (fRef_in == 0.0) ? f_min : fRef_in;

    REAL8 f_max_prime = ComputeIMRPhenomHMfmax(Mf_CUT_HM, f_min, f_max, M);


    /* Compute PhenomD amp and phase coefficients*/

    IMRPhenomDAmplitudeCoefficients *pAmp;
    pAmp = XLALMalloc(sizeof(IMRPhenomDAmplitudeCoefficients));
    ComputeIMRPhenomDAmplitudeCoefficients(pAmp, eta, chi1z, chi2z, PhenomDQuantities.finspin);

    // printf("pAmp->chi1 = %f\n", pAmp->chi1);

    if (!pAmp) XLAL_ERROR(XLAL_EFUNC);
    AmpInsPrefactors amp_prefactors;
    errcode = init_amp_ins_prefactors(&amp_prefactors, pAmp);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_amp_ins_prefactors() failed.");
    const LALSimInspiralTestGRParam *extraParams = NULL;
    IMRPhenomDPhaseCoefficients *pPhi;
    pPhi = XLALMalloc(sizeof(IMRPhenomDPhaseCoefficients));
    ComputeIMRPhenomDPhaseCoefficients(pPhi, eta, chi1z, chi2z, PhenomDQuantities.finspin, extraParams);
    if (!pPhi) XLAL_ERROR(XLAL_EFUNC);
    PNPhasingSeries *pn = NULL;
    XLALSimInspiralTaylorF2AlignedPhasing(&pn, m1Msun, m2Msun, chi1z, chi2z, 1.0, 1.0, LAL_SIM_INSPIRAL_SPIN_ORDER_35PN, extraParams);
    if (!pn) XLAL_ERROR(XLAL_EFUNC);


    // printf("pPhi->beta1 = %f\n", pPhi->beta1);

    // Subtract 3PN spin-spin term below as this is in LAL's TaylorF2 implementation
    // (LALSimInspiralPNCoefficients.c -> XLALSimInspiralPNPhasing_F2), but
    REAL8 testGRcor=1.0;
    if (extraParams!=NULL)
    {
  	  if (XLALSimInspiralTestGRParamExists(extraParams,"dchi6"))  testGRcor += XLALSimInspiralGetTestGRParam(extraParams,"dchi6");
    }

    // was not available when PhenomD was tuned.
    pn->v[6] -= (Subtract3PNSpinSpin(m1Msun, m2Msun, M, eta, chi1z, chi2z) * pn->v[0])* testGRcor;

    PhiInsPrefactors phi_prefactors;
    errcode = init_phi_ins_prefactors(&phi_prefactors, pPhi, pn);
    XLAL_CHECK(XLAL_SUCCESS == errcode, errcode, "init_phi_ins_prefactors failed");
    /* NOTE: There seems to be a problem here, with the phi_prefactors as in IMRPhenomHMMultiModehlm they work fine here but in this
    function they do not. */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */
    /*FIXME */

    /* Compute the amplitude pre-factor */
    const REAL8 amp0 = M * LAL_MRSUN_SI * M * LAL_MTSUN_SI / distance;

    LIGOTimeGPS ligotimegps_zero = LIGOTIMEGPSZERO; // = {0, 0}

    /* Coalesce at t=0 */
    REAL8 InvDeltaF = 1./deltaF;
    // shift by overall length in time
    XLAL_CHECK ( XLALGPSAdd(&ligotimegps_zero, -InvDeltaF), XLAL_EFUNC, "Failed to shift coalescence time to t=0, tried to apply shift of -1.0/deltaF with deltaF=%g.", deltaF);

    /* compute array sizes */
    size_t n = NextPow2(f_max_prime * InvDeltaF) + 1;
    /* range that will have actual non-zero waveform values generated */
    size_t ind_min = (size_t) (f_min * InvDeltaF);
    size_t ind_max = (size_t) (f_max_prime * InvDeltaF);
    XLAL_CHECK ( !(ind_max>n) && !(ind_min>ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=hptilde->data>length=%zu.", ind_min, ind_max, n);


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
     /* NOTE: We could compute the t0 here as the time of the peak of the
        22 mode and make that the assumtion.*/


     /*
      * NOTE: I'm not sure what Mf should be used for the reference time... I think it should be the scaled one. And it should come from the amplitude
      */

     const REAL8 t0 = Computet0(eta, chi1z, chi2z, PhenomDQuantities.finspin);

     double Rholm = PhenomDQuantities.Rholm[ell][mm];
     double Taulm = PhenomDQuantities.Taulm[ell][mm];

     // Compute coefficients to make phase C^1 continuous (phase and first derivative)
     ComputeIMRPhenDPhaseConnectionCoefficients(pPhi, pn, &phi_prefactors, Rholm, Taulm);
     // printf("pPhi->C1Int = %f\n", pPhi->C1Int);
     /* compute phenomHM pre computations */
     /* NOTE: Need to make this an input and NOT part of the frequency loop! */
     HMPhasePreComp z;
     errcode = XLALSimIMRPhenomHMPhasePreComp(&z, ell, mm, eta, chi1z, chi2z, &PhenomDQuantities);
     if (errcode != XLAL_SUCCESS){
         XLALPrintError("XLAL Error - XLALSimIMRPhenomHMPhasePreComp failed\n");
         XLAL_ERROR(XLAL_EDOM);
     }

     double HMphaseRef = XLALSimIMRPhenomHMPhase( MfRef, ell, mm, &z, pn, pPhi, &phi_prefactors, Rholm, Taulm );
     // printf("HMphaseRef = %f\n",HMphaseRef);
     /* compute reference phase at reference frequency */

     // factor of m spherical harmonic mode b/c phi0 is orbital phase
     /* NOTE: Does HMphaseRef already have the mm scaling? as it's m*(phi0 + phiref) */
     /* NOTE: Because phenomD and phenomHM use the function XLALSimInspiralTaylorF2AlignedPhasing
     to generate the inspiral SPA TF2 phase it does NOT contain the -LAL_PI / 4. phase shift.
     if it did there would be an extra -(2.0-m) * LAL_PI/8.0 term here.*/
     REAL8 phi_precalc = mm*phi0 + HMphaseRef;
     //FP: malloc and then free?
     /* compute amplitude ratio correction to take 22 mode in to (ell, mm) mode amplitude */
     double MfAtScale_22_amp = XLALSimIMRPhenomHMFreqDomainMap( MfAtScale_wf_amp, ell, mm, &PhenomDQuantities, AmpFlagTrue );
     UsefulMfPowers powers_of_MfAtScale_22_amp;
     errcode = init_useful_mf_powers(&powers_of_MfAtScale_22_amp, MfAtScale_22_amp);
     XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "init_useful_mf_powers failed for MfAtScale_22_amp");
     UsefulPowers downsized_powers_of_MfAtScale_22_amp;
     errcode = downsize_useful_mf_powers(&downsized_powers_of_MfAtScale_22_amp, &powers_of_MfAtScale_22_amp);
     XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "downsized_init_useful_powers failed");

     /* NOTE: Do I need this bit? */
     /* XLALUnitMultiply(hlm->sampleUnits, hlm->sampleUnits, &lalSecondUnit); */

     /* LOOP OVER FREQUENCY */
     /* Now generate the waveform for a single (l,m) mode */
     REAL8 M_sec_dF = M_sec * deltaF;
     COMPLEX16 It0 = I*t0; 
     #pragma omp parallel for
     for (size_t i = ind_min; i < ind_max; i++)
     {
        REAL8 Mf = i * M_sec_dF; /* geometric frequency */
        /* now we can compute the hlm */
        /* TODO: fix phase and time offsets */
        /* TODO: inclusion of reference frequency */
        // phi -= t0*(Mf-MfRef) + phi_precalc;
        // ((*htilde)->data->data)[i] = amp0 * amp * cexp(-I * phi);
        /* construct hlm at single frequency point and return */
        // (hlm->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(eta, chi1z, chi2z, ell, mm, Mf, MfRef, phi0, &z);

        ((*hlmtilde)->data->data)[i] = amp0 * IMRPhenomHMSingleModehlm(ell, mm, Mf, &z, pAmp, &amp_prefactors, pn, pPhi, &phi_prefactors, Rholm, Taulm, phi_precalc, &PhenomDQuantities, &powers_of_MfAtScale_22_amp, &downsized_powers_of_MfAtScale_22_amp);

        /* NOTE: The frequency used in the time shift term is the fourier variable of the gravitational wave frequency. i.e., Not rescaled. */
        ((*hlmtilde)->data->data)[i] *= cexp(-It0*(Mf-MfRef)*(2.0-mm) );
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
