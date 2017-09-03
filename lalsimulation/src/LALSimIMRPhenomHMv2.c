/*
 *  Copyright (C) 2017 Sebastian Khan
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

/**
 * \author Sebastian Khan
 *
 * \file
 *
 * \brief PhenomHM model
 *
 * Inspiral-merger and ringdown phenomenological, frequecny domain
 * waveform model for binary black holes systems.
 * Models not only the dominant (l,|m|) = (2,2) modes
 * but also some of the sub-domant modes too.
 */


#include <lal/LALSimIMR.h>
#include <lal/SphericalHarmonics.h>
#include "LALSimIMRPhenomHMv2.h"
#include "LALSimIMRPhenomInternalUtils.h"
#include "LALSimIMRPhenomUtils.h"
#include "LALSimRingdownCW.h"
#include "LALSimIMRPhenomD.h"

/*
 * Phase shift due to leading order complex amplitude
 * [L.Blancet, arXiv:1310.1528 (Sec. 9.5)]
 * "Spherical hrmonic modes for numerical relativity"
 */
/* List of phase shifts: the index is the azimuthal number m */
static const double cShift[7] = {0.0,
                                 LAL_PI_2 /* i shift */,
                                 0.0,
                                 -LAL_PI_2 /* -i shift */,
                                 LAL_PI /* 1 shift */,
                                 LAL_PI_2 /* -1 shift */,
                                 0.0};

/* This allows us to reuse internal IMRPhenomD functions without making those functions XLAL */
/* DIDN"T WANT TO HAVE TO DO THIS!
 * Need to rewrite parts of PhenomD code to make this nicer. */
// #include "LALSimIMRPhenomD_internals.c"

/**
 *
 */
int PhenomHM_init_useful_mf_powers(PhenomHMUsefulMfPowers *p, REAL8 number)
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

int PhenomHM_init_useful_powers(PhenomHMUsefulPowers *p, REAL8 number)
{
	XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");
	XLAL_CHECK(number >= 0 , XLAL_EDOM, "number must be non-negative");

	// consider changing pow(x,1/6.0) to cbrt(x) and sqrt(x) - might be faster
	double sixth = pow(number, 1.0/6.0);
	p->third = sixth * sixth;
	//p->third = cbrt(number);
	p->two_thirds = p->third * p->third;
	p->four_thirds = number * p->third;
	p->five_thirds = p->four_thirds * p->third;
	p->two = number * number;
	p->seven_thirds = p->third * p->two;
	p->eight_thirds = p->two_thirds * p->two;
    p->inv = 1./number;
	double m_sixth = 1.0/sixth;
    p->m_seven_sixths = p->inv * m_sixth;
    p->m_third = m_sixth * m_sixth;
    p->m_two_thirds = p->m_third * p->m_third;
    p->m_five_thirds = p->inv * p->m_two_thirds;

	return XLAL_SUCCESS;
}


/**
 * returns the real and imag parts of the complex ringdown frequency
 * for the (l,m) mode.
 */
int IMRPhenomHMGetRingdownFrequency(
    REAL8 *fringdown,
    REAL8 *fdamp,
    UINT4 ell,
    INT4 mm,
    REAL8 finalmass,
    REAL8 finalspin
)
{
    const REAL8 inv2Pi = 0.5/LAL_PI;
    complex double ZZ;
    ZZ = SimRingdownCW_CW07102016( SimRingdownCW_KAPPA(finalspin, ell, mm), ell, mm, 0 );
    const REAL8 Mf_RD_tmp = inv2Pi * creal(ZZ); /* GW ringdown frequency, converted from angular frequency */
    *fringdown = Mf_RD_tmp / finalmass; /* scale by predicted final mass */
    /* lm mode ringdown damping time (imaginary part of ringdown), geometric units */
    const REAL8 f_DAMP_tmp = inv2Pi * cimag(ZZ); /* this is the 1./tau in the complex QNM */
    *fdamp = f_DAMP_tmp / finalmass; /* scale by predicted final mass */

    return XLAL_SUCCESS;
}

/**
 * Precompute a bunch of PhenomHM related quantities and store them filling in a
 * PhenomHMStorage variable
 */
static int init_PhenomHM_Storage(
    PhenomHMStorage* p,
    const REAL8 m1_SI,
    const REAL8 m2_SI,
    const REAL8 chi1z,
    const REAL8 chi2z,
    REAL8Sequence *freqs,
    const REAL8 deltaF,
    const REAL8 f_ref
)
{
    int retcode;
    XLAL_CHECK(0 != p, XLAL_EFAULT, "p is NULL");

    /*TODO: Check which mass is larger and swap masses and spins if needed
     * to enforce m1>=m2
     */


    p->m1 = m1_SI / LAL_MSUN_SI;
    p->m2 = m2_SI / LAL_MSUN_SI;
    p->Mtot = p->m1 + p->m2;
    p->eta = p->m1 * p->m2 / (p->Mtot*p->Mtot);
    p->chi1z = chi1z;
    p->chi2z = chi2z;

    if (p->eta > 0.25)
        PhenomInternal_nudge(&(p->eta), 0.25, 1e-6);
    if (p->eta > 0.25 || p->eta < 0.0)
        XLAL_ERROR(XLAL_EDOM, "Unphysical eta. Must be between 0. and 0.25\n");


    /*XLAL_PRINT_INFO("before swap: m1 = %f, m2 = %f, chi1z = %f, chi2z = %f\n",
            p->m1, p->m2, p->chi1z, p->chi2z); */

    retcode = PhenomInternal_AlignedSpinEnforcePrimaryIsm1(
        &(p->m1),
        &(p->m2),
        &(p->chi1z),
        &(p->chi2z)
        );
    XLAL_CHECK(
        XLAL_SUCCESS == retcode,
        XLAL_EFUNC,
        "PhenomInternal_AlignedSpinEnforcePrimaryIsm1 failed");

    /*XLAL_PRINT_INFO("after swap: m1 = %f, m2 = %f, chi1z = %f, chi2z = %f\n",
            p->m1, p->m2, p->chi1z, p->chi2z); */

    /* sanity checks on frequencies */
    /* determine how to populate frequency sequence */
    /* if len(freqs) == 2 and deltaF > 0. then
     * f_min = freqs[0]
     * f_max = freqs[1]
     * else if len(freqs) != 2 and deltaF <= 0. then
     * user has given an arbitrary set of frequencies to evaluate the model at.
     */
    p->freqs = freqs;
    p->deltaF = deltaF;
    if ( (p->freqs->length == 2 ) && ( p->deltaF > 0. ) )
    { /* This case we use regularly spaced frequencies */
        p->freq_is_uniform = 1;
        p->f_min = p->freqs->data[0];
        p->f_max = p->freqs->data[1];

        /* If p->f_max == 0. Then we default to the ending frequency
         * for PhenomHM
         * TODO: Check the ending frequency
         * TODO: Want to implement variable ending frequency for each mode.
         */
        if ( p->f_max == 0. )
        {
            p->f_max = XLALSimPhenomUtilsMftoHz(
                        PHENOMHM_DEFAULT_MF_MAX, p->Mtot
                    );
        }
    }
    else if ( ( p->freqs->length != 2 ) && ( p->deltaF <= 0. ) )
    { /* This case we possibly use irregularly spaced frequencies */
        /* Check that the frequencies are always increasing */
        p->freq_is_uniform = 0;
        for(UINT4 i=0; i < p->freqs->length - 1; i++){
            XLAL_CHECK(
                p->freqs->data[i] - p->freqs->data[i+1] < 0.,
                XLAL_EFUNC,
                "custom frequencies must be increasing."
            );
        }

        XLAL_PRINT_INFO("Using custom frequency input.\n");
        p->f_min = p->freqs->data[0];
        p->f_max = p->freqs->data[ p->freqs->length - 1 ]; /* Last element */
    }
    else
    { /* Throw an informative error. */
        XLAL_PRINT_ERROR("Input sequence of frequencies and deltaF is not \
compatible.\nSpecify a f_min and f_max by using a REAL8Sequence of length = 2 \
along with a deltaF > 0.\
\nIf you want to supply an arbitrary list of frequencies to evaluate the with \
then supply those frequencies using a REAL8Sequence and also set deltaF <= 0.");
    }

    /* Fix default behaviour for f_ref */
    /* If f_ref = 0. then set f_ref = f_min */
    p->f_ref = f_ref;
    if (p->f_ref == 0.)
    {
        p->f_ref = p->f_min;
    }


    p->finmass = IMRPhenomDFinalMass(p->m1, p->m2, p->chi1z, p->chi2z);
    p->finspin = XLALSimIMRPhenomDFinalSpin(p->m1, p->m2, p->chi1z, p->chi2z); /* dimensionless final spin */
    if (p->finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");


    /* populate the ringdown frequency array */
    /* If you want to model a new mode then you have to add it here. */
    /* TODO: Fix this to run over a loop of allowed modes. */
    /* (l,m) = (2,2) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[2][2],
        &p->PhenomHMfdamp[2][2],
        2, 2,
        p->finmass, p->finspin
    );

    /* (l,m) = (2,1) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[2][1],
        &p->PhenomHMfdamp[2][1],
        2, 1,
        p->finmass, p->finspin
    );

    /* (l,m) = (3,3) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[3][3],
        &p->PhenomHMfdamp[3][3],
        3, 3,
        p->finmass, p->finspin
    );

    /* (l,m) = (3,2) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[3][2],
        &p->PhenomHMfdamp[3][2],
        3, 2,
        p->finmass, p->finspin
    );

    /* (l,m) = (4,4) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[4][4],
        &p->PhenomHMfdamp[4][4],
        4, 4,
        p->finmass, p->finspin
    );

    /* (l,m) = (4,3) */
    IMRPhenomHMGetRingdownFrequency(
        &p->PhenomHMfring[4][3],
        &p->PhenomHMfdamp[4][3],
        4, 3,
        p->finmass, p->finspin
    );

    p->Mf_RD_22 = p->PhenomHMfring[2][2];
    p->Mf_DM_22 = p->PhenomHMfdamp[2][2];

    /* (l,m) = (2,2) */
    int ell, mm;
    ell = 2;
    mm = 2;
    p->Rholm[ell][mm] = p->Mf_RD_22/p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm]/p->Mf_DM_22;
    /* (l,m) = (2,1) */
    ell = 2;
    mm = 1;
    p->Rholm[ell][mm] = p->Mf_RD_22/p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm]/p->Mf_DM_22;
    /* (l,m) = (3,3) */
    ell = 3;
    mm = 3;
    p->Rholm[ell][mm] = p->Mf_RD_22/p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm]/p->Mf_DM_22;
    /* (l,m) = (3,2) */
    ell = 3;
    mm = 2;
    p->Rholm[ell][mm] = p->Mf_RD_22/p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm]/p->Mf_DM_22;
    /* (l,m) = (4,4) */
    ell = 4;
    mm = 4;
    p->Rholm[ell][mm] = p->Mf_RD_22/p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm]/p->Mf_DM_22;
    /* (l,m) = (4,3) */
    ell = 4;
    mm = 3;
    p->Rholm[ell][mm] = p->Mf_RD_22/p->PhenomHMfring[ell][mm];
    p->Taulm[ell][mm] = p->PhenomHMfdamp[ell][mm]/p->Mf_DM_22;


    /* populate Blm_prefactor array */
    /* see equ.(4) from 1708.00404 */
    /* A bunch of useful powers used in XLALSimIMRPhenomHMPNAmplitudeLeadingOrder */
    /* pow_Mf_wf_prefactor are the coefficients from the leading order amplitude terms for each l,m mode we consider.  */
    /* If A_lm(f) = alpha_lm * f^klm then pow_Mf_wf_prefactor = alpha_lm */
    /* note The factors of pi's normally multiplying with the frequency are factored into these terms. */
    /* note delta == sqrt(1. - 4.*eta) */
    /* note delta2 == 1. - 4.*eta */
    // REAL8 sqrteta = sqrt(p->eta);
    REAL8 Seta = sqrt( 1.0 - 4.0 * p->eta );
    REAL8 delta = Seta;
    REAL8 delta2 = 1.0 - 4.0 * p->eta ;
    p->Blm_prefactor[2][2] = 1.0;
    p->Blm_prefactor[2][1] = delta * pow(LAL_PI, 1.0/3.0) /  3.0;
    p->Blm_prefactor[3][3] = (3.0/4.0) * sqrt(15.0/14.0) * pow(LAL_PI, 1.0/3.0) * delta;
    p->Blm_prefactor[3][2] = sqrt(5.0/63.0) * pow(LAL_PI, 2.0/3.0) * (delta2 + p->eta);
    p->Blm_prefactor[4][4] = sqrt(320.0/567.0) * pow(LAL_PI, 2.0/3.0) * (delta2 + p->eta);
    p->Blm_prefactor[4][3] = sqrt(81.0/1120.0) * LAL_PI * (delta2 + 2*p->eta) * delta;

    return XLAL_SUCCESS;
};

/**
 * mathematica function postfRDflm
 * domain mapping function - ringdown
 */
double IMRPhenomHMTrd(
    REAL8 Mf,
    REAL8 Mf_RD_22,
    REAL8 Mf_RD_lm,
    const INT4 AmpFlag,
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage* pHM
)
{
    double ans = 0.0;
    if ( AmpFlag==1 ) {
        /* For amplitude */
        ans = Mf-Mf_RD_lm+Mf_RD_22; /*Used for the Amplitude as an approx fix for post merger powerlaw slope */
    } else {
        /* For phase */
        REAL8 Rholm = pHM->Rholm[ell][mm];
        ans = Rholm * Mf;          /* Used for the Phase */
    }

    return ans;
}

/**
 * mathematica function Ti
 * domain mapping function - inspiral
 */
double IMRPhenomHMTi(REAL8 Mf, const INT4 mm){
    return 2.0 * Mf / mm;
}

/**
 * helper function for IMRPhenomHMFreqDomainMap
 */
int IMRPhenomHMSlopeAmAndBm(
    double *Am,
    double *Bm,
    const INT4 mm,
    REAL8 fi,
    REAL8 fr,
    REAL8 Mf_RD_22,
    REAL8 Mf_RD_lm,
    const INT4 AmpFlag,
    const INT4 ell,
    PhenomHMStorage* pHM
)
{
    REAL8 Trd = IMRPhenomHMTrd(fr, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, mm, pHM);
    REAL8 Ti = IMRPhenomHMTi(fi, mm);

    //Am = ( Trd[fr]-Ti[fi] )/( fr - fi );
    *Am = ( Trd-Ti )/( fr - fi );

    //Bm = Ti[fi] - fi*Am;
    *Bm = Ti - fi*(*Am);

    return XLAL_SUCCESS;
}

/**
 * helper function for IMRPhenomHMFreqDomainMap
 */
int IMRPhenomHMMapParams(
    REAL8 *a,
    REAL8 *b,
    REAL8 flm,
    REAL8 fi,
    REAL8 fr,
    REAL8 Ai,
    REAL8 Bi,
    REAL8 Am,
    REAL8 Bm,
    REAL8 Ar,
    REAL8 Br
)
{
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


/**
 * helper function for IMRPhenomHMFreqDomainMap
 */
int IMRPhenomHMFreqDomainMapParams(
    REAL8 *a,/**< [Out]  */
    REAL8 *b,/**< [Out]  */
    REAL8 *fi,/**< [Out]  */
    REAL8 *fr,/**< [Out]  */
    REAL8 *f1,/**< [Out]  */
    const REAL8 flm, /**< input waveform frequency */
    const INT4 ell, /**< spherical harmonics ell mode */
    const INT4 mm, /**< spherical harmonics m mode */
    PhenomHMStorage *pHM, /**< Stores quantities in order to calculate them only once */
    const int AmpFlag /**< is ==1 then computes for amplitude, if ==0 then computes for phase */
)
{

    /*check output points are NULL*/
    XLAL_CHECK(a != NULL, XLAL_EFAULT);
    XLAL_CHECK(b != NULL, XLAL_EFAULT);
    XLAL_CHECK(fi != NULL, XLAL_EFAULT);
    XLAL_CHECK(fr != NULL, XLAL_EFAULT);
    XLAL_CHECK(f1 != NULL, XLAL_EFAULT);

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

    REAL8 Mf_RD_22 = pHM->Mf_RD_22;
    REAL8 Mf_RD_lm = pHM->PhenomHMfring[ell][mm];

    // Define a ratio of QNM frequencies to be used for scaling various quantities
    REAL8 Rholm = pHM->Rholm[ell][mm];

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
    IMRPhenomHMSlopeAmAndBm(&Am, &Bm, mm, *fi, *fr, Mf_RD_22, Mf_RD_lm, AmpFlag, ell, pHM);

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
    int ret = IMRPhenomHMMapParams(a, b, flm, *fi, *fr, Ai, Bi, Am, Bm, Ar, Br);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - IMRPhenomHMMapParams failed in IMRPhenomHMFreqDomainMapParams (1)\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    // *f2lm = 0.5 * ( *fr + *fi );

    return XLAL_SUCCESS;
}

/**
 * IMRPhenomHMFreqDomainMap
 * Input waveform frequency in Geometric units (Mflm)
 * and computes what frequency this corresponds
 * to scaled to the 22 mode.
 */
double IMRPhenomHMFreqDomainMap(
    REAL8 Mflm,
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage* pHM,
    const int AmpFlag
)
{

    /* Mflm here has the same meaning as Mf_wf in XLALSimIMRPhenomHMFreqDomainMapHM (old deleted function). */
    REAL8 a = 0.;
    REAL8 b = 0.;
    /* Following variables not used in this funciton but are returned in IMRPhenomHMFreqDomainMapParams */
    REAL8 fi = 0.;
    REAL8 fr = 0.;
    REAL8 f1 = 0.;
    int ret = IMRPhenomHMFreqDomainMapParams(&a, &b, &fi, &fr, &f1, Mflm, ell, mm, pHM, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - IMRPhenomHMFreqDomainMapParams failed in IMRPhenomHMFreqDomainMap\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    REAL8 Mf22 = a * Mflm + b;
    return Mf22;
}

/*
 * For a given frequency and ell and m spherical harmonic mode
 * return the frequency scaled to give the leading order PN
 * amplitude for the given ell and m modes.
 * Calcuated from mathematica function: FrequencyPower[f, {ell, m}] / FrequencyPower[f, {2, 2}]
 * FrequencyPower function just returns the leading order PN term in the amplitude.
 */
 /**
   * If A_lm(f) = alpha_lm * f^klm then this function
   * returns f^klm / f^k22
   */
double IMRPhenomHMPNFrequencyScale(
    PhenomHMUsefulPowers *p,
    REAL8 Mf,
    INT4 ell,
    INT4 mm
)
{

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

/**
  * If A_lm(f) = alpha_lm * f^klm then this function
  * returns f^(klm) where klm is the appropriate exponent for th l,m mode
  */
double IMRPhenomHMPNAmplitudeLeadingOrderFpow(
    INT4 ell,
    INT4 mm,
    REAL8 Mf
)
{
    /* Initialise answer */
    REAL8 ans = 0.0;

    PhenomHMUsefulMfPowers powers_of_Mf;
    int errcode = XLAL_SUCCESS;
    errcode = PhenomHM_init_useful_mf_powers(&powers_of_Mf, Mf);
    XLAL_CHECK(errcode == XLAL_SUCCESS, errcode, "PhenomHM_init_useful_mf_powers failed for Mf");

    //FP: some of these can be computed directly here rather than for each mm and ll
    if ( ell==2 ) {
        if ( mm==2 ) {
          ans = powers_of_Mf.m_seven_sixths;
        } else { //mm==1
          ans = powers_of_Mf.m_five_sixths;
        }
    } else if ( ell==3 ) {
        if ( mm==3 ) {
          ans = powers_of_Mf.m_five_sixths;
        }
        else { //mm==2
          ans = powers_of_Mf.m_sqrt;
        }
    } else if ( ell==4 ) {
        if ( mm==4 ) {
          ans = powers_of_Mf.m_sqrt;
        }
        else { //mm==3
          ans = powers_of_Mf.m_sixth;
        }
    } else if ( ell==5 ) {
        if ( mm==5 ) {
          ans = powers_of_Mf.m_sixth;
        }
        else { //mm==4
          ans = powers_of_Mf.sixth;
        }
    } else if ( ell==6 ) {
        if ( mm==6 ) {
          ans = powers_of_Mf.sixth;
        }
        else { //mm==5
          ans = powers_of_Mf.sqrt;
        }
    } else {
        XLALPrintError("XLAL Error - requested ell = %i and m = %i mode not available, check documentation for available modes\n", ell, mm);
        XLAL_ERROR(XLAL_EDOM);
    }

    return ans;
}

int IMRPhenomHMPhasePreComp(
    HMPhasePreComp *q, /**< [out] */
    const INT4 ell,
    const INT4 mm,
    PhenomHMStorage *pHM,
    UNUSED LALDict *extraParams
)
{
    REAL8 ai = 0.0;
    REAL8 bi = 0.0;
    REAL8 am = 0.0;
    REAL8 bm = 0.0;
    REAL8 ar = 0.0;
    REAL8 br = 0.0;
    REAL8 fi = 0.0;
    REAL8 f1 = 0.0;
    REAL8 fr = 0.0;

    const INT4 AmpFlag = 0;

    /* NOTE: As long as Mfshit + f2lm isn't >= fr then the value of the shift is arbitrary. */
    const REAL8 Mfshift = 0.0001;

    int ret = IMRPhenomHMFreqDomainMapParams(&ai, &bi, &fi, &fr, &f1, Mfshift, ell, mm, pHM, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - IMRPhenomHMFreqDomainMapParams failed in IMRPhenomHMFreqDomainMapParams - inspiral\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    q->ai = ai;
    q->bi = bi;

    ret = IMRPhenomHMFreqDomainMapParams(&am, &bm, &fi, &fr, &f1, fi+Mfshift, ell, mm, pHM, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - IMRPhenomHMFreqDomainMapParams failed in IMRPhenomHMFreqDomainMapParams - intermediate\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    q->am = am;
    q->bm = bm;

    ret = IMRPhenomHMFreqDomainMapParams(&ar, &br, &fi, &fr, &f1, fr+Mfshift, ell, mm, pHM, AmpFlag);
    if (ret != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - IMRPhenomHMFreqDomainMapParams failed in IMRPhenomHMFreqDomainMapParams - merger-ringdown\n");
        XLAL_ERROR(XLAL_EDOM);
    }

    q->ar = ar;
    q->br = br;

    q->fi = fi;
    q->fr = fr;

    // printf("q->ai = %f\n", q->ai);
    // printf("q->bi = %f\n", q->bi);
    // printf("q->am = %f\n", q->am);
    // printf("q->bm = %f\n", q->bm);
    // printf("q->ar = %f\n", q->ar);
    // printf("q->br = %f\n", q->br);

    double Rholm = pHM->Rholm[ell][mm];
    double Taulm = pHM->Taulm[ell][mm];

    REAL8 PhDBMf = am*fi + bm;
    q->PhDBconst = IMRPhenomDPhase_OneFrequency(PhDBMf, pHM, Rholm, Taulm, extraParams)/am;

    REAL8 PhDCMf = ar*fr + br;
    q->PhDCconst = IMRPhenomDPhase_OneFrequency(PhDCMf, pHM, Rholm, Taulm, extraParams)/ar;

    REAL8 PhDBAMf = ai*fi + bi;
    q->PhDBAterm = IMRPhenomDPhase_OneFrequency(PhDBAMf, pHM, Rholm, Taulm, extraParams)/ai;
    return XLAL_SUCCESS;

}


/**
 * h = h+ - i*hx = Sum hlm * Ylm - FIXME: check sign
 * Returns h+ and hx in the frequency domain
 * evaluated at a user defined set of frequencies.
 * Loops over IMRPhenomHMOneFrequency
 *
 * This function can be called in the usual sense
 * where you supply a f_min, f_max and deltaF.
 * This is the case when deltaF > 0.
 * If f_max = 0. then the default ending frequnecy is used.
 * or you can also supply a custom set of discrete
 * frequency points with which to evaluate the waveform.
 * To do this you must call this function with
 * deltaF <= 0.
 *
 */
UNUSED int XLALSimIMRPhenomHM(
    UNUSED COMPLEX16FrequencySeries **hptilde,        /**< [out] Frequency-domain waveform h+ */
    UNUSED COMPLEX16FrequencySeries **hctilde,        /**< [out] Frequency-domain waveform hx */
    UNUSED REAL8Sequence *freqs,                   /**< Frequency points at which to evaluate the waveform (Hz) */
    UNUSED REAL8 m1_SI,                                /**< mass of companion 1 (kg) */
    UNUSED REAL8 m2_SI,                                /**< mass of companion 2 (kg) */
    UNUSED REAL8 chi1z,                                  /**< z-component of the dimensionless spin of object 1 w.r.t. Lhat = (0,0,1) */
    UNUSED REAL8 chi2z,                                  /**< z-component of the dimensionless spin of object 2 w.r.t. Lhat = (0,0,1) */
    UNUSED const REAL8 distance,                       /**< distance of source (m) */
    UNUSED const REAL8 inclination,                    /**< inclination of source (rad) */
    UNUSED const REAL8 phiRef,                         /**< reference orbital phase (rad) */
    UNUSED const REAL8 deltaF,                         /**< Sampling frequency (Hz). To use arbitrary frequency points set deltaF <= 0. */
    UNUSED REAL8 f_ref,                          /**< Reference frequency */
    UNUSED LALDict *extraParams /**<linked list containing the extra testing GR parameters */
)
{
    /* define and init return code for this function */
    int retcode;
    /* sanity checks on input parameters: check pointers, etc. */

    /* Check inputs for sanity */
    XLAL_CHECK(NULL != hptilde, XLAL_EFAULT);
    XLAL_CHECK(NULL != hctilde, XLAL_EFAULT);
    XLAL_CHECK(*hptilde == NULL, XLAL_EFAULT);
    XLAL_CHECK(*hctilde == NULL, XLAL_EFAULT);
    XLAL_CHECK(m1_SI > 0, XLAL_EDOM, "m1 must be positive.\n");
    XLAL_CHECK(m2_SI > 0, XLAL_EDOM, "m2 must be positive.\n");
    XLAL_CHECK(fabs(chi1z) <= 1.0, XLAL_EDOM, "Aligned spin chi1z=%g \
must be <= 1 in magnitude!\n", chi1z);
    XLAL_CHECK(fabs(chi2z) <= 1.0, XLAL_EDOM, "Aligned spin chi2z=%g \
must be <= 1 in magnitude!\n", chi2z);
    XLAL_CHECK(distance > 0, XLAL_EDOM, "distance must be positive.\n");
    XLAL_CHECK(f_ref >= 0, XLAL_EDOM, "Reference frequency must be \
positive.\n"); /* FIXME: check this one */

    /* setup ModeArray */
    if (extraParams==NULL)
      extraParams=XLALCreateDict();
    LALValue* ModeArray = XLALSimInspiralWaveformParamsLookupModeArray(extraParams);
    if ( ModeArray == NULL )
    {   /* Default behaviour */
        /* TODO: Move this into a function */
        XLAL_PRINT_INFO("Using default modes for PhenomHM.\n");
        ModeArray = XLALSimInspiralCreateModeArray();
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, 1);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 2, -1);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 3);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -3);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 3, 2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 3, -2);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 4);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -4);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 4, 3);
        XLALSimInspiralModeArrayActivateMode(ModeArray, 4, -3);
        // XLALSimInspiralModeArrayPrintModes(ModeArray);
        /* Don't forget to insert ModeArray back into extraParams. */
        XLALSimInspiralWaveformParamsInsertModeArray(extraParams, ModeArray);
    }
    else
    {
        XLAL_PRINT_INFO("Using custom modes for PhenomHM.\n");
        /* TODO: loop over mode array and throw an error if
         there is a mode that is not included in the model */
        /* Have a variable that contains the list of modes currently
        implemented and check against this. */
    }

     /* setup PhenomHM model storage struct / structs */
     /* Compute quantities/parameters related to PhenomD only once and store them */
     PhenomHMStorage PhenomHMQuantities;
     retcode = 0;
     retcode = init_PhenomHM_Storage(
                                    &PhenomHMQuantities,
                                    m1_SI,
                                    m2_SI,
                                    chi1z,
                                    chi2z,
                                    freqs,
                                    deltaF,
                                    f_ref
                                );
     XLAL_CHECK(XLAL_SUCCESS == retcode, XLAL_EFUNC, "init_PhenomHM_Storage \
failed");



     /* main: evaluate model at given frequencies */
     retcode = 0;
     retcode = IMRPhenomHMCore(
                            &PhenomHMQuantities,
                            extraParams
                        );
     XLAL_CHECK(retcode == XLAL_SUCCESS,
         XLAL_EFUNC, "IMRPhenomHMCore failed in XLALSimIMRPhenomHM.");

     /* cleanup */
     /* XLALDestroy and XLALFree any pointers. */

    return XLAL_SUCCESS;
}

/**
 * internal function that returns h+ and hx.
 * Inside this function the my bulk of the work is done
 * like the loop over frequencies.s
 */
int IMRPhenomHMCore(
    UNUSED PhenomHMStorage *PhenomHMQuantities,
    UNUSED LALDict *extraParams
)
{
    int retcode;
    // printf("PhenomHMQuantities.m1 = %f\n", PhenomHMQuantities->m1);
    /* precompute all frequency independent terms here. */

    // for i, f in enumerate(PhenomHMQuantities.freqs):
    //     hptilde[i], hctilde[i] = IMRPhenomHMOneFrequency(f, PhenomHMQuantities)

    /* evaluate all hlm modes */
    SphHarmFrequencySeries **hlms=XLALMalloc(sizeof(SphHarmFrequencySeries));
    *hlms=NULL;
    retcode = 0;
    retcode = IMRPhenomHMEvaluatehlmModes(hlms,
        PhenomHMQuantities, extraParams);
    XLAL_CHECK(XLAL_SUCCESS == retcode,
        XLAL_EFUNC, "IMRPhenomHMEvaluatehlmModes failed");

    return XLAL_SUCCESS;
}


/**
 * Function to compute the hlm modes.
 * Note this is not static so that IMRPhenomPv3HM
 * can also use this function
 */
int IMRPhenomHMEvaluatehlmModes(
    UNUSED SphHarmFrequencySeries **hlms,
    UNUSED PhenomHMStorage *pHM,
    UNUSED LALDict *extraParams
)
{
    int retcode;

    /* setup frequency sequency */

    REAL8Sequence *phases = NULL;
    REAL8Sequence *freqs = NULL; /* freqs is in Hz */
    REAL8Sequence *freqs_geom = NULL; /* freqs is in geometric units */

    size_t npts;
    /* these used to compute the waveform where we want non-zero values. */
    size_t ind_min;
    size_t ind_max;

    /* Two possibilities */
    if (pHM->freq_is_uniform==1)
    { /* 1. uniformly spaced */
        printf("freq_is_uniform = True\n");

        /* we only need to evaluate the phase from
         * f_min to f_max with a spacing of deltaF
         */
         npts = PhenomInternal_NextPow2(pHM->f_max / pHM->deltaF) + 1;
         ind_min = (size_t) (pHM->f_min / pHM->deltaF);
         ind_max = (size_t) (pHM->f_max / pHM->deltaF) + 1; /* TODO: double check this ''+1' I added at the end... I did it so that the loop made sense but I think this makes more sense. */
         XLAL_CHECK ( (ind_max<=npts) && (ind_min<=ind_max), XLAL_EDOM, "minimum freq index %zu and maximum freq index %zu do not fulfill 0<=ind_min<=ind_max<=freqs->length=%zu.", ind_min, ind_max, npts);

         freqs = XLALCreateREAL8Sequence( npts );
         phases = XLALCreateREAL8Sequence( npts );

        for( size_t i=0; i < npts; i++ )
        { /* populate the frequency unitformly from zero - this is the standard
             convention we use when generating waveforms in LAL. */
             freqs->data[i] = i * pHM->deltaF; /* This is in Hz */
             phases->data[i] = 0; /* initalise all phases to zero. */
        }

    }
    else if (pHM->freq_is_uniform==0)
    { /* 2. arbitrarily space */
        printf("freq_is_uniform = False\n");
        freqs = pHM->freqs; /* This is in Hz */
        npts = freqs->length;
        ind_min = 0;
        ind_max = npts;
        phases = XLALCreateREAL8Sequence( freqs->length );
        for( size_t i=0; i <npts; i++ )
         {
             phases->data[i] = 0; /* initalise all phases to zero. */
         }
    }
    else
    {
        XLAL_ERROR(XLAL_EDOM, "freq_is_uniform = %i and should be either 0 or 1.", pHM->freq_is_uniform);
    }

    /* PhenomD functions take geometric frequencies */
    freqs_geom = XLALCreateREAL8Sequence( npts );
    for( size_t i=0; i < npts; i++ )
     {
         freqs_geom->data[i] = XLALSimPhenomUtilsHztoMf(freqs->data[i], pHM->Mtot); /* initalise all phases to zero. */
     }

    // retcode = 0;
    // retcode = IMRPhenomDPhaseFrequencySequence(
    //     phases,
    //     freqs_geom,
    //     ind_min, ind_max,
    //     pHM->m1, pHM->m2,
    //     pHM->chi1z, pHM->chi2z,
    //     1., 1.,
    //     extraParams
    // );
    // XLAL_CHECK(XLAL_SUCCESS == retcode,
    //     XLAL_EFUNC, "IMRPhenomDPhaseFrequencySequence failed");

    /* TODO: modify ind_min and ind_max
     * based on the (l,m) number to try and
     */

    retcode = 0;
    retcode = IMRPhenomHMPhase(
        phases,
        freqs_geom,
        pHM,
        2, 2,
        ind_min,
        ind_max,
        extraParams
    );
    XLAL_CHECK(XLAL_SUCCESS == retcode,
        XLAL_EFUNC, "IMRPhenomHMPhase failed");


    for(UINT4 i=0; i < phases->length; i++)
    {
        printf("freqs->data[%i] = %f, freqs_geom->data[%i] = %f, phases->data[%i] = %.8f\n", i, freqs->data[i], i, freqs_geom->data[i], i, phases->data[i]);
    }

    /* cleanup */
    XLALDestroyREAL8Sequence(freqs_geom);

    if (pHM->freq_is_uniform==1)
    { /* 1. uniformly spaced */
        XLALDestroyREAL8Sequence(phases);
        XLALDestroyREAL8Sequence(freqs);
    }
    else if (pHM->freq_is_uniform==0)
    { /* 2. arbitrarily space */
        XLALDestroyREAL8Sequence(phases);
    }
    else
    {
        XLAL_ERROR(XLAL_EDOM, "freq_is_uniform = %i and should be either 0 or 1.", pHM->freq_is_uniform);
    }

    return XLAL_SUCCESS;
}



/**
 * returns IMRPhenomHM amplitude evaluated at a set of input frequencies
 * for the l,m mode
 */
int IMRPhenomHMAmplitude(
    UNUSED REAL8Sequence *amps, /**< [out] */
    UNUSED REAL8Sequence *freqs_geom,
    UNUSED PhenomHMStorage *pHM,
    UNUSED UINT4 ell,
    UNUSED INT4 mm,
    UNUSED size_t ind_min,
    UNUSED size_t ind_max,
    UNUSED LALDict *extraParams
)
{
    int retcode;

    /* scale input frequencies according to PhenomHM model */
    REAL8Sequence *freqs_amp = XLALCreateREAL8Sequence( freqs_geom->length );
    for(UINT4 i=0; i<freqs_amp->length; i++)
    {
        freqs_amp->data[i] = IMRPhenomHMFreqDomainMap(
            freqs_geom->data[i], ell, mm, pHM, AmpFlagTrue
        );
    }

    retcode = 0;
    retcode = IMRPhenomDAmpFrequencySequence(
        amps,
        freqs_amp,
        ind_min, ind_max,
        pHM->m1, pHM->m2,
        pHM->chi1z, pHM->chi2z
    );
    XLAL_CHECK(XLAL_SUCCESS == retcode,
        XLAL_EFUNC, "IMRPhenomDAmpFrequencySequence failed");

    /* coefficients of leading order PN amplitude */
    double Blm_prefac = pHM->Blm_prefactor[ell][mm];

    /* (m/2)^klm NOTE in the paper this is (2/m)^(-klm) i.e. inverted. */
    double m_over_2_pow_klm = IMRPhenomHMPNAmplitudeLeadingOrderFpow(ell, mm, mm/2.0);

    int status_in_for = XLAL_SUCCESS;
    for(UINT4 i=0; i<freqs_amp->length; i++)
    {
        PhenomHMUsefulPowers powers_of_freq_amp;
        status_in_for = PhenomHM_init_useful_powers(
            &powers_of_freq_amp, freqs_amp->data[i]
        );
        if (XLAL_SUCCESS != status_in_for)
        {
          XLALPrintError("PhenomHM_init_useful_powers failed for Mf, status_in_for=%d", status_in_for);
          retcode = status_in_for;
        }
        /* ratio of frequency term in leadering order PN */
        /* at the scaled frequency */
        double f_frac = IMRPhenomHMPNFrequencyScale(
            &powers_of_freq_amp, freqs_amp->data[i], ell, mm
        );
        double Blm = Blm_prefac * f_frac;
        double betalm = Blm * m_over_2_pow_klm;
        amps->data[i] *= betalm;
    }

    /* cleanup */
    // XLALDestroyREAL8Sequence(freqs_amp);


    return XLAL_SUCCESS;
}
/**
 * returns IMRPhenomHM phase evaluated at a set of input frequencies
 * for the l,m mode
 */

int IMRPhenomHMPhase(
    UNUSED REAL8Sequence *phases, /**< [out] */
    UNUSED REAL8Sequence *freqs_geom,
    UNUSED PhenomHMStorage *pHM,
    UNUSED UINT4 ell,
    UNUSED INT4 mm,
    UNUSED size_t ind_min,
    UNUSED size_t ind_max,
    UNUSED LALDict *extraParams
)
{
    int retcode;

    // retcode = 0;
    // retcode = IMRPhenomDPhaseFrequencySequence(
    //     phases,
    //     freqs_geom,
    //     ind_min, ind_max,
    //     pHM->m1, pHM->m2,
    //     pHM->chi1z, pHM->chi2z,
    //     1., 1.,
    //     extraParams
    // );
    // XLAL_CHECK(XLAL_SUCCESS == retcode,
    //     XLAL_EFUNC, "IMRPhenomDPhaseFrequencySequence failed");

    // HMPhasePreComp q = NULL;
    HMPhasePreComp q;
    retcode = 0;
    retcode = IMRPhenomHMPhasePreComp(&q, ell, mm, pHM, extraParams);
    if (retcode != XLAL_SUCCESS){
        XLALPrintError("XLAL Error - IMRPhenomHMPhasePreComp failed\n");
        XLAL_ERROR(XLAL_EDOM);
    }
    REAL8 Rholm = pHM->Rholm[ell][mm];
    REAL8 Taulm = pHM->Taulm[ell][mm];

    REAL8 Mf_wf = 0.0;
    REAL8 Mf = 0.0;
    REAL8 Mfr = 0.0;
    REAL8 tmpphaseC = 0.0;
    for(UINT4 i=ind_min; i<ind_max; i++)
    {
        phases->data[i] += cShift[mm];
        Mf_wf = freqs_geom->data[i];
        // This if ladder is in the mathematica function HMPhase. PhenomHMDev.nb
        if ( !(Mf_wf > q.fi) )
        { /* in mathematica -> IMRPhenDPhaseA */
            Mf = q.ai * Mf_wf + q.bi;
            phases->data[i] += IMRPhenomDPhase_OneFrequency(Mf, pHM, Rholm, Taulm, extraParams) / q.ai;
        }
        else if ( !(Mf_wf > q.fr) )
        { /* in mathematica -> IMRPhenDPhaseB */
            Mf = q.am*Mf_wf + q.bm;
            phases->data[i] += IMRPhenomDPhase_OneFrequency(Mf, pHM, Rholm, Taulm, extraParams) / q.am - q.PhDBconst + q.PhDBAterm;
        }
        else if ( !(Mf_wf > q.fr) )
        { /* in mathematica -> IMRPhenDPhaseC */
            Mfr = q.am*q.fr + q.bm;
            tmpphaseC = IMRPhenomDPhase_OneFrequency(Mfr, pHM, Rholm, Taulm, extraParams) / q.am - q.PhDBconst + q.PhDBAterm;
            Mf = q.ar*Mf_wf + q.br;
            phases->data[i] += IMRPhenomDPhase_OneFrequency(Mf, pHM, Rholm, Taulm, extraParams) / q.ar - q.PhDCconst + tmpphaseC;
        }
        else
        {
            XLALPrintError("XLAL_ERROR - should not get here - in function IMRPhenomHMPhase");
            XLAL_ERROR(XLAL_EDOM);
        }
    }

    return XLAL_SUCCESS;
}

/**
 * helper function that returns the PhenomD phase
 * generalised with Rholm and Taulm.
 * A wrapper of IMRPhenomDPhaseFrequencySequence
 * for convenience.
 */
double IMRPhenomDPhase_OneFrequency(
    REAL8 Mf,
    PhenomHMStorage *pHM,
    REAL8 Rholm,
    REAL8 Taulm,
    LALDict *extraParams
)
{
    REAL8Sequence *one_freq = XLALCreateREAL8Sequence(1);
    REAL8Sequence *phase = XLALCreateREAL8Sequence(1);
    one_freq->data[0] = Mf;

    size_t ind_min = 0;
    size_t ind_max = 1;

    int retcode = 0;
    retcode = IMRPhenomDPhaseFrequencySequence(
        phase,
        one_freq,
        ind_min, ind_max,
        pHM->m1, pHM->m2,
        pHM->chi1z, pHM->chi2z,
        Rholm,
        Taulm,
        extraParams
    );
    XLAL_CHECK(XLAL_SUCCESS == retcode,
        XLAL_EFUNC, "IMRPhenomDPhaseFrequencySequence failed");

    REAL8 out = phase->data[0];

    XLALDestroyREAL8Sequence(phase);
    XLALDestroyREAL8Sequence(one_freq);

    return out;
}



/**
 * returns h+ and hx in the frequency domain at one frequency. Optimised inputs
 */
// int IMRPhenomHMOneFrequency(
//     REAL8 fHz /** < Single frequency at which to evaluate the waveform (Hz) */
// )

/**
 * Returns hlm at a set of given frequencies
 * * Loops over IMRPhenomHMSingleModehlmOneFrequency
 */
// int XLALSimIMRPhenomHMSingleModehlmFrequencySequence

/**
 * returns hlm at one frequency. Optimised inputs
 */
// int IMRPhenomHMSingleModehlmOneFrequency
