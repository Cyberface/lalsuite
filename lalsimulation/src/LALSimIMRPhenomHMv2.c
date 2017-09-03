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

#include "LALSimIMRPhenomD.h"

/* This allows us to reuse internal IMRPhenomD functions without making those functions XLAL */
/* DIDN"T WANT TO HAVE TO DO THIS!
 * Need to rewrite parts of PhenomD code to make this nicer. */
// #include "LALSimIMRPhenomD_internals.c"


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


    // p->Inv1MinusEradRational0815 = 1.0/(1.0-EradRational0815(p->eta, p->chi1z, p->chi2z));
    // p->finspin = XLALSimIMRPhenomDFinalSpin(p->m1, p->m2, p->chi1z, p->chi2z); /* dimensionless final spin */
    // if (p->finspin > 1.0) XLAL_ERROR(XLAL_EDOM, "PhenomD fring function: final spin > 1.0 not supported\n");

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
         ind_max = (size_t) (pHM->f_max / pHM->deltaF);
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

    retcode = 0;
    retcode = IMRPhenomDPhaseFrequencySequence(
        phases,
        freqs_geom,
        ind_min, ind_max,
        pHM->m1, pHM->m2,
        pHM->chi1z, pHM->chi2z,
        1., 1.,
        extraParams
    );
    XLAL_CHECK(XLAL_SUCCESS == retcode,
        XLAL_EFUNC, "IMRPhenomDPhaseFrequencySequence failed");

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
        // freqs_amp->data[i] = IMRPhenomHMFreqDomainMap(
            // freqs_geom->data[i], ell, mm, PhenomDQuantities, AmpFlagTrue
        // );
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

    // /* coefficients of leading order PN amplitude */
    // double Blm_prefac = PhenomDQuantities->Blm_prefactor[ell][mm];
    //
    // /* ratio of frequency term in leadering order PN */
    // /* at the scaled frequency */
    // double f_frac = XLALSimIMRPhenomHMPNFrequencyScale( &powers_of_Mf_22, Mf_22, ell, mm);
    //
    // double Blm = Blm_prefac * f_frac;
    //
    // /* (m/2)^klm NOTE in the paper this is (2/m)^(-klm) i.e. inverted. */
    // double m_over_2_pow_klm = XLALSimIMRPhenomHMPNAmplitudeLeadingOrderFpow(ell, mm, mm/2.0);
    //
    // double betalm = Blm * m_over_2_pow_klm;
    //
    // double HMamp = PhenDamp * betalm;




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

    retcode = 0;
    retcode = IMRPhenomDPhaseFrequencySequence(
        phases,
        freqs_geom,
        ind_min, ind_max,
        pHM->m1, pHM->m2,
        pHM->chi1z, pHM->chi2z,
        1., 1.,
        extraParams
    );
    XLAL_CHECK(XLAL_SUCCESS == retcode,
        XLAL_EFUNC, "IMRPhenomDPhaseFrequencySequence failed");

    return XLAL_SUCCESS;
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
