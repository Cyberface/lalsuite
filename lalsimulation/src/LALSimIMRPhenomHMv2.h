#ifndef _LALSIM_IMR_PHENOMHMv2_H
#define _LALSIM_IMR_PHENOMHMv2_H


#include <lal/LALDatatypes.h>
#include <lal/Sequence.h>
#include <lal/LALDict.h>
#include <lal/LALConstants.h>
#include <lal/XLALError.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

/**
 * Dimensionless frequency (Mf) at which define the end of the waveform
 */
#define PHENOMHM_DEFAULT_MF_MAX 0.5


/**
 * Maximum number of (l,m) mode paris PhenomHM models.
 * Only count positive 'm' modes.
 * Used to set default mode array
 */
#define NMODES_MAX 6

/**
 * Highest ell multipole PhenomHM models + 1.
 * Used to set array sizes
 */
#define L_MAX_PLUS_1 5

/* Activates amplitude part of the model */
#define AmpFlagTrue 1
#define AmpFlagFalse 0

/**
 * Structure storing pre-determined quantities
 * complying to the conventions of the PhenomHM model.
 * convensions such as m1>=m2
 */
typedef struct tagPhenomHMStorage
{
    REAL8 m1; /**< mass of larger body in solar masses */
    REAL8 m2; /**< mass of lighter body in solar masses */
    REAL8 Mtot;  /**< total mass in solar masses */
    REAL8 eta;  /**< symmetric mass-ratio */
    REAL8 chi1z; /**< dimensionless aligned component spin of larger body */
    REAL8 chi2z; /**< dimensionless aligned component spin of lighter body */
    REAL8Sequence *freqs;
    REAL8 deltaF;
    REAL8 f_min;
    REAL8 f_max;
    REAL8 f_ref;
    UINT4 freq_is_uniform; /**< If = 1 then assume uniform spaced, If = 0 then assume arbitrarily spaced. */
    // REAL8 Inv1MinusEradRational0815;
    // REAL8 finspin;
    // REAL8 Mf_RD_22;
    // REAL8 Mf_DM_22;
    // REAL8 PhenomHMfring[L_MAX_PLUS_1][L_MAX_PLUS_1];
    // REAL8 PhenomHMfdamp[L_MAX_PLUS_1][L_MAX_PLUS_1];
    // REAL8 Rholm[L_MAX_PLUS_1][L_MAX_PLUS_1];
    // REAL8 Taulm[L_MAX_PLUS_1][L_MAX_PLUS_1];
    REAL8 Blm_prefactor[L_MAX_PLUS_1][L_MAX_PLUS_1];
}
PhenomHMStorage;

static int init_PhenomHM_Storage(
    PhenomHMStorage* p, /**< [out] PhenomHMStorage struct */
    const REAL8 m1_SI,
    const REAL8 m2_SI,
    const REAL8 chi1z,
    const REAL8 chi2z,
    REAL8Sequence *freqs,
    const REAL8 deltaF,
    const REAL8 f_ref
);



int IMRPhenomHMCore(
    PhenomHMStorage *PhenomHMQuantities,
    LALDict *extraParams
);


int IMRPhenomHMEvaluatehlmModes(
    SphHarmFrequencySeries **hlms,
    PhenomHMStorage *PhenomHMQuantities,
    LALDict *extraParams
);

int IMRPhenomHMAmplitude(
    UNUSED REAL8Sequence *amps,
    UNUSED REAL8Sequence *freqs_geom,
    UNUSED PhenomHMStorage *pHM,
    UNUSED UINT4 ell,
    UNUSED INT4 mm,
    size_t ind_min,
    size_t ind_max,
    UNUSED LALDict *extraParams
);

int IMRPhenomHMPhase(
    UNUSED REAL8Sequence *phases,
    UNUSED REAL8Sequence *freqs_geom,
    UNUSED PhenomHMStorage *pHM,
    UNUSED UINT4 ell,
    UNUSED INT4 mm,
    size_t ind_min,
    size_t ind_max,
    UNUSED LALDict *extraParams
);


#ifdef __cplusplus
}
#endif

#endif /* _LALSIM_IMR_PHENOMHMv2_H */
