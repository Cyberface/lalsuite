#ifndef _LALSIM_IMR_PHENOMP_H
#define _LALSIM_IMR_PHENOMP_H

/*
 *  Copyright (C) 2013,2014,2015,2016,2017 Michael Puerrer, Alejandro Bohe, Sebastian Khan
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

#include <lal/LALStdlib.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConstants.h>

#include "LALSimIMRPhenomC_internals.h"
#include "LALSimIMRPhenomD_internals.h"

#include <lal/FrequencySeries.h>
#include <lal/LALSimInspiral.h>


/* CONSTANTS */

/**
 * Tolerance used below which numbers are treated as zero for the calculation of atan2
 */
#define MAX_TOL_ATAN 1.0e-15


/* ************************** PhenomP internal function prototypes *****************************/
/* atan2 wrapper that returns 0 when both magnitudes of x and y are below tol, otherwise it returns
   atan2(x, y) */
static REAL8 atan2tol(REAL8 x, REAL8 y, REAL8 tol);

/* PhenomC parameters for modified ringdown: Uses final spin formula of Barausse & Rezzolla, Astrophys.J.Lett.704:L40-L44, 2009 */
static BBHPhenomCParams *ComputeIMRPhenomCParamsRDmod(
  const REAL8 m1,   /**< Mass of companion 1 (solar masses) */
  const REAL8 m2,   /**< Mass of companion 2 (solar masses) */
  const REAL8 chi,  /**< Reduced aligned spin of the binary chi = (m1*chi1 + m2*chi2)/M */
  const REAL8 chip,  /**< Dimensionless spin in the orbital plane */
  LALDict *extraParams /**< linked list containing the extra testing GR parameters */
);

typedef struct tagNNLOanglecoeffs {
    REAL8 alphacoeff1; /* Coefficient of omega^(-1)   in alphaNNLO */
    REAL8 alphacoeff2; /* Coefficient of omega^(-2/3) in alphaNNLO */
    REAL8 alphacoeff3; /* Coefficient of omega^(-1/3) in alphaNNLO */
    REAL8 alphacoeff4; /* Coefficient of log(omega)   in alphaNNLO */
    REAL8 alphacoeff5; /* Coefficient of omega^(1/3)  in alphaNNLO */

    REAL8 epsiloncoeff1; /* Coefficient of omega^(-1)   in epsilonNNLO */
    REAL8 epsiloncoeff2; /* Coefficient of omega^(-2/3) in epsilonNNLO */
    REAL8 epsiloncoeff3; /* Coefficient of omega^(-1/3) in epsilonNNLO */
    REAL8 epsiloncoeff4; /* Coefficient of log(omega)   in epsilonNNLO */
    REAL8 epsiloncoeff5; /* Coefficient of omega^(1/3)  in epsilonNNLO */
} NNLOanglecoeffs;

static void ComputeNNLOanglecoeffs(
  NNLOanglecoeffs *angcoeffs,  /**< Output: Structure to store results */
  const REAL8 q,               /**< Mass-ratio (convention q>1) */
  const REAL8 chil,            /**< Dimensionless aligned spin of the largest BH */
  const REAL8 chip             /**< Dimensionless spin component in the orbital plane */
);

typedef struct tagSpinWeightedSphericalHarmonic_l2 {
  COMPLEX16 Y2m2, Y2m1, Y20, Y21, Y22;
} SpinWeightedSphericalHarmonic_l2;

/* Internal core function to calculate PhenomP polarizations for a sequence of frequences. */
static int PhenomPCore(
  COMPLEX16FrequencySeries **hptilde,   /**< Output: Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde,   /**< Output: Frequency-domain waveform hx */
  const REAL8 chi1_l_in,                /**< Dimensionless aligned spin on companion 1 */
  const REAL8 chi2_l_in,                /**< Dimensionless aligned spin on companion 2 */
  const REAL8 chip,                     /**< Effective spin in the orbital plane */
  const REAL8 thetaJ,                   /**< Angle between J0 and line of sight (z-direction) */
  const REAL8 m1_SI_in,                 /**< Mass of companion 1 (kg) */
  const REAL8 m2_SI_in,                 /**< Mass of companion 2 (kg) */
  const REAL8 distance,                 /**< Distance of source (m) */
  const REAL8 alpha0,                   /**< Initial value of alpha angle (azimuthal precession angle) */
  const REAL8 phic,                     /**< Orbital phase at the peak of the underlying non precessing model (rad) */
  const REAL8 f_ref,                    /**< Reference frequency */
  const REAL8Sequence *freqs,           /**< Frequency points at which to evaluate the waveform (Hz) */
  double deltaF,                        /**< Sampling frequency (Hz).
   * If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  IMRPhenomP_version_type IMRPhenomP_version, /**< IMRPhenomPv1 uses IMRPhenomC, IMRPhenomPv2 uses IMRPhenomD, IMRPhenomPv3 is same as IMRPhenomPv2 but uses the precession angles from arXiv 1703.03967. */
  LALDict *extraParams /**< linked list containing the extra testing GR parameters */
);

/* Internal core function to calculate PhenomP polarizations for a single frequency. */
static int PhenomPCoreOneFrequency(
  const REAL8 fHz,                        /**< Frequency (Hz) */
  const REAL8 eta,                        /**< Symmetric mass ratio */
  const REAL8 chi1_l,                     /**< Dimensionless aligned spin on companion 1 */
  const REAL8 chi2_l,                     /**< Dimensionless aligned spin on companion 2 */
  const REAL8 chip,                       /**< Dimensionless spin in the orbital plane */
  const REAL8 distance,                   /**< Distance of source (m) */
  const REAL8 M,                          /**< Total mass (Solar masses) */
  const REAL8 phic,                       /**< Orbital phase at the peak of the underlying non precessing model (rad) */
  IMRPhenomDAmplitudeCoefficients *pAmp,  /**< Internal IMRPhenomD amplitude coefficients */
  IMRPhenomDPhaseCoefficients *pPhi,      /**< Internal IMRPhenomD phase coefficients */
  BBHPhenomCParams *PCparams,             /**< Internal PhenomC parameters */
  PNPhasingSeries *PNparams,              /**< PN inspiral phase coefficients */
  NNLOanglecoeffs *angcoeffs,             /**< Struct with PN coeffs for the NNLO angles */
  SpinWeightedSphericalHarmonic_l2 *Y2m,  /**< Struct of l=2 spherical harmonics of spin weight -2 */
  const REAL8 alphaoffset,                /**< f_ref dependent offset for alpha angle (azimuthal precession angle) */
  const REAL8 epsilonoffset,              /**< f_ref dependent offset for epsilon angle */
  COMPLEX16 *hp,                          /**< Output: tilde h_+ */
  COMPLEX16 *hc,                          /**< Output: tilde h_+ */
  REAL8 *phasing,                         /**< Output: overall phasing */
  const UINT4 IMRPhenomP_version,         /**< Version number: 1 uses IMRPhenomC, 2 uses IMRPhenomD */
  AmpInsPrefactors *amp_prefactors,       /**< pre-calculated (cached for saving runtime) coefficients for amplitude. See LALSimIMRPhenomD_internals.c*/
  PhiInsPrefactors *phi_prefactors        /**< pre-calculated (cached for saving runtime) coefficients for phase. See LALSimIMRPhenomD_internals.*/
);

/* Simple 2PN version of L, without any spin terms expressed as a function of v */
static REAL8 L2PNR(
  const REAL8 v,   /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 eta  /**< Symmetric mass-ratio */
);

static REAL8 L2PNR_v1(
  const REAL8 v,   /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 eta  /**< Symmetric mass-ratio */
);

static void WignerdCoefficients(
  REAL8 *cos_beta_half,   /**< Output: cos(beta/2) */
  REAL8 *sin_beta_half,   /**< Output: sin(beta/2) */
  const REAL8 v,          /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 SL,         /**< Dimensionfull aligned spin */
  const REAL8 eta,        /**< Symmetric mass-ratio */
  const REAL8 Sp          /**< Dimensionfull spin component in the orbital plane */
);

static void WignerdCoefficients_SmallAngleApproximation(
  REAL8 *cos_beta_half, /**< Output: cos(beta/2) */
  REAL8 *sin_beta_half, /**< Output: sin(beta/2) */
  const REAL8 v,        /**< Cubic root of (Pi * Frequency (geometric)) */
  const REAL8 SL,       /**< Dimensionfull aligned spin */
  const REAL8 eta,      /**< Symmetric mass-ratio */
  const REAL8 Sp        /**< Dimensionfull spin component in the orbital plane */
);

static void CheckMaxOpeningAngle(
  const REAL8 m1,     /**< Mass of companion 1 (solar masses) */
  const REAL8 m2,     /**< Mass of companion 2 (solar masses) */
  const REAL8 chi1_l, /**< Aligned spin of BH 1 */
  const REAL8 chi2_l, /**< Aligned spin of BH 2 */
  const REAL8 chip    /**< Dimensionless spin in the orbital plane */
);

static REAL8 FinalSpinIMRPhenomD_all_in_plane_spin_on_larger_BH(
  const REAL8 m1,     /**< Mass of companion 1 (solar masses) */
  const REAL8 m2,     /**< Mass of companion 2 (solar masses) */
  const REAL8 chi1_l, /**< Aligned spin of BH 1 */
  const REAL8 chi2_l, /**< Aligned spin of BH 2 */
  const REAL8 chip    /**< Dimensionless spin in the orbital plane */
);

static REAL8 FinalSpinBarausse2009_all_spin_on_larger_BH(
  const REAL8 nu,     /**< Symmetric mass-ratio */
  const REAL8 chi,    /**< Effective aligned spin of the binary:  chi = (m1*chi1 + m2*chi2)/M  */
  const REAL8 chip    /**< Dimensionless spin in the orbital plane */
);

static REAL8 FinalSpinBarausse2009(  /* Barausse & Rezzolla, Astrophys.J.Lett.704:L40-L44, 2009, arXiv:0904.2577 */
  const REAL8 nu,               /**< Symmetric mass-ratio */
  const REAL8 a1,               /**< |a_1| norm of dimensionless spin vector for BH 1 */
  const REAL8 a2,               /**< |a_2| norm of dimensionless spin vector for BH 2 */
  const REAL8 cos_alpha,        /**< cos(alpha) = \\hat a_1 . \\hat a_2 (Eq. 7) */
  const REAL8 cos_beta_tilde,   /**< cos(\\tilde beta)  = \\hat a_1 . \\hat L (Eq. 9) */
  const REAL8 cos_gamma_tilde   /**< cos(\\tilde gamma) = \\hat a_2 . \\hat L (Eq. 9)*/
);

static bool approximately_equal(REAL8 x, REAL8 y, REAL8 epsilon);
static void nudge(REAL8 *x, REAL8 X, REAL8 epsilon);

/* BEGIN IMRPhenomPv3 */

/* CONSTANTS */

/* default and constant value places lhat = (0,0,1) */
#define LHAT_COS_THETA 1.0 /* Cosine of Polar angle of orbital angular momentum */
#define LHAT_PHI 0.0 /* Azimuthal angle of orbital angular momentum */

/**
 * function to convert from 3d cartesian components to polar angles and vector magnitude.
 * https://en.wikipedia.org/wiki/Spherical_coordinate_system#Cartesian_coordinates
 */
static void ComputeIMRPhenomPv3CartesianToPolar(REAL8 *polar, REAL8 *azimuthal, REAL8 *magnitude, REAL8 x, REAL8 y, REAL8 z);

/**
 * Structure storing initial and derived variables for IMRPhenomPv3
 */
typedef struct tagPhenomPv3Storage
{
    /* input parameters */
    REAL8 m1_SI; /**< mass of primary in SI (kg) */
    REAL8 m2_SI; /**< mass of secondary in SI (kg) */
    REAL8 chi1x; /**< x-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    REAL8 chi1y; /**< y-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    REAL8 chi1z; /**< z-component of dimensionless spin on primary w.r.t. Lhat = (0,0,1) */
    REAL8 chi2x; /**< x-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    REAL8 chi2y; /**< y-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    REAL8 chi2z; /**< z-component of dimensionless spin on secondary w.r.t. Lhat = (0,0,1) */
    REAL8 distance_SI; /**< distance to source in SI (m) */
    REAL8 inclination; /**< inclination - used to compute the angle thetaJN (rad) */
    REAL8 phiRef; /**< */
    REAL8 deltaF; /**< frequency spacing (Hz) */
    REAL8 f_min; /**< starting GW frequency (Hz) */
    REAL8 f_max; /**< ending GW frequency (Hz) */
    REAL8 f_ref; /**< reference GW frequency (Hz) */
    /* derived parameters */
    REAL8 m1_Msun; /**< mass of primary in solar masses */
    REAL8 m2_Msun; /**< mass of secondary in solar masses */
    REAL8 Mtot_SI; /**< total mass in SI (kg) */
    REAL8 Mtot_Msun; /**< total mass in solar masses */
    REAL8 eta; /**< Symmetric mass ratio*/
    REAL8 q; /* with m1>=m2 so q>=1 */
    REAL8 Msec; /**< Total mass in seconds */
    REAL8 piM; /**< LAL_PI * Msec */
    REAL8 f_ref_Orb_Hz; /**< Reference orbital frequency (Hz) [It's the reference GW frequency converted to orbital frequency] */
    /* variables used when rotating input parameters (LAL frame) into PhenomP intrinsic parameters  */
    REAL8 chip; /**< effective precessing parameter */
    REAL8 thetaJN; /**< Angle between J0 and line of sight (z-direction) */
    REAL8 alpha0; /**< Initial value of alpha angle (azimuthal precession angle) */
    REAL8 phi_aligned; /**< Initial phase to feed the underlying aligned-spin model */
    REAL8 zeta_polariz; /**< Angle to rotate the polarizations */
    /* compute spins in polar coordinates */
    REAL8 chi1_mag; /**< dimensionless spin magnitude on primary */
    REAL8 chi1_theta; /**< polar angle w.r.t. Lhat = (0,0,1) on primary */
    REAL8 chi1_phi; /**< azimuthal angle w.r.t. Lhat = (0,0,1) on primary */
    REAL8 chi2_mag; /**< dimensionless spin magnitude on secondary */
    REAL8 chi2_theta; /**< polar angle w.r.t. Lhat = (0,0,1) on secondary */
    REAL8 chi2_phi; /**< azimuthal angle w.r.t. Lhat = (0,0,1) on secondary */
    /* Precession angles at reference frequency */
    REAL8 alphaRef; /**< azimuthal precession angle at f_ref */
    REAL8 epsilonRef; /**< epsilon precession angle at f_ref */
    REAL8 betaRef; /**< beta (opening angle) precession angle at f_ref */
    /* spherical harmonics */
    SpinWeightedSphericalHarmonic_l2 Y2m;
    /* PhenomD parameters */
    REAL8 finspin; /**< dimensionless final spin */
} PhenomPv3Storage;

/**
 * function to swap masses and spins to enforece m1 >= m2
 */
static int PhenomPv3EnforcePrimaryIsm1(REAL8 *m1, REAL8 *m2, REAL8 *chi1x, REAL8 *chi1y, REAL8 *chi1z, REAL8 *chi2x, REAL8 *chi2y, REAL8 *chi2z);

/**
 * must be called before the first usage of *p
 */
static int init_PhenomPv3_Storage(PhenomPv3Storage *p, REAL8 m1_SI, REAL8 m2_SI, REAL8 S1x, REAL8 S1y, REAL8 S1z, REAL8 S2x, REAL8 S2y, REAL8 S2z, const REAL8 distance, const REAL8 inclination, const REAL8 phiRef, const REAL8 deltaF, const REAL8 f_min, const REAL8 f_max, const REAL8 f_ref);

/* Internal core function to calculate PhenomP version 3 polarizations.
 * This function computes all quantities that are independent of frequency
 * to be passed into PhenomPv3CoreOneFrequency
 */
static int PhenomPv3Core(
  COMPLEX16FrequencySeries **hptilde,        /**< [out] Frequency-domain waveform h+ */
  COMPLEX16FrequencySeries **hctilde,        /**< [out] Frequency-domain waveform hx */
  PhenomPv3Storage *pv3,       /**< PhenomPv3Storage Struct for storing internal variables */
  const REAL8Sequence *freqs_in,             /**< Frequency points at which to evaluate the waveform (Hz) */
  double deltaF,                             /**< Sampling frequency (Hz).
   * If deltaF > 0, the frequency points given in freqs are uniformly spaced with
   * spacing deltaF. Otherwise, the frequency points are spaced non-uniformly.
   * Then we will use deltaF = 0 to create the frequency series we return. */
  LALDict *extraParams /**<linked list containing the extra testing GR parameters */
  );

//   static int PhenomPv3CoreOneFrequency(
//     COMPLEX16 *hp,                              /**< [out] plus polarization \f$\tilde h_+\f$ */
//     COMPLEX16 *hc,                              /**< [out] cross polarization \f$\tilde h_x\f$ */
//     REAL8 *phasing,                             /**< [out] overall phasing */
//     const REAL8 fHz,                            /**< Frequency (Hz) */
//     PhenomPv3Storage *pv3,                      /**< PhenomPv3Storage Struct for storing internal variables */
//     IMRPhenomDAmplitudeCoefficients *pAmp,      /**< Internal IMRPhenomD amplitude coefficients */
//     IMRPhenomDPhaseCoefficients *pPhi,          /**< Internal IMRPhenomD phase coefficients */
//     PNPhasingSeries *PNparams,                  /**< PN inspiral phase coefficients */
//     AmpInsPrefactors *amp_prefactors,           /**< pre-calculated (cached for saving runtime) coefficients for amplitude. See LALSimIMRPhenomD_internals.c*/
//     PhiInsPrefactors *phi_prefactors            /**< pre-calculated (cached for saving runtime) coefficients for phase. See LALSimIMRPhenomD_internals.*/
// );

/* END IMRPhenomPv3 */

#endif	// of #ifndef _LALSIM_IMR_PHENOMP_H
