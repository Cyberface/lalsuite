/*
 * Copyright (C) 2017 Katerina Chatziioannou, Sebastian Khan
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

#include "LALSimInspiralFDPrecAngles_internals.c"
#include "LALSimIMR.h"

/* *********************************************************************************/
/* XLAL function that does everything.                                             */
/* *********************************************************************************/

int XLALComputeAngles(
    REAL8Sequence *phiz_of_f, /**< [out] azimuthal angle of L around J */
    REAL8Sequence *zeta_of_f, /**< [out] Third euler angle to describe L w.r.t. J  */
    REAL8Sequence *costhetaL_of_f, /**< [out] Cosine of polar angle between L and J */
    const REAL8Sequence *f, /**< list of input Orbtial frequencies (Hz) */
    const double m1, /**< Primary mass in SI (kg) */
    const double m2, /**< Secondary mass in SI (kg) */
    const double mul, /**< Cosine of Polar angle of orbital angular momentum */
    const double phl, /**< Azimuthal angle of orbital angular momentum  */
    const double mu1, /**< Cosine of Polar angle of primary spin w.r.t. orbital angular momentum */
    const double ph1, /**< Azimuthal angle of primary spin  */
    const double ch1, /**< Dimensionless spin magnitude of primary spin */
    const double mu2, /**< Cosine of Polar angle of secondary spin w.r.t. orbital angular momentum */
    const double ph2, /**< Azimuthal angle of secondary spin  */
    const double ch2, /**< Dimensionless spin magnitude of secondary spin */
    const double f_0 /**< Reference orbital frequency (Hz) */
){
    sysq system  = InitializeSystem(m1,m2,mul,phl,mu1,ph1,ch1,mu2,ph2,ch2,f_0);

    double xi;
    const double twopiGM_over_cthree = LAL_TWOPI*LAL_G_SI*(m1+m2)/LAL_C_SI/LAL_C_SI/LAL_C_SI;
    /* twopiGM_over_cthree is the same as -> LAL_TWOPI * LAL_MTSUN_SI * (m1+m2) / LAL_MSUN_SI */
    
    vector angles;

    for(UINT4 i = 0; i < (*f).length; i++){
         xi = pow(((*f).data[i])*twopiGM_over_cthree, system.onethird);
        angles = compute_phiz_zeta_costhetaL(xi,&system);
        (*phiz_of_f).data[i] = angles.x;
        (*zeta_of_f).data[i] = angles.y;
        (*costhetaL_of_f).data[i] = angles.z;
    }
    return XLAL_SUCCESS;
}

/* *********************************************************************************/
/* XLAL function that does everything at 3PN.                                             */
/* *********************************************************************************/

int XLALComputeAngles3PN(
    REAL8Sequence *phiz_of_f, /**< [out] azimuthal angle of L around J */
    REAL8Sequence *zeta_of_f, /**< [out] Third euler angle to describe L w.r.t. J  */
    REAL8Sequence *costhetaL_of_f, /**< [out] Cosine of polar angle between L and J */
    const REAL8Sequence *f, /**< list of input Orbtial frequencies (Hz) */
    const double m1, /**< Primary mass in SI (kg) */
    const double m2, /**< Secondary mass in SI (kg) */
    const double mul, /**< Cosine of Polar angle of orbital angular momentum */
    const double phl, /**< Azimuthal angle of orbital angular momentum  */
    const double mu1, /**< Cosine of Polar angle of primary spin w.r.t. orbital angular momentum */
    const double ph1, /**< Azimuthal angle of primary spin  */
    const double ch1, /**< Dimensionless spin magnitude of primary spin */
    const double mu2, /**< Cosine of Polar angle of secondary spin w.r.t. orbital angular momentum */
    const double ph2, /**< Azimuthal angle of secondary spin  */
    const double ch2, /**< Dimensionless spin magnitude of secondary spin */
    const double f_0 /**< Reference orbital frequency (Hz) */
){
    sysq system  = InitializeSystem(m1,m2,mul,phl,mu1,ph1,ch1,mu2,ph2,ch2,f_0);

    double xi;
    const double twopiGM_over_cthree = LAL_TWOPI*LAL_G_SI*(m1+m2)/LAL_C_SI/LAL_C_SI/LAL_C_SI;
    /* twopiGM_over_cthree is the same as -> LAL_TWOPI * LAL_MTSUN_SI * (m1+m2) / LAL_MSUN_SI */

    vector angles;

    for(UINT4 i = 0; i < (*f).length; i++){
        xi = pow(((*f).data[i])*twopiGM_over_cthree, system.onethird);
        angles = compute_phiz_zeta_costhetaL3PN(xi,&system);
        (*phiz_of_f).data[i] = angles.x;
        (*zeta_of_f).data[i] = angles.y;
        (*costhetaL_of_f).data[i] = angles.z;
    }
    return XLAL_SUCCESS;
}
