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

/* *********************************************************************************/
/* XLAL function that does everything.                                             */
/* *********************************************************************************/

void XLALComputeAngles(
    REAL8Sequence *phiz_of_f,
    REAL8Sequence *zeta_of_f,
    REAL8Sequence *costhetaL_of_f,
    REAL8Sequence *costhetaL3PN_of_f,
    const REAL8Sequence *f,
    const double m1,
    const double m2,
    const double mul,
    const double phl,
    const double mu1,
    const double ph1,
    const double ch1,
    const double mu2,
    const double ph2,
    const double ch2,
    const double f_0
){
    sysq system  = InitializeSystem(m1,m2,mul,phl,mu1,ph1,ch1,mu2,ph2,ch2,f_0);
    
    double xi;
    const double twopiGM_over_cthree = LAL_TWOPI*LAL_G_SI*(m1+m2)/LAL_C_SI/LAL_C_SI/LAL_C_SI;
    
    for(UINT4 i = 0; i < (*f).length; i++){
         xi = pow(((*f).data[i])*twopiGM_over_cthree, system.onethird);
        (*phiz_of_f).data[i] = phiz_of_xi(xi,&system);
        (*zeta_of_f).data[i] = zeta_of_xi(xi,&system);
        (*costhetaL3PN_of_f).data[i] = costhetaL_3PN(xi,&system);
        (*costhetaL_of_f).data[i] = costhetaL(xi,&system);
    }
}
