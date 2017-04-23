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

#ifndef _LALSIM_INS_FD_PREC_ANGLES
#define _LALSIM_INS_FD_PREC_ANGLES

#include <lal/LALConstants.h>

typedef struct tagvector
{
    double x;
    double y;
    double z;
} vector;

typedef struct tagsystemquantites
{
    int flag;
    double onethird;
    double constants_u[6];
    double constants_phiz[6];
    double constants_zeta[6];
    double constants_L[6];
    double phiz_0, zeta_0, constant_of_S;
    double c_1, Ssqave, sqrtSsqave, Seff, c1_2, nu_2, nu_4;
    double S1_norm_square, S2_norm_square, S1_normalized_2, S2_normalized_2;
    double dot1, dot2, dot12;
    double deltam_over_M, GMsquare_over_c, Mfour_in_L, nu, q;
} sysq;

void XLALComputeAngles(REAL8Sequence *phiz_of_f, REAL8Sequence *zeta_of_f, REAL8Sequence *costhetaL_of_f, REAL8Sequence *costhetaL3PN_of_f, const REAL8Sequence *f, const double m1, const double m2, const double mul, const double phl, const double mu1, const double ph1, const double ch1, const double mu2, const double ph2, double ch2, const double f_0);

#endif /* _LALSIM_INS_FD_PREC_ANGLES */
