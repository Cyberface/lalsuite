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

typedef struct tagvector
{
    double x;
    double y;
    double z;
} vector;

void XLALInitializePrecession(void);
void XLALFreePrecession(void);
void XLALInitializeSystem(double m1, double m2, double costhetaL, double phiL, double costheta1, double phi1, double chi1, double costheta2, double phi2, double chi2, double f0);
double XLALcosthetaL(double f);
double XLALcosthetaL_3PN(double f);
double XLALphiz_of_xi(double f);
double XLALzeta_of_xi(double f);

double onethird = 0.;
double onesixth = 0.;
double twotothe_onethird = 0.;
double sqrt_three = 0.;
double pi = 0.;
double pi_over_two = 0.;
double pi_over_four = 0.;
double pisquare = 0.;
double pi_four = 0.;
double log_2 = 0.;
double log_2_2 = 0.;
double log_3 = 0.;
double log_3_2 = 0.;
double log_5 = 0.;
double log_7 = 0.;
double log_3_over_2 = 0.;
double twopi = 0.;
double gamma_Euler = 0.;
double Light_vel = 0.;
double Newton_c = 0.;
double SOLAR_MASS = 0.;
double G_over_c = 0.;
double G_over_csquare = 0.;
double G_over_cthree = 0.;
double cthree_over_G = 0.;
double csix_over_Gsquare = 0.;
double cseven_over_Gthree = 0.;
double ceight_over_Gfour = 0.;

double domegadt_csts_nonspin[17] = {0.};
double domegadt_csts_spinorbit[18] = {0.};
double domegadt_csts_spinspin[4] = {0.};
double L_csts_nonspin[9] = {0.};
double L_csts_spinorbit[6] = {0.};

double m1;
double m2;
double mul;
double phl;
double mu1;
double ph1;
double ch1;
double mu2;
double ph2;
double ch2;
double f_0;
double th1;
double th2;
double thl;

double M;
double mu;
double nu;
double nu_2;
double Msquare;
double Mthree;
double m1square;
double m2square;
double q;
double deltam;
double deltam_over_M;
double one_over_m1_square;
double one_over_m2_square;
double GMsquare_over_c;
double twopiGM_over_cthree;

double S1_norm;
double S2_norm;
double S1_norm_square;
double S2_norm_square;
double Seff;
double Ssqave;

double xi_0;

vector L_0;
vector Lhat_0;
vector S1_0;
vector S2_0;

double c_1;

double constants_u[3] = {0.};
double constants_phiz[6] = {0.};
double constants_zeta[6] = {0.};
double constants_L[5] = {0.};

double constant_of_S;
double phiz_0;
double zeta_0;

#endif /* _LALSIM_INS_FD_PREC_ANGLES */
