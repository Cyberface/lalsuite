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

#include "LALSimInspiralFDPrecAngles_internals.h"

static sysq InitializeSystem(const double m1, const double m2, const double mul, const double phl, const double mu1, const double ph1, const double ch1, const double mu2, const double ph2, const double ch2, const double f_0)
{
    sysq system = {0,0.,{0.},{0.},{0.},{0.},0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    
    system.onethird = 1./3.;
    
    const double domegadt_csts_nonspin[17] = {96./5.,-1486./35.,-264./5.,384.*LAL_PI/5.,34103./945.,13661./105.,944./15.,LAL_PI*(-4159./35.),LAL_PI*(-2268./5.),(16447322263./7276500. + LAL_PI*LAL_PI*512./5. - LAL_LN2*109568./175. -LAL_GAMMA*54784./175.),(-56198689./11340. + LAL_PI*LAL_PI*902./5.),1623./140.,-1121./27.,-54784./525.,-LAL_PI*883./42.,LAL_PI*71735./63.,LAL_PI*73196./63.};
    const double domegadt_csts_spinorbit[18] = {-904./5.,-120.,-62638./105.,4636./5.,-6472./35.,3372./5.,-LAL_PI*720.,-LAL_PI*2416./5.,-208520./63.,796069./105.,-100019./45.,-1195759./945.,514046./105.,-8709./5.,-LAL_PI*307708./105.,LAL_PI*44011./7.,-LAL_PI*7992./7.,LAL_PI*151449./35.};
    const double domegadt_csts_spinspin[4] = {-494./5.,-1442./5.,-233./5.,-719./5.};
    const double L_csts_nonspin[9] = {3./2.,1./6.,27./8.,-19./8.,1./24.,135./16.,-6889/144.+ 41./24.*LAL_PI*LAL_PI,31./24.,7./1296.};
    const double L_csts_spinorbit[6] = {-14./6.,-3./2.,-11./2.,133./72.,-33./8.,7./4.};
    
    const double G_over_c = LAL_G_SI/LAL_C_SI;
    const double M = m1 + m2;
    system.q = m2/m1;
    const double q = system.q;
    system.nu = q/(1+q)/(1+q);
    const double nu = system.nu;
    system.nu_2 = nu*nu;
    const double nu_2 = system.nu_2;
    const double piGM_over_cthree = LAL_PI*G_over_c/LAL_C_SI/LAL_C_SI*(m1+m2);
    system.deltam_over_M = (1-q)/(1+q);
    const double deltam_over_M = system.deltam_over_M;
    
    const double S1_norm = fabs(ch1)*nu/q;//in geometric units and with M=1
    const double S2_norm = fabs(ch2)*nu*q;

    const double xi_0 = pow(piGM_over_cthree*f_0, system.onethird);
    const double xi0_2 = xi_0*xi_0;
    const vector Lhat_0 = CreateSphere(1.,acos(mul),phl);
    const vector S1_0 = CreateSphere(S1_norm,acos(mu1),ph1);//in geometric units and with M=1
    const vector S2_0 = CreateSphere(S2_norm,acos(mu2),ph2);
    const vector L_0 = ScalarProd(nu/xi_0,Lhat_0);
    
    system.dot1 = DotProd(S1_0,Lhat_0);
    system.dot2 = DotProd(S2_0,Lhat_0);
    system.dot12 = DotProd(S1_0,S2_0);

    system.Seff = ((system.dot1)/m1 +(system.dot2)/m2)*M;
    
    const vector S_0 = Sum(S1_0,S2_0);
    const vector J_0 = Sum(L_0,S_0);
    
    const double S_0_norm = Norm(S_0);
    const double J_0_norm = Norm(J_0);
    const double L_0_norm = Norm(L_0);
    system.S1_norm_2 = S1_norm*S1_norm;
    system.S2_norm_2 = S2_norm*S2_norm;
    
    const vector roots = Roots(L_0_norm,J_0_norm,&system);

    const double q_2 = q*q;
    const double one_m_q_sq = (1.-q)*(1.-q);
    const double one_m_q_4 = one_m_q_sq*one_m_q_sq;
    const double one_p_q_sq = (1.+q)*(1.+q);
    
    const double Save_square = 0.5*(roots.z+roots.y);
    const double S1_norm_2 = system.S1_norm_2;
    const double S2_norm_2 = system.S2_norm_2;
    const double Seff = system.Seff;
    const double Seff_2 = Seff*Seff;

    //computed with initial spin couplings, they're not exactly accurate for generic precession, but the correction should be 4PN
    //these constants are used in TaylorT1 where domega/dt is expressed as a polynomial
    const double a0 = nu*domegadt_csts_nonspin[0];
    const double a2 = nu*(domegadt_csts_nonspin[1] + nu*(domegadt_csts_nonspin[2]));
    const double a3 = nu*(domegadt_csts_nonspin[3] + beta(domegadt_csts_spinorbit[0], domegadt_csts_spinorbit[1], &system));
    const double a4 = nu*(domegadt_csts_nonspin[4] + nu*(domegadt_csts_nonspin[5] + nu*(domegadt_csts_nonspin[6])) + sigma(domegadt_csts_spinspin[0], domegadt_csts_spinspin[1], &system) + tau(domegadt_csts_spinspin[2], domegadt_csts_spinspin[3], &system));
    const double a5 = nu*(domegadt_csts_nonspin[7] + nu*(domegadt_csts_nonspin[8]) + beta((domegadt_csts_spinorbit[2] + nu*(domegadt_csts_spinorbit[3])), (domegadt_csts_spinorbit[4] + nu*(domegadt_csts_spinorbit[5])), &system));
    
    const double a0_2 = a0*a0;
    const double a0_3 = a0_2*a0;
    const double a2_2 = a2*a2;

    //these constants are used in TaylorT2 where domega/dt is expressed as an inverse polynomial
    const double c0 = 1./a0;
    const double c2 = -a2/a0_2;
    const double c3 = -a3/a0_2;
    const double c4 = (a2_2 - a0*a4)/a0_3;
    const double c5 = (2.*a2*a3 - a0*a5)/a0_3;

    system.nu_4 = nu_2*nu_2;
    const double nu_4 = system.nu_4;
    
    system.c_1 =0.5*(J_0_norm*J_0_norm - L_0_norm*L_0_norm - Save_square)/L_0_norm*nu;
    double c_1 = system.c_1;
    system.c_1_over_nu = system.c_1/nu;
    double c_1_over_nu = system.c_1_over_nu;
    system.c1_2 = c_1*c_1;
    double c1_2 = system.c1_2;
    double c1_over_nu_2 = c_1_over_nu*c_1_over_nu;
    
    const double Delta = sqrt((4.*c1_over_nu_2*one_p_q_sq - 8.*c_1_over_nu*q*(1.+q)*Seff - 4.*((1.-q_2)*(1.-q_2)*S1_norm_2-q_2*Seff_2))*(4.*c1_over_nu_2*q_2*one_p_q_sq - 8.*c_1_over_nu*q_2*(1.+q)*Seff - 4.*((1.-q_2)*(1.-q_2)*S2_norm_2-q_2*Seff_2)));
    
    system.constants_u[0] = -c0;
    system.constants_u[1] = 6.*Seff*nu/deltam_over_M/deltam_over_M - 3.*c_1_over_nu/deltam_over_M/deltam_over_M;
    system.constants_u[2] = 3.*c2/c0 + 0.75*one_p_q_sq/one_m_q_4*(-20.*c1_over_nu_2*q_2*one_p_q_sq + 2.*(1.-q_2)*(1.-q_2)*(q*(2.+q)*S1_norm_2 + (1.+2.*q)*S2_norm_2 -2.*q*Save_square) + 2.*q_2*(7.+6.*q+7.*q_2)*2.*c_1_over_nu*Seff - 2.*q_2*(3.+4.*q+3.*q_2)*Seff_2 + q*Delta);
    
    system.Ssqave = 0.5*(roots.z+roots.y);
    system.sqrtSsqave = sqrt(system.Ssqave);

    
    const double Rm = roots.z - roots.y;
    const double Rm_2 = Rm*Rm;
    const double cp = roots.z*nu_2 - c1_2;
    const double cm = cp-Rm*nu_2;
    const double S0m = S1_norm_2 - S2_norm_2;
    const double cpcm = fabs(cp*cm);
    const double sqrt_cpcm = sqrt(cpcm);
    
    const double A1t = 0.5+0.75/nu;//xi^6
    const double A2t = -0.75*Seff/nu;//xi^7
    const double A1ave = (cp-sqrt_cpcm)/nu_2 ;
    const double Bave = -0.5*Rm*sqrt_cpcm/nu_2 - cp/nu_4*(sqrt_cpcm-cp);
    
    const double aw = (-3.*(1. + q)/q*(2.*(1. + q)*nu_2*Seff*c_1 - (1. + q)*c1_2 + (1. - q)*nu_2*S0m));
    const double cw = 3./32./nu*Rm_2;
    const double dw = 4.*cp - 4.*A1ave*nu_2;
    const double hw = -2*(2*A1ave - Rm)*c_1;
    const double fw = Rm*A1ave-Bave-0.25*Rm_2;
    
    const double ad = aw/dw;
    const double hd = hw/dw;
    const double cd = cw/dw;
    const double fd = fw/dw;
    
    const double hd_2 = hd*hd;
    const double hd_4 = hd_2*hd_2;
    const double adfd = ad*fd;
    const double adfdhd = adfd*hd;
    const double adfdhd_2 = adfd*hd_2;
    const double adhd = ad*hd;
    const double adhd_2 = ad*hd_2;
    const double adhd_3 = ad*hd_2*hd;
    const double adhd_4 = ad*hd_4;
    const double cdfd = cd*fd;
    const double cdhd = cd*hd;
    const double cdhd_2 = cd*hd_2;
    
    double Omegaz0 = A1t + ad;
    double Omegaz1 = A2t - ad*Seff - adhd;
    double Omegaz2 = cd - adfd + adhd_2 + adhd*Seff;
    double Omegaz3 = (adfd - cd - adhd_2)*Seff + 2*adfdhd - adhd_3 - cdhd;
    double Omegaz4 = -(2*adfdhd - adhd_3 - cdhd)*Seff + adfd*fd - cdfd + cdhd_2 - 3*adfdhd_2 + adhd_4;
    double Omegaz5 = -(adfd*fd - cdfd + cdhd_2 - 3*adfdhd_2 + adhd_4)*Seff + hd*(2*cdfd - 3*adfd*fd - cdhd_2 + 4*adfdhd_2 - adhd_4);
    
    system.constants_phiz[0] = 3.*Omegaz0*c0;
    system.constants_phiz[1] = 3.*Omegaz1*c0;
    system.constants_phiz[2] = 3.*(Omegaz2*c0 + Omegaz0*c2);
    system.constants_phiz[3] = 3.*(Omegaz0*c3 + Omegaz1*c2 + Omegaz3*c0);
    system.constants_phiz[4] = 3.*(Omegaz0*c4 + Omegaz1*c3 + Omegaz2*c2 + Omegaz4*c0);
    system.constants_phiz[5] = 3.*(c5*Omegaz0 + c4*Omegaz1 + c3*Omegaz2 + c2*Omegaz3 + c0*Omegaz5);

    const double gw = 3./16./nu_2/nu*Rm_2*(c_1 - nu_2*Seff);
    const double gd = gw/dw;
    
    Omegaz5 += Omegaz4*c_1/nu_2 - fd*gd + gd*hd_2 + gd*hd*Seff;
    Omegaz4 += Omegaz3*c_1/nu_2 - gd*hd - gd*Seff;
    Omegaz3 += Omegaz2*c_1/nu_2 + gd;
    Omegaz2 += Omegaz1*c_1/nu_2;
    Omegaz1 += Omegaz0*c_1/nu_2;
    
    system.constants_zeta[0] = 3.*Omegaz0*c0;
    system.constants_zeta[1] = 3.*Omegaz1*c0;
    system.constants_zeta[2] = 3.*(Omegaz2*c0 + Omegaz0*c2);
    system.constants_zeta[3] = 3.*(Omegaz0*c3 + Omegaz1*c2 + Omegaz3*c0);
    system.constants_zeta[4] = 3.*(Omegaz0*c4 + Omegaz1*c3 + Omegaz2*c2 + Omegaz4*c0);
    system.constants_zeta[5] = 3.*(Omegaz0*c5 + Omegaz1*c4 + Omegaz2*c3 + Omegaz3*c2 + Omegaz5*c0);
    
    const double m = sqrt((roots.y - roots.z)/(roots.x - roots.z));
    
    const double B = (S_0_norm*S_0_norm-roots.z)/(roots.y-roots.z);
    const double volumeellement = DotProd(CrossProd(L_0,S1_0),S2_0);
    const int sign_num = (volumeellement > 0) - (volumeellement < 0);
    
    if(S1_norm ==0 || S2_norm ==0) system.constant_of_S = 0;
    else{
        if(B < 0. || B > 1.) {
            if(B > 1 && B-1. < 0.00001) system.constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(1.)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0,xi0_2,&system);
            if(B < 0 && B > -0.00001) system.constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(0.)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0,xi0_2,&system);
        }
        else system.constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(B)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0,xi0_2,&system);
    }
    
    system.constants_L[0] = (L_csts_nonspin[0] + nu*L_csts_nonspin[1]);
    system.constants_L[1] = beta(L_csts_spinorbit[0], L_csts_spinorbit[1], &system);
    system.constants_L[2] = (L_csts_nonspin[2] + nu*L_csts_nonspin[3] + nu*nu*L_csts_nonspin[4]);
    system.constants_L[3] = beta((L_csts_spinorbit[2]+L_csts_spinorbit[3]*nu), (L_csts_spinorbit[4]+L_csts_spinorbit[5]*nu), &system);
    system.constants_L[4] = (L_csts_nonspin[5]+L_csts_nonspin[6]*nu +L_csts_nonspin[7]*nu*nu+L_csts_nonspin[8]*nu*nu*nu);
    
    system.flag = 0;
    if(S1_norm == 0 || S2_norm ==0 || (fabs(DotProd(S1_0,Lhat_0)/S1_norm) == 1 && fabs(DotProd(S2_0,Lhat_0)/S2_norm) == 1)) system.flag = 1;
    
    vector MScorrections = {0.,0.,0.};
    if(system.flag == 0){
        MScorrections = computeMScorrections(xi_0,xi0_2,L_0_norm,J_0_norm,roots,&system);
    }
    
    system.phiz_0 = 0.;
    system.phiz_0 = - phiz_of_xi(xi_0,xi0_2,J_0_norm,&system) - MScorrections.x;
    system.zeta_0 = 0.;
    system.zeta_0 = - zeta_of_xi(xi_0,xi0_2,&system) - MScorrections.y;
    
    
    return system;
}

/* *********************************************************************************/
/* Internal function that computes phiz, zeta, and costhetaL at 3PN                */
/* *********************************************************************************/

static vector compute_phiz_zeta_costhetaL3PN(const double xi, const sysq *system)
{
    vector angles;
    const double xi_2 = xi*xi;
    const double L_norm = ((*system).nu)/xi;
    const double L_norm3PN = L_norm_3PN_of_xi(xi,xi_2,L_norm,system);
    const double J_norm3PN = J_norm_of_xi(L_norm3PN,system);
    const double J_norm = J_norm_of_xi(L_norm,system);
    const vector roots = Roots(L_norm,J_norm,system);
    const double S_norm = S_norm_of_xi(xi,xi_2,roots,system);
    
    vector MScorrections = {0.,0.,0.};
    if((*system).flag == 0){
        MScorrections = computeMScorrections(xi,xi_2,L_norm,J_norm,roots,system);
    }
    
    angles.x = phiz_of_xi(xi,xi_2,J_norm,system) + MScorrections.x;
    angles.y = zeta_of_xi(xi,xi_2,system) + MScorrections.y;
    angles.z = costhetaL(J_norm3PN,L_norm3PN,S_norm);//costhetaL 3PN
    
    return angles;
}

/* *********************************************************************************/
/* Internal function that computes phiz, zeta, and costhetaL at Newtonian order    */
/* *********************************************************************************/

static vector compute_phiz_zeta_costhetaL(const double xi, const sysq *system)
{
    vector angles;
    const double xi_2 = xi*xi;
    const double L_norm = ((*system).nu)/xi;
    const double J_norm = J_norm_of_xi(L_norm,system);
    const vector roots = Roots(L_norm,J_norm,system);
    const double S_norm = S_norm_of_xi(xi,xi_2,roots,system);
    
    vector MScorrections = {0.,0.,0.};
    if((*system).flag == 0){
        MScorrections = computeMScorrections(xi,xi_2,L_norm,J_norm,roots,system);
    }
    
    angles.x = phiz_of_xi(xi,xi_2,J_norm,system) + MScorrections.x;
    angles.y = zeta_of_xi(xi,xi_2,system) + MScorrections.y;
    angles.z = costhetaL(J_norm,L_norm,S_norm);//costhetaL Newtonian
    
    return angles;
}

/* *********************************************************************************/
/* Internal function that computes the roots of Eq. 22 in arxiv:1703.03967         */
/* *********************************************************************************/

static vector Roots(const double L_norm, const double J_norm, const sysq *system)
{
    vector out;
    vector coeffs = BCDcoeff(L_norm, J_norm, system);
    double A1, A2, A3;
    const double B_2 = coeffs.x*coeffs.x;
    const double B_3 = B_2*coeffs.x;
    const double B_C =  coeffs.x*coeffs.y;
    
    const double p = coeffs.y - ((*system).onethird)*B_2;
    const double qc = 2./27.*B_3 - ((*system).onethird)*B_C + coeffs.z;
    const double sqrtarg = sqrt(-3./p);
    double acosarg = 1.5*qc/p*sqrtarg;
    if((acosarg + 1) < 0.00000000001) acosarg = -1;
    
    if(p < 0 && acosarg <= 1 && acosarg >= -1.){
        out.x = 2./sqrtarg*cos(((*system).onethird)*acos(acosarg)) - ((*system).onethird)*coeffs.x;
        out.y = 2./sqrtarg*cos(((*system).onethird)*acos(acosarg) - LAL_TWOPI*((*system).onethird)) - ((*system).onethird)*coeffs.x;
        out.z = 2./sqrtarg*cos(((*system).onethird)*acos(acosarg) - 2.*LAL_TWOPI*((*system).onethird)) - ((*system).onethird)*coeffs.x;
        
        A3 = fmax(fmax(out.x,out.y),out.z);
        A1 = fmin(fmin(out.x,out.y),out.z);
        if((A3 - out.z) > 0 && (A1 - out.z) < 0) A2 = out.z;
        else if((A3 - out.x) > 0 && (A1 - out.x) < 0) A2 = out.x;
        else A2 = out.y;
        
        out.x = A1;
        out.y = A2;
        out.z = A3;
    }
    else{
        printf("complex roots, %f, %f\n",p, acosarg);
        out.x = 0;
        out.y = 0;
        out.z = 0;
    }
    
    return out;
}

/* *********************************************************************************/
/* Internal function that computes the coefficients of Eq. 22 in arxiv:1703.03967  */
/* *********************************************************************************/

static vector BCDcoeff(const double L_norm, const double J_norm, const sysq *system)
{
    vector out;
    const double J_norm_2 = J_norm*J_norm;
    const double L_norm_2 = L_norm*L_norm;
    
    out.x = (L_norm_2+((*system).S1_norm_2))*((*system).q) + (L_norm_2+((*system).S2_norm_2))/((*system).q) + 2.*L_norm*((*system).Seff) - 2.*J_norm_2 - ((*system).S1_norm_2) - ((*system).S2_norm_2);
    
    out.y = (L_norm_2 - J_norm_2)*(L_norm_2 - J_norm_2) + 2.*((*system).Seff)*L_norm*(L_norm_2 - J_norm_2) - 2.*L_norm_2*(((*system).S1_norm_2) - ((*system).S2_norm_2)*((*system).q))*(1./((*system).q)-1) + 4.*((*system).nu)*((*system).Seff)*((*system).Seff)*L_norm_2 - 2.*((*system).Seff)*L_norm*(((*system).S1_norm_2)-((*system).S2_norm_2))*((*system).deltam_over_M) + 2.*J_norm_2*(((*system).S1_norm_2)*((*system).q) - ((*system).S2_norm_2))*(1./((*system).q)-1);
    
    out.z = -(L_norm_2 - J_norm_2)*(L_norm_2 - J_norm_2)*(((*system).S1_norm_2)*((*system).q) - ((*system).S2_norm_2))*(1./((*system).q)-1) - 2.*((*system).Seff)*((*system).deltam_over_M)*(((*system).S1_norm_2) - ((*system).S2_norm_2))*L_norm*(L_norm_2 - J_norm_2) + (((*system).S1_norm_2) - ((*system).S2_norm_2))*(((*system).S1_norm_2) - ((*system).S2_norm_2))*L_norm_2*((*system).deltam_over_M)*((*system).deltam_over_M)/((*system).nu);
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of J to Newtonian order            */
/* *********************************************************************************/

static double J_norm_of_xi(const double L_norm, const sysq *system)
{
    return sqrt(L_norm*L_norm + 2.*L_norm*((*system).c_1_over_nu) + (*system).Ssqave);
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of S divided by GMsquare_over_c    */
/* *********************************************************************************/

static double S_norm_of_xi(const double xi, const double xi_2, const vector roots, const sysq *system)
{
    double sn, cn, dn;
    
    const double m = (roots.y - roots.z)/(roots.x - roots.z);
    const double u = u_of_xi(xi,xi_2,system)+(*system).constant_of_S;
    
    gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
    
    if((*system).S1_norm_2 == 0. || (*system).S2_norm_2 == 0. || ((*system).dot1 == sqrt((*system).S1_norm_2)  && (*system).dot2 == sqrt((*system).S2_norm_2))) sn = 0.;
    
    const double S_norm_square_bar = roots.z + (roots.y - roots.z)*sn*sn;
    return sqrt(S_norm_square_bar);
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of L divided by GMsquare_over_c to */
/* 3PN order                                                                       */
/* *********************************************************************************/

static double L_norm_3PN_of_xi(const double xi, const double xi_2, const double L_norm, const sysq *system)
{
    //from 0605140 and Blanchet LRR
    return L_norm*(1. + xi_2*((*system).constants_L[0] + xi*(*system).constants_L[1] + xi_2*((*system).constants_L[2] + xi*(*system).constants_L[3] + xi_2*((*system).constants_L[4]))));
}

/* *********************************************************************************/
/* Internal function that returns a coefficient. I don't really remember           */
/* what it is right now.                                                           */
/* *********************************************************************************/

static vector c(const double xi, const double xi_2, const double J_norm, const vector roots, const sysq *system)
{
    const double xi_3 = xi_2*xi;
    const double xi_4 = xi_3*xi;
    const double xi_6 = xi_3*xi_3;
    
    const double J_norm_2 = J_norm*J_norm;
    
    vector out;
    
    out.x = -0.75*((J_norm_2-roots.z)*(J_norm_2-roots.z)*xi_4/((*system).nu) - 4.*((*system).nu)*((*system).Seff)*(J_norm_2-roots.z)*xi_3-2.*(J_norm_2-roots.z+2*(((*system).S1_norm_2)-((*system).S2_norm_2))*((*system).deltam_over_M))*((*system).nu)*xi_2+(4.*((*system).Seff)*xi+1)*((*system).nu)*((*system).nu_2))*J_norm*xi_2*(((*system).Seff)*xi-1.);
    out.y = 1.5*(roots.y-roots.z)*J_norm*((J_norm_2-roots.z)/((*system).nu)*xi_2-2.*((*system).nu)*((*system).Seff)*xi-((*system).nu))*(((*system).Seff)*xi-1.)*xi_4;
    out.z = -0.75*J_norm*(((*system).Seff)*xi-1.)*(roots.z-roots.y)*(roots.z-roots.y)*xi_6/((*system).nu);
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns a coefficient. I don't really remember           */
/* what it is right now.                                                           */
/* *********************************************************************************/

static vector d(const double L_norm, const double J_norm, const vector roots)
{
    vector out;
    
    out.x = -((L_norm+J_norm)*(L_norm+J_norm)-roots.z)*((L_norm-J_norm)*(L_norm-J_norm)-roots.z);
    out.y = 2.*(roots.y-roots.z)*(J_norm*J_norm+L_norm*L_norm-roots.z);
    out.z = -(roots.z-roots.y)*(roots.z-roots.y);
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns the cosine of the angle between L and J          */
/* *********************************************************************************/

static double costhetaL(const double J_norm, const double L_norm, const double S_norm)
{
    return 0.5*(J_norm*J_norm + L_norm*L_norm - S_norm*S_norm)/L_norm/J_norm;
}


/* *********************************************************************************/
/* Internal function that returns the phase of the magnitude of S                  */
/* *********************************************************************************/

static double u_of_xi(const double xi, const double xi_2, const sysq *system)
{
    return 0.75*((*system).deltam_over_M)*(*system).constants_u[0]/xi_2/xi*(1. + xi*((*system).constants_u[1] + xi*((*system).constants_u[2])));
}

/* *********************************************************************************/
/* Internal function that returns the derivative of the phase of the magnitude of S*/
/* *********************************************************************************/

static double psidot(const double xi, const double xi_2, const vector roots, const sysq *system)
{
    const double xi_3 = xi_2*xi;
    const double xi_6 = xi_3*xi_3;
    
    return 0.75/sqrt((*system).nu)*(1.-(*system).Seff*xi)*xi_6*sqrt(roots.z-roots.x);
}


/* *********************************************************************************/
/* Internal function that returns phiz                                             */
/* *********************************************************************************/

static double phiz_of_xi(const double xi, const double xi_2, const double J_norm, const sysq *system)
{
    const double log1 = log(fabs(((*system).c_1) + J_norm*((*system).nu)+((*system).nu_2)/xi));
    const double log2 = log(fabs(((*system).c_1) + J_norm*((*system).sqrtSsqave)*xi + ((*system).Ssqave)*xi));
    
    const double phiz0 = J_norm*(0.5*((*system).c1_2)/((*system).nu_4) - 0.5*((*system).onethird)*((*system).c_1)/xi/((*system).nu_2) - ((*system).onethird)*((*system).Ssqave)/((*system).nu_2) - ((*system).onethird)/xi_2) - 0.5*((*system).c_1)*(((*system).c1_2)/((*system).nu_4) - ((*system).Ssqave)/((*system).nu_2))/((*system).nu)*log1;
    const double phiz1 = -J_norm*0.5*(((*system).c_1)/((*system).nu_2) + 1./xi) + 0.5*(((*system).c1_2)/((*system).nu_2) - ((*system).Ssqave))/((*system).nu)*log1 ;
    const double phiz2 = -J_norm + ((*system).sqrtSsqave)*log2 - ((*system).c_1)/((*system).nu)*log1;
    const double phiz3 = J_norm*xi -((*system).nu)*log1 + ((*system).c_1)/((*system).sqrtSsqave)*log2;
    const double phiz4 = J_norm*0.5*xi*(((*system).c_1)/((*system).Ssqave) + xi) - 0.5*(((*system).c1_2)/((*system).Ssqave) - ((*system).nu_2))/((*system).sqrtSsqave)*log2;
    const double phiz5 = J_norm*xi*(-0.5*((*system).c1_2)/((*system).Ssqave)/((*system).Ssqave) + 0.5*((*system).onethird)*((*system).c_1)*xi/((*system).Ssqave) + ((*system).onethird)*(xi_2 + ((*system).nu_2)/((*system).Ssqave))) + 0.5*((*system).c_1)*(((*system).c1_2)/((*system).Ssqave) - ((*system).nu_2))/((*system).Ssqave)/((*system).sqrtSsqave)*log2;
    
    return phiz0*(*system).constants_phiz[0] + phiz1*(*system).constants_phiz[1] + phiz2*(*system).constants_phiz[2] + phiz3*(*system).constants_phiz[3] + phiz4*(*system).constants_phiz[4] + phiz5*(*system).constants_phiz[5] + (*system).phiz_0;
}

/* *********************************************************************************/
/* Internal function that returns zeta                                             */
/* *********************************************************************************/

static double zeta_of_xi(const double xi, const double xi_2, const sysq *system)
{
    const double logxi = log(xi);
    const double xi_3 = xi_2*xi;
    
    return ((*system).nu)*(-(*system).onethird*(*system).constants_zeta[0]/xi_3 - 0.5*(*system).constants_zeta[1]/xi_2 - (*system).constants_zeta[2]/xi + (*system).constants_zeta[3]*logxi + (*system).constants_zeta[4]*xi  + 0.5*(*system).constants_zeta[5]*xi_2) + (*system).zeta_0;
}

/* *********************************************************************************/
/* Internal function that computes the MS corrections for phiz and zeta            */
/* *********************************************************************************/

static vector computeMScorrections (const double xi, const double xi_2, const double L_norm, const double J_norm, const vector roots, const sysq *system)
{
    vector out;
    const vector c_return = c(xi,xi_2,J_norm,roots,system);
    const vector d_return = d(L_norm,J_norm,roots);
    const double c0 = c_return.x;
    const double c2 = c_return.y;
    const double c4 = c_return.z;
    const double d0 = d_return.x;
    const double d2 = d_return.y;
    const double d4 = d_return.z;
    
    const double sqt = sqrt(fabs(d2*d2-4.*d0*d4));
    
    const double Aa = 0.5*(J_norm/L_norm + L_norm/J_norm - roots.z/J_norm/L_norm);
    const double Bb = 0.5*(roots.z - roots.y)/L_norm/J_norm;
    
    const double twod0_d2 = 2*d0+d2;
    const double twod0d4 = 2*d0*d4;
    const double d2_twod4 = d2+2.*d4;
    const double d2_d4 =d2+d4;
    const double d0_d2_d4 = d0+d2+d4;
    const double n3 = (2.*d0_d2_d4/(twod0_d2+sqt));
    const double n4 = ((twod0_d2+sqt)/(2.*d0));
    const double sqtn3 = sqrt(fabs(n3));
    const double sqtn4 = sqrt(fabs(n4));
    const double den = 2.*d0*d0_d2_d4*sqt;
    
    const double u = u_of_xi(xi,xi_2,system)+(*system).constant_of_S;
    const double tanu = tan(u);
    const double atanu = atan(tan(u));
    const double psdot = psidot(xi,xi_2,roots,system);
    
    double L3, L4;
    if(n3 == 1) L3 = 0;
    else L3 = -fabs((c4*d0*(twod0_d2+sqt)-c2*d0*(d2_twod4-sqt)-c0*(twod0d4-d2_d4*(d2-sqt)))/(den))*(sqtn3/(n3-1.)*(atanu - atan(sqtn3*tanu)))/psdot;
    
    if(n4 == 1) L4 = 0;
    else L4 = -fabs((-c4*d0*(twod0_d2-sqt)+c2*d0*(d2_twod4+sqt)-c0*(-twod0d4+d2_d4*(d2+sqt))))/(den)*(sqtn4/(n4-1.)*(atanu - atan(sqtn4*tanu)))/psdot;
    
    out.x = L3+L4;
    if (out.x != out.x) out.x=0;
    
    out.y = Aa*out.x + 2.*Bb*d0*(L3/(sqt-d2)-L4/(sqt+d2));
    if (out.y != out.y) out.y=0;
    
    out.z = 0;
    
    return out;
}


/* *********************************************************************************/
/* Internal function that computes the spin-orbit couplings                        */
/* *********************************************************************************/

static double beta(const double a, const double b, const sysq *system)
{
    return ((((*system).dot1)*(a + b*((*system).q))) + (((*system).dot2)*(a + b/((*system).q))));
}

/* *********************************************************************************/
/* Internal function that computes the spin-spin couplings                         */
/* *********************************************************************************/

static double sigma(const double a, const double b, const sysq *system)
{
    return (a*((*system).dot12) - b*((*system).dot1)*((*system).dot2))/((*system).nu);
}

/* *********************************************************************************/
/* Internal function that computes the spin-spin couplings                         */
/* *********************************************************************************/

static double tau(const double a, const double b, const sysq *system)
{
    return (((*system).q)*((*system).S1_norm_2*a - b*((*system).dot1)*((*system).dot1)) + ((*system).S2_norm_2*a - b*((*system).dot2)*((*system).dot2))/((*system).q))/((*system).nu);
}

/* *********************************************************************************/
/* Internal function that returns the dot product of two vectors                   */
/* *********************************************************************************/

static double DotProd(const vector vec1, const vector vec2)
{
    return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
}

/* *********************************************************************************/
/* Internal function that returns the norm of a vector                             */
/* *********************************************************************************/

static double Norm(const vector vec)
{
    return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

/* *********************************************************************************/
/* Internal function that calculates a vector fro its spherical components         */
/* *********************************************************************************/

static vector CreateSphere(const double r, const double th, const double ph)
{
    vector out;
    const double fact = r*sin(th);
    out.x = fact*cos(ph);
    out.y = fact*sin(ph);
    out.z = r*cos(th);
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns the scalar product of a vector with a scalar     */
/* *********************************************************************************/

static vector ScalarProd(double c, vector vec)
{
    vector out;
    out.x = c*vec.x;
    out.y = c*vec.y;
    out.z = c*vec.z;
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns the sum of two vectors                           */
/* *********************************************************************************/

static vector Sum(vector vec1, vector vec2)
{
    vector out;
    out.x = vec1.x+vec2.x;
    out.y = vec1.y+vec2.y;
    out.z = vec1.z+vec2.z;
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns the cross product of two vectors                 */
/* *********************************************************************************/


static vector CrossProd(vector vec1, vector vec2)
{
    vector out;
    out.x = vec1.y*vec2.z-vec1.z*vec2.y;
    out.y = vec1.z*vec2.x-vec1.x*vec2.z;
    out.z = vec1.x*vec2.y-vec1.y*vec2.x;
    
    return out;
}



