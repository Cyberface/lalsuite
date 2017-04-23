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
    sysq system = {0,0.,{0.},{0.},{0.},{0.},0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0., 0.};
    
    system.onethird = 1./3.;
    
    const double domegadt_csts_nonspin[17] = {96./5.,-1486./35.,-264./5.,384.*LAL_PI/5.,34103./945.,13661./105.,944./15.,LAL_PI*(-4159./35.),LAL_PI*(-2268./5.),(16447322263./7276500. + LAL_PI*LAL_PI*512./5. - LAL_LN2*109568./175. -LAL_GAMMA*54784./175.),(-56198689./11340. + LAL_PI*LAL_PI*902./5.),1623./140.,-1121./27.,-54784./525.,-LAL_PI*883./42.,LAL_PI*71735./63.,LAL_PI*73196./63.};
    const double domegadt_csts_spinorbit[18] = {-904./5.,-120.,-62638./105.,4636./5.,-6472./35.,3372./5.,-LAL_PI*720.,-LAL_PI*2416./5.,-208520./63.,796069./105.,-100019./45.,-1195759./945.,514046./105.,-8709./5.,-LAL_PI*307708./105.,LAL_PI*44011./7.,-LAL_PI*7992./7.,LAL_PI*151449./35.};
    const double domegadt_csts_spinspin[4] = {-494./5.,-1442./5.,-233./5.,-719./5.};
    const double L_csts_nonspin[9] = {3./2.,1./6.,27./8.,-19./8.,1./24.,135./16.,-6889/144.+ 41./24.*LAL_PI*LAL_PI,31./24.,7./1296.};
    const double L_csts_spinorbit[6] = {-14./6.,-3./2.,-11./2.,133./72.,-33./8.,7./4.};
    
//    if(ch1 < 0.){
//        ch1 = -ch1;
//        mu1 = -mu1;
//        if(ph1 >= LAL_PI) ph1 -= LAL_PI;
//        else ph1 += LAL_PI;
//    }
//    if(ch2 < 0.){
//        ch2 = -ch2;
//        mu2 = -mu2;
//        if(ph2 >= LAL_PI) ph2 -= LAL_PI;
//        else ph2 += LAL_PI;
//    }
    
    const double G_over_c = LAL_G_SI/LAL_C_SI;
    const double M = m1 + m2;
    const double Msquare = M*M;
    system.nu = m1*m2/Msquare;
    const double nu = system.nu;
    system.nu_2 = nu*nu;
    const double nu_2 = system.nu_2;
    system.q = m2/m1;
    const double q = system.q;
    const double m1square = m1*m1;
    const double m2square = m2*m2;
    const double M_4_nu = Msquare*Msquare*nu;
    const double piGM_over_cthree = LAL_PI*G_over_c/LAL_C_SI/LAL_C_SI*M;
    system.deltam_over_M = (m1 - m2)/M;
    const double deltam_over_M = system.deltam_over_M;
    system.GMsquare_over_c = G_over_c*Msquare;
    const double GMsquare_over_c = system.GMsquare_over_c;
    system.Mfour_in_L = GMsquare_over_c*GMsquare_over_c;
    
    const double S1_norm = fabs(ch1)*m1square*G_over_c;
    const double S2_norm = fabs(ch2)*m2square*G_over_c;
    system.S1_norm_square = S1_norm*S1_norm;
    system.S2_norm_square = S2_norm*S2_norm;

    const double xi_0 = pow(piGM_over_cthree*f_0, system.onethird);
    const vector Lhat_0 = CreateSphere(1.,acos(mul),phl);
    const vector S1_0 = CreateSphere(S1_norm,acos(mu1),ph1);
    const vector S2_0 = CreateSphere(S2_norm,acos(mu2),ph2);
    const vector L_0 = ScalarProd(GMsquare_over_c*nu/xi_0,Lhat_0);
    
    system.dot1 = DotProd(S1_0,Lhat_0);
    system.dot2 = DotProd(S2_0,Lhat_0);
    system.dot12 = DotProd(S1_0,S2_0);

    system.Seff = ((system.dot1)/m1 +(system.dot2)/m2)/M/G_over_c;
    
    const vector S_0 = Sum(S1_0,S2_0);
    const vector J_0 = Sum(L_0,S_0);
    
    const double S_0_norm = Norm(S_0);
    const double J_0_norm = Norm(J_0);
    const double L_0_norm = Norm(L_0);
    const double S1_normalized = S1_norm/GMsquare_over_c;
    const double S2_normalized = S2_norm/GMsquare_over_c;
    system.S1_normalized_2 = S1_normalized*S1_normalized;
    system.S2_normalized_2 = S2_normalized*S2_normalized;
    
    const vector roots = Roots(xi_0,J_0_norm,&system);
    double A1, A2, A3;

    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;

    const double q_2 = q*q;
    const double one_m_q_sq = (1.-q)*(1.-q);
    const double one_m_q_4 = one_m_q_sq*one_m_q_sq;
    const double one_p_q_sq = (1.+q)*(1.+q);
    
    const double Save_square = 0.5*(A3+A2)*GMsquare_over_c*GMsquare_over_c;
    const double Save_square_normalized = 0.5*(A3+A2);
    const double S1_normalized_2 = system.S1_normalized_2;
    const double S2_normalized_2 = system.S2_normalized_2;
    const double Seff = system.Seff;
    const double Seff_2 = Seff*Seff;

    system.c_1 = 0.5*(J_0_norm*J_0_norm - L_0_norm*L_0_norm - Save_square)/L_0_norm/GMsquare_over_c*nu;

    const double factor_a = nu;
    
    //computed with initial spin couplings, they're not exactly accurate for generic precession, but the correction should be 4PN
    //these constants are used in TaylorT1 where domega/dt is expressed as a polynomial
    const double a0 = factor_a*domegadt_csts_nonspin[0];
    const double a2 = factor_a*(domegadt_csts_nonspin[1] + nu*(domegadt_csts_nonspin[2]));
    const double a3 = factor_a*(domegadt_csts_nonspin[3] + beta(domegadt_csts_spinorbit[0]/Msquare, domegadt_csts_spinorbit[1]/Msquare, &system));
    const double a4 = factor_a*(domegadt_csts_nonspin[4] + nu*(domegadt_csts_nonspin[5] + nu*(domegadt_csts_nonspin[6])) + sigma(domegadt_csts_spinspin[0]/M_4_nu, domegadt_csts_spinspin[1]/M_4_nu, &system) + tau(domegadt_csts_spinspin[2]/M_4_nu, domegadt_csts_spinspin[3]/M_4_nu, &system));
    const double a5 = factor_a*(domegadt_csts_nonspin[7] + nu*(domegadt_csts_nonspin[8]) + beta((domegadt_csts_spinorbit[2] + nu*(domegadt_csts_spinorbit[3]))/Msquare, (domegadt_csts_spinorbit[4] + nu*(domegadt_csts_spinorbit[5]))/Msquare, &system));
    
    const double a0_2 = a0*a0;
    const double a0_3 = a0_2*a0;
    const double a2_2 = a2*a2;

    //these constants are used in TaylorT2 where domega/dt is expressed as an inverse polynomial
    const double c0 = 1./a0;
    const double c2 = -a2/a0_2;
    const double c3 = -a3/a0_2;
    const double c4 = (a2_2 - a0*a4)/a0_3;
    const double c5 = (2.*a2*a3 - a0*a5)/a0_3;

    system.c_1 /= nu;
    double c_1 = system.c_1;
    system.nu_4 = nu_2*nu_2;
    const double nu_4 = system.nu_4;
    system.c1_2 = c_1*c_1;
    double c1_2 = system.c1_2;
    const double Delta = sqrt((4.*c1_2*one_p_q_sq - 8.*c_1*q*(1.+q)*Seff - 4.*((1.-q_2)*(1.-q_2)*S1_normalized_2-q_2*Seff_2))*(4.*c1_2*q_2*one_p_q_sq - 8.*c_1*q_2*(1.+q)*Seff - 4.*((1.-q_2)*(1.-q_2)*S2_normalized_2-q_2*Seff_2)));
    
    system.constants_u[0] = -c0;
    system.constants_u[1] = 6.*Seff*nu/deltam_over_M/deltam_over_M - 3.*c_1/deltam_over_M/deltam_over_M;
    system.constants_u[2] = 3.*c2/c0 + 0.75*one_p_q_sq/one_m_q_4*(-20.*c1_2*q_2*one_p_q_sq + 2.*(1.-q_2)*(1.-q_2)*(q*(2.+q)*S1_normalized_2 + (1.+2.*q)*S2_normalized_2 -2.*q*Save_square_normalized) + 2.*q_2*(7.+6.*q+7.*q_2)*2.*c_1*Seff - 2.*q_2*(3.+4.*q+3.*q_2)*Seff_2 + q*Delta);
    
    system.Ssqave = 0.5*(A3+A2);
    system.sqrtSsqave = sqrt(system.Ssqave);
    system.c_1 = 0.5*(J_0_norm*J_0_norm - L_0_norm*L_0_norm - Save_square)/L_0_norm/GMsquare_over_c*nu;
    c_1 = system.c_1;
    system.c1_2 = c_1*c_1;
    c1_2 = system.c1_2;
    
    const double Rm = A3 - A2;
    const double Rm_2 = Rm*Rm;
    const double cp = A3*nu_2 - c1_2;
    const double cm = cp-Rm*nu_2;
    const double S0m = S1_normalized_2 - S2_normalized_2;
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
    
    double Omegaz0 = A1t + ad;
    double Omegaz1 = A2t - ad*Seff - ad*hd;
    double Omegaz2 = cd - ad*fd + ad*hd*hd + ad*hd*Seff;
    double Omegaz3 = (ad*fd - cd - ad*hd*hd)*Seff + 2*ad*fd*hd - ad*hd*hd*hd - cd*hd;
    double Omegaz4 = -(2*ad*fd*hd - ad*hd*hd*hd - cd*hd)*Seff + ad*fd*fd - cd*fd + cd*hd*hd - 3*ad*fd*hd*hd + ad*hd*hd*hd*hd;
    double Omegaz5 = -(ad*fd*fd - cd*fd + cd*hd*hd - 3*ad*fd*hd*hd + ad*hd*hd*hd*hd)*Seff + hd*(2*cd*fd - 3*ad*fd*fd - cd*hd*hd + 4*ad*fd*hd*hd - ad*hd*hd*hd*hd);
    
    system.constants_phiz[0] = 3.*Omegaz0*c0;
    system.constants_phiz[1] = 3.*Omegaz1*c0;
    system.constants_phiz[2] = 3.*(Omegaz2*c0 + Omegaz0*c2);
    system.constants_phiz[3] = 3.*(Omegaz0*c3 + Omegaz1*c2 + Omegaz3*c0);
    system.constants_phiz[4] = 3.*(Omegaz0*c4 + Omegaz1*c3 + Omegaz2*c2 + Omegaz4*c0);
    system.constants_phiz[5] = 3.*(c5*Omegaz0 + c4*Omegaz1 + c3*Omegaz2 + c2*Omegaz3 + c0*Omegaz5);

    const double gw = 3./16./nu_2/nu*Rm_2*(c_1 - nu_2*Seff);
    const double gd = gw/dw;
    
    Omegaz5 += Omegaz4*c_1/nu_2 - fd*gd + gd*hd*hd + gd*hd*Seff;
    Omegaz4 += Omegaz3*c_1/nu_2 - gd*hd - gd*Seff;
    Omegaz3 += Omegaz2*c_1/nu_2 + gd;
    Omegaz2 += Omegaz1*c_1/nu_2;
    Omegaz1 += Omegaz0*c_1/nu_2;
    
    system.constants_zeta[0] = 3.*Omegaz0*c0;
    system.constants_zeta[1] = 3.*Omegaz1*c0;
    system.constants_zeta[2] = 3.*(Omegaz2*c0 + Omegaz0*c2);
    system.constants_zeta[3] = 3.*(Omegaz0*c3 + Omegaz1*c2 + Omegaz3*c0);
    system.constants_zeta[4] = 3.*(Omegaz0*c4 + Omegaz1*c3 + Omegaz2*c2 + Omegaz4*c0);
    system.constants_zeta[5] = 3.*(c5*Omegaz0 + c4*Omegaz1 + c3*Omegaz2 + c2*Omegaz3 + c0*Omegaz5);
    
    const double m = sqrt((A2 - A3)/(A1 - A3));
    
    const double B = (S_0_norm*S_0_norm/GMsquare_over_c/GMsquare_over_c-A3)/(A2-A3);
    const double volumeellement = DotProd(CrossProd(L_0,S1_0),S2_0);
    const int sign_num = (volumeellement > 0) - (volumeellement < 0);
    
    if(S1_norm ==0 || S2_norm ==0) system.constant_of_S = 0;
    else{
        if(B < 0. || B > 1.) {
            if(B > 1 && B-1. < 0.00001) system.constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(1.)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0,&system);
            if(B < 0 && B > -0.00001) system.constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(0.)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0,&system);
        }
        else system.constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(B)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0,&system);
    }
    
    system.constants_L[0] = (L_csts_nonspin[0] + nu*L_csts_nonspin[1]);
    system.constants_L[1] = beta(L_csts_spinorbit[0]/Msquare, L_csts_spinorbit[1]/Msquare, &system);
    system.constants_L[2] = (L_csts_nonspin[2] + nu*L_csts_nonspin[3] + nu*nu*L_csts_nonspin[4]);
    system.constants_L[3] = beta((L_csts_spinorbit[2]+L_csts_spinorbit[3]*nu)/Msquare, (L_csts_spinorbit[4]+L_csts_spinorbit[5]*nu)/Msquare, &system);
    system.constants_L[4] = (L_csts_nonspin[5]+L_csts_nonspin[6]*nu +L_csts_nonspin[7]*nu*nu+L_csts_nonspin[8]*nu*nu*nu);
    
    system.flag = 0;
    if(S1_norm == 0 || S2_norm ==0 || (fabs(DotProd(S1_0,Lhat_0)/S1_norm) == 1 && fabs(DotProd(S2_0,Lhat_0)/S2_norm) == 1)) system.flag = 1;
    
    system.phiz_0 = 0.;
    system.phiz_0 = - phiz_of_xi(xi_0,&system);
    system.zeta_0 = 0.;
    system.zeta_0 = - zeta_of_xi(xi_0,&system);
    
    
    return system;
}

/* *********************************************************************************/
/* Internal function that computes the roots of Eq. 22 in arxiv:1703.03967         */
/* *********************************************************************************/

static vector Roots(const double xi, const double J_norm, const sysq *system)
{
    vector out;
    vector coeffs = BCDcoeff(xi, J_norm,system);
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

static vector BCDcoeff(const double xi, const double J_norm, const sysq *system)
{
    vector out;
    const double L_normalized = ((*system).nu)/xi;
    const double J_norm_2 = J_norm*J_norm;
    const double L_square_normalized = L_normalized*L_normalized;
    const double J_square_normalized = J_norm_2/((*system).Mfour_in_L);
    
    out.x = (L_square_normalized+((*system).S1_normalized_2))*((*system).q) + (L_square_normalized+((*system).S2_normalized_2))/((*system).q) + 2.*L_normalized*((*system).Seff) - 2.*J_square_normalized - ((*system).S1_normalized_2) - ((*system).S2_normalized_2);
    
    out.y = (L_square_normalized - J_square_normalized)*(L_square_normalized - J_square_normalized) + 2.*((*system).Seff)*L_normalized*(L_square_normalized - J_square_normalized) - 2.*L_square_normalized*(((*system).S1_normalized_2) - ((*system).S2_normalized_2)*((*system).q))*(1./((*system).q)-1) + 4.*((*system).nu)*((*system).Seff)*((*system).Seff)*L_square_normalized - 2.*((*system).Seff)*L_normalized*(((*system).S1_normalized_2)-((*system).S2_normalized_2))*((*system).deltam_over_M) + 2.*J_square_normalized*(((*system).S1_normalized_2)*((*system).q) - ((*system).S2_normalized_2))*(1./((*system).q)-1);
    
    out.z = -(L_square_normalized - J_square_normalized)*(L_square_normalized - J_square_normalized)*(((*system).S1_normalized_2)*((*system).q) - ((*system).S2_normalized_2))*(1./((*system).q)-1) - 2.*((*system).Seff)*((*system).deltam_over_M)*(((*system).S1_normalized_2) - ((*system).S2_normalized_2))*L_normalized*(L_square_normalized - J_square_normalized) + (((*system).S1_normalized_2) - ((*system).S2_normalized_2))*(((*system).S1_normalized_2) - ((*system).S2_normalized_2))*L_square_normalized*((*system).deltam_over_M)*((*system).deltam_over_M)/((*system).nu);
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of J to Newtonian order            */
/* *********************************************************************************/

static double J_norm_of_xi(const double xi, const sysq *system)
{
    const double L_norm = ((*system).nu)/xi;
    const double constant_1=(*system).c_1/((*system).nu);
    
    return sqrt(L_norm*L_norm + 2.*L_norm*constant_1 + (*system).Ssqave)*((*system).GMsquare_over_c);
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of S                               */
/* *********************************************************************************/

static double S_norm_of_xi(const double xi, const sysq *system)
{
    const double J_norm = J_norm_of_xi(xi,system);
    const vector roots = Roots(xi,J_norm,system);
    double A1, A2, A3;
    double sn, cn, dn;
    
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;
    
    const double m = (A2 - A3)/(A1 - A3);
    const double u = u_of_xi(xi,system)+(*system).constant_of_S;
    
    gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
    
    if((*system).S1_norm_square == 0. || (*system).S2_norm_square == 0.) sn = 0.;
    
    const double S_norm_square_bar = A3 + (A2 - A3)*sn*sn;
    return sqrt(S_norm_square_bar)*((*system).GMsquare_over_c);
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of J to 3PN order                  */
/* *********************************************************************************/

static double J_norm_3PN_of_xi(const double xi, const sysq *system)
{
    const double L_norm = L_norm_3PN_of_xi(xi,system)/((*system).GMsquare_over_c);
    const double constant_1=(*system).c_1/((*system).nu);
    
    return sqrt(L_norm*L_norm + 2.*L_norm*constant_1 + (*system).Ssqave)*((*system).GMsquare_over_c);
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of L to 3PN order                  */
/* *********************************************************************************/

static double L_norm_3PN_of_xi(const double xi, const sysq *system)
{
    const double L_norm_0PN = ((*system).GMsquare_over_c)*((*system).nu)/xi;
    const double xi_2 = xi*xi;
    
    //from 0605140 and Blanchet LRR
    return L_norm_0PN*(1. + xi_2*((*system).constants_L[0] + xi*(*system).constants_L[1] + xi_2*((*system).constants_L[2] + xi*(*system).constants_L[3] + xi_2*((*system).constants_L[4]))));
}

/* *********************************************************************************/
/* Internal function that return a coefficient. I don't really remember            */
/* what it is right now.                                                           */
/* *********************************************************************************/

static vector c(const double xi, const sysq *system)
{
    const double xi_2 = xi*xi;
    const double xi_3 = xi_2*xi;
    const double xi_4 = xi_3*xi;
    const double xi_6 = xi_3*xi_3;
    
    const double J_norm = J_norm_of_xi(xi,system);
    const vector roots = Roots(xi,J_norm,system);
    double A1, A2, A3;
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;
    
    const double J_normalized = J_norm/((*system).GMsquare_over_c);
    const double J_normalized_2 = J_normalized*J_normalized;
    
    vector out;
    
    out.x = -0.75*((J_normalized_2-A3)*(J_normalized_2-A3)*xi_4/((*system).nu) - 4.*((*system).nu)*((*system).Seff)*(J_normalized_2-A3)*xi_3-2.*(J_normalized_2-A3+2*(((*system).S1_normalized_2)-((*system).S2_normalized_2))*((*system).deltam_over_M))*((*system).nu)*xi_2+(4.*((*system).Seff)*xi+1)*((*system).nu)*((*system).nu_2))*J_normalized*xi_2*(((*system).Seff)*xi-1.);
    out.y = 1.5*(A2-A3)*J_normalized*((J_normalized*J_normalized-A3)/((*system).nu)*xi_2-2.*((*system).nu)*((*system).Seff)*xi-((*system).nu))*(((*system).Seff)*xi-1.)*xi_4;
    out.z = -0.75*J_normalized*(((*system).Seff)*xi-1.)*(A3-A2)*(A3-A2)*xi_6/((*system).nu);
    
    return out;
}

/* *********************************************************************************/
/* Internal function that return a coefficient. I don't really remember            */
/* what it is right now.                                                           */
/* *********************************************************************************/

static vector d(const double xi, const sysq *system)
{
    const double J_norm = J_norm_of_xi(xi,system);
    const vector roots = Roots(xi,J_norm,system);
    double A1, A2, A3;
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;
    
    const double L_normalized = ((*system).nu)/xi;
    const double J_normalized = J_norm/((*system).GMsquare_over_c);
    vector out;
    
    out.x = -((L_normalized+J_normalized)*(L_normalized+J_normalized)-A3)*((L_normalized-J_normalized)*(L_normalized-J_normalized)-A3);
    out.y = 2.*(A2-A3)*(J_normalized*J_normalized+L_normalized*L_normalized-A3);
    out.z = -(A3-A2)*(A3-A2);
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns the cosine of the angle between the Newtonian L  */
/* and J                                                                           */
/* *********************************************************************************/

static double costhetaL(const double xi, const sysq *system)
{
    const double L_norm = ((*system).GMsquare_over_c)*((*system).nu)/xi;
    const double J_norm = J_norm_of_xi(xi,system);
    const double S_norm = S_norm_of_xi(xi,system);
    
    return 0.5*(J_norm*J_norm + L_norm*L_norm - S_norm*S_norm)/L_norm/J_norm;
}

/* *********************************************************************************/
/* Internal function that returns the cosine of the angle between the 3PN L and J  */
/* *********************************************************************************/

static double costhetaL_3PN(const double xi, const sysq *system)
{
    const double L_norm = L_norm_3PN_of_xi(xi,system);
    const double J_norm = J_norm_3PN_of_xi(xi,system);
    const double S_norm = S_norm_of_xi(xi,system);
    
    return 0.5*(J_norm*J_norm + L_norm*L_norm - S_norm*S_norm)/L_norm/J_norm;
}

/* *********************************************************************************/
/* Internal function that returns the phase of the magnitude of S                  */
/* *********************************************************************************/

static double u_of_xi(const double xi, const sysq *system)
{
    const double xi_2 = xi*xi;
    return 0.75*((*system).deltam_over_M)*(*system).constants_u[0]/xi_2/xi*(1. + xi*((*system).constants_u[1] + xi*((*system).constants_u[2])));
}

/* *********************************************************************************/
/* Internal function that returns the derivative of the phase of the magnitude of S*/
/* *********************************************************************************/

static double psidot(const double xi,  const sysq *system)
{
    const double J_norm = J_norm_of_xi(xi,system);
    const vector roots = Roots(xi,J_norm,system);
    
    const double A3 = fmax(fmax(roots.x,roots.y),roots.z);
    const double A1 = fmin(fmin(roots.x,roots.y),roots.z);
    
    const double xi_2 = xi*xi;
    const double xi_3 = xi_2*xi;
    const double xi_6 = xi_3*xi_3;
    
    return 0.75/sqrt((*system).nu)*(1.-(*system).Seff*xi)*xi_6*sqrt(A3-A1);
}

/* *********************************************************************************/
/* Internal function that returns phiz                                             */
/* *********************************************************************************/

static double phiz_of_xi(const double xi, const sysq *system)
{
    const double J_norm = J_norm_of_xi(xi,system);
    const double xi_2 = xi*xi;
    const double J_norm_normalized = J_norm/((*system).GMsquare_over_c);
    const double log1 = log(fabs(((*system).c_1) + J_norm_normalized*((*system).nu)+((*system).nu_2)/xi));
    const double log2 = log(fabs(((*system).c_1) + J_norm_normalized*((*system).sqrtSsqave)*xi + ((*system).Ssqave)*xi));
    
    const double phiz0 = J_norm_normalized*(0.5*((*system).c1_2)/((*system).nu_4) - 0.5*((*system).onethird)*((*system).c_1)/xi/((*system).nu_2) - ((*system).onethird)*((*system).Ssqave)/((*system).nu_2) - ((*system).onethird)/xi_2) - 0.5*((*system).c_1)*(((*system).c1_2)/((*system).nu_4) - ((*system).Ssqave)/((*system).nu_2))/((*system).nu)*log1;
    const double phiz1 = -J_norm_normalized*0.5*(((*system).c_1)/((*system).nu_2) + 1./xi) + 0.5*(((*system).c1_2)/((*system).nu_2) - ((*system).Ssqave))/((*system).nu)*log1 ;
    const double phiz2 = -J_norm_normalized + ((*system).sqrtSsqave)*log2 - ((*system).c_1)/((*system).nu)*log1;
    const double phiz3 = J_norm_normalized*xi -((*system).nu)*log1 + ((*system).c_1)/((*system).sqrtSsqave)*log2;
    const double phiz4 = J_norm_normalized*0.5*xi*(((*system).c_1)/((*system).Ssqave) + xi) - 0.5*(((*system).c1_2)/((*system).Ssqave) - ((*system).nu_2))/((*system).sqrtSsqave)*log2;
    const double phiz5 = J_norm_normalized*xi*(-0.5*((*system).c1_2)/((*system).Ssqave)/((*system).Ssqave) + 0.5*((*system).onethird)*((*system).c_1)*xi/((*system).Ssqave) + ((*system).onethird)*(xi_2 + ((*system).nu_2)/((*system).Ssqave))) + 0.5*((*system).c_1)*(((*system).c1_2)/((*system).Ssqave) - ((*system).nu_2))/((*system).Ssqave)/((*system).sqrtSsqave)*log2;
    
    double phMS;
    if((*system).flag == 0) phMS = phiz_MS_corr(xi,system);
    else phMS = 0;
    
    return phiz0*(*system).constants_phiz[0] + phiz1*(*system).constants_phiz[1] + phiz2*(*system).constants_phiz[2] + phiz3*(*system).constants_phiz[3] + phiz4*(*system).constants_phiz[4] + phiz5*(*system).constants_phiz[5] + (*system).phiz_0 + phMS;
}

/* *********************************************************************************/
/* Internal function that returns zeta                                             */
/* *********************************************************************************/

static double zeta_of_xi(const double xi, const sysq *system)
{
    const double logxi = log(xi);
    const double xi_2 = xi*xi;
    const double xi_3 = xi_2*xi;
    
    double zMS;
    if((*system).flag == 0) zMS = zeta_MS_corr(xi,system);
    else zMS = 0;
    
    return ((*system).nu)*(-(*system).onethird*(*system).constants_zeta[0]/xi_3 - 0.5*(*system).constants_zeta[1]/xi_2 - (*system).constants_zeta[2]/xi + (*system).constants_zeta[3]*logxi + (*system).constants_zeta[4]*xi  + 0.5*(*system).constants_zeta[5]*xi_2) + (*system).zeta_0 + zMS;
}

/* *********************************************************************************/
/* Internal function that returns the MS corrections of phiz                       */
/* *********************************************************************************/

static double phiz_MS_corr(const double xi, const sysq *system)
{
    const double c0 = c(xi,system).x;
    const double c2 = c(xi,system).y;
    const double c4 = c(xi,system).z;
    const double d0 = d(xi,system).x;
    const double d2 = d(xi,system).y;
    const double d4 = d(xi,system).z;
    
    const double sqt = sqrt(fabs(d2*d2-4.*d0*d4));
    
    const double n3 = (2.*(d0+d2+d4)/(2.*d0+d2+sqt));
    const double n4 = ((2*d0+d2+sqt)/(2.*d0));
    const double sqtn3 = sqrt(fabs(n3));
    const double sqtn4 = sqrt(fabs(n4));
    
    const double u = u_of_xi(xi,system)+(*system).constant_of_S;
    
    double L3, L4;
    if(n3 == 1) L3 = 0;
    else L3 = -fabs((c4*d0*(2.*d0+d2+sqt)-c2*d0*(d2+2.*d4-sqt)-c0*(2.*d0*d4-(d2+d4)*(d2-sqt)))/(2.*d0*(d0+d2+d4)*sqt))*(sqtn3/(n3-1.)*(atan(tan(u)) - atan(sqtn3*tan(u))))/psidot(xi,system);
    if(n4 == 1) L4 = 0;
    else L4 = -fabs((-c4*d0*(2.*d0+d2-sqt)+c2*d0*(d2+2.*d4+sqt)-c0*(-2.*d0*d4+(d2+d4)*(d2+sqt))))/(2.*d0*(d0+d2+d4)*sqt)*(sqtn4/(n4-1.)*(atan(tan(u)) - atan(sqtn4*tan(u))))/psidot(xi,system);
    
    double corr = L3+L4;
    if (corr != corr) corr=0;
    return corr;
}

/* *********************************************************************************/
/* Internal function that returns the MS corrections of zeta                       */
/* *********************************************************************************/

static double zeta_MS_corr(const double xi, const sysq *system)
{
    const double J_norm = J_norm_of_xi(xi,system);
    const double c0 = c(xi,system).x;
    const double c2 = c(xi,system).y;
    const double c4 = c(xi,system).z;
    const double d0 = d(xi,system).x;
    const double d2 = d(xi,system).y;
    const double d4 = d(xi,system).z;
    
    const vector roots = Roots(xi,J_norm,system);
    double A1, A2, A3;
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;
    
    const double J_normalized = J_norm/((*system).GMsquare_over_c);
    const double L_normalized = ((*system).nu)/xi;
    
    const double Aa = 0.5*(J_normalized/L_normalized + L_normalized/J_normalized - A3/J_normalized/L_normalized);
    const double Bb = 0.5*(A3 - A2)/L_normalized/J_normalized;
    const double sqt = sqrt(d2*d2-4.*d0*d4);
    const double n3 = 2.*(d0+d2+d4)/(2.*d0+d2+sqt);
    const double n4 = (2*d0+d2+sqt)/(2.*d0);
    const double sqtn3 = sqrt(fabs(n3));
    const double sqtn4 = sqrt(fabs(n4));
    
    const double u = u_of_xi(xi,system)+(*system).constant_of_S;
    
    double L3, L4;
    if(n3 == 1) L3 = 0;
    else L3 = -fabs((c4*d0*(2.*d0+d2+sqt)-c2*d0*(d2+2.*d4-sqt)-c0*(2.*d0*d4-(d2+d4)*(d2-sqt)))/(2.*d0*(d0+d2+d4)*sqt))*(sqtn3/(n3-1.)*(atan(tan(u)) - atan(sqtn3*tan(u))))/psidot(xi,system);
    if(n4 == 1) L4 = 0;
    else L4 = -fabs((-c4*d0*(2.*d0+d2-sqt)+c2*d0*(d2+2.*d4+sqt)-c0*(-2.*d0*d4+(d2+d4)*(d2+sqt))))/(2.*d0*(d0+d2+d4)*sqt)*(sqtn4/(n4-1.)*(atan(tan(u)) - atan(sqtn4*tan(u))))/psidot(xi,system);
    
    double corr = Aa*phiz_MS_corr(xi,system) + 2.*Bb*d0*(L3/(sqt-d2)-L4/(sqt+d2));
    if (corr != corr) corr=0;
    return corr;
}

/* *********************************************************************************/
/* Internal function that computes the spin-orbit couplings                        */
/* *********************************************************************************/

static double beta(const double a, const double b, const sysq *system)
{
    double coeff = LAL_G_SI/LAL_C_SI;
    return ((((*system).dot1)*(a + b*((*system).q))) + (((*system).dot2)*(a + b/((*system).q))))/coeff;
}

/* *********************************************************************************/
/* Internal function that computes the spin-spin couplings                         */
/* *********************************************************************************/

static double sigma(const double a, const double b, const sysq *system)
{
    double coeff = LAL_G_SI/LAL_C_SI*LAL_G_SI/LAL_C_SI;
    return (a*((*system).dot12) - b*((*system).dot1)*((*system).dot2))/coeff;
}

/* *********************************************************************************/
/* Internal function that computes the spin-spin couplings                         */
/* *********************************************************************************/

static double tau(const double a, const double b, const sysq *system)
{
    double coeff = LAL_G_SI/LAL_C_SI*LAL_G_SI/LAL_C_SI;
    return (((*system).q)*((*system).S1_norm_square*a - b*((*system).dot1)*((*system).dot1)) + ((*system).S2_norm_square*a - b*((*system).dot2)*((*system).dot2))/((*system).q))/coeff;
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



