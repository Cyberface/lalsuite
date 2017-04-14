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

/* **********************************************************************************/
/* Internal function that initalizes all constants needed for the precession angles.*/
/* It needs to be called only once at the beginning.                                */
/* **********************************************************************************/

static void InitializePrecession()
{
    onethird = 1./3.;
    onesixth = 1./6.;
    twotothe_onethird = pow(2.,onethird);
    sqrt_three = sqrt(3.);
    pi = 3.141592653589793238463;
    pi_over_two = 0.5*pi;
    pi_over_four = 0.25*pi;
    pisquare = pi*pi;
    pi_four = pisquare*pisquare;
    log_2 = log(2.);
    log_2_2 = log_2*log_2;
    log_3 = log(3.);
    log_3_2 = log_3*log_3;
    log_5 = log(5.);
    log_7 = log(7.);
    log_3_over_2 = log(1.5);
    twopi = 2.*pi;
    gamma_Euler = 0.57721566490153286060;
    Light_vel = 299792458.;
    Newton_c = 6.67384e-11;
    SOLAR_MASS = 1.9891e30;
    G_over_c = Newton_c/Light_vel;
    G_over_csquare = G_over_c/Light_vel;
    G_over_cthree = G_over_csquare/Light_vel;
    cthree_over_G = 1./G_over_cthree;
    csix_over_Gsquare = cthree_over_G*cthree_over_G;
    cseven_over_Gthree = csix_over_Gsquare/G_over_c;
    ceight_over_Gfour = cseven_over_Gthree/G_over_c;
    
    double domegadt_csts_nonspin_internal[17] = {
        csix_over_Gsquare*96./5., // 0 0
        -csix_over_Gsquare*1486./35., // 1 2
        -csix_over_Gsquare*264./5.,
        csix_over_Gsquare*384.*pi/5., // 3 3
        csix_over_Gsquare*34103./945., // 4 4
        csix_over_Gsquare*13661./105.,
        csix_over_Gsquare*944./15.,
        csix_over_Gsquare*pi*(-4159./35.), // 7 5
        csix_over_Gsquare*pi*(-2268./5.),
        csix_over_Gsquare*(16447322263./7276500. + pisquare*512./5. - log_2*109568./175. - gamma_Euler*54784./175.), // 9 6
        csix_over_Gsquare*(-56198689./11340. + pisquare*902./5.),
        csix_over_Gsquare*1623./140.,
        -csix_over_Gsquare*1121./27.,
        -csix_over_Gsquare*54784./525., // 13 6l
        -csix_over_Gsquare*pi*883./42., // 14 7
        csix_over_Gsquare*pi*71735./63.,
        csix_over_Gsquare*pi*73196./63.,
    };
    
    double domegadt_csts_spinorbit_internal[18] = {
        -cseven_over_Gthree*904./5.,
        -cseven_over_Gthree*120.,
        -cseven_over_Gthree*62638./105.,
        cseven_over_Gthree*4636./5.,
        -cseven_over_Gthree*6472./35.,
        cseven_over_Gthree*3372./5.,
        -cseven_over_Gthree*pi*720.,
        -cseven_over_Gthree*pi*2416./5.,
        -cseven_over_Gthree*208520./63.,
        cseven_over_Gthree*796069./105.,
        -cseven_over_Gthree*100019./45.,
        -cseven_over_Gthree*1195759./945.,
        cseven_over_Gthree*514046./105.,
        -cseven_over_Gthree*8709./5.,
        -cseven_over_Gthree*pi*307708./105.,
        cseven_over_Gthree*pi*44011./7.,
        -cseven_over_Gthree*pi*7992./7.,
        cseven_over_Gthree*pi*151449./35.
    };
    
    
    double domegadt_csts_spinspin_internal[4] = {
        -ceight_over_Gfour*494./5.,
        -ceight_over_Gfour*1442./5.,
        -ceight_over_Gfour*233./5.,
        -ceight_over_Gfour*719./5.
    };
    
    double L_csts_nonspin_internal[9] = {
        3./2.,
        1./6.,
        27./8.,
        -19./8.,
        1./24.,
        135./16.,
        -6889/144.+ 41./24.*pisquare,
        31./24.,
        7./1296.
    };
    
    double L_csts_spinorbit_internal[6] = {
        -14./6.,
        -3./2.,
        -11./2.,
        133./72.,
        -33./8.,
        7./4.
    };
    
    
    memcpy(domegadt_csts_spinorbit, domegadt_csts_spinorbit_internal, 18*sizeof(double));
    memcpy(domegadt_csts_nonspin, domegadt_csts_nonspin_internal, 17*sizeof(double));
    memcpy(domegadt_csts_spinspin, domegadt_csts_spinspin_internal, 4*sizeof(double));
    memcpy(L_csts_spinorbit, L_csts_spinorbit_internal, 6*sizeof(double));
    memcpy(L_csts_nonspin, L_csts_nonspin_internal, 9*sizeof(double));
    
}

/* *********************************************************************************/
/* XLAL function that initalizes all constants that depend on the particular       */
/* system. It needs to be called once when the system parameters have been defined.*/
/* *********************************************************************************/

static void InitializeSystem()
{
    if(ch1 < 0.){
        ch1 = -ch1;
        mu1 = -mu1;
        if(ph1 >= pi) ph1 -= pi;
        else ph1 += pi;
    }
    if(ch2 < 0.){
        ch2 = -ch2;
        mu2 = -mu2;
        if(ph2 >= pi) ph2 -= pi;
        else ph2 += pi;
    }

    th1 = acos(mu1);
    th2 = acos(mu2);
    thl = acos(mul);
    
    M = m1 + m2;
    mu = m1 * m2 / M;
    nu = mu / M;
    nu_2 = nu*nu;
    Msquare = M*M;
    Mthree = Msquare*M;
    m1square = m1*m1;
    m2square = m2*m2;
    q = m2/m1;
    deltam = (m1 - m2);
    deltam_over_M = deltam/M;
    one_over_m1_square = 1. / m1square;
    one_over_m2_square = 1. / m2square;
    double GM_over_cthree = G_over_cthree*M;
    double GM_over_csquare = G_over_csquare*M;
    GMsquare_over_c = GM_over_csquare*M*Light_vel;
    double piGM_over_cthree = pi*GM_over_cthree;
    twopiGM_over_cthree = 2*piGM_over_cthree;
    
    S1_norm = fabs(ch1) * m1square * G_over_c;
    S2_norm = fabs(ch2) * m2square * G_over_c;
    S1_norm_square = S1_norm*S1_norm;
    S2_norm_square = S2_norm*S2_norm;

    xi_0 = pow(piGM_over_cthree*f_0, onethird);
    Lhat_0 = CreateSphere(1.,thl,phl);
    S1_0 = CreateSphere(S1_norm,th1,ph1);
    S2_0 = CreateSphere(S2_norm,th2,ph2);
    L_0 = ScalarProd(GMsquare_over_c*nu/xi_0,Lhat_0);

    Seff = (DotProd(S1_0,Lhat_0)/m1 + DotProd(S2_0,Lhat_0)/m2)/M/G_over_c;
    
    vector S_0 = Sum(S1_0,S2_0);
    vector J_0 = Sum(L_0,S_0);
    
    double S_0_norm = Norm(S_0);
    double J_0_norm = Norm(J_0);
    double L_0_norm = Norm(L_0);
    
    vector roots = Roots(xi_0,J_0_norm);
    double A1, A2, A3;

    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;

    double q_2 = q*q;
    double one_m_q_sq = (1.-q)*(1.-q);
    double one_m_q_4 = one_m_q_sq*one_m_q_sq;
    double one_p_q_sq = (1.+q)*(1.+q);
    double S1_normalized = S1_norm/GMsquare_over_c;
    double S2_normalized = S2_norm/GMsquare_over_c;
    double Save_square = 0.5*(A3+A2)*GMsquare_over_c*GMsquare_over_c;
    double Save_square_normalized = 0.5*(A3+A2);
    double S1_normalized_2 = S1_normalized*S1_normalized;
    double S2_normalized_2 = S2_normalized*S2_normalized;
    double Seff_2 = Seff*Seff;

    c_1 = 0.5*(J_0_norm*J_0_norm - L_0_norm*L_0_norm - Save_square)/L_0_norm/GMsquare_over_c*nu;

    double factor_a = nu/csix_over_Gsquare;
    double a0, a2, a3, a4, a5, a6, a6log, a7;
    double c0, c2, c3, c4, c5, c6, c6log, c7;
    
    
    //computed with initial spin couplings, they're not exactly accurate for generic precession, but the correction should be 4PN
    //these constants are used in TaylorT1 where domega/dt is expressed as a polynomial
    a0 = factor_a*domegadt_csts_nonspin[0];
    a2 = factor_a*(domegadt_csts_nonspin[1] + nu*(domegadt_csts_nonspin[2]));
    a3 = factor_a*(domegadt_csts_nonspin[3] + beta(domegadt_csts_spinorbit[0]/Msquare, domegadt_csts_spinorbit[1]*nu));
    a4 = factor_a*(domegadt_csts_nonspin[4] + nu*(domegadt_csts_nonspin[5] + nu*(domegadt_csts_nonspin[6])) + sigma(domegadt_csts_spinspin[0]/(Msquare*Msquare*nu), domegadt_csts_spinspin[1]/(Msquare*Msquare*nu)) + tau(domegadt_csts_spinspin[2]/Msquare, domegadt_csts_spinspin[3]/Msquare));
    a5 = factor_a*(domegadt_csts_nonspin[7] + nu*(domegadt_csts_nonspin[8]) + beta((domegadt_csts_spinorbit[2] + nu*(domegadt_csts_spinorbit[3]))/Msquare, (domegadt_csts_spinorbit[4] + nu*(domegadt_csts_spinorbit[5]))*nu));
    a6 = factor_a*(domegadt_csts_nonspin[9] + nu*(domegadt_csts_nonspin[10] + nu*(domegadt_csts_nonspin[11] + nu*(domegadt_csts_nonspin[12]))) + beta(domegadt_csts_spinorbit[6]/Msquare, domegadt_csts_spinorbit[7]*nu));
    a6log = factor_a*domegadt_csts_nonspin[13];
    a7 = factor_a*(domegadt_csts_nonspin[14] + nu*(domegadt_csts_nonspin[15] + nu*(domegadt_csts_nonspin[16])) + beta((domegadt_csts_spinorbit[8] + nu*(domegadt_csts_spinorbit[9] + nu*(domegadt_csts_spinorbit[10])))/Msquare, (domegadt_csts_spinorbit[11] + nu*(domegadt_csts_spinorbit[12] + nu*(domegadt_csts_spinorbit[13])))*nu));
    
    double a0_2 = a0*a0;
    double a0_3 = a0_2*a0;
    double a0_4 = a0_3*a0;
    double a2_2 = a2*a2;
    double a2_3 = a2_2*a2;
    double a3_2 = a3*a3;

    //these constants are used in TaylorT2 where domega/dt is expressed as an inverse polynomial
    c0 = 1./a0;
    c2 = -a2/a0_2;
    c3 = -a3/a0_2;
    c4 = (a2_2 - a0*a4)/a0_3;
    c5 = (2.*a2*a3 - a0*a5)/a0_3;
    c6 = (-a2_3 + a0*a3_2 + 2.*a0*a2*a4 - a0_2*a6)/a0_4;
    c6log = -a6log/a0_2;
    c7 = (-3.*a2_2*a3 + 2.*a0*a3*a4 + 2.*a0*a2*a5 - a0_2*a7)/a0_4;

    c_1 /= nu;
    double nu_4 = nu_2*nu_2;
    double c1_2 = c_1*c_1;
    double c1_3 = c1_2*c_1;
    double c1_4 = c1_3*c_1;
    double Delta = sqrt((4.*c1_2*one_p_q_sq - 8.*c_1*q*(1.+q)*Seff - 4.*((1.-q_2)*(1.-q_2)*S1_normalized_2-q_2*Seff_2))*(4.*c1_2*q_2*one_p_q_sq - 8.*c_1*q_2*(1.+q)*Seff - 4.*((1.-q_2)*(1.-q_2)*S2_normalized_2-q_2*Seff_2)));
    
    constants_u[0] = -c0;
    constants_u[1] = 6.*Seff*nu/deltam_over_M/deltam_over_M - 3.*c_1/deltam_over_M/deltam_over_M;
    constants_u[2] = 3.*c2/c0 + 0.75*one_p_q_sq/one_m_q_4*(-20.*c1_2*q_2*one_p_q_sq + 2.*(1.-q_2)*(1.-q_2)*(q*(2.+q)*S1_normalized_2 + (1.+2.*q)*S2_normalized_2 -2.*q*Save_square_normalized) + 2.*q_2*(7.+6.*q+7.*q_2)*2.*c_1*Seff - 2.*q_2*(3.+4.*q+3.*q_2)*Seff_2 + q*Delta);
    
    Ssqave = 0.5*(A3+A2);
    c_1 = 0.5*(J_0_norm*J_0_norm - L_0_norm*L_0_norm - Save_square)/L_0_norm/GMsquare_over_c*nu;
    c1_2 = c_1*c_1;
    c1_3 = c1_2*c_1;
    c1_4 = c1_3*c_1;
    
    double Rm = A3 - A2;
    double Rm_2 = Rm*Rm;
    double cp = A3*nu_2 - c1_2;
    double cm = cp-Rm*nu_2;
    double S0m = S1_normalized_2 - S2_normalized_2;
    double cpcm = fabs(cp*cm);
    double sqrt_cpcm = sqrt(cpcm);
    
    double A1t = 0.5+0.75/nu;//xi^6
    double A2t = -0.75*Seff/nu;//xi^7
    double A1ave = (cp-sqrt_cpcm)/nu_2 ;
    double Bave = -0.5*Rm*sqrt_cpcm/nu_2 - cp/nu_4*(sqrt_cpcm-cp);
    
    double aw = (-3.*(1. + q)/q*(2.*(1. + q)*nu_2*Seff*c_1 - (1. + q)*c1_2 + (1. - q)*nu_2*S0m));
    double cw = 3./32./nu*Rm_2;
    double dw = 4.*cp - 4.*A1ave*nu_2;
    double hw = -2*(2*A1ave - Rm)*c_1;
    double fw = Rm*A1ave-Bave-0.25*Rm_2;
    
    double ad = aw/dw;
    double hd = hw/dw;
    double cd = cw/dw;
    double fd = fw/dw;
    
    double Omegaz0 = A1t + ad;
    double Omegaz1 = A2t - ad*Seff - ad*hd;
    double Omegaz2 = cd - ad*fd + ad*hd*hd + ad*hd*Seff;
    double Omegaz3 = (ad*fd - cd - ad*hd*hd)*Seff + 2*ad*fd*hd - ad*hd*hd*hd - cd*hd;
    double Omegaz4 = -(2*ad*fd*hd - ad*hd*hd*hd - cd*hd)*Seff + ad*fd*fd - cd*fd + cd*hd*hd - 3*ad*fd*hd*hd + ad*hd*hd*hd*hd;
    double Omegaz5 = -(ad*fd*fd - cd*fd + cd*hd*hd - 3*ad*fd*hd*hd + ad*hd*hd*hd*hd)*Seff + hd*(2*cd*fd - 3*ad*fd*fd - cd*hd*hd + 4*ad*fd*hd*hd - ad*hd*hd*hd*hd);
    
    constants_phiz[0] = 3.*Omegaz0*c0;
    constants_phiz[1] = 3.*Omegaz1*c0;
    constants_phiz[2] = 3.*(Omegaz2*c0 + Omegaz0*c2);
    constants_phiz[3] = 3.*(Omegaz0*c3 + Omegaz1*c2 + Omegaz3*c0);
    constants_phiz[4] = 3.*(Omegaz0*c4 + Omegaz1*c3 + Omegaz2*c2 + Omegaz4*c0);
    constants_phiz[5] = 3.*(c5*Omegaz0 + c4*Omegaz1 + c3*Omegaz2 + c2*Omegaz3 + c0*Omegaz5);

    double gw = 3./16./nu_2/nu*Rm_2*(c_1 - nu_2*Seff);
    double gd = gw/dw;
    
    Omegaz5 += Omegaz4*c_1/nu_2 - fd*gd + gd*hd*hd + gd*hd*Seff;
    Omegaz4 += Omegaz3*c_1/nu_2 - gd*hd - gd*Seff;
    Omegaz3 += Omegaz2*c_1/nu_2 + gd;
    Omegaz2 += Omegaz1*c_1/nu_2;
    Omegaz1 += Omegaz0*c_1/nu_2;
    
    constants_zeta[0] = 3.*Omegaz0*c0;
    constants_zeta[1] = 3.*Omegaz1*c0;
    constants_zeta[2] = 3.*(Omegaz2*c0 + Omegaz0*c2);
    constants_zeta[3] = 3.*(Omegaz0*c3 + Omegaz1*c2 + Omegaz3*c0);
    constants_zeta[4] = 3.*(Omegaz0*c4 + Omegaz1*c3 + Omegaz2*c2 + Omegaz4*c0);
    constants_zeta[5] = 3.*(c5*Omegaz0 + c4*Omegaz1 + c3*Omegaz2 + c2*Omegaz3 + c0*Omegaz5);
    
    double m = sqrt((A2 - A3)/(A1 - A3));
    
    double B = (S_0_norm*S_0_norm/GMsquare_over_c/GMsquare_over_c-A3)/(A2-A3);
    double volumeellement = DotProd(CrossProd(L_0,S1_0),S2_0);
    double sign_num = (volumeellement > 0) - (volumeellement < 0);
    
    if(S1_norm ==0 || S2_norm ==0) constant_of_S = 0;
    else{
        if(B < 0. || B > 1.) {
            //std::cout<<std::setprecision(30)<<B<<"\n";
            if(B > 1 && B-1. < 0.00001) constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(1.)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0);
            if(B < 0 && B > -0.00001) constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(0.)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0);
        }
        else constant_of_S = gsl_sf_ellint_F(asin(sign_num*sqrt(B)), m, GSL_PREC_DOUBLE) - u_of_xi(xi_0);
    }
    
    phiz_0 = 0.;
    phiz_0 = - phiz_of_xi(xi_0);
    zeta_0 = 0.;
    zeta_0 = - zeta_of_xi(xi_0);
    
    constants_L[0] = (L_csts_nonspin[0] + nu*L_csts_nonspin[1]);
    constants_L[1] = beta(L_csts_spinorbit[0]/GMsquare_over_c, L_csts_spinorbit[1]*nu/G_over_c);
    constants_L[2] = (L_csts_nonspin[2] + nu*L_csts_nonspin[3] + nu*nu*L_csts_nonspin[4]);
    constants_L[3] = beta((L_csts_spinorbit[2]+L_csts_spinorbit[3]*nu)/GMsquare_over_c, (L_csts_spinorbit[4]+L_csts_spinorbit[5]*nu)*nu/G_over_c);
    constants_L[4] = (L_csts_nonspin[5]+L_csts_nonspin[6]*nu +L_csts_nonspin[7]*nu*nu+L_csts_nonspin[8]*nu*nu*nu);
}

/* *********************************************************************************/
/* Internal function that computes the roots of Eq. 22 in arxiv:1703.03967         */
/* *********************************************************************************/

static vector Roots(double xi, double J_norm)
{
    vector out;
    vector coeffs = BCDcoeff(xi, J_norm);
    double B_2 = coeffs.x*coeffs.x;
    double B_3 = B_2*coeffs.x;
    double B_C =  coeffs.x*coeffs.y;
    
    double p = coeffs.y - onethird*B_2;
    double qc = 2./27.*B_3 - onethird*B_C + coeffs.z;
    double sqrtarg = sqrt(-3./p);
    double acosarg = 1.5*qc/p*sqrtarg;
    if((acosarg + 1) < 0.00000000001) acosarg = -1;
    
    if(p < 0 && acosarg <= 1 && acosarg >= -1.){
        out.x = 2./sqrtarg*cos(onethird*acos(acosarg)) - onethird*coeffs.x;
        out.y = 2./sqrtarg*cos(onethird*acos(acosarg) - twopi*onethird) - onethird*coeffs.x;
        out.z = 2./sqrtarg*cos(onethird*acos(acosarg) - 2.*twopi*onethird) - onethird*coeffs.x;
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

static vector BCDcoeff(double xi, double J_norm)
{
    vector out;
    double L_norm = GMsquare_over_c*nu/xi;
    double L_norm_2 = L_norm*L_norm;
    double J_norm_2 = J_norm*J_norm;
    double Mfour_in_L = GMsquare_over_c*GMsquare_over_c;
    double S1_square_normalized = S1_norm_square/Mfour_in_L;
    double S2_square_normalized = S2_norm_square/Mfour_in_L;
    double L_square_normalized = L_norm_2/Mfour_in_L;
    double J_square_normalized = J_norm_2/Mfour_in_L;
    double L_normalized = L_norm/GMsquare_over_c;
    
    out.x = (L_square_normalized+S1_square_normalized)*q + (L_square_normalized+S2_square_normalized)/q + 2.*L_normalized*Seff - 2.*J_square_normalized - S1_square_normalized - S2_square_normalized;
    out.y = (L_square_normalized - J_square_normalized)*(L_square_normalized - J_square_normalized) + 2.*Seff*L_normalized*(L_square_normalized - J_square_normalized) - 2.*L_square_normalized*(S1_square_normalized - S2_square_normalized*q)*deltam/m2 + 4.*nu*Seff*Seff*L_square_normalized - 2.*Seff*L_normalized*(S1_square_normalized-S2_square_normalized)*deltam_over_M + 2.*J_square_normalized*(S1_square_normalized*q - S2_square_normalized)*deltam/m2;
    out.z = -(L_square_normalized - J_square_normalized)*(L_square_normalized - J_square_normalized)*(S1_square_normalized*q - S2_square_normalized)*deltam/m2 - 2.*Seff*deltam_over_M*(S1_square_normalized - S2_square_normalized)*L_normalized*(L_square_normalized - J_square_normalized) + (S1_square_normalized - S2_square_normalized)*(S1_square_normalized - S2_square_normalized)*L_square_normalized*deltam_over_M*deltam_over_M/nu;
    
    return out;
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of J to Newtonian order            */
/* *********************************************************************************/

static double J_norm_of_xi(double xi)
{
    double L_norm = nu/xi;
    double constant_1=c_1/nu;
    
    return sqrt(L_norm*L_norm + 2.*L_norm*constant_1 + Ssqave)*GMsquare_over_c;
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of S                               */
/* *********************************************************************************/

static double S_norm_of_xi(double xi)
{
    double J_norm = J_norm_of_xi(xi);
    vector roots = Roots(xi,J_norm);
    double A1, A2, A3;
    double sn, cn, dn;
    double u, m;
    
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;
    
    m = (A2 - A3)/(A1 - A3);
    u = u_of_xi(xi)+constant_of_S;
    
    gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
    
    if(S1_norm == 0. || S2_norm == 0.) sn = 0.;
    
    double S_norm_square_bar = A3 + (A2 - A3)*sn*sn;
    return sqrt(S_norm_square_bar)*GMsquare_over_c;
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of J to 3PN order                  */
/* *********************************************************************************/

static double J_norm_3PN_of_xi(double xi)
{
    double L_norm = L_norm_3PN_of_xi(xi)/GMsquare_over_c;
    double constant_1=c_1/nu;
    
    return sqrt(L_norm*L_norm + 2.*L_norm*constant_1 + Ssqave)*GMsquare_over_c;
}

/* *********************************************************************************/
/* Internal function that returns the magnitude of L to 3PN order                  */
/* *********************************************************************************/

static double L_norm_3PN_of_xi(double xi)
{
    double L_norm_0PN = GMsquare_over_c*nu/xi;
    double xi_2 = xi*xi;
    
    //from 0605140 and Blanchet LRR
    return L_norm_0PN*(1. + xi_2*(constants_L[0] + xi*constants_L[1] + xi_2*(constants_L[2] + xi*constants_L[3] + xi_2*(constants_L[4]))));
}

/* *********************************************************************************/
/* Internal function that return a coefficient. I don't really remember            */
/* what it is right now.                                                           */
/* *********************************************************************************/

static vector c(double xi)
{
    double xi_2 = xi*xi;
    double xi_3 = xi_2*xi;
    double xi_4 = xi_3*xi;
    double xi_6 = xi_3*xi_3;
    
    double J_norm = J_norm_of_xi(xi);
    vector roots = Roots(xi,J_norm);
    double A1, A2, A3;
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;
    
    double J_normalized = J_norm/GMsquare_over_c;
    double J_normalized_2 = J_normalized*J_normalized;
    double Mfour_in_L = GMsquare_over_c*GMsquare_over_c;
    double S1_square_normalized = S1_norm_square/Mfour_in_L;
    double S2_square_normalized = S2_norm_square/Mfour_in_L;
    
    vector out;
    
    out.x = -0.75*((J_normalized_2-A3)*(J_normalized_2-A3)*xi_4/nu - 4.*nu*Seff*(J_normalized_2-A3)*xi_3-2.*(J_normalized_2-A3+2*(S1_square_normalized-S2_square_normalized)*deltam_over_M)*nu*xi_2+(4.*Seff*xi+1)*nu*nu*nu)*J_normalized*xi_2*(Seff*xi-1.);
    out.y = 1.5*(A2-A3)*J_normalized*((J_normalized*J_normalized-A3)/nu*xi_2-2.*nu*Seff*xi-nu)*(Seff*xi-1.)*xi_4;
    out.z = -0.75*J_normalized*(Seff*xi-1.)*(A3-A2)*(A3-A2)*xi_6/nu;
    
    return out;
}

/* *********************************************************************************/
/* Internal function that return a coefficient. I don't really remember            */
/* what it is right now.                                                           */
/* *********************************************************************************/

static vector d(double xi)
{
    double J_norm = J_norm_of_xi(xi);
    vector roots = Roots(xi,J_norm);
    double A1, A2, A3;
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;
    
    double L_normalized = nu/xi;
    double J_normalized = J_norm/GMsquare_over_c;
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

static double costhetaL(double xi)
{
    double L_norm = GMsquare_over_c*nu/xi;
    double J_norm = J_norm_of_xi(xi);
    double S_norm = S_norm_of_xi(xi);
    
    return 0.5*(J_norm*J_norm + L_norm*L_norm - S_norm*S_norm)/L_norm/J_norm;
}

/* *********************************************************************************/
/* Internal function that returns the cosine of the angle between the 3PN L and J  */
/* *********************************************************************************/

static double costhetaL_3PN(double xi)
{
    double L_norm = L_norm_3PN_of_xi(xi);
    double J_norm = J_norm_3PN_of_xi(xi);
    double S_norm = S_norm_of_xi(xi);
    
    return 0.5*(J_norm*J_norm + L_norm*L_norm - S_norm*S_norm)/L_norm/J_norm;
}

/* *********************************************************************************/
/* Internal function that returns the phase of the magnitude of S                  */
/* *********************************************************************************/

static double u_of_xi(double xi)
{
    double xi_square = xi*xi;
    return 0.75*deltam_over_M*constants_u[0]/xi_square/xi*(1. + xi*(constants_u[1] + xi*(constants_u[2])));
}

/* *********************************************************************************/
/* Internal function that returns the derivative of the phase of the magnitude of S*/
/* *********************************************************************************/

static double psidot(double xi)
{
    double J_norm = J_norm_of_xi(xi);
    vector roots = Roots(xi,J_norm);
    double A1, A3;
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    
    double xi_2 = xi*xi;
    double xi_3 = xi_2*xi;
    double xi_6 = xi_3*xi_3;
    
    return 0.75/sqrt(nu)*(1.-Seff*xi)*xi_6*sqrt(A3-A1);
}

/* *********************************************************************************/
/* Internal function that returns phiz                                             */
/* *********************************************************************************/

static double phiz_of_xi(double xi)
{
    double J_norm = J_norm_of_xi(xi);
    double xi_2 = xi*xi;
    double J_norm_normalized = J_norm/GMsquare_over_c;
    double log1 = log(fabs(c_1 + J_norm_normalized*nu+nu_2/xi));
    double log2 = log(fabs(c_1 + J_norm_normalized*sqrt(Ssqave)*xi + Ssqave*xi));
    double c1_2 = c_1*c_1;
    double nu_4 = nu_2*nu_2;
    
    double phiz0 = J_norm_normalized*(0.5*c1_2/nu_4 - onesixth*c_1/xi/nu_2 - onethird*Ssqave/nu_2 - onethird/xi_2) - 0.5*c_1*(c1_2/nu_4 - Ssqave/nu_2)/nu*log1;
    double phiz1 = -J_norm_normalized*0.5*(c_1/nu_2 + 1./xi) + 0.5*(c1_2/nu_2 - Ssqave)/nu*log1 ;
    double phiz2 = -J_norm_normalized + sqrt(Ssqave)*log2 - c_1/nu*log1;
    double phiz3 = J_norm_normalized*xi -nu*log1 + c_1/sqrt(Ssqave)*log2;
    double phiz4 = J_norm_normalized*0.5*xi*(c_1/Ssqave + xi) - 0.5*(c1_2/Ssqave - nu_2)/sqrt(Ssqave)*log2;
    double phiz5 = J_norm_normalized*xi*(-0.5*c1_2/Ssqave/Ssqave + onesixth*c_1*xi/Ssqave + onethird*(xi_2 + nu_2/Ssqave)) + 0.5*c_1*(c1_2/Ssqave - nu_2)/Ssqave/sqrt(Ssqave)*log2;
    
    return phiz0*constants_phiz[0] + phiz1*constants_phiz[1] + phiz2*constants_phiz[2] + phiz3*constants_phiz[3] + phiz4*constants_phiz[4] + phiz5*constants_phiz[5] + phiz_0 + phiz_MS_corr(xi);
}

/* *********************************************************************************/
/* Internal function that returns zeta                                             */
/* *********************************************************************************/

static double zeta_of_xi(double xi)
{
    double logxi = log(xi);
    double xi_2 = xi*xi;
    double xi_3 = xi_2*xi;
    
    return nu*(-onethird*constants_zeta[0]/xi_3 - 0.5*constants_zeta[1]/xi_2 - constants_zeta[2]/xi + constants_zeta[3]*logxi + constants_zeta[4]*xi  + 0.5*constants_zeta[5]*xi_2) + zeta_0 + zeta_MS_corr(xi);
}

/* *********************************************************************************/
/* Internal function that returns the MS corrections of phiz                       */
/* *********************************************************************************/

static double phiz_MS_corr(double xi)
{
    double c0 = c(xi).x;
    double c2 = c(xi).y;
    double c4 = c(xi).z;
    double d0 = d(xi).x;
    double d2 = d(xi).y;
    double d4 = d(xi).z;
    
    double sqt = sqrt(fabs(d2*d2-4.*d0*d4));
    
    double n3 = (2.*(d0+d2+d4)/(2.*d0+d2+sqt));
    double n4 = ((2*d0+d2+sqt)/(2.*d0));
    double sqtn3 = sqrt(fabs(n3));
    double sqtn4 = sqrt(fabs(n4));
    
    double u = u_of_xi(xi)+constant_of_S;
    
    double L3, L4;
    if(n3 == 1) L3 = 0;
    else L3 = -fabs((c4*d0*(2.*d0+d2+sqt)-c2*d0*(d2+2.*d4-sqt)-c0*(2.*d0*d4-(d2+d4)*(d2-sqt)))/(2.*d0*(d0+d2+d4)*sqt))*(sqtn3/(n3-1.)*(atan(tan(u)) - atan(sqtn3*tan(u))))/psidot(xi);
    if(n4 == 1) L4 = 0;
    else L4 = -fabs((-c4*d0*(2.*d0+d2-sqt)+c2*d0*(d2+2.*d4+sqt)-c0*(-2.*d0*d4+(d2+d4)*(d2+sqt))))/(2.*d0*(d0+d2+d4)*sqt)*(sqtn4/(n4-1.)*(atan(tan(u)) - atan(sqtn4*tan(u))))/psidot(xi);
    
    double corr;
    if(S1_norm == 0 || S2_norm ==0 || (fabs(DotProd(S1_0,Lhat_0)/S1_norm) == 1 && fabs(DotProd(S2_0,Lhat_0)/S2_norm) == 1)) corr = 0;
    else corr = L3+L4;
    if (corr != corr) corr=0;
    return corr;
}

/* *********************************************************************************/
/* Internal function that returns the MS corrections of zeta                       */
/* *********************************************************************************/

static double zeta_MS_corr(double xi)
{
    double J_norm = J_norm_of_xi(xi);
    double c0 = c(xi).x;
    double c2 = c(xi).y;
    double c4 = c(xi).z;
    double d0 = d(xi).x;
    double d2 = d(xi).y;
    double d4 = d(xi).z;
    
    vector roots = Roots(xi,J_norm);
    double A1, A2, A3;
    
    A3 = fmax(fmax(roots.x,roots.y),roots.z);
    A1 = fmin(fmin(roots.x,roots.y),roots.z);
    if((A3 - roots.z) > 0 && (A1 - roots.z) < 0) A2 = roots.z;
    else if((A3 - roots.x) > 0 && (A1 - roots.x) < 0) A2 = roots.x;
    else A2 = roots.y;
    
    double J_normalized = J_norm/GMsquare_over_c;
    double L_normalized = nu/xi;
    
    double Aa = 0.5*(J_normalized/L_normalized + L_normalized/J_normalized - A3/J_normalized/L_normalized);
    double Bb = 0.5*(A3 - A2)/L_normalized/J_normalized;
    double sqt = sqrt(d2*d2-4.*d0*d4);
    double n3 = 2.*(d0+d2+d4)/(2.*d0+d2+sqt);
    double n4 = (2*d0+d2+sqt)/(2.*d0);
    double sqtn3 = sqrt(fabs(n3));
    double sqtn4 = sqrt(fabs(n4));
    
    double u = u_of_xi(xi)+constant_of_S;
    
    double L3, L4;
    if(n3 == 1) L3 = 0;
    else L3 = -fabs((c4*d0*(2.*d0+d2+sqt)-c2*d0*(d2+2.*d4-sqt)-c0*(2.*d0*d4-(d2+d4)*(d2-sqt)))/(2.*d0*(d0+d2+d4)*sqt))*(sqtn3/(n3-1.)*(atan(tan(u)) - atan(sqtn3*tan(u))))/psidot(xi);
    if(n4 == 1) L4 = 0;
    else L4 = -fabs((-c4*d0*(2.*d0+d2-sqt)+c2*d0*(d2+2.*d4+sqt)-c0*(-2.*d0*d4+(d2+d4)*(d2+sqt))))/(2.*d0*(d0+d2+d4)*sqt)*(sqtn4/(n4-1.)*(atan(tan(u)) - atan(sqtn4*tan(u))))/psidot(xi);
    
    double corr;
    if(S1_norm == 0 || S2_norm ==0 || (fabs(DotProd(S1_0,Lhat_0)/S1_norm) == 1 && fabs(DotProd(S2_0,Lhat_0)/S2_norm) == 1)) corr = 0;
    else corr = Aa*phiz_MS_corr(xi) + 2.*Bb*d0*(L3/(sqt-d2)-L4/(sqt+d2));
    if (corr != corr) corr=0;
    return corr;
}


/* *********************************************************************************/
/* Internal function that computes the spin-orbit couplings                        */
/* *********************************************************************************/

static double beta(double a, double b)
{
    return (DotProd(S1_0,Lhat_0)*(a + b*one_over_m1_square)) + (DotProd(S2_0,Lhat_0)*(a + b*one_over_m2_square));
}

/* *********************************************************************************/
/* Internal function that computes the spin-spin couplings                         */
/* *********************************************************************************/

static double sigma(double a, double b)
{
    return a*DotProd(S1_0,S2_0) - b*DotProd(Lhat_0,S1_0)*DotProd(Lhat_0,S2_0);
}

/* *********************************************************************************/
/* Internal function that computes the spin-spin couplings                         */
/* *********************************************************************************/

static double tau(double a, double b)
{
    double cos1 = DotProd(Lhat_0,S1_0);
    double cos2 = DotProd(Lhat_0,S2_0);
    return one_over_m1_square*(S1_norm_square*a - b*cos1*cos1) + one_over_m2_square*(S2_norm_square*a - b*cos2*cos2);
}


/* *********************************************************************************/
/* Internal function that frees the memory allocated in XALInitializeAngles        */
/* *********************************************************************************/

static void FreePrecession()
{
//    free(domegadt_csts_nonspin);
//    free(domegadt_csts_spinorbit);
//    free(domegadt_csts_spinspin);
//    free(L_csts_nonspin);
//    free(L_csts_spinorbit);
}


/* *********************************************************************************/
/* Internal function that returns the dot product of two vectors                   */
/* *********************************************************************************/

static double DotProd(vector vec1, vector vec2)
{
    return vec1.x*vec2.x + vec1.y*vec2.y + vec1.z*vec2.z;
}

/* *********************************************************************************/
/* Internal function that returns the norm of a vector                             */
/* *********************************************************************************/

static double Norm(vector vec)
{
    return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

/* *********************************************************************************/
/* Internal function that calculates a vector fro its spherical components         */
/* *********************************************************************************/

static vector CreateSphere(double r, double th, double ph)
{
    vector out;
    double fact = r*sin(th);
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



