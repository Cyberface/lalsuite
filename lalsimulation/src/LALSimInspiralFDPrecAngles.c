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
/* XLAL function that initalizes all constants needed for the precession angles. It*/
/* needs to be called only once at the beginning.                                  */
/* *********************************************************************************/

void XLALInitializePrecession()
{
    InitializePrecession();
}


/* *********************************************************************************/
/* XLAL function that initalizes all constants that depend on the particular       */
/* system. It needs to be called once when the system parameters have been defined.*/
/* *********************************************************************************/

void XLALInitializeSystem(
    double m1_in,                  /**< mass of body 1 in SI */
    double m2_in,                  /**< mass of body 2 in SI */
    double costhetaL_in,           /**< declination of the orbital angular momentum */
    double phiL_in,                /**< right ascension of the orbital angular momentum  */
    double costheta1_in,           /**< declination of the spin S1 */
    double phi1_in,                /**< right ascension of the spin S1 */
    double chi1_in,                /**< magnitude of the spin S1 (between 0 and 1) */
    double costheta2_in,           /**< declination of the spin S2 */
    double phi2_in,                /**< right ascension of the spin S2 */
    double chi2_in,                /**< magnitude of the spin S2 (between 0 and 1) */
    double f0_in                   /**< reference frequency (Hz) */
){
    m1 = m1_in;
    m2 = m2_in;
    mul = costhetaL_in;
    phl = phiL_in;
    mu1 = costheta1_in;
    ph1 = phi1_in;
    ch1 = chi1_in;
    mu2 = costheta2_in;
    ph2 = phi2_in;
    ch2 = chi2_in;
    f_0 = f0_in;
    
    InitializeSystem();
}

/* *********************************************************************************/
/* XLAL function that returns the cosine of the angle between the Newtonian L and J*/
/* *********************************************************************************/

double XLALcosthetaL(double f)
{
    double xi = pow(f*twopiGM_over_cthree, onethird);
    return costhetaL(xi);
}

/* *********************************************************************************/
/* XLAL function that returns the cosine of the angle between the 3PN L and J      */
/* *********************************************************************************/

double XLALcosthetaL_3PN(double f)
{
    double xi = pow(f*twopiGM_over_cthree, onethird);
    return costhetaL_3PN(xi);
}

/* *********************************************************************************/
/* XLAL function that returns phiz                                                 */
/* *********************************************************************************/

double XLALphiz_of_xi(double f)
{
    double xi = pow(f*twopiGM_over_cthree, onethird);
    return phiz_of_xi(xi);
}

/* *********************************************************************************/
/* XLAL function that returns zeta                                                 */
/* *********************************************************************************/

double XLALzeta_of_xi(double f)
{
    double xi = pow(f*twopiGM_over_cthree, onethird);
    return zeta_of_xi(xi);
}

/* *********************************************************************************/
/* XLAL function that frees the memory allocated in XALInitializeAngles            */
/* *********************************************************************************/

void XLALFreePrecession()
{
    FreePrecession();
}
