/*
 *  Copyright (C) 2017 Sebastian Khan
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

/**
 * \author Sebastian Khan
 *
 * \file
 *
 * \brief External (SWIG'd) Auxiliary functions for phenomenological model development
 *
 * Helper functions for phenomenological waveform models
 * Can be used through python SWIG wrapping
 * NOTE: The convention for naming functions in there is to use
 * the prefix 'XLALSimPhenom_'
 */


#include <lal/LALSimIMRPhenomUtils.h>
#include "LALSimIMRPhenomInternalUtils.h"

// /**
//  * Example how to write an external XLAL phenom function
//  */
// void XLALSimPhenomUtilsTest(){
//     printf("Hello! I am the XLALSimPhenomUtilsTest function\n");
// }

/**
 * Tests the PhenomInternal_nudge function
 */
void XLALSimPhenomUtilsTestNudge(){
    REAL8 eta = 0.25 + 1e-5;
    REAL8 target = 0.25;
    REAL8 epsilon = 1e-3;

    if (eta > target)
    {
        printf("eta = %.15f is greater than target = %.15f.\n\
Nudging eta by epsilon = %.15f\n", eta, target, epsilon);
    }

    printf("before nudge: eta = %.15f\n", eta);
    PhenomInternal_nudge(&eta, target, epsilon);
    printf("after nudge: eta = %.15f\n", eta);
}

/**
 * Convert from geometric frequency to frequency in Hz
 */
double XLALSimPhenomUtilsMftoHz(
    REAL8 Mf, /**< Geometric frequency */
    REAL8 Mtot_Msun /**< Total mass in solar masses */
)
{
    return Mf / ( LAL_MTSUN_SI * Mtot_Msun );
}

/**
 * Convert from frequency in Hz to geometric frequency
 */
double XLALSimPhenomUtilsHztoMf(
    REAL8 fHz, /**< Frequency in Hz */
    REAL8 Mtot_Msun  /**< Total mass in solar masses */
)
{
    return fHz * ( LAL_MTSUN_SI * Mtot_Msun );
}
