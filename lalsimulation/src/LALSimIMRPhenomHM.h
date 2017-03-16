/*
 * Copyright (C) 2017 Sebastian Khan, Francesco Pannarale
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

// #define NMODES_MAX 8
#define NMODES_MAX 5
// #define NMODES_MAX 4

// #include <LALSimIMRPhenomHM.c>

/**
  * Useful powers in Mf: 1/6, 1/3, 2/3, 4/3, 5/3, 2, 7/3, 8/3, -7/6, -5/6, -1/2, -1/6, 1/2
  * calculated using only one invocation of 'pow', the rest are just multiplications and divisions
  */
typedef struct tagUsefulMfPowers
{
    REAL8 sixth;
    REAL8 third;
    REAL8 two_thirds;
    REAL8 four_thirds;
    REAL8 five_thirds;
    REAL8 two;
    REAL8 seven_thirds;
    REAL8 eight_thirds;
    REAL8 m_seven_sixths;
    REAL8 m_five_sixths;
    REAL8 m_sqrt;
    REAL8 m_sixth;
    REAL8 sqrt;
} UsefulMfPowers;

/**
 * must be called before the first usage of *p
 */
int init_useful_mf_powers(UsefulMfPowers * p, REAL8 number);

/**
 * Structure storing (2,2) mode quantities and other parameters/properties
 * that are needed to compute all other modes
 */
typedef struct tagPhenomDStorage
{
    REAL8 Inv1MinusEradRational0815;
    REAL8 finspin;
    REAL8 Mf_RD_22;
    REAL8 Mf_DM_22;
    REAL8 PhenomHMfring[NMODES_MAX][NMODES_MAX];
    REAL8 PhenomHMfdamp[NMODES_MAX][NMODES_MAX];
    REAL8 Rholm[NMODES_MAX][NMODES_MAX];
    REAL8 Taulm[NMODES_MAX][NMODES_MAX];
    REAL8 pow_Mf_wf_prefactor[NMODES_MAX][NMODES_MAX];
} PhenomDStorage;

/**
 * must be called before the first usage of *p
 */
int init_PhenomD_Storage(PhenomDStorage* p, const REAL8 m1, const REAL8 m2, const REAL8 chi1z, const REAL8 chi2z);

/**
  * Structure holding Higher Mode Phase pre-computations
  */
typedef struct tagHMPhasePreComp {
 double ai;
 double bi;
 double a2lm;
 double b2lm;
 double ar;
 double br;
 double fi;
 double f2lm;
 double fr;
 double PhDBconst;
 double PhDCconst;
 double PhDDconst;
 double PhDBAterm;
} HMPhasePreComp;



int XLALSimIMRPhenomHMPhasePreComp(struct tagHMPhasePreComp *q, const INT4 ell, const INT4 mm, const REAL8 eta, const REAL8 chi1z, const REAL8 chi2z, const double finspin);
