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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "gsl/gsl_sf_elljac.h"
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_mode.h>
#include <lal/LALError.h>

#include "LALSimInspiralFDPrecAngles.h"

static void InitializePrecession(void);
static void FreePrecession(void);
static void InitializeSystem(void);

static double DotProd(vector vec1, vector vec2);
static double Norm(vector vec1);
static vector CreateSphere(double r, double th, double ph);
static vector ScalarProd(double c, vector vec);
static vector Sum(vector vec1, vector vec2);
static vector CrossProd(vector vec1, vector vec2);

static vector Roots(double xi, double J_norm);
static vector BCDcoeff(double xi, double J_norm);

static double beta(double a, double b);
static double sigma(double a, double b);
static double tau(double a, double b);

static double J_norm_of_xi(double xi);
static double S_norm_of_xi(double xi);
static double J_norm_3PN_of_xi(double xi);
static double L_norm_3PN_of_xi(double xi);

static vector c(double xi);
static vector d(double xi);

static double costhetaL(double xi);
static double costhetaL_3PN(double xi);

static double u_of_xi(double xi);
static double psidot(double xi);

static double phiz_MS_corr(double xi);
static double zeta_MS_corr(double xi);
static double phiz_of_xi(double xi);
static double zeta_of_xi(double xi);

