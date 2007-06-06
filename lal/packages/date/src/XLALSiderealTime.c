/*
 * $Id$
 *
 * Copyright (C) 2007  Jolien Creighton, and Kipp Cannon
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 */


/** \file
 * \ingroup std
 * \author Creighton, J., and Cannon, K.
 * \brief XLAL routines for computing the sidereal time.
 */


#include <math.h>
#include <lal/Date.h>

#include <lal/LALRCSID.h>
NRCSID (XLALSIDEREALTIMEC,"$Id$");


/*
 * Returns the Greenwich Sidereal Time IN RADIANS corresponding to a
 * specified GPS time.  Aparent sidereal time is computed by providing the
 * equation of equinoxes in units of seconds.  For mean sidereal time, set
 * this parameter to 0.
 *
 * The output of this code is in radians in the range [0,2pi).
 *
 * Inspired by the function sidereal_time() in the NOVAS-C library, version
 * 2.0.1, which is dated December 10th, 1999, and carries the following
 * references:
 *
 * Aoki, et al. (1982) Astronomy and Astrophysics 105, 359-361.
 * Kaplan, G. H. "NOVAS: Naval Observatory Vector Astrometry
 *   Subroutines"; USNO internal document dated 20 Oct 1988;
 *   revised 15 Mar 1990.
 *
 * See http://aa.usno.navy.mil/software/novas for more information.
 *
 * Note:  rather than maintaining this code separately, it would be a good
 * idea for LAL to simply link to the NOVAS-C library directly.  Something
 * to do when we have some spare time.
 */


REAL8 XLALGreenwichSiderealTime(
	const LIGOTimeGPS *gpstime,
	REAL8 equation_of_equinoxes
)
{
	static const char *func = "XLALGreenwichSiderealTime";
	struct tm utc;
	double julian_day;
	double t_hi, t_lo;
	double t;
	double sidereal_time;

	/*
	 * Convert GPS seconds to UTC.  This is where we pick up knowledge
	 * of leap seconds which are required for the mapping of atomic
	 * time scales to celestial time scales.  We deal only with integer
	 * seconds.
	 */

	if(!XLALGPSToUTC(&utc, gpstime->gpsSeconds))
		XLAL_ERROR_REAL8(func, XLAL_EFUNC);

	/*
	 * And now to Julian day number.  Again, only accurate to integer
	 * seconds.
	 */

	julian_day = XLALJulianDay(&utc);
	if(XLAL_IS_REAL8_FAIL_NAN(julian_day))
		XLAL_ERROR_REAL8(func, XLAL_EFUNC);

	/*
	 * Convert Julian day number to the number of centuries since the
	 * Julian epoch (1 century = 36525.0 days).  Here, we incorporate
	 * the fractional part of the seconds.  For precision, we keep
	 * track of the most significant and least significant parts of the
	 * time separately.  The original code in NOVAS-C determined t_hi
	 * and t_lo from Julian days, with t_hi receiving the integer part
	 * and t_lo the fractional part.  Because LAL's Julian day routine
	 * is accurate to the second, here the hi/lo split is most
	 * naturally done at the integer seconds boundary.  Note that the
	 * "hi" and "lo" components have the same units and so the split
	 * can be done anywhere.
	 */

	t_hi = (julian_day - XLAL_EPOCH_J2000_0_JD) / 36525.0;
	t_lo = gpstime->gpsNanoSeconds / (1e9 * 36525.0 * 86400.0);

	/*
	 * Compute sidereal time in sidereal seconds.  (magic)
	 */

	t = t_hi + t_lo;

	sidereal_time = equation_of_equinoxes + (-6.2e-6 * t + 0.093104) * t * t + 67310.54841;
	sidereal_time += 8640184.812866 * t_lo;
	sidereal_time += 3155760000.0 * t_lo;
	sidereal_time += 8640184.812866 * t_hi;
	sidereal_time += 3155760000.0 * t_hi;

	/*
	 * Return radians (2 pi radians in 1 sidereal day = 86400 sidereal
	 * seconds).
	 */

	return sidereal_time * LAL_PI / 43200.0;
}


REAL8 XLALGreenwichMeanSiderealTime(
	const LIGOTimeGPS *gpstime
)
{
	return XLALGreenwichSiderealTime(gpstime, 0.0);
}


LIGOTimeGPS *XLALGreenwichMeanSiderealTimeToGPS(
	REAL8 gmst,
	LIGOTimeGPS *gps
)
{
	static const char func[] = "XLALGreenwichMeanSiderealTimeToGPS";
	const double gps_seconds_per_sidereal_radian = 13713.44;  /* approx */
	const double precision = 1e-14;
	int iterations = 10;
	double error;

	gps->gpsSeconds = gps->gpsNanoSeconds = 0;

	do {
		error = gmst - XLALGreenwichMeanSiderealTime(gps);
		if(fabs(error / gmst) <= precision)
			return gps;
		XLALAddFloatToGPS(gps, error * gps_seconds_per_sidereal_radian);
	} while(--iterations);
	XLAL_ERROR_NULL(func, XLAL_EMAXITER);
}


LIGOTimeGPS *XLALGreenwichSiderealTimeToGPS(
	REAL8 gmst,
	REAL8 equation_of_equinoxes,
	LIGOTimeGPS *gps
)
{
	return XLALGreenwichMeanSiderealTimeToGPS(gmst - equation_of_equinoxes, gps);
}
