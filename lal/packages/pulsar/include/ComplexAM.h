/*
*  Copyright (C) 2007 John Whelan
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

/*** <lalVerbatim file="ComplexAMHV">
Author: Whelan, J.T.
$Id$
*** </lalVerbatim> */

/* <lalLaTeX>
\section{Header \texttt{ComplexAM.h}}
\label{s:ComplexAM.h}
Computes filter components for amplitude demodulation.

\subsection*{Synposis}
\begin{verbatim}
#include <lal/ComplexAM.h>
\end{verbatim}

\noindent  structures and prototypes associated with complex AM coefficients

</lalLaTeX> */


#ifndef _COMPLEXAM_H
#define _COMPLEXAM_H

#include <math.h>
#include <lal/DetResponse.h>
#include <lal/DetectorSite.h>
#include <lal/ComputeFstat.h>
#include "LALBarycenter.h"

#ifdef __cplusplus
extern "C" {
#endif

NRCSID (COMPLEXAMH, "$Id: ComplexAM.h");


/*----- Error-codes -----*/
#define COMPLEXAMC_ENULL 		1
#define COMPLEXAMC_ENONULL 		2
#define COMPLEXAMC_EINPUT   		3
#define COMPLEXAMC_EMEM   		4
#define COMPLEXAMC_EXLAL		5
#define COMPLEXAMC_EIEEE		6
#define COMPLEXAMC_ERAALISA		7

#define COMPLEXAMC_MSGENULL 		"Arguments contained an unexpected null pointer"
#define COMPLEXAMC_MSGENONULL 	"Output pointer is non-NULL"
#define COMPLEXAMC_MSGEINPUT   	"Invalid input"
#define COMPLEXAMC_MSGEMEM   	"Out of memory. Bad."
#define COMPLEXAMC_MSGEXLAL	"XLAL function call failed"
#define COMPLEXAMC_MSGEIEEE	"Floating point failure"
#define COMPLEXAMC_MSGERAALISA	"RAA response only available for LISA"

  /* <lalLaTeX>
\subsection*{Structures}

</lalLaTeX> */

/** The 'detector tensor' for a GW-detector: symmetric 3x3 matrix, storing only the upper triangle.
 * The coordinate-system is SSB-fixed Cartesian coordinates, in particular EQUATORIAL coords for
 * Earth-based detectors and ECLIPTIC coords for LISA.
 */
typedef struct
{
  SymmTensor3 re;	/**< tensor holding real-parts of all components */
  SymmTensor3 im;	/**< tensor holding imaginary-parts of all components */
} CmplxDetectorTensor;

  /* <lalLaTeX>
\newpage\input{ComplexAMHV}
% \newpage\input{ComplexAMC}
</lalLaTeX> */

/*---------- exported prototypes [API] ----------*/

void
LALGetCmplxAMCoeffs( LALStatus *,
		     CmplxAMCoeffs *coeffs,
		     const DetectorStateSeries *DetectorStates,
		     PulsarDopplerParams doppler);

void
LALGetMultiCmplxAMCoeffs( LALStatus *,
			  MultiCmplxAMCoeffs **multiAMcoef,
			  const MultiDetectorStateSeries *multiDetStates,
			  PulsarDopplerParams doppler );

int
XLALWeighMultiCmplxAMCoeffs ( MultiCmplxAMCoeffs *multiAMcoef, const MultiNoiseWeights *multiWeights );

/* destructors */
void XLALDestroyMultiCmplxAMCoeffs ( MultiCmplxAMCoeffs *multiAMcoef );

#ifdef __cplusplus
}
#endif

#endif /* _COMPLEXAM_H */
