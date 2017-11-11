# Copyright (C) 2017 Sebastian Khan
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

# \defgroup lalsimulation_py_waveform Waveform
# \ingroup lalsimulation_python
"""This module provides code for generating LAL waveforms in python
"""

from .. import git_version

import lalsimulation as lalsim

__author__ = "Sebastian Khan <sebastian.khan@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date

# just an example for when there are functions
# __all__ = ['CacheEntry', 'lalcache_from_gluecache']


def test():
    print("Hello! I am the function test() from waveform.py")


td_required_args = [
    "m1",
    "m2",
    "s1x",
    "s1y",
    "s1z",
    "s2x",
    "s2y",
    "s2z",
    "distance",
    "inclination",
    "phiRef",
    "longAscNodes",
    "eccentricity",
    "meanPerAno",
    "deltaT",
    "f_min",
    "f_ref",
    "params",
    "approximant"
]

fd_required_args = [
    "m1",
    "m2",
    "S1x",
    "S1y",
    "S1z",
    "S2x",
    "S2y",
    "S2z",
    "distance",
    "inclination",
    "phiRef",
    "longAscNodes",
    "eccentricity",
    "meanPerAno",
    "deltaF",
    "f_min",
    "f_max",
    "f_ref",
    "LALpars",
    "approximant"
]

def WrapperSimInspiralChooseTDWaveform(**kwargs):
    """wrapper for XLALSimInspiralChooseTDWaveform

    Example
    -------

    >>> import lal
    >>> import lalsimulation as lalsim
    >>> import lalsimulation.waveform as waveform
    >>> hp, hc = waveform.WrapperSimInspiralChooseTDWaveform(m1=50*lal.MSUN_SI,m2=40.*lal.MSUN_SI, s1x=0., s1y=0., s1z=0., s2x=0., s2y=0., s2z=0., distance=1e6*lal.PC_SI, inclination=0., phiRef=0., longAscNodes=0., eccentricity=0., meanPerAno=0., deltaT=1./1024, f_min=10, f_ref=10., params=None, approximant=lalsim.IMRPhenomD)

    """
    for arg in td_required_args:
        if arg not in kwargs:
            raise ValueError("Please provide " + str(arg))

    hp, hc = lalsim.SimInspiralChooseTDWaveform(**kwargs)
    return hp, hc


def WrapperSimInspiralChooseFDWaveform(**kwargs):
    """wrapper for XLALSimInspiralChooseFDWaveform

    Example
    -------

    >>> import lal
    >>> import lalsimulation as lalsim
    >>> import lalsimulation.waveform as waveform
    >>> hptilde, hctilde = waveform.WrapperSimInspiralChooseFDWaveform(m1=50*lal.MSUN_SI,m2=40.*lal.MSUN_SI, S1x=0., S1y=0., S1z=0., S2x=0., S2y=0., S2z=0., distance=1e6*lal.PC_SI, inclination=0., phiRef=0., longAscNodes=0., eccentricity=0., meanPerAno=0., deltaF=1./128, f_min=10, f_max=0, f_ref=10., LALpars=None, approximant=lalsim.IMRPhenomD)

    """
    for arg in fd_required_args:
        if arg not in kwargs:
            raise ValueError("Please provide " + str(arg))

    hptilde, hctilde = lalsim.SimInspiralChooseFDWaveform(**kwargs)
    return hptilde, hctilde
