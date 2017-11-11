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

## \defgroup lalsimulation_py_parameters Parameters
## \ingroup lalsimulation_python


#
# =============================================================================
#
#                                   Preamble
#
# =============================================================================
#

"""This module defines the parameters used to generate gravitaional waveforms
based upon the pycbc parameters.py.
see https://github.com/ligo-cbc/pycbc/blob/master/pycbc/waveform/parameters.py
"""

from .. import git_version
__author__ = "Sebastian Khan <sebastian.khan@ligo.org>"
__version__ = git_version.id
__date__ = git_version.date


#
# =============================================================================
#
#                                  Base definitions
#
# =============================================================================
#

class Parameter(str):
    """A class that stores information about a parameter. This is done by
    sub-classing string, adding additional attributes.
    ****STOLEN FROM PYCBC/WAVEFORM/PARAMETERS.PY****
    """

    def __new__(cls, name, dtype=None, default=None, label=None,
                description="No description."):
        obj = str.__new__(cls, name)
        obj.name = name
        obj.dtype = dtype
        obj.default = default
        obj.label = label
        obj.description = description
        return obj

    def docstr(self, prefix='', include_label=True):
        """Returns a string summarizing the parameter. Format is:
        <prefix>``name`` : {``default``, ``dtype``}
        <prefix>   ``description`` Label: ``label``.
        """
        outstr = "%s%s : {%s, %s}\n" % (prefix, self.name, str(self.default),
                                        str(self.dtype).replace("<type '", '').replace("'>", '')) + \
            "%s    %s" % (prefix, self.description)
        if include_label:
            outstr += " Label: %s" % (self.label)
        return outstr


#
# =============================================================================
#
#                        Parameter definitions
#
# =============================================================================
#


#
#   CBC intrinsic parameters
#
mass1 = Parameter("mass1",
                  dtype=float, default=None, label=r"$m_1~(\mathrm{M}_\odot)$",
                  description="The mass of the first component object in the "
                  "binary (in solar masses).")
mass2 = Parameter("mass2",
                  dtype=float, default=None, label=r"$m_2~(\mathrm{M}_\odot)$",
                  description="The mass of the second component object in the "
                  "binary (in solar masses).")
spin1x = Parameter("spin1x",
                   dtype=float, default=0., label=r"$\chi_{1x}$",
                   description="The x component of the first binary component's "
                   "dimensionless spin.")
spin1y = Parameter("spin1y",
                   dtype=float, default=0., label=r"$\chi_{1y}$",
                   description="The y component of the first binary component's "
                   "dimensionless spin.")
spin1z = Parameter("spin1z",
                   dtype=float, default=0., label=r"$\chi_{1z}$",
                   description="The z component of the first binary component's "
                   "dimensionless spin.")
spin2x = Parameter("spin2x",
                   dtype=float, default=0., label=r"$\chi_{2x}$",
                   description="The x component of the second binary component's "
                   "dimensionless spin.")
spin2y = Parameter("spin2y",
                   dtype=float, default=0., label=r"$\chi_{2y}$",
                   description="The y component of the second binary component's "
                   "dimensionless spin.")
spin2z = Parameter("spin2z",
                   dtype=float, default=0., label=r"$\chi_{2z}$",
                   description="The z component of the second binary component's "
                   "dimensionless spin.")
eccentricity = Parameter("eccentricity",
                         dtype=float, default=0., label=r"$e$",
                         description="Eccentricity.")
