#=====================================================================#
# Author: Willem van Opstal
# Created: 2024-12-19
# Last modified: 2024-12-19
#=====================================================================#

#=====================================================================#
# This document is part of:
#
# NCPTrans - Simple horizontal coordinate system and vertical
# reference system transformations on the Dutch Continental Shelf.
#
# Copyright (C) 2024  Willem van Opstal <willemvanopstal@home.nl>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#=====================================================================#

#=====================================================================#
# Transformations between RD-NAP and ETRS89 are based on the document:
# 'RDNAPTRANStm2018: Coordinate transformation to and from Stelsel van
# de Rijksdriehoeksmeting and Normaal Amsterdams Peil' version 27 June
# 2022 by NSGI, available for download at https://nsgi.nl. 
#
# Although this software meets the accuracy requirements set by prior
# document and the NSGI validation service,
#
#    THIS SOFTWARE IS NOT APPROVED TO USE THE NAME RDNAPTRANS2018
#
# While no guarantees are provided in the first place, one should
# be extremely cautious using the results from this software in
# critical or official circumstances. Especially for usage beyond the
# last modification date.
#=====================================================================#

from math import sqrt, cos, sin, tan, radians, log, atan
import os

from importlib import reload
import helpers
reload(helpers)

RESOURCES = helpers.get_script_path()
DB_FILE = "rdnaptrans_grids.sqlite"

ETRS_H0 = 43.0
BESSEL_H0 = 0.0

GRS80_SEMI_MAJOR = 6378137.0
GRS80_FLATTENING = 1 / 298.257222101
GRS80_SQ_ECCENTRICITY = GRS80_FLATTENING * ( 2.0 - GRS80_FLATTENING )
GRS80_ECCENTRICITY = sqrt(GRS80_SQ_ECCENTRICITY)

BESSEL_SEMI_MAJOR = 6377397.155
BESSEL_FLATTENING = 1 / 299.1528128
BESSEL_SQ_ECCENTRICITY = BESSEL_FLATTENING * ( 2.0 - BESSEL_FLATTENING )
BESSEL_ECCENTRICITY = sqrt(BESSEL_SQ_ECCENTRICITY)

ETRSBESSEL_TX = -565.7346
ETRSBESSEL_TY = -50.4058
ETRSBESSEL_TZ = -465.2895
ETRSBESSEL_A = -1.91513e-6
ETRSBESSEL_B = 1.60365e-6
ETRSBESSEL_Y = -9.09546e-6
ETRSBESSEL_D = -4.07242e-6
ETRSBESSEL_S = 1.0 + ETRSBESSEL_D

ETRSBESSEL = {
    'Tx': ETRSBESSEL_TX,
    'Ty': ETRSBESSEL_TY,
    'Tz': ETRSBESSEL_TZ,
    'Alpha': ETRSBESSEL_A,
    'Beta': ETRSBESSEL_B,
    'Gamma': ETRSBESSEL_Y,
    'Delta': ETRSBESSEL_D,
    'Scale': ETRSBESSEL_S,
    'R11': cos(ETRSBESSEL_Y) * cos(ETRSBESSEL_B),
    'R12': cos(ETRSBESSEL_Y) * sin(ETRSBESSEL_B) * sin(ETRSBESSEL_A) + sin(ETRSBESSEL_Y) * cos(ETRSBESSEL_A),
    'R13': -1 * cos(ETRSBESSEL_Y) * sin(ETRSBESSEL_B) * cos(ETRSBESSEL_A) + sin(ETRSBESSEL_Y) * sin(ETRSBESSEL_A),
    'R21': -1 * sin(ETRSBESSEL_Y) * cos(ETRSBESSEL_B),
    'R22': -1 * sin(ETRSBESSEL_Y) * sin(ETRSBESSEL_B) * sin(ETRSBESSEL_A) + cos(ETRSBESSEL_Y) * cos(ETRSBESSEL_A),
    'R23': sin(ETRSBESSEL_Y) * sin(ETRSBESSEL_B) * cos(ETRSBESSEL_A) +  cos(ETRSBESSEL_Y) * sin(ETRSBESSEL_A),
    'R31': sin(ETRSBESSEL_B),
    'R32': -1 * cos(ETRSBESSEL_B) * sin(ETRSBESSEL_A),
    'R33': cos(ETRSBESSEL_B) * cos(ETRSBESSEL_A)
}

CTG_ITERATION_THRESHOLD = 2e-11 #rad
CTG_ROUNDING = 4 #digits

RDCORR2018_PHI_MIN = 50.0
RDCORR2018_PHI_MAX = 56.0
RDCORR2018_LMBDA_MIN = 2.0
RDCORR2018_LMBDA_MAX = 8.0
RDCORR2018_PHI_DELTA = 0.0125
RDCORR2018_LMBDA_DELTA = 0.02

RDCORR2018_ITERATION_THRESHOLD = 0.000000001
RDCORR2018_OUT_OF_BOUNDS_CORRECTION = 0.0

RDCORR2018 = {
    'TYPE': 'ASCII',
    'DB': os.path.join(RESOURCES, DB_FILE),
    'TABLE': 'rdcorr2018',
    'N0': 1,
    'PHI_MIN': RDCORR2018_PHI_MIN,
    'PHI_MAX': RDCORR2018_PHI_MAX,
    'LMBDA_MIN': RDCORR2018_LMBDA_MIN,
    'LMBDA_MAX': RDCORR2018_LMBDA_MAX,
    'PHI_DELTA': RDCORR2018_PHI_DELTA,
    'LMBDA_DELTA': RDCORR2018_LMBDA_DELTA,
    'PHI_N': 1 + (RDCORR2018_PHI_MAX - RDCORR2018_PHI_MIN) / RDCORR2018_PHI_DELTA,
    'LMBDA_N': 1 + (RDCORR2018_LMBDA_MAX - RDCORR2018_LMBDA_MIN) / RDCORR2018_LMBDA_DELTA,
    'C0': RDCORR2018_OUT_OF_BOUNDS_CORRECTION,
    'ITERATION_THRESHOLD': RDCORR2018_ITERATION_THRESHOLD
}

RD_PROJECTION_BESSEL_LAT_DEG = 52.15616056
RD_PROJECTION_BESSEL_LON_DEG = 5.38763889
RD_PROJECTION_SCALE = 0.9999079
RD_PROJECTION_FEASTING = 155000.0
RD_PROJECTION_FNORTHING = 463000.0

RD_PROJECTION_RN = helpers.second_curvature(radians(RD_PROJECTION_BESSEL_LAT_DEG), BESSEL_SEMI_MAJOR, BESSEL_SQ_ECCENTRICITY)
RD_PROJECTION_RM = helpers.first_curvature(RD_PROJECTION_RN, radians(RD_PROJECTION_BESSEL_LAT_DEG), BESSEL_SEMI_MAJOR, BESSEL_SQ_ECCENTRICITY)
RD_PROJECTION_R = sqrt(RD_PROJECTION_RM * RD_PROJECTION_RN)
RD_PROJECTION_SPHERE_LAT_RAD = atan( sqrt(RD_PROJECTION_RM) / sqrt(RD_PROJECTION_RN) * tan(radians(RD_PROJECTION_BESSEL_LAT_DEG)) )
RD_PROJECTION_SPHERE_LON_RAD = radians(RD_PROJECTION_BESSEL_LON_DEG)
RD_PROJECTION_W0 = log( tan((RD_PROJECTION_SPHERE_LAT_RAD + radians(90.0)) / 2) )
RD_PROJECTION_Q0 = log( tan((radians(RD_PROJECTION_BESSEL_LAT_DEG) + radians(90.0)) / 2) ) - (BESSEL_ECCENTRICITY / 2.0) * log(( 1 + BESSEL_ECCENTRICITY * sin(radians(RD_PROJECTION_BESSEL_LAT_DEG))) / ( 1 - BESSEL_ECCENTRICITY * sin(radians(RD_PROJECTION_BESSEL_LAT_DEG))))
RD_PROJECTION_N = sqrt(1 + (BESSEL_SQ_ECCENTRICITY * cos(radians(RD_PROJECTION_BESSEL_LAT_DEG))**4 ) / (1 - BESSEL_SQ_ECCENTRICITY))
RD_PROJECTION_M = RD_PROJECTION_W0 - RD_PROJECTION_N * RD_PROJECTION_Q0

RD_SPHERICAL = {
    'E': BESSEL_ECCENTRICITY,
    'M': RD_PROJECTION_M,
    'N': RD_PROJECTION_N,
    'R': RD_PROJECTION_R,
    'L0_SPHERICAL': RD_PROJECTION_SPHERE_LON_RAD,
    'L0_BESSEL': radians(RD_PROJECTION_BESSEL_LON_DEG),
    'P0_SPHERICAL': RD_PROJECTION_SPHERE_LAT_RAD,
    'K': RD_PROJECTION_SCALE,
    'X0': RD_PROJECTION_FEASTING,
    'Y0': RD_PROJECTION_FNORTHING
}

NLGEO2018_PHI_MIN = 50.0
NLGEO2018_PHI_MAX = 56.0
NLGEO2018_LMBDA_MIN = 2.0
NLGEO2018_LMBDA_MAX = 8.0
NLGEO2018_PHI_DELTA = 0.0125
NLGEO2018_LMBDA_DELTA = 0.02

NLGEO2018_OUT_OF_BOUNDS_VALUE = None

NLGEO2018 = {
    'TYPE': 'ASCII',
    'DB': os.path.join(RESOURCES, DB_FILE),
    'TABLE': 'nlgeo2018',
    'N0': 1,
    'PHI_MIN': NLGEO2018_PHI_MIN,
    'PHI_MAX': NLGEO2018_PHI_MAX,
    'LMBDA_MIN': NLGEO2018_LMBDA_MIN,
    'LMBDA_MAX': NLGEO2018_LMBDA_MAX,
    'PHI_DELTA': NLGEO2018_PHI_DELTA,
    'LMBDA_DELTA': NLGEO2018_LMBDA_DELTA,
    'PHI_N': 1 + (NLGEO2018_PHI_MAX - NLGEO2018_PHI_MIN) / NLGEO2018_PHI_DELTA,
    'LMBDA_N': 1 + (NLGEO2018_LMBDA_MAX - NLGEO2018_LMBDA_MIN) / NLGEO2018_LMBDA_DELTA,
    'C0': NLGEO2018_OUT_OF_BOUNDS_VALUE
}

RD_INV_PROJECTION_THRESHOLD = 2e-11

BESSELETRS_TX = 565.7381
BESSELETRS_TY = 50.4018
BESSELETRS_TZ = 465.2904
BESSELETRS_A = 1.91514e-6
BESSELETRS_B = -1.60363e-6
BESSELETRS_Y = 9.09546e-6
BESSELETRS_D = 4.07244e-6
BESSELETRS_S = 1.0 + BESSELETRS_D

BESSELETRS = {
    'Tx': BESSELETRS_TX,
    'Ty': BESSELETRS_TY,
    'Tz': BESSELETRS_TZ,
    'Alpha': BESSELETRS_A,
    'Beta': BESSELETRS_B,
    'Gamma': BESSELETRS_Y,
    'Delta': BESSELETRS_D,
    'Scale': BESSELETRS_S,
    'R11': cos(BESSELETRS_Y) * cos(BESSELETRS_B),
    'R12': cos(BESSELETRS_Y) * sin(BESSELETRS_B) * sin(BESSELETRS_A) + sin(BESSELETRS_Y) * cos(BESSELETRS_A),
    'R13': -1 * cos(BESSELETRS_Y) * sin(BESSELETRS_B) * cos(BESSELETRS_A) + sin(BESSELETRS_Y) * sin(BESSELETRS_A),
    'R21': -1 * sin(BESSELETRS_Y) * cos(BESSELETRS_B),
    'R22': -1 * sin(BESSELETRS_Y) * sin(BESSELETRS_B) * sin(BESSELETRS_A) + cos(BESSELETRS_Y) * cos(BESSELETRS_A),
    'R23': sin(BESSELETRS_Y) * sin(BESSELETRS_B) * cos(BESSELETRS_A) +  cos(BESSELETRS_Y) * sin(BESSELETRS_A),
    'R31': sin(BESSELETRS_B),
    'R32': -1 * cos(BESSELETRS_B) * sin(BESSELETRS_A),
    'R33': cos(BESSELETRS_B) * cos(BESSELETRS_A)
}

NLLAT2018_PHI_MIN = 51.0
NLLAT2018_PHI_MAX = 56.0
NLLAT2018_LMBDA_MIN = 2.5
NLLAT2018_LMBDA_MAX = 7.5
NLLAT2018_PHI_DELTA = 0.00625
NLLAT2018_LMBDA_DELTA = 0.01

NLLAT2018_OUT_OF_BOUNDS_VALUE = None

NLLAT2018 = {
    'TYPE': 'ASCII',
    'DB': os.path.join(RESOURCES, DB_FILE),
    'TABLE': 'nllat2018',
    'N0': 0,
    'PHI_MIN': NLLAT2018_PHI_MIN,
    'PHI_MAX': NLLAT2018_PHI_MAX,
    'LMBDA_MIN': NLLAT2018_LMBDA_MIN,
    'LMBDA_MAX': NLLAT2018_LMBDA_MAX,
    'PHI_DELTA': NLLAT2018_PHI_DELTA,
    'LMBDA_DELTA': NLLAT2018_LMBDA_DELTA,
    'PHI_N': 1 + (NLLAT2018_PHI_MAX - NLLAT2018_PHI_MIN) / NLLAT2018_PHI_DELTA,
    'LMBDA_N': 1 + (NLLAT2018_LMBDA_MAX - NLLAT2018_LMBDA_MIN) / NLLAT2018_LMBDA_DELTA,
    'C0': NLLAT2018_OUT_OF_BOUNDS_VALUE
}