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

import os
from NCPTrans import *
from constants import RESOURCES

xy_tolerance = 0.0010
h_tolerance = 0.0010
latlon_tolerance = 0.000000010
replace_none = 'NaN'

def self_validation():
    
    with open(os.path.join(RESOURCES, 'validation', 'Z001_ETRS89andRDNAP.txt')) as inf:
        for line in inf.readlines()[1:]:
            sline = line.strip().split('\t')
            
            pt = {
                'id': int(sline[0]),
                'etrs_lat': float(sline[1]),
                'etrs_lon': float(sline[2]),
                'etrs_h': float(sline[3]),
                'rd_x': float(sline[4]),
                'rd_y': float(sline[5]),
                'nap': None if sline[6] == '*' else float(sline[6])
            }
                
            lat, lon, h = rdnap_to_etrsh(pt['rd_x'], pt['rd_y'], pt['nap'])
            if pt['etrs_lat'] - lat >= latlon_tolerance:
                print(f"{pt['id']}: lat failed test")
            if pt['etrs_lon'] - lon >= latlon_tolerance:
                print(f"{pt['id']}: lon failed test")
            if pt['nap'] is not None and pt['etrs_h'] - h >= h_tolerance:
                print(f"{pt['id']}: height failed test")
                
            x, y, nap = etrsh_to_rdnap(pt['etrs_lat'], pt['etrs_lon'], pt['etrs_h'])
            if pt['rd_x'] - x >= xy_tolerance:
                print(f"{pt['id']}: lat failed test")
            if pt['rd_y'] - y >= xy_tolerance:
                print(f"{pt['id']}: lon failed test")
            if pt['nap'] is not None and pt['nap'] - nap >= h_tolerance:
                print(f"{pt['id']}: height failed test")
        
def validate_rdnap_to_etrsh():

    with open(os.path.join(constants.RESOURCES, 'validation', '002_RDNAP.txt')) as inf:
        with open(os.path.join(constants.RESOURCES, 'validation', '002_RDNAP_to_ETRSH.txt'), 'w') as onf:
            onf.write('point_id\tlatitude\tlongitude\theight\n')
            for line in inf.readlines()[1:]:
                s = line.strip().split(' ')
                s = [v for v in s if v]
                lat, lon, h = rdnap_to_etrsh(float(s[1]), float(s[2]), float(s[3]))
                
                onf.write(f'{s[0]}\t{lat}\t{lon}\t{h if h else replace_none}\n')   
                
def validate_etrsh_to_rdnap():

    with open(os.path.join(constants.RESOURCES, 'validation', '002_ETRS89.txt')) as inf:
        with open(os.path.join(constants.RESOURCES, 'validation', '002_ETRS89_to_RDNAP.txt'), 'w') as onf:
            onf.write('point_id\tx_coordinate\ty_coordinate\theight\n')
            for line in inf.readlines()[1:]:
                s = line.strip().split(' ')
                s = [v for v in s if v]
                x, y, h = etrsh_to_rdnap(float(s[1]), float(s[2]), float(s[3]))
                
                onf.write(f'{s[0]}\t{x}\t{y}\t{h if h else replace_none}\n')
    
    
