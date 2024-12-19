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

import math
import sqlite3
import numpy as np
import os
import sys

from importlib import reload
import constants
reload(constants)

def get_script_path():
    return os.path.dirname(os.path.realpath(__file__))

def radius_curvature(phi, semi_major, squared_eccentricity):
    # R
    second_curvature_value = second_curvature(phi, semi_major, squared_eccentricity)
    return math.sqrt(first_curvature(second_curvature_value, phi, semi_major, squared_eccentricity) * second_curvature_value)
    
def first_curvature(second_curvature_value, phi, semi_major, squared_eccentricity):
    # Rm
    return (second_curvature_value * (1.0 - squared_eccentricity)) / ( 1.0 - squared_eccentricity * math.sin(phi)**2 )
    
def second_curvature(phi, semi_major, squared_eccentricity):
    # Rn
    return semi_major / math.sqrt( 1.0 - squared_eccentricity * math.sin(phi)**2 )
    
def is_within_bounds(phi, lmbda, grid):
    if grid['PHI_MIN'] <= phi <= grid['PHI_MAX'] and \
       grid['LMBDA_MIN'] <= lmbda <= grid['LMBDA_MAX']:
        return True
    return False
    
def get_grid_indices(phi_norm, lmbda_norm, grid):
    if grid['TYPE'] == 'ASCII':
        i_nw = math.ceil(phi_norm) * grid['LMBDA_N'] + math.floor(lmbda_norm) + grid['N0']
        i_ne = math.ceil(phi_norm) * grid['LMBDA_N'] + math.ceil(lmbda_norm) + grid['N0']
        i_sw = math.floor(phi_norm) * grid['LMBDA_N'] + math.floor(lmbda_norm) + grid['N0']
        i_se = math.floor(phi_norm) * grid['LMBDA_N'] + math.ceil(lmbda_norm) + grid['N0']
    
    return i_nw, i_ne, i_sw, i_se
    
#======================================================#

def insert_correction_grid_into_table(db, table, text_file):

    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()

        sql = f'''
            CREATE TABLE IF NOT EXISTS {table} (
                line INTEGER PRIMARY KEY,
                lat REAL,
                lon REAL,
                lat_corr REAL,
                lon_corr REAL
            )
        '''
        cursor.execute(sql)
        conn.commit()

        sql = f'''
            DELETE FROM {table};
        '''
        cursor.execute(sql)
        conn.commit()

        sql = f'''
            INSERT INTO {table} (line, lat, lon, lat_corr, lon_corr)
            VALUES (?, ?, ?, ?, ?);
        '''

        with open(text_file) as inf:
            lines = inf.readlines()

            for l, line in enumerate(lines):
                if not l:
                    continue
                sline = line.strip().split('\t')
                fline = [float(s) for s in sline]

                cursor.execute(sql, (l, fline[0], fline[1], fline[2], fline[3]))
        conn.commit()

def insert_geoid_grid_into_table(db, table, text_file):

    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()

        sql = f'''
            CREATE TABLE IF NOT EXISTS {table} (
                line INTEGER PRIMARY KEY,
                lat REAL,
                lon REAL,
                geoid REAL
            )
        '''
        cursor.execute(sql)
        conn.commit()

        sql = f'''
            DELETE FROM {table};
        '''
        cursor.execute(sql)
        conn.commit()

        sql = f'''
            INSERT INTO {table} (line, lat, lon, geoid)
            VALUES (?, ?, ?, ?);
        '''

        with open(text_file) as inf:
            lines = inf.readlines()

            for l, line in enumerate(lines):
                if not l:
                    continue
                sline = line.strip().split('\t')
                fline = [float(s) for s in sline]

                cursor.execute(sql, (l, fline[0], fline[1], fline[2]))
        conn.commit()
        
def insert_lat_grid_into_table(db, table, text_file):

    with sqlite3.connect(db) as conn:
        cursor = conn.cursor()

        sql = f'''
            CREATE TABLE IF NOT EXISTS {table} (
                line INTEGER PRIMARY KEY,
                lat REAL,
                lon REAL,
                geoid_ellipsoid REAL,
                lat_ellipsoid REAL,
                geoid_lat REAL
            )
        '''
        cursor.execute(sql)
        conn.commit()

        sql = f'''
            DELETE FROM {table};
        '''
        cursor.execute(sql)
        conn.commit()

        sql = f'''
            INSERT INTO {table} (line, lat, lon, geoid_ellipsoid, lat_ellipsoid, geoid_lat)
            VALUES (?, ?, ?, ?, ?, ?);
        '''

        pts = {}
        with open(text_file) as inf:
            for line in inf.readlines()[5:]:
                sline = line.strip().split(' ')
                fline = [float(s) for s in sline if s]
                pts[(fline[0], fline[1])] = (fline[2], fline[3], fline[4])
        
        i = 0
        for p in np.linspace(constants.NLLAT2018['PHI_MIN'], constants.NLLAT2018['PHI_MAX'], int((constants.NLLAT2018['PHI_MAX'] - constants.NLLAT2018['PHI_MIN']) / constants.NLLAT2018['PHI_DELTA'] + 1)):
            for l in np.linspace(constants.NLLAT2018['LMBDA_MIN'], constants.NLLAT2018['LMBDA_MAX'], int((constants.NLLAT2018['LMBDA_MAX'] - constants.NLLAT2018['LMBDA_MIN']) / constants.NLLAT2018['LMBDA_DELTA'] + 1)):
                p = round(p, 5)
                l = round(l, 5)
                if (p, l) in pts:
                    vals = pts[(p,l)]
                    cursor.execute(sql, (i, p, l, vals[0], vals[1], vals[2]))
                i += 1
        conn.commit()
        
def prepare_grid_db():
    db_file = os.path.join(constants.RESOURCES, constants.DB_FILE)
    
    rdcorr2018_file = os.path.join(constants.RESOURCES, 'rdcorr2018.txt')
    rdcorr2018_table = constants.RDCORR2018['TABLE']
    nlgeo2018_file = os.path.join(constants.RESOURCES, 'nlgeo2018.txt')
    nlgeo2018_table = constants.NLGEO2018['TABLE']
    nllat2018_file = os.path.join(constants.RESOURCES, 'NLLAT2018.txt')
    nllat2018_table = constants.NLLAT2018['TABLE']

    insert_correction_grid_into_table(db_file, rdcorr2018_table, rdcorr2018_file)
    insert_geoid_grid_into_table(db_file, nlgeo2018_table, nlgeo2018_file)
    insert_lat_grid_into_table(db_file, nllat2018_table, nllat2018_file)