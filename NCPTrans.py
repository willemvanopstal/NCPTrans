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
import os
import sqlite3

from importlib import reload
import constants
import helpers

reload(constants)
reload(helpers)

class Transformer():

    cache = {}

    def __init__(self):
        
        self.etrs_lat = None
        self.etrs_lon = None
        self.etrs_lat_rad = None
        self.etrs_lon_rad = None
        self.etrs_iterations = 0
        self.etrs_error = None
        self.etrs_x = None
        self.etrs_y = None
        self.etrs_z = None
        
        self.pseudo_bessel_x = None
        self.pseudo_bessel_y = None
        self.pseudo_bessel_z = None
        self.pseudo_bessel_lat_deg = None
        self.pseudo_bessel_lon_deg = None
        self.pseudo_bessel_lat_rad = None
        self.pseudo_bessel_lon_rad = None
        self.pseudo_bessel_iterations = 0
        self.pseudo_bessel_error = None
        self.pseudo_bessel_corrections = None
        
        self.rd_correction_out_of_bounds = False
        self.bessel_lat_deg = None
        self.bessel_lon_deg = None
        self.bessel_lat_rad = None
        self.bessel_lon_rad = None
        self.bessel_iterations = 0
        self.bessel_corrections = None
        self.bessel_error = None
        
        self.geoid_ellipsoid_separation = None
        self.lat_ellipsoid_separation = None
        self.geoid_lat_separation = None
        self.ellipsoidal_height = None
        self.nap_height = None
        self.lat_height = None
        
    def __str__(self):
        
        return
    #==============================================================================================#
    
    def load_rd_xy(self, x, y):
        self.rd_x = x
        self.rd_y = y
    
    def load_etrs_latlon(self, lat, lon):
        self.etrs_lat = lat
        self.etrs_lon = lon
        self.etrs_lat_rad = math.radians(lat)
        self.etrs_lon_rad = math.radians(lon)
    
    def load_height(self, h, ref):
        if ref == 'ETRS89':
            self.ellipsoidal_height = h
        elif ref == 'NAP':
            self.nap_height = h
        elif ref == 'LAT':
            self.lat_height = h
    
    def compute_etrs_latlon_from_rd_xy(self):
        self.bessel_lat_rad, self.bessel_lon_rad, self.bessel_iterations, self.bessel_error = self.inverse_rd_projection(self.rd_x, self.rd_y, constants.RD_SPHERICAL)
        self.bessel_lat_deg = math.degrees(self.bessel_lat_rad)
        self.bessel_lon_deg = math.degrees(self.bessel_lon_rad)
        
        self.pseudo_bessel_lat_deg, self.pseudo_bessel_lon_deg, self.pseudo_bessel_corrections = self.pseudo_bessel_correction(
            self.bessel_lat_deg, 
            self.bessel_lon_deg,
            constants.RDCORR2018
        )
        self.pseudo_bessel_lat_rad = math.radians(self.pseudo_bessel_lat_deg)
        self.pseudo_bessel_lon_rad = math.radians(self.pseudo_bessel_lon_deg)
        
        self.pseudo_bessel_x, self.pseudo_bessel_y, self.pseudo_bessel_z = self.geographic_to_cartesian(
            self.pseudo_bessel_lat_rad,
            self.pseudo_bessel_lon_rad,
            constants.BESSEL_H0,
            constants.BESSEL_SEMI_MAJOR,
            constants.BESSEL_SQ_ECCENTRICITY
        )
        
        self.etrs_x, self.etrs_y, self.etrs_z = self.similarity_transformation(
            self.pseudo_bessel_x,
            self.pseudo_bessel_y,
            self.pseudo_bessel_z,
            constants.BESSELETRS
        )
        
        self.etrs_lat_rad, self.etrs_lon_rad, self.etrs_iterations, self.etrs_error = self.cartesian_to_geographic(
            self.etrs_x,
            self.etrs_y,
            self.etrs_z,
            constants.GRS80_SEMI_MAJOR,
            constants.GRS80_SQ_ECCENTRICITY
        )
        self.etrs_lat = math.degrees(self.etrs_lat_rad)
        self.etrs_lon = math.degrees(self.etrs_lon_rad)
    
    def compute_rd_xy_from_etrs_latlon(self):
        self.etrs_x, self.etrs_y, self.etrs_z = self.geographic_to_cartesian(
            self.etrs_lat_rad,
            self.etrs_lon_rad,
            constants.ETRS_H0,
            constants.GRS80_SEMI_MAJOR,
            constants.GRS80_SQ_ECCENTRICITY
        )
        
        self.pseudo_bessel_x, self.pseudo_bessel_y, self.pseudo_bessel_z = self.similarity_transformation(
            self.etrs_x,
            self.etrs_y,
            self.etrs_z,
            constants.ETRSBESSEL
        )
        
        self.pseudo_bessel_lat_rad, self.pseudo_bessel_lon_rad, self.pseudo_bessel_iterations, self.pseudo_bessel_error = self.cartesian_to_geographic(
            self.pseudo_bessel_x,
            self.pseudo_bessel_y,
            self.pseudo_bessel_z,
            constants.BESSEL_SEMI_MAJOR,
            constants.BESSEL_SQ_ECCENTRICITY
        )
        self.pseudo_bessel_lat_deg = math.degrees(self.pseudo_bessel_lat_rad)
        self.pseudo_bessel_lon_deg = math.degrees(self.pseudo_bessel_lon_rad)
        
        self.bessel_lat_deg, self.bessel_lon_deg, self.bessel_iterations, self.bessel_corrections = self.bessel_correction(
            self.pseudo_bessel_lat_deg,
            self.pseudo_bessel_lon_deg,
            constants.RDCORR2018,
            constants.RDCORR2018_ITERATION_THRESHOLD
        )
        
        self.bessel_lat_rad = math.radians(self.bessel_lat_deg)
        self.bessel_lon_rad = math.radians(self.bessel_lon_deg)
        
        if self.bessel_iterations == 0:
            self.rd_correction_out_of_bounds = True

        self.rd_x, self.rd_y = self.rd_projection(
            self.bessel_lat_rad,
            self.bessel_lon_rad,
            constants.RD_SPHERICAL
        )
    
    def compute_separations(self, include_lat=True):
        self.geoid_ellipsoid_separation = self.get_geoid_ellipsoid_separation(self.etrs_lat, self.etrs_lon, constants.NLGEO2018)
        
        if include_lat:
            self.lat_ellipsoid_separation = self.get_lat_ellipsoid_separation(self.etrs_lat, self.etrs_lon, constants.NLLAT2018)
            if self.geoid_ellipsoid_separation is not None and self.lat_ellipsoid_separation is not None:
                self.geoid_lat_separation = self.geoid_ellipsoid_separation - self.lat_ellipsoid_separation
                
            if self.lat_height is not None and self.ellipsoidal_height is None and self.lat_ellipsoid_separation is not None:
                self.ellipsoidal_height = self.lat_height + self.lat_ellipsoid_separation
        
        if self.nap_height is not None and self.ellipsoidal_height is None and self.geoid_ellipsoid_separation is not None:
            self.ellipsoidal_height = self.nap_height + self.geoid_ellipsoid_separation
        elif self.ellipsoidal_height is not None and self.nap_height is None and self.geoid_ellipsoid_separation is not None:
            self.nap_height = self.ellipsoidal_height - self.geoid_ellipsoid_separation
            
        if include_lat and self.lat_height is None and self.lat_ellipsoid_separation is not None and self.ellipsoidal_height is not None:
            self.lat_height = self.ellipsoidal_height - self.lat_ellipsoid_separation
        
    #==============================================================================================#
    
    def get_corrections_at_index(self, index, grid):
        
        if grid['TABLE'] in Transformer.cache.keys():
            if index in Transformer.cache[grid['TABLE']].keys():
                return Transformer.cache[grid['TABLE']][index]
        else:
            Transformer.cache[grid['TABLE']] = dict()
        
        with sqlite3.connect(os.path.join(constants.RESOURCES, constants.DB_FILE)) as conn:
            sql = f'''
                SELECT lat_corr, lon_corr
                FROM {grid["TABLE"]}
                WHERE line = ?;
            '''
            cursor = conn.cursor()
            cursor.execute(sql, (index,))
            results = cursor.fetchall()[0]
            Transformer.cache[grid['TABLE']][index] = results
        return results
        
    def get_val_at_line(self, index, grid, col='geoid'):
        
        if grid['TABLE'] in Transformer.cache.keys():
            if index in Transformer.cache[grid['TABLE']].keys():
                return Transformer.cache[grid['TABLE']][index]
        else:
            Transformer.cache[grid['TABLE']] = dict()
        
        with sqlite3.connect(os.path.join(constants.RESOURCES, constants.DB_FILE)) as conn:
            sql = f'''
                SELECT {col}
                FROM {grid["TABLE"]}
                WHERE line = ?;
            '''
            cursor = conn.cursor()
            cursor.execute(sql, (index,))
            try:
                results = cursor.fetchall()[0][0]
            except:
                results = None
            Transformer.cache[grid['TABLE']][index] = results
            
        return results

    def get_geoid_ellipsoid_separation(self, phi, lmbda, grid):
        if not helpers.is_within_bounds(phi, lmbda, grid):
            return grid['C0']
        
        phi_norm = (phi - grid['PHI_MIN']) / grid['PHI_DELTA']
        lmbda_norm = (lmbda - grid['LMBDA_MIN']) / grid['LMBDA_DELTA']
            
        i_nw, i_ne, i_sw, i_se = helpers.get_grid_indices(phi_norm, lmbda_norm, grid)
        g_nw = self.get_val_at_line(i_nw, grid, 'geoid')
        g_ne = self.get_val_at_line(i_ne, grid, 'geoid')
        g_sw = self.get_val_at_line(i_sw, grid, 'geoid')
        g_se = self.get_val_at_line(i_se, grid, 'geoid')
                
        N = (phi_norm - math.floor(phi_norm)) * (g_nw * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_ne * (lmbda_norm - math.floor(lmbda_norm))) + (math.floor(phi_norm) + 1 - phi_norm) * (g_sw * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_se * (lmbda_norm - math.floor(lmbda_norm)))
            
        return N
    
    def get_lat_ellipsoid_separation(self, phi, lmbda, grid):
        if not helpers.is_within_bounds(phi, lmbda, grid):
            return grid['C0']
        
        phi_norm = (phi - grid['PHI_MIN']) / grid['PHI_DELTA']
        lmbda_norm = (lmbda - grid['LMBDA_MIN']) / grid['LMBDA_DELTA']
            
        i_nw, i_ne, i_sw, i_se = helpers.get_grid_indices(phi_norm, lmbda_norm, grid)
        g_nw = self.get_val_at_line(i_nw, grid, 'lat_ellipsoid')
        g_ne = self.get_val_at_line(i_ne, grid, 'lat_ellipsoid')
        g_sw = self.get_val_at_line(i_sw, grid, 'lat_ellipsoid')
        g_se = self.get_val_at_line(i_se, grid, 'lat_ellipsoid')
        
        if g_nw and g_ne and g_sw and g_se:
            N = (phi_norm - math.floor(phi_norm)) * (g_nw * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_ne * (lmbda_norm - math.floor(lmbda_norm))) + (math.floor(phi_norm) + 1 - phi_norm) * (g_sw * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_se * (lmbda_norm - math.floor(lmbda_norm)))
            return N
        
        return None
    
    #==============================================================================================#

    def geographic_to_cartesian(self, phi, lmbda, h, a, e2):
        Rn = helpers.second_curvature(phi, a, e2)
        cos_phi = math.cos(phi)
        
        X = (Rn + h) * cos_phi * math.cos(lmbda)
        Y = (Rn + h) * cos_phi * math.sin(lmbda)
        Z = (Rn * (1.0 - e2) + h) * math.sin(phi)
        
        return X, Y, Z
        
    def similarity_transformation(self, X, Y, Z, params):
        X2 = params['Scale'] * ( params['R11'] * X + params['R12'] * Y + params['R13'] * Z ) + params['Tx']
        Y2 = params['Scale'] * ( params['R21'] * X + params['R22'] * Y + params['R23'] * Z ) + params['Ty']
        Z2 = params['Scale'] * ( params['R31'] * X + params['R32'] * Y + params['R33'] * Z ) + params['Tz']
        
        return X2, Y2, Z2
    
    def cartesian_to_geographic(self, X, Y, Z, a, e2, threshold=constants.CTG_ITERATION_THRESHOLD, zero_rounding=constants.CTG_ROUNDING):
        Xr = round(X, zero_rounding)
        Yr = round(Y, zero_rounding)
        Zr = round(Z, zero_rounding)
        
        if Xr == 0.0 and Yr == 0.0:
            if Zr >= 0.0:
                phi = math.radians(90.0)
            else:
                phi = math.radians(-90.0)
            lmbda = 0.0
            
            return phi, lmbda, 0, 0.0
            
        else:
            phi_i = 0.0
            phi_prev = phi_i + 2*threshold
            i = 0
            while abs(phi_i - phi_prev) >= threshold:
                phi_prev = phi_i
                Rn = helpers.second_curvature(phi_i, a, e2)
                phi_i = math.atan( (Z + e2 * Rn * math.sin(phi_i)) / math.sqrt(X**2 + Y**2) )
                i += 1
                
            if Xr > 0.0:
                lmbda = math.atan(Y/X)
            elif Xr < 0.0 and Yr >= 0.0:
                lmbda = math.atan(Y/X) + math.radians(180.0)
            elif Xr < 0.0 and Yr <= 0.0:
                lmbda = math.atan(Y/X) - math.radians(180.0)
            elif Xr == 0.0 and Yr > 0.0:
                lmbda = math.radians(90.0)
            elif Xr == 0.0 and Yr < 0.0:
                lmbda = math.radians(-90.0)
                
            return phi_i, lmbda, i, abs(phi_i - phi_prev)
        
    def bessel_correction(self, phi, lmbda, grid, threshold=constants.RDCORR2018_ITERATION_THRESHOLD):    
        # from etrs to rd
        if not helpers.is_within_bounds(phi, lmbda, grid):
            return phi, lmbda, 0, (grid['C0'], grid['C0']) 
    
        i = 0
        phi_i = phi
        lmbda_i = lmbda
        phi_prev = phi_i + 10*threshold
        lmbda_prev = lmbda_i + 10*threshold
        
        while abs(phi_i - phi_prev) >= threshold and abs(lmbda_i - lmbda_prev) >= threshold:
            phi_prev = phi_i
            lmbda_prev = lmbda_prev
        
            phi_norm = (phi_i - grid['PHI_MIN']) / grid['PHI_DELTA']
            lmbda_norm = (lmbda_i - grid['LMBDA_MIN']) / grid['LMBDA_DELTA']
            
            i_nw, i_ne, i_sw, i_se = helpers.get_grid_indices(phi_norm, lmbda_norm, grid)
            g_nw = self.get_corrections_at_index(i_nw, grid)
            g_ne = self.get_corrections_at_index(i_ne, grid)
            g_sw = self.get_corrections_at_index(i_sw, grid)
            g_se = self.get_corrections_at_index(i_se, grid)
            
            c_phi = (phi_norm - math.floor(phi_norm)) * (g_nw[0] * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_ne[0] * (lmbda_norm - math.floor(lmbda_norm))) \
                    + (math.floor(phi_norm) + 1 - phi_norm) * (g_sw[0] * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_se[0] * (lmbda_norm - math.floor(lmbda_norm)))
            c_lmbda = (phi_norm - math.floor(phi_norm)) * (g_nw[1] * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_ne[1] * (lmbda_norm - math.floor(lmbda_norm))) \
                    + (math.floor(phi_norm) + 1 - phi_norm) * (g_sw[1] * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_se[1] * (lmbda_norm - math.floor(lmbda_norm)))
            
            phi_i = phi - c_phi
            lmbda_i = lmbda - c_lmbda     
            i += 1
        
        return phi_i, lmbda_i, i, (c_phi, c_lmbda)
    
    def pseudo_bessel_correction(self, phi, lmbda, grid):
        # from rd to etrs
        if not helpers.is_within_bounds(phi, lmbda, grid):
            return phi, lmbda, (grid['C0'], grid['C0'])
            
        phi_norm = (phi - grid['PHI_MIN']) / grid['PHI_DELTA']
        lmbda_norm = (lmbda - grid['LMBDA_MIN']) / grid['LMBDA_DELTA']
        i_nw, i_ne, i_sw, i_se = helpers.get_grid_indices(phi_norm, lmbda_norm, grid)
        g_nw = self.get_corrections_at_index(i_nw, grid)
        g_ne = self.get_corrections_at_index(i_ne, grid)
        g_sw = self.get_corrections_at_index(i_sw, grid)
        g_se = self.get_corrections_at_index(i_se, grid)
        
        c_phi = (phi_norm - math.floor(phi_norm)) * (g_nw[0] * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_ne[0] * (lmbda_norm - math.floor(lmbda_norm))) \
                + (math.floor(phi_norm) + 1 - phi_norm) * (g_sw[0] * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_se[0] * (lmbda_norm - math.floor(lmbda_norm)))
        c_lmbda = (phi_norm - math.floor(phi_norm)) * (g_nw[1] * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_ne[1] * (lmbda_norm - math.floor(lmbda_norm))) \
                + (math.floor(phi_norm) + 1 - phi_norm) * (g_sw[1] * (math.floor(lmbda_norm) + 1 - lmbda_norm) + g_se[1] * (lmbda_norm - math.floor(lmbda_norm)))
                
        phi += c_phi
        lmbda += c_lmbda
        
        return phi, lmbda, (c_phi, c_lmbda)
    
    def rd_projection(self, phi, lmbda, params, zero_rounding=9):
        
        if math.degrees(phi) == 90.0 or math.degrees(phi) == -90.0:
            phi_spherical = phi
        else:
            q = math.log(math.tan( (phi + math.radians(90.0)) / 2.0 )) - (params['E'] / 2) * math.log((1 + params['E']*math.sin(phi))/(1 - params['E']*math.sin(phi)))
            w = params['N'] * q + params['M']
            phi_spherical = 2 * math.atan(math.exp(w)) - math.radians(90.0)
            
        lmbda_spherical = params['L0_SPHERICAL'] + params['N'] * (lmbda - params['L0_BESSEL'])
                
        if round(phi_spherical, zero_rounding) == params['P0_SPHERICAL'] and round(lmbda_spherical, zero_rounding) == params['L0_SPHERICAL']:
            x = params['X0']
            y = params['Y0']
        elif round(phi_spherical, zero_rounding) == -1*params['P0_SPHERICAL'] and round(lmbda_spherical, zero_rounding) == params['L0_SPHERICAL'] - math.radians(180.0):
            x = None
            y = None
        else:
            sinpsi = math.sqrt((math.sin((phi_spherical - params['P0_SPHERICAL'])/2))**2 + math.sin((lmbda_spherical - params['L0_SPHERICAL'])/2)**2 * math.cos(phi_spherical) * math.cos(params['P0_SPHERICAL']))
            cospsi = math.sqrt(1 - sinpsi**2)
            tanpsi = sinpsi / cospsi
            
            sina = (math.sin(lmbda_spherical - params['L0_SPHERICAL']) * math.cos(phi_spherical)) / (2 * sinpsi * cospsi)
            cosa = ( math.sin(phi_spherical)  - math.sin(params['P0_SPHERICAL']) + 2*math.sin(params['P0_SPHERICAL'])*sinpsi**2 ) / (2 * math.cos(params['P0_SPHERICAL']) * sinpsi * cospsi )
            r = 2 * params['K'] * params['R'] * tanpsi
            
            x = r * sina + params['X0']
            y = r * cosa + params['Y0']
       
        return x, y
    
    def inverse_rd_projection(self, x, y, params, threshold=constants.RD_INV_PROJECTION_THRESHOLD):
        
        r = math.sqrt( (x - params['X0'])**2 + (y - params['Y0'])**2 )
        sina = (x - params['X0']) / r
        cosa = (y - params['Y0']) / r
        psi = 2 * math.atan(r / (2 * params['K'] * params['R']))
        
        if round(x, 4) == params['X0'] and round(y, 4) == params['Y0']:
            Xn = math.cos(params['P0_SPHERICAL'])
            Yn = 0.0
            Zn = math.sin(params['P0_SPHERICAL'])
        else:
            Xn = math.cos(params['P0_SPHERICAL']) * math.cos(psi) - cosa * math.sin(params['P0_SPHERICAL']) * math.sin(psi)
            Yn = sina * math.sin(psi)
            Zn = cosa * math.cos(params['P0_SPHERICAL']) * math.sin(psi) + math.sin(params['P0_SPHERICAL']) * math.cos(psi)
            
        phi_spherical = math.asin(Zn)
        if round(Xn, 11) > 0.0:
            lmbda_spherical = params['L0_SPHERICAL'] + math.atan(Yn / Xn)
        elif round(Xn, 11) < 0.0 and round(x, 4) >= params['X0']:
            lmbda_spherical = params['L0_SPHERICAL'] + math.atan(Yn / Xn) + math.radians(180.0)
        elif round(Xn, 11) < 0.0 and round(x, 4) < params['X0']:
            lmbda_spherical = params['L0_SPHERICAL'] + math.atan(Yn / Xn) - math.radians(180.0)
        elif round(Xn, 11) == 0.0 and round(x, 4) > params['X0']:
            lmbda_spherical = params['L0_SPHERICAL'] + math.radians(90.0)
        elif round(Xn, 11) == 0.0 and round(x, 4) < params['X0']:
            lmbda_spherical = params['L0_SPHERICAL'] - math.radians(90.0)
        elif round(Xn, 11) == 0.0 and round(x, 4) == params['X0']:
            lmbda_spherical = params['L0_SPHERICAL']
        
        # second part of projection
        if phi_spherical == math.radians(90.0) or phi_spherical == math.radians(-90.0):
            phi_i = phi_spherical
        else:
            w = math.log( math.tan( (phi_spherical + math.radians(90)) / 2.0 ) )
            q = ( w - params['M'] ) / params['N']
            phi_i = 0.0
            phi_prev = phi_i + 10*threshold
            
            i = 0
            while abs(phi_i - phi_prev) > threshold:
                phi_prev = phi_i
                phi_i = 2 * math.atan(math.exp(q + (params['E']/2.0) * math.log((1 + params['E'] * math.sin(phi_i)) / (1 - params['E'] * math.sin(phi_i))))) - math.radians(90.0)   
                i += 1
                
        lmbda_n = (lmbda_spherical - params['L0_SPHERICAL']) / params['N'] + params['L0_BESSEL']
        lmbda = lmbda_n + math.radians(360.0) * math.floor((math.radians(180.0) - lmbda_n) / math.radians(360.0))        
        
        return phi_i, lmbda, i, abs(phi_i - phi_prev)
    
    #==============================================================================================#

def transform(x, y, h=None, src_xy=None, trg_xy=None, src_h=None, trg_h=None, full=False):
    
    # VALIDATION
    if src_xy is None:
        raise ValueError("src_xy is always needed!")
    if h is not None and src_h is None:
        raise ValueError("if h is supplied, define a src_h as well!")
    
    t = Transformer()
    r = [None, None, None]
    
    # LOAD
    if src_xy == 'RD':
        t.load_rd_xy(x, y)
    elif src_xy == 'ETRS89':
        t.load_etrs_latlon(y, x)
        
    if h is not None and src_h is not None:
        t.load_height(h, src_h)

    # TRANSFORM
    if trg_xy == 'ETRS89' and src_xy == 'RD':
        t.compute_etrs_latlon_from_rd_xy()
    elif trg_xy == 'RD' and src_xy == 'ETRS89':
        t.compute_rd_xy_from_etrs_latlon()
        
    if trg_h and trg_h != src_h:
        if t.etrs_lat is None:
            # we need ETRSlatlon for all heights, if no trg_xy is defined, 
            t.compute_etrs_latlon_from_rd_xy()
        
        if (trg_h is not None and 'LAT' in trg_h) or (src_h is not None and 'LAT' in src_h):
            t.compute_separations(include_lat=True)
        else:
            t.compute_separations(include_lat=False)
        
    # EXTRACT
    if trg_xy == 'RD':
        r[0], r[1] = t.rd_x, t.rd_y
    elif trg_xy == 'ETRS89':
        r[0], r[1] = t.etrs_lon, t.etrs_lat
        
    if trg_h == 'ETRS89':
        r[2] = t.ellipsoidal_height
    elif trg_h == 'NAP':
        r[2] = t.nap_height
    elif trg_h == 'LAT':
        r[2] = t.lat_height
    elif trg_h == 'NAP_LAT_SEPARATION':
        r[2] = t.geoid_lat_separation
    else:
        r[2] = h
    
    return r[0], r[1], r[2]
    
def rd_to_etrs(x, y):
    t = Transformer()
    t.load_rd_xy(x, y)
    t.compute_etrs_latlon_from_rd_xy()
    return t.etrs_lat, t.etrs_lon
    
def etrs_to_rd(lat, lon):
    t = Transformer()
    t.load_etrs_latlon(lat, lon)
    t.compute_rd_xy_from_etrs_latlon()
    return t.rd_x, t.rd_y
    
def rdnap_to_etrsh(x, y, h):
    t = Transformer()
    t.load_rd_xy(x, y)
    t.load_height(h, 'NAP')
    t.compute_etrs_latlon_from_rd_xy()
    t.compute_separations(include_lat=False)
    return t.etrs_lat, t.etrs_lon, t.ellipsoidal_height

def rdnap_to_etrslat(x, y, h):
    t = Transformer()
    t.load_rd_xy(x, y)
    t.load_height(h, 'NAP')
    t.compute_etrs_latlon_from_rd_xy()
    t.compute_separations(include_lat=True)
    return t.etrs_lat, t.etrs_lon, t.lat_height
    
def etrsh_to_rdnap(lat, lon, h):
    t = Transformer()
    t.load_etrs_latlon(lat, lon)
    t.load_height(h, 'ETRS89')
    t.compute_rd_xy_from_etrs_latlon()
    t.compute_separations(include_lat=False)
    return t.rd_x, t.rd_y, t.nap_height