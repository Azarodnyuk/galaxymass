# -*- coding: utf-8 -*-

import numpy as np
from astropy.cosmology import WMAP9 as cosmo
from astropy import constants as const



def to_xyz_coordinates(data,name_dec_col,name_ra_col):
    
    """The function gets spherical angles ra and dec and gives 
       for each point on a identical sphere its x,y,z coordinates"""
    
    dec = data[name_dec_col]*np.pi/180
    ra = data[name_ra_col]*np.pi/180
    
    x = np.cos(dec)*np.cos(ra)
    y = np.cos(dec)*np.sin(ra)
    z = np.sin(dec)
    
    return np.array([x,y,z])

def linear_angle(list_xyz_0, list_xyz):
    
    """For each pair of points on an idetical sphere the function 
          allows us to obtain the flat angle between them"""
    
    cos_a = sum(list_xyz_0*list_xyz)
    
    return np.arccos(cos_a)

def data_query_id(data, clust_id, name_column_id):
    
    """The function allows us to obtain the query 
         of our data for particular cluster ID"""
    
    C = const.c.to('km/s')
    
    data_query = data[data[name_column_id] == clust_id][['RAJ2000_gal',   'DEJ2000_gal',   'z_gal', 
                                                         'RAJ2000_group', 'DEJ2000_group', 'z_group', 
                                                         'iGrID', 'DL']]
    

    xyz_0 = to_xyz_coordinates(data_query,'DEJ2000_group','RAJ2000_group')
    xyz_gal = to_xyz_coordinates(data_query,'DEJ2000_gal','RAJ2000_gal')

    a = linear_angle(xyz_0, xyz_gal)

    data_query['r_pr'] = cosmo.angular_diameter_distance(data_query.z_group)*a
    data_query['v'] = C*(data_query['z_gal'] - data_query['z_group'])
    
    return data_query