# import sys
# import os
# import numpy as np
# import MulensModel as mm
# from glob import glob
# import astropy.units as u
# from astropy.coordinates import SkyCoord
# from MulensModel.utils import Utils
#
# #This path should point to the MulensModel
# local_path='/Users/adrianhernandez/JPL_2021/'
# sys.path.append(local_path + 'MulensModel_master/source')
#
# #Import data
# #This path should point to the folder that contains your data
# dir_casu = '/Users/adrianhernandez/JPL_2021/ukirt_casu/'
# #Glob has the ability to iterate throuhg lots of text files but specify which file
# #I want it to grab by stating *C58505*.txt
# casu_data = glob(f"{dir_casu}/*C58505*.txt")
#
# #Load data
# ukirt_casu_H = np.loadtxt(casu_data[0], usecols=range(3))
# ukirt_casu_K = np.loadtxt(casu_data[2], usecols=range(3))
#
# #Arrays for H and K bands.
# H_data = mm.MulensData(file_name=casu_data[0])
# K_data = mm.MulensData(file_name=casu_data[2])
#
#
# start = 2458658.
# stop = 2458670.
# # magnification_methods = [start, 'VBBL', stop]
#
# magnification_methods = [2457820, 'point_source_point_lens', start, 'VBBL', stop,
#     'point_source_point_lens', 2458750.]
#
# #Flag data related to the planet
# flag_planet = (K_data.time > start) & (K_data.time < stop) | np.isnan(K_data.err_mag)
#
# # Exclude those data from the fitting (for no
# K_data.bad = flag_planet
#
# time_flag = (K_data.time[np.invert(flag_planet)])
# mag_flag = (K_data.mag[np.invert(flag_planet)])
# err_flag = (K_data.err_mag[np.invert(flag_planet)])
