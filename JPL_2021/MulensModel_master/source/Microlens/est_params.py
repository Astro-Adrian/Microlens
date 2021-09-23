import sys
import os

local_path='/Users/adrianhernandez/JPL_2021/'
sys.path.append(local_path + 'MulensModel_master/source')
#
import MulensModel as mm
import numpy as np
import matplotlib.pyplot as pl
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
#
from glob import glob
dir_casu = '/Users/adrianhernandez/JPL_2021/ukirt_casu/'
casu_data = glob(f"{dir_casu}/*C58505*.txt")
#
# # fln = casu_data[0][44:77]
# # print(fln)
ukirt_casu_H = np.loadtxt(casu_data[0], usecols=range(3)) #H band data
ukirt_casu_K = np.loadtxt(casu_data[2], usecols=range(3)) # K band data
H_data = mm.MulensData(file_name=casu_data[0]) #used only for the mulens code
K_data = mm.MulensData(file_name=casu_data[2]) #used only for the mulens code

#
ex_start = 2458658.
ex_stop = 2458690.
start = 2458658.
stop = 2458690.

# start = 2458660.
# stop = 2458690.
# # magnification_methods = [start, 'VBBL', stop]

magnification_methods = [2457820, 'point_source_point_lens', start, 'VBBL', stop,
    'point_source_point_lens', 2458750.]

  # Flag data related to the planet
flag_planet = (K_data.time > ex_start) & (K_data.time < ex_stop) | np.isnan(K_data.err_mag)

# Exclude those data from the fitting (for no
K_data.bad = flag_planet

time_flag = (K_data.time[np.invert(flag_planet)])
mag_flag = (K_data.mag[np.invert(flag_planet)])
err_flag = (K_data.err_mag[np.invert(flag_planet)])


#Estimate for t0
index_t_0 = np.argmin(K_data.mag)
t_o = K_data.time[index_t_0]

#Estimate for u0
baseline_mag = np.min([K_data.mag[0],K_data.mag[-1]]) # A crude estimate
A_max =10.**((K_data.mag[index_t_0] - baseline_mag) / -2.5)
u_o = 1. / A_max # True in the high-magnification limit

#Estimate for tE
t_1 = np.interp( baseline_mag - 0.3, K_data.mag[index_t_0:0:-1], K_data.time[index_t_0:0:-1])
t_e = np.abs((t_o - t_1) / np.sqrt(1. - u_o**2))


# Approximate time of the planetary perturbation
t_planet = (stop + start) / 2.

# Position of the source at the time of the planetary perturbation
# tau = ((t_planet - (t_o+2450000)) /t_e)
# u_planet = np.sqrt(u_o**2 + tau**2)
# Position of the lens images at the time of the planetary perturbation
# --> Estimate of the planet location
# s_minus = 0.5 * (np.sqrt(u_planet**2 + 4.) - u_planet)
# s_plus = 0.5 * (np.sqrt(u_planet**2 + 4.) + u_planet)
#
# alpha_est = np.rad2deg(-np.arctan2(u_o, tau))
# print(alpha_est)


# pspl_model= mm.Model({'t_0': t_o, 'u_0':u_o, 't_E':t_e})
# pspl_event = mm.Event(datasets=[da.K_data,da.H_data], model=pspl_model)
#
# pl.figure(figsize=(15,10))
# pspl_event.plot_model(t_range = [2457900,2458800],subtract_2450000=True,color='black')
# pspl_event.plot_data(subtract_2450000=True,show_bad=True)
# pl.xlim(8620,8690)
# pl.show()
