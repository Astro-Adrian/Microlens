#!/usr/bin/env python
# coding: utf-8
class color:
   PURPLE = '\033[1;35;48m'
   CYAN = '\033[1;36;48m'
   BOLD = '\033[1;37;48m'
   BLUE = '\033[1;34;48m'
   GREEN = '\033[1;32;48m'
   YELLOW = '\033[1;33;48m'
   RED = '\033[1;31;48m'
   BLACK = '\033[1;30;48m'
   UNDERLINE = '\033[4;37;48m'
   END = '\033[1;37;0m'

#import matplotlib.pyplot as pl
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
import MulensModel as mm
from MulensModel.utils import Utils
import est_params as dat

t0_est = dat.t_o
u0_est = dat.u_o
tE_est = dat.t_e
print(color.BLUE +'t0 estimate = '+ color.END,t0_est)
print(color.BLUE +'u0 estimate = '+ color.END,u0_est)
print(color.BLUE +'tE estimate = '+ color.END,tE_est)


print(color.PURPLE +'Beginning MCMC'+ color.END)
import mcmcFit as mcfit
filename = 'test'
u0,t0,tE,Ftot,fb,u0err,t0err,tEerr,Ftoterr,fberr = mcfit.mcmcFit(filename,dat.time_flag -2450000, dat.mag_flag, dat.err_flag,u0_est,t0_est-2450000,tE_est)

# mcmc_model= mm.Model({'t_0': t0+2450000, 'u_0':u0, 't_E':tE})
# mcmc_event = mm.Event(datasets=[dat.K_data,dat.H_data], model=mcmc_model)
# pl.figure(figsize=(15,10))
# mcmc_event.plot_model(t_range = [2457900,2458800],subtract_2450000=True,color='black')
# mcmc_event.plot_data(subtract_2450000=True,show_bad=True)
# pl.xlim(8600,8800)
# pl.show()

#Estimate Planetary parameters
print(color.PURPLE +'Estimating the planetary parameter alpha'+ color.END)
tau = ((dat.t_planet - (t0+2450000)) /tE)
alpha_est = np.rad2deg(np.arctan2(u0, tau))
print(color.BLUE +'alpha planet = '+ color.END,alpha_est)

#unflag planet data
dat.K_data.bad = np.isnan(dat.K_data.err_mag)

print(color.PURPLE +'Beginning grid fitting'+ color.END)
import grid_code as gc
gc.grid_fit(t0+2450000,u0,tE,alpha_est)
