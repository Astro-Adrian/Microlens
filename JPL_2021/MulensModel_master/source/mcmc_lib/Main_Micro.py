#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from grid_fitter import load_data as dat


# In[ ]:


t0_est = dat.to_test
u0_est = dat.uo_test
tE_est = dat.tE_test 
print(t0_est)
print(u0_est)
print(tE_est)


# In[ ]:


import mcmcFit as mcfit
filename = 'test'
u0,t0,tE,Ftot,fb,u0err,t0err,tEerr,Ftoterr,fberr = mcfit.mcmcFit(filename,dat.time_flag -2450000, dat.mag_flag, dat.err_flag,u0_est,t0_est-2450000,tE_est)

print(t0,u0,tE)


# In[ ]:


import MulensModel as mm
from MulensModel.utils import Utils
import matplotlib.pyplot as pl

mcmc_model= mm.Model({'t_0': t0+2450000, 'u_0':u0, 't_E':tE})
mcmc_event = mm.Event(datasets=[dat.K_data,dat.H_data], model=mcmc_model)

# pl.figure(figsize=(15,10)) 
# mcmc_event.plot_model(t_range = [2457900,2458800],subtract_2450000=True,color='black')
# mcmc_event.plot_data(subtract_2450000=True,show_bad=True)
# pl.xlim(8600,8800)
# pl.show()


# In[ ]:


import numpy as np
dat.K_data.bad = np.isnan(dat.K_data.err_mag)


# In[ ]:


from grid_fitter import grid_code


# In[ ]:


grid_code.grid_fit(t0+2450000,u0,tE,-10)


# In[ ]:





# In[ ]:





# In[ ]:




