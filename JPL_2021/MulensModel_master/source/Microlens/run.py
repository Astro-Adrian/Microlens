from mcmcFit import mcmcFit
import numpy as np
import time
start_time = time.time()

from glob import glob
im_dir = '/Users/adrianhernandez/JPL_2021/ogle_files/'
ogle_file = glob(f"{im_dir}/OGLE_*.dat")
#print(ogle_file[3]) #change indicies to update the file you are interested with

ogle_data = np.loadtxt(ogle_file[7])

name = 'OGLE_190651'

#print ogle_data

ogle_time = ogle_data[:,0]-2450000 #Time
reg_ogle_time = ogle_data[:,0]
ogle_mag = ogle_data[:,1]
ogle_error = ogle_data[:,2]
data_points = len(ogle_data)

initial_u0 = 0.78
initial_t0 = 2458610.06
initial_tE = 5.52

u0,t0,tE,Ftot,fb,u0err,t0err,tEerr,Ftoterr,fberr = mcmcFit(name,ogle_time,ogle_mag,ogle_error,initial_u0,initial_t0,initial_tE)

print('u0, error',u0,u0err)
print('t0, error',t0,t0err)
print('tE, error',tE,tEerr)
print('Ftot, error',Ftot,Ftoterr)
print('Fb, error',fb,fberr)

print("--- %s seconds ---" % (time.time() - start_time))
