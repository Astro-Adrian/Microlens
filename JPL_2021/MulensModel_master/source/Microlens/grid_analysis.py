import matplotlib.pyplot as pl
import numpy as np
from glob import glob
from matplotlib import cm
dir = '/Users/adrianhernandez/JPL_2021/MulensModel_master/source/Microlens/'
data = glob(f"{dir}/*grid_50*.txt")

print(data)

grid_data = np.loadtxt(data[0])


import matplotlib.mlab as mlab
import matplotlib.colors

chi2 = np.array(grid_data[:,0])
s = np.array(grid_data[:,1])
q = np.array(grid_data[:,2])
t0 = 8653.49134
u0 = 0.047196
tE = 37.9749
planet_df_title = 'best_8_model_params.csv'
pspl_df_title = 'pspl_values.csv'
pspl_params = np.column_stack([t0+2450000,u0,tE])

colors = ['magenta', 'green', 'cyan','yellow','blue', 'white','orange']
cmap = cm.get_cmap(name='viridis', lut=None)
index_sorted = np.argsort(chi2)
least_chi = np.where(np.min(chi2)==chi2)
norm= matplotlib.colors.Normalize(vmax=np.max(chi2), vmin=np.min(chi2))
fig, ax = pl.subplots(figsize=(10,8))
tcf = ax.tricontourf(s,q, chi2,cmap = cmap, norm = norm, linestyles = 'dotted')
#sm = pl.cm.ScalarMappable(norm=norm, cmap = tcf.cmap)
sm = pl.cm.ScalarMappable(norm=norm, cmap = cmap)
sm.set_array([])
fig.colorbar(sm,ticks=tcf.levels)
ax.set_title('Refined Grid Search for s > 0')
pl.scatter(grid_data[least_chi,1],grid_data[least_chi,2],marker='o', color='red')
for j, index in enumerate(index_sorted[1:8]):
    pl.scatter(grid_data[index,1],grid_data[index,2],marker='.', color=colors[j - 1])
pl.xlabel('log s')
pl.ylabel('log q')
pl.show()

import pandas as pd

best_8 = index_sorted[0:8]
best_models = grid_data[best_8]
planet_df = pd.DataFrame(best_models, columns = ['chi2','s','q','alpha','rho'])
pspl_df = pd.DataFrame(pspl_params, columns = ['t0','u0','tE'])
planet_df.to_csv(planet_df_title,index=False)
pspl_df.to_csv(pspl_df_title, index=False)
