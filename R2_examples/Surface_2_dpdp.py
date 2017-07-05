# -*- coding: utf-8 -*-
"""
Created on Thu May 25 08:06:54 2017

Example of a pyres implementation of "Surface electorde array 2" from the 
R2_v3.1_readme.pdf

@author: kbefus
"""
from __future__ import print_function
import os
import numpy as np

# Setup pyplot for saving figure to pdf
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

from pyres import mesh_tools, mesh_utils
from pyres import r2_tools
from pyres import pyres_utils, plot_utils

#%%
# Define directory where example scripts are located
ex_dir = r'C:\ER' 
work_dir = os.path.join(ex_dir,r'R2_examples\Surface_2\dpdp\Forward')
inv_work_dir = os.path.join(ex_dir,r'R2_examples\Surface_2\dpdp\Inverse')

if not os.path.isdir(work_dir):
    os.makedirs(work_dir)

# Define directory with R2
exe_dir = os.path.join(ex_dir,r'R2\Distribition_3.1\Progs')

# -------- Create mesh -------------
n_electrodes = 25.
e_spacing = 2. # meters is unit of distance for this example
ndivx,ndivy = 8.,8. # mesh divisions between electrodes
dx_foreground = e_spacing/ndivx # 0.25 m
dy_foreground = e_spacing/ndivy # 0.25 m  
xbuff,ybuff = 2.,10. # meters to extend foreground beyond electrode x,y locations
dx_background = [0.5,1., 2., 5., 10., 20., 50., 100.]
dy_background = [0.5,1., 2., 5., 10., 20., 50., 100.]

# Define mesh nodes
electrode_x = e_spacing*np.arange(n_electrodes) # first electrode at x=0, last electrode at x=48
xx_mesh = np.arange(electrode_x[0]-xbuff,electrode_x[-1]+xbuff+dx_foreground,dx_foreground)
xx_mesh = np.unique(np.hstack([(xx_mesh[0]-np.cumsum(dx_background))[::-1],
                               xx_mesh,
                               xx_mesh[-1]+np.cumsum(dx_background),
                               electrode_x]))

yy_mesh = np.arange(0,ybuff+dy_foreground,dy_foreground)
yy_mesh = np.hstack([yy_mesh,yy_mesh[-1]+np.cumsum(dy_background)])

# Define topography
topo_func = lambda m,x,b: m*x+b
topog = topo_func(4.8/electrode_x[-1],xx_mesh,0.)

# Define mesh creation input
node_dict = {'xx':xx_mesh,'yy':yy_mesh,'topog':topog}

# Define electrode rows and columns
e_rows = np.ones_like(electrode_x) # all in top row
e_cols = np.array([(xx_mesh==ix).nonzero()[0][0] for ix in electrode_x])+1 # find columns of electrode nodes
electrode_array = np.column_stack([np.arange(e_cols.shape[0])+1,e_cols,e_rows]).astype(int)
# mesh_ax = plot_utils.plot_quad_mesh([xx_mesh,yy_mesh],label_nodes=False) # uncomment to plot mesh

# ----------------- Define target -----------------
main_res = 100 # ohm m
target_res = 10 # ohm m
target_xypos = [[14,1],[16,1],
                [16,4],[14,4]] # do not want a closed polygon
region_dict = {'mesh_xy':[xx_mesh,yy_mesh],'target_xypos':target_xypos,
               'target_res':target_res,'background_res':main_res}
region_elems = mesh_utils.quad_mesh_region(**region_dict)


# ----------------- Setup R2 inputs ----------------
job_type = 0 # 0=forward, 1=inverse model
r2_dict = {'survey_name':'Surface_2_dpdp',
           'survey_dir': work_dir, 'job_type':job_type,
           'exe_dir':exe_dir,
           'exe_name':'R2.exe'}
iR2 = r2_tools.R2(**r2_dict)

mesh_type = 4 # 4= simple quadrilateral element mesh
output_xy = [[electrode_x[0],0.],[electrode_x[-1],topog[xx_mesh==48.]],
             [electrode_x[-1],-8.],[electrode_x[0],-8.],
             [electrode_x[0],0.]] # closed polygon outlining main domain
r2_in_dict = {'electrode_array':electrode_array,
              'output_domain':output_xy,'mesh_type':mesh_type,
              'reg_elems':region_elems,
              'res_matrix':0,
              'singular_type':0,'node_dict':node_dict}

# Make synthetic electrode combinations for forward model
synth_dict = {'n_electrodes':n_electrodes,'array_type':'dpdp',
              'max_dipole':1,'max_separation':7}
protocol_data = pyres_utils.make_synthetic_survey(**synth_dict)
protocol_data = protocol_data[:,np.array([0,3,4,1,2])] # reorganize to match example
protocol_dict = {'meas_data':protocol_data}

# ------------- Run R2 forward model ----------------
run_r2_dict = {'protocol_dict':protocol_dict, 'r2_in_dict':r2_in_dict,
                'run_bool':True}

iR2.run_all(**run_r2_dict)


#%% Invert data

if not os.path.isdir(inv_work_dir):
    os.makedirs(inv_work_dir)

inv_r2_dict = {}
inv_r2_dict.update(r2_dict) # take inputs from forward model
inv_r2_dict.update({'job_type':1,
                    'survey_dir': inv_work_dir}) # run inverse model and save to new directory

iR2_inv = r2_tools.R2(**inv_r2_dict)

# Assign inversion options when different from r2_tools.inv_defaults
inv_options = {'a_wgt':1e-3,'b_wgt':2e-2,
               'inverse_type':1,'patch_size_xy':[4,4],
               'rho_min':-1e3,'rho_max':1e3,}

inv_r2_in_dict = {}
inv_r2_in_dict.update(r2_in_dict)
inv_r2_in_dict.update({'reg_elems':[[1,int((xx_mesh.shape[0]-1)*(yy_mesh.shape[0]-1)),main_res]],
                        'res_matrix':1,
                        'inv_dict':inv_options})

# Load forward model data
meas_fwd_data = pyres_utils.load_fwd_output(work_dir=work_dir)

inv_protocol_dict = {'meas_data':meas_fwd_data}
run_inv_r2_dict = {'protocol_dict':inv_protocol_dict, 'r2_in_dict':inv_r2_in_dict,
                   'run_bool':True}

iR2_inv.run_all(**run_inv_r2_dict)

#%%
# Plot inverse model results
plot_dict = {'work_dir':inv_work_dir,
             'topog_xy':np.column_stack([xx_mesh,topog]),
             'plt_opts':{'vmin':10**1.3,'vmax':10.**2},
             'invert_y':False,
             'cmap':'rainbow','keep_log':True}
fig,ax,[X,Y,ER] = plot_utils.plot_res(**plot_dict)
ax.plot(electrode_x,topo_func(4.8/electrode_x[-1],electrode_x,0.),'k.-')

# Plot target
txy = np.array(target_xypos)
txy = np.vstack([txy,txy[0,:]])
topo_y = pyres_utils.extrap(txy[:,0],xx_mesh,topog)
ax.plot(txy[:,0],topo_y-txy[:,1],'k-')

#%% Compare performance of inversion

# Define normalized misfit equation
errfunc = lambda obs,mod: (obs-mod)/obs

obs_array = np.log10(main_res)*np.ones_like(ER)
target_bool = (X>=txy[0,0]) & (X<=txy[1,0]) &\
              (Y<=(topo_y[0]-txy[0,1])) & (Y>=(topo_y[2]-txy[2,1]))
obs_array[target_bool] = np.log10(target_res)

err_array = errfunc(10**obs_array,10**ER)
# get labels for target area
label_xy = np.column_stack([X[target_bool],Y[target_bool]])
target_err = err_array[target_bool]

ratio=1
im2=ax.pcolormesh(X,Y,err_array,vmin=-ratio,vmax=ratio,cmap='RdBu')

#for err,(ex,ey) in zip(target_err,label_xy):
#    ax.text(ex,ey,'{}'.format(err))

plt.colorbar(im2,extend='both')
                          
                          
