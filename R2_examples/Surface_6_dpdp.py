# -*- coding: utf-8 -*-
"""
Created on Thu May 25 08:06:54 2017

Example of a pyres implementation of "Surface electorde array 6" from the 
R2_v3.1_readme.pdf

@author: kbefus
"""
from __future__ import print_function
import os
import numpy as np

from pyres import mesh_tools, mesh_utils
from pyres import r2_tools
from pyres import pyres_utils, plot_utils

#%%
# Define directory where example scripts are located
ex_dir = r'C:\ER' 
work_dir = os.path.join(ex_dir,r'R2_examples\Surface_1\dpdp\Forward')
inv_main_work_dir = os.path.join(ex_dir,r'R2_examples\Surface_6\dpdp\case_{0:02.0f}\Inverse')

# Define directory with R2
exe_dir = os.path.join(ex_dir,r'R2\Distribition_3.1\Progs')

alpha_s = [10, 50, 10]
rho_back = [100, 100, 50]


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

# Define mesh creation input
node_dict = {'xx':xx_mesh,'yy':yy_mesh,'topog':None}

# Define electrode rows and columns
e_rows = np.ones_like(electrode_x) # all in top row
e_cols = np.array([(xx_mesh==ix).nonzero()[0][0] for ix in electrode_x])+1 # find columns of electrode nodes
electrode_array = np.column_stack([np.arange(e_cols.shape[0])+1,e_cols,e_rows]).astype(int)

out_inv_dirs = []

# Loop through cases
for icase,(temp_alpha_s,temp_rho_back) in enumerate(zip(alpha_s,rho_back)):
    inv_work_dir = inv_main_work_dir.format(icase+1)
    
    if not os.path.isdir(inv_work_dir):
        os.makedirs(inv_work_dir)

    # Save Inverse directories for DOI below
    out_inv_dirs.append(inv_work_dir)

    main_res = temp_rho_back

    # Only invert data, use Forward model from Surface 1 (requires Surface 1 to have been run)
    inv_r2_dict = {'survey_name':'Surface_6_dpdp_case{0:02.0f}'.format(icase+1),
                   'survey_dir': inv_work_dir, 'job_type':1,
                   'exe_dir':exe_dir,
                   'exe_name':'R2.exe'}
    
    iR2_inv = r2_tools.R2(**inv_r2_dict)
    
    # Assign inversion options when different from r2_tools.inv_defaults
    inv_options = {'a_wgt':1e-3,'b_wgt':2e-2,
                   'inverse_type':1,'patch_size_xy':[4,4],
                   'rho_min':-1e3,'rho_max':1e3,
                   'reg_mode':1,'alpha_s':temp_alpha_s} # assign reg_mode=1 and alpha_s
    
    mesh_type = 4 # 4= simple quadrilateral element mesh
    output_xy = [[electrode_x[0],0.],[electrode_x[-1],0.],
         [electrode_x[-1],-8.],[electrode_x[0],-8.],
         [electrode_x[0],0.]] # closed polygon outlining main domain
                 
    inv_r2_in_dict = {'electrode_array':electrode_array,
                      'output_domain':output_xy,'mesh_type':mesh_type,
                      'reg_elems':[[1,int((xx_mesh.shape[0]-1)*(yy_mesh.shape[0]-1)),main_res]],
                      'res_matrix':1,'inv_dict':inv_options,
                      'singular_type':0,'node_dict':node_dict}
    
    # Load forward model data
    meas_fwd_data = pyres_utils.load_fwd_output(work_dir=work_dir)
    
    inv_protocol_dict = {'meas_data':meas_fwd_data}
    run_inv_r2_dict = {'protocol_dict':inv_protocol_dict, 'r2_in_dict':inv_r2_in_dict,
                       'run_bool':True}
    
    iR2_inv.run_all(**run_inv_r2_dict)

    # Plot inverse model results
    plot_dict = {'work_dir':inv_work_dir,
                 'plt_opts':{'vmin':10**1.3,'vmax':10.**2,'aspect':2},
                 'invert_y':False,'nxny':[5e2,1e2],
                 'cmap':'rainbow','keep_log':True}
    fig,ax,_ = plot_utils.plot_res(**plot_dict)

#%% DOI analysis

doi_dirs = out_inv_dirs[::2]
doi_ref_res = rho_back[::2]
doi_in_dict = {'model_dirs':doi_dirs,'ref_res_array':doi_ref_res}
xy_array,doi_array = pyres_utils.doi_analysis(**doi_in_dict)


# Arrange data into gridded cells for plotting
grid_dict = {'inv_data':np.column_stack([xy_array,doi_array]),
             'inv_col':2,'nxny':[200,200]}
X,Y,DOI = plot_utils.grid_inv_data(**grid_dict)

# Plot contours
plt = plot_utils.plt
sc1=plt.contourf(X,Y,DOI,np.linspace(0,1,10),extend='both')
cbar=plt.colorbar(sc1,extend='both')
cbar.ax.set_ylabel('DOI [-]')


