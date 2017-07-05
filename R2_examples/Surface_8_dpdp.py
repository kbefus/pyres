# -*- coding: utf-8 -*-
"""
Created on Thu May 25 08:06:54 2017

Example of a pyres implementation of "Surface electorde array 8" from the 
R2_v3.1_readme.pdf

@author: kbefus
"""
from __future__ import print_function
import os
from shutil import copyfile
import numpy as np

from pyres import mesh_tools, mesh_utils
from pyres import r2_tools
from pyres import pyres_utils, plot_utils

#%%
# Define directory where example scripts are located
ex_dir = r'C:\ER' 
work_dir = os.path.join(ex_dir,r'R2_examples\Surface_8\dpdp\Forward')
inv_work_dir = os.path.join(ex_dir,r'R2_examples\Surface_8\dpdp\Inverse')
if not os.path.isdir(work_dir):
    os.makedirs(work_dir)

# Define directory with R2
exe_dir = os.path.join(ex_dir,r'R2\Distribition_3.1\Progs')

# Define gmsh dir
gmsh_dir = os.path.join(ex_dir,r'gmsh-2.16.0-Windows')

# -------- Create mesh -------------
n_electrodes = 25.
e_spacing = 2. # meters is unit of distance for this example
ndivx,ndivy = 4.,4. # mesh divisions between electrodes
xbuff,ybuff = 2.,10. # meters to extend foreground beyond electrode x,y locations


clen = e_spacing/ndivx
bg_size_scaler = 50.
foreground_dxdy = [xbuff,ybuff] # buffer around foreground domain
background_dxdy = np.array([bg_size_scaler,bg_size_scaler]) + np.array(foreground_dxdy) # buffer of background domain
background_clen = clen * bg_size_scaler
boundary_dict = {'fdxdy':foreground_dxdy,
                 'bdxdy':background_dxdy,
                 'bclen':background_clen}

mr2_dict = {'name':'Surface_8',
            'nelectrodes_per_line':n_electrodes,
            'electrode_spacing':e_spacing}
mr2 = mesh_tools.meshR2(**mr2_dict)

write_dict = {'out_fname':os.path.join(work_dir,"{}.geo".format(mr2_dict['name'])),
              'boundary_dict':boundary_dict}
geom_inputs = {'clen':clen}
gmsh_dict = {'gmshdir':gmsh_dir}
# ----------------- Define target -----------------
main_res = 100 # ohm m
target_res = 10 # ohm m
target_xypos = [[[14,-1,0],[16,-1,0],
                [16,-4,0],[14,-4,0],
                ]] # don't want a closed polygon
region_dict = {'region_xyzpts':target_xypos,'boundary_dict':boundary_dict,
               'clen':clen,'outside_foreground':False}

# Make the triangular mesh with target region
make_dict = {'geom_dict':geom_inputs,'write_dict':write_dict,
             'run_gmsh':True,'gmsh_dict':gmsh_dict,
             'region_dict':region_dict}
mr2.make_mesh(**make_dict)
#%%
# Convert the .msh to R2 .dat format
to_dat_dict={'param_zones':{1:{'param':False,'zone':False},
                            2:{'param':False,'zone':False},
                            3:{'param':False,'zone':False}}} # change parame and zone values to true to constrain regions for inversion
mr2.msh_to_dat(**to_dat_dict)
#%%
# ----------------- Setup R2 inputs ----------------
job_type = 0 # 0=forward, 1=inverse model
r2_dict = {'survey_name':'Surface_1_dpdp',
           'survey_dir': work_dir, 'job_type':job_type,
           'exe_dir':exe_dir,
           'exe_name':'R2.exe'}
iR2 = r2_tools.R2(**r2_dict)

mesh_type = 3 # 3= triangular element mesh
electrode_x = mr2.line_strexyz.squeeze()[:,2]
output_xy = [[electrode_x[0],0.],[electrode_x[-1],0.],
             [electrode_x[-1],-8.],[electrode_x[0],-8.],
             [electrode_x[0],0.]] # closed polygon outlining main domain
             
output_xy = None # for plotting with triangular mesh

# Define two regions using mesh creation information
target_region = np.unique(mr2.mesh_dict['zones'])[0] # target regions made first if in foreground
reg_dict = {'mesh_dict':mr2.mesh_dict,'target_zones':target_region,
            'main_res':main_res,'target_res':target_res}
reg_elems = mesh_utils.make_tri_region_elements(**reg_dict)
              
r2_in_dict = {'electrode_array':mr2.electrode_array,
              'output_domain':output_xy,'mesh_type':mesh_type,
              'reg_elems':reg_elems,
              'res_matrix':0,
              'singular_type':0}

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

# Copy mesh.dat to Inverse folder
fwd_meshdat = os.path.join(work_dir,'mesh.dat')
inv_meshdat = os.path.join(inv_work_dir,'mesh.dat')
copyfile(fwd_meshdat,inv_meshdat)

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
inv_r2_in_dict.update({'reg_elems':[[1,len(mr2.mesh_dict['elements']),main_res]],
                        'res_matrix':1,
                        'inv_dict':inv_options})

# Load forward model data
meas_fwd_data = pyres_utils.load_fwd_output(work_dir=work_dir)

inv_protocol_dict = {'meas_data':meas_fwd_data}
run_inv_r2_dict = {'protocol_dict':inv_protocol_dict, 'r2_in_dict':inv_r2_in_dict,
                   'run_bool':True}

iR2_inv.run_all(**run_inv_r2_dict)

#%% ######### Plot gmsh mesh ################

#m_fig,m_ax,m_tri_obj = plot_utils.plot_tri_mesh(mr2.mesh_dict)

#%%
# Plot inverse model results
#plot_dict = {'work_dir':inv_work_dir,'method':'linear',
#             'topog_xy':np.column_stack([electrode_x,np.zeros_like(electrode_x)]),
#             'plt_opts':{'vmin':10**1.3,'vmax':10.**2},
#             'invert_y':False,
#             'cmap':'rainbow','keep_log':True}
#fig,ax,_ = plot_utils.plot_res(**plot_dict)
#ax.plot(electrode_x,np.zeros_like(electrode_x),'ko-')
#
## Plot target
#txy = np.array(target_xypos).squeeze()
#txy = np.row_stack([txy,txy[0,:]])
#ax.plot(txy[:,0],txy[:,1],'k-')

#%%

plot_dict = {'work_dir':inv_work_dir,
             'mesh_dict':mr2.mesh_dict,
             'xylims':[electrode_x[0],electrode_x[-1],-10,0],
             'plt_opts':{'vmin':10**1.3,'vmax':10.**2,'cmap':'rainbow'},
             'keep_log':True}
fig,ax = plot_utils.plot_tri_res(**plot_dict)