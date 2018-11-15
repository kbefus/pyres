# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 15:38:23 2017

@author: kbefus
"""

from __future__ import print_function
import os
import numpy as np

from pyres import mesh_tools
from pyres import r2_tools
from pyres import pyres_utils, plot_utils

#%% ########## Load ER data ################
# Define directory where example scripts are located
ex_dir = r'C:\ER' 
work_dir = os.path.join(ex_dir,r'R2_examples\ExFieldData')
stg_file = 'PER3A.stg'
stg_fname = os.path.join(work_dir,stg_file)

# Define directory with R2
exe_dir = os.path.join(ex_dir,r'R2\Distribition_3.1\Progs')

load_dict = {'fname':stg_fname,
             'data_fmt':'AGI'}
             
meas_data, electrode_data = pyres_utils.load_er_data(**load_dict)  


electrode_data[:,3] = electrode_data[:,4] # force y to be elevation
nlines = np.unique(electrode_data[:,0]).shape[0] # assumes each line is a new string
nelectrodes_per_line = int(electrode_data.shape[0]/nlines)
line_strexyz = electrode_data.reshape((nlines,nelectrodes_per_line,-1))
#%% ############## Create triangular mesh ######################
# Initialize meshR2 object with basic line information
mr2 = mesh_tools.meshR2(nlines=1,name=stg_file,nelectrodes_per_line=nelectrodes_per_line)

# Setup near and far field domains
dx = np.diff(electrode_data[:,2])[0] # x-spacing of electrodes
clen = dx/6. # Foreground mesh characteristic length scale
bg_size = 50. # One-sided size of the background domain, m
foreground_dxdy = [clen*20,clen*50] # buffer around foreground domain
background_dxdy = np.array([bg_size,bg_size]) + np.array(foreground_dxdy) # buffer of background domain
background_clen = clen * bg_size*10. # Background mesh characteristic length scale
boundary_dict = {'fdxdy':foreground_dxdy,
                 'bdxdy':background_dxdy,
                 'bclen':background_clen}
              
write_dict = {'out_fname':os.path.join(work_dir,"{}.geo".format(os.path.splitext(stg_file)[0])),
              'boundary_dict':boundary_dict}
geom_inputs = {'line_strexyz':line_strexyz,'clen':clen}
gmsh_dict = {'gmshdir':r'D:\Research\Software\gmsh-2.16.0-Windows'}
make_dict = {'geom_dict':geom_inputs,'write_dict':write_dict,
             'run_gmsh':True,'gmsh_dict':gmsh_dict}
mr2.make_mesh(**make_dict)

topo_fname = os.path.join(work_dir,'PER3_topo.trn')
topo_dict = {'topo_file':topo_fname,'topo_load_dict':{'cols2read':[0,1],
                                                      'nheaders':3,
                                                      'delimiter':','},
             'write_gmsh_bool':True,'method':'constant'}

to_dat_dict = {'param_zones':{1:{'param':False,
                                 'zone':True}},
               'topo_correct':True,
               'topo_dict':topo_dict}
mr2.msh_to_dat(**to_dat_dict)
#%% ############ Prepare R2 files ##############
r2_dict = {'survey_name':os.path.basename(stg_fname).split('.')[0],
           'survey_dir': work_dir, 'job_type':1,
           'exe_dir':exe_dir,
           'exe_name':'R2.exe'}
iR2 = r2_tools.R2(**r2_dict)

# Define inversion output domain as foreground area
foreground_nodes = np.hstack([mr2.boundaries['foreground'][ikey][0] for ikey in mr2.boundaries['foreground']['order']])
foreground_nodes = np.hstack([foreground_nodes,foreground_nodes[0]])
foreground_xy = np.array([mr2.mesh_dict['nodes'][inode-1][:2] for inode in foreground_nodes])


mesh_type = 3
inv_options = {'a_wgt':1e-2,'error_mod':2,
               'inverse_type':1}
r2_in_dict = {'electrode_array':mr2.electrode_array,
              'output_domain':foreground_xy,'mesh_type':mesh_type,
              'reg_elems':[[1,len(mr2.mesh_dict['elements']),1e3]],
              'inv_dict':inv_options,
              'singular_type':0}

meas_temp = meas_data.copy() # keep meas_data in memory, and unchanged                           
cols2d = [0,2,4,6,8,9]# ,10] include last column for error
meas_temp = meas_temp[:,cols2d]
meas_temp[:,5] = -meas_temp[:,5] # correct polarity, or use [4,2,8,6] for p+, p-

# Zero bad datapoints
plot_utils.plt.close('all')
# Use very generous data culling - won't remove any data
_,_,_,bad_inds,ax_raw = plot_utils.plot_raw(meas_data=meas_temp,npt_neighbor=7,
                                     plot_bool=True,median_thresh=.75)
meas_temp[bad_inds,-1] = 0.
protocol_dict = {'meas_data':meas_temp}
run_r2_dict = {'protocol_dict':protocol_dict, 'r2_in_dict':r2_in_dict,
                'run_bool':True}

iR2.run_all(**run_r2_dict)

#%% ################## Plot inverse model results ##################
plot_dict = {'work_dir':work_dir,
             'method':'linear',
             'topog_xy':mr2.topo_xyz,
             'plt_opts':{'vmin':1.5,'vmax':200.,'ticks':[1.5,5,10,50,100,200.]},
             'invert_y':False,
             'cmap':'rainbow','keep_log':False}
fig,ax,[X,Y,ER] = plot_utils.plot_res(**plot_dict)
c1=ax.contour(X,Y,ER,[5,10,50,100], colors='k')
plot_utils.plt.clabel(c1, inline=0, fontsize=10,fmt='%3.1f')
#ax.set_ylim([-6,2])
#ax.set_xlim([0,68.75])

#%% ################## Plot model sensitivity results ##################
plot_dict = {'fname':os.path.join(work_dir,'f001_sen.dat'),
             'method':'linear',
             'topog_xy':mr2.topo_xyz,
             'plt_opts':{'vmin':1e-6,'vmax':1e1},
             'invert_y':False,
             'cmap':'rainbow','keep_log':True}
fig2,ax2,[X,Y,sens] = plot_utils.plot_res(**plot_dict)

c1=ax2.contour(X,Y,sens,[-5,-4,-3,-2,-1], colors='k')
plot_utils.plt.clabel(c1, inline=0, fontsize=10,fmt='%1.0f')
#ax.set_ylim([-6,2])
ax2.set_xlim([0,68.75])

#%% ################## Plot ER with sensitivity transparency ##################

# Normalize sensitivity matrix to 0-1 alpha (transparency) values
min_sens = -5
full_sens = -2
sens_norm = 1-(full_sens-sens)/(full_sens-min_sens)
sens_norm[sens_norm<0] = 0
sens_norm[sens_norm>1] = 1
sens_norm = np.ma.filled(sens_norm,0.)

ax.cla()
c2 = ax.pcolormesh(X,Y,np.ma.filled(ER,np.nan),vmin=1.5,vmax=200.,cmap='rainbow')
fig.canvas.draw()

new_alphas = []
new_facecolors = []
for iface,ialpha in zip(c2.get_facecolors(),sens_norm[1:,1:].flatten()):
    iface[3]=ialpha # set alpha value of RGBA tuple
fig.canvas.draw()
plot_utils.plt.show()





