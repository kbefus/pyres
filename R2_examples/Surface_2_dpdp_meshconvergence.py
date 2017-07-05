# -*- coding: utf-8 -*-
"""
Created on Thu May 25 08:06:54 2017

Example of a pyres implementation of "Surface electorde array 2" from the 
R2_v3.1_readme.pdf

@author: kbefus
"""
from __future__ import print_function
import os,time
import numpy as np

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
import matplotlib.pyplot as plt

from pyres import mesh_tools, mesh_utils
from pyres import r2_tools
from pyres import pyres_utils, plot_utils

#%%

def make_mesh(elocs=None,ndivxy=None,
              dx_background = [0.5,1., 2., 5., 10., 20., 50., 100.],
              dy_background = [0.5,1., 2., 5., 10., 20., 50., 100.],
              xybuff = [2.,10.],ndec=4):
    
    xbuff,ybuff = xybuff
    e_spacing = np.min(np.abs(np.diff(elocs)))
    dx_foreground = e_spacing/ndivxy[0] # 0.25 m
    dy_foreground = e_spacing/ndivxy[1] # 0.25 m
    xx_mesh = np.arange(elocs[0]-xbuff,elocs[-1]+xbuff+dx_foreground,dx_foreground)
    xx_mesh = np.unique(np.round(np.hstack([(xx_mesh[0]-np.cumsum(dx_background))[::-1],
                                   xx_mesh,
                                   xx_mesh[-1]+np.cumsum(dx_background),
                                   elocs]),decimals=ndec))
    
    yy_mesh = np.arange(0,ybuff+dy_foreground,dy_foreground)
    yy_mesh = np.hstack([yy_mesh,yy_mesh[-1]+np.cumsum(dy_background)])
    
    return xx_mesh,yy_mesh

errfunc = lambda obs,mod: (obs-mod)/obs
topo_func = lambda m,x,b: m*x+b

def make_divisble(array=None,ndiv=None,step_size=None):
    
    isuccess=False
    while not isuccess:
        if np.fmod(array.shape[0]-1,ndiv)!=0:
            array = np.hstack([array,array[-1]+step_size])
        else:
            isuccess=True
    return array
        

#%%
# Define directory where example scripts are located
ex_dir = r'C:\ER' 
work_dir_main = os.path.join(ex_dir,r'R2_examples\Surface_2_meshconvg\dpdp\Forward')

if not os.path.isdir(work_dir_main):
    os.makedirs(work_dir_main)

# Define directory with R2
exe_dir = os.path.join(ex_dir,r'R2\Distribition_3.1\Progs')

# -------- Create mesh -------------
n_electrodes = 25.
e_spacing = 2. # meters is unit of distance for this example
ndivx,ndivy = 8.,8. # mesh divisions between electrodes
xbuff,ybuff = 2.,10. # meters to extend foreground beyond electrode x,y locations
electrode_x = e_spacing*np.arange(n_electrodes) # first electrode at x=0, last electrode at x=48
mesh_divs = np.arange(1,16.,1)
mesh_divxy = np.column_stack([mesh_divs,mesh_divs])

# Target properties
main_res = 100 # ohm m
target_res = 10 # ohm m
target_xypos = [[14,1],[16,1],
                [16,4],[14,4]] # do not want a closed polygon
# Convert target location to array
txy = np.array(target_xypos)

all_outputs=[]
# Loop through forward mesh division sizes
for fwddivxy in mesh_divxy[1::2]:# only even size mesh divisions
    work_dir = os.path.join(work_dir_main,'mesh_fwd_div_x{0:.0f}y{1:.0f}'.format(*fwddivxy))
    if not os.path.isdir(work_dir):
        os.makedirs(work_dir)
        
    # Define mesh nodes
    ndivx,ndivy = fwddivxy
    make_mesh_dict = {'elocs':electrode_x,'ndivxy':[ndivx,ndivy],
                      'xybuff':[xbuff,ybuff]}
    xx_mesh,yy_mesh = make_mesh(**make_mesh_dict)
    
    # Define topography
    topog = topo_func(4.8/electrode_x[-1],xx_mesh,0.)
    
    # Define mesh creation input
    node_dict = {'xx':xx_mesh,'yy':yy_mesh,'topog':topog}
    
    # Define electrode rows and columns
    e_rows = np.ones_like(electrode_x) # all in top row
    e_cols = np.array([(xx_mesh==ix).nonzero()[0][0] for ix in electrode_x])+1 # find columns of electrode nodes
    electrode_array = np.column_stack([np.arange(e_cols.shape[0])+1,e_cols,e_rows]).astype(int)
    # mesh_ax = plot_utils.plot_quad_mesh([xx_mesh,yy_mesh],label_nodes=False) # uncomment to plot mesh
    
    # ----------------- Define target -----------------

    region_dict = {'mesh_xy':[xx_mesh,yy_mesh],'target_xypos':txy,
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
    output_xy = None
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
    
    
    ################ Invert data ########################

    # Load forward model data
    meas_fwd_data = pyres_utils.load_fwd_output(work_dir=work_dir)
        
    error_outputs = []
    icatch_warnings = []
    npatch = 4
    model_elapsed=[]
    force_run=False
    # Loop through mesh division sizes
    for indivxy in mesh_divxy:
        model_start = time.time()
        inv_work_dir = os.path.join(work_dir,'inv_mesh_div_x{0:.0f}y{1:.0f}'.format(*indivxy))
        
        if not os.path.isdir(inv_work_dir):
            os.makedirs(inv_work_dir)
        
        # Make new inverse mesh
        # Define mesh creation input
        make_mesh_dict_inv = {}
        make_mesh_dict_inv.update(make_mesh_dict)
        make_mesh_dict_inv.update({'ndivxy':indivxy})
        xx_mesh,yy_mesh = make_mesh(**make_mesh_dict_inv)    
        
        # Force xx and yy_mesh to be divisible by npatch
        xx_mesh = make_divisble(xx_mesh,npatch,5)
        yy_mesh = make_divisble(yy_mesh,npatch,5)
        
        topog = topo_func(4.8/electrode_x[-1],xx_mesh,0.)
        
        if not os.path.isfile(os.path.join(inv_work_dir,'f001_res.dat')) and not force_run:
        
            inv_r2_dict = {}
            inv_r2_dict.update(r2_dict) # take inputs from forward model
            inv_r2_dict.update({'job_type':1,
                                'survey_dir': inv_work_dir}) # run inverse model and save to new directory
            
            iR2_inv = r2_tools.R2(**inv_r2_dict)
    
            # Assign inversion options when different from r2_tools.inv_defaults
            inv_options = {'a_wgt':1e-3,'b_wgt':2e-2,
                           'inverse_type':1,'patch_size_xy':[npatch,npatch],
                           'rho_min':-1e3,'rho_max':1e3}
            
            new_node_dict = {'xx':xx_mesh,'yy':yy_mesh,'topog':topog}
            
            e_cols = np.array([(xx_mesh==ix).nonzero()[0][0] for ix in electrode_x])+1 # find columns of electrode nodes
            electrode_array = np.column_stack([np.arange(e_cols.shape[0])+1,e_cols,e_rows]).astype(int)
            
            inv_r2_in_dict = {}
            inv_r2_in_dict.update(r2_in_dict)
            inv_r2_in_dict.update({'reg_elems':[[1,int((xx_mesh.shape[0]-1)*(yy_mesh.shape[0]-1)),main_res]],
                                    'res_matrix':1,
                                    'inv_dict':inv_options,
                                    'node_dict':new_node_dict,
                                    'electrode_array':electrode_array})
            
            inv_protocol_dict = {'meas_data':meas_fwd_data}
            run_inv_r2_dict = {'protocol_dict':inv_protocol_dict, 'r2_in_dict':inv_r2_in_dict,
                               'run_bool':True}
            
            iR2_inv.run_all(**run_inv_r2_dict)
            if 'FATAL error' in iR2_inv.run_output:
                break
            
            catch_warnings=[itemp.strip() for itemp in iR2_inv.run_output if 'WARNING' in itemp]
            if len(catch_warnings)>0:
                iwarn = catch_warnings[-1]
            else:
                iwarn = None
                
            icatch_warnings.append(iwarn)
            
            # Collect final inversion RMS
            rms1 = [float(i.split()[-1]) for i in iR2_inv.run_output if 'Final RMS' in i][-1]
            convg_status = 1
            if iwarn is not None:
                if 'not converged' in iwarn:
                    convg_status = -1

        else:
            output = pyres_utils.load_r2out(inv_work_dir)
            catch_warnings=[itemp.strip() for itemp in output if 'WARNING' in itemp]
            if len(catch_warnings)>0:
                iwarn = catch_warnings[-1]
            else:
                iwarn = None
                
            icatch_warnings.append(iwarn)
            
            rms1 = [float(i.split()[-1]) for i in output if 'Final RMS' in i][-1]
            convg_status = 1
            if iwarn is not None:
                if 'not converged' in iwarn:
                    convg_status = -1
            
        # Post process runs
        inv_data = pyres_utils.load_inv_output(work_dir=inv_work_dir)
        X,Y,ER = pyres_utils.grid_inv_data(inv_data,inv_col=2) # use non-log-transformed data
        
        obs_array = main_res*np.ones_like(ER)
        
        topo_y = pyres_utils.extrap(txy[:,0],xx_mesh,topog)
        target_bool = (X>=txy[0,0]) & (X<=txy[1,0]) &\
                      (Y<=(topo_y[0]-txy[0,1])) & (Y>=(topo_y[2]-txy[2,1]))
        obs_array[target_bool] = target_res
        
        err_array = errfunc(obs_array,ER)
        target_err = err_array[target_bool]
        # Integrate error over target weighted by areas
        # Use area weighting of errors, only for quad meshes
        x2,y2 = np.meshgrid(xx_mesh,yy_mesh)
        xdiff = np.hstack([np.diff(x2,axis=1),np.zeros((x2.shape[0],1))])
        ydiff = np.vstack([np.diff(y2,axis=0),np.zeros((1,y2.shape[1]))])
        area_array = xdiff[:-1,:-1]*ydiff[:-1,:-1]
        
        target_areas = area_array[target_bool]
        areal_rms = np.sqrt((np.sum(target_areas*target_err)/np.sum(target_areas))**2.)
        domain_rms = np.sqrt((np.sum(area_array*err_array)/np.sum(area_array))**2.)
        

        error_outputs.append(np.hstack([main_res,target_res,fwddivxy,indivxy,convg_status,
                                        rms1,areal_rms,domain_rms,np.min(target_err),np.median(target_err),
                                        np.max(target_err),target_err.size]))
        model_elapsed.append(time.time()-model_start)
    all_outputs.append(np.vstack(error_outputs))
    
total_time = np.sum(model_elapsed)

plt.close('all')
all_err_data = np.dstack(all_outputs) # 3d array
all_err_data[all_err_data==-99.]=np.nan
#%%
rowY,colX = np.meshgrid(mesh_divs[1::2],mesh_divs)
fig,ax = plt.subplots()
p1=ax.imshow(np.abs(all_err_data[:,9,:].squeeze()),
                 vmin=0,vmax=5,cmap='RdBu',interpolation='none',
                 extent=[mesh_divs[1]-1,mesh_divs[-2]+1,mesh_divs[-1]+.5,mesh_divs[0]-.5])
#ax.plot(rowY.T,all_err_data[:,9,:].squeeze().T)

plt.colorbar(p1)
    
#%% 3d plot
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(rowY,colX, np.log10(np.abs(all_err_data[:,9,:].squeeze())), cmap='coolwarm',
                       linewidth=0.1, antialiased=False,rstride=1,cstride=1)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()
#%% Fwd vs invs mesh effects
plt.close('all')
from matplotlib.ticker import MultipleLocator
minorLocator = MultipleLocator(1)
fig,ax = plt.subplots()
colors = plt.cm.Accent(np.linspace(0,1,all_err_data.shape[-1]))
# Fwd mesh
for ifwd in range(all_err_data.shape[-1]):
    ydata = all_err_data[:,8,ifwd].squeeze().copy()
#    ydata = (ydata-np.min(ydata))/np.mean(ydata)
#    ax[0].fill_between(mesh_divs,all_err_data[:,8,ifwd].squeeze(),
#                      y2=all_err_data[:,10,ifwd].squeeze(),
#                      facecolor=colors[ifwd],alpha=0.15,edgecolor='k')
    ax.plot(mesh_divs,ydata,'o-',color=colors[ifwd],label='target, fwd {0:1.0f}'.format(mesh_divs[1::2][ifwd]),lw=2)
    ax.plot(mesh_divs,all_err_data[:,9,ifwd].squeeze(),'s-',lw=2,color=colors[ifwd],label='domain, fwd {0:1.0f}'.format(mesh_divs[1::2][ifwd]))
ax.set_yscale('log')
ax.set_ylabel('Integrated misfit')
ax.set_xlabel('Inverse model mesh divisions')
ax.xaxis.set_minor_locator(minorLocator)
ax.set_xlim([1,15])
plt.legend()

nnodes = all_err_data[:,-1,0].squeeze()
ax2 = ax.twinx()
ax2.plot(mesh_divs,nnodes,'ko-',lw=2)
ax2.set_xlim([1,15])
ax2.set_ylabel('Number of elements')
               
