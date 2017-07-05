# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 15:38:23 2017

@author: kbefus
"""

from __future__ import print_function
import os
from shutil import copyfile
import numpy as np

from pyres import mesh_tools
from pyres import r2_tools
from pyres import pyres_utils, plot_utils

#%% ########## Load ER data ################
# Define directory where example scripts are located
ex_dir = r'C:\ER' 
work_dir = os.path.join(ex_dir,r'R2_examples\ExFieldData\SensitivityTest')
stg_file = 'PER3A.stg'
stg_fname = os.path.join(os.path.split(work_dir)[0],stg_file)

if not os.path.isdir(work_dir):
    os.makedirs(work_dir)

# Define directory with R2
exe_dir = os.path.join(ex_dir,r'R2\Distribition_3.1\Progs')

load_dict = {'fname':stg_fname,
             'data_fmt':'AGI'}
             
meas_data, electrode_data = pyres_utils.load_er_data(**load_dict)  


electrode_data[:,3] = electrode_data[:,4] # force y to be elevation
nlines = np.unique(electrode_data[:,0]).shape[0] # assumes each line is a new string
nelectrodes_per_line = electrode_data.shape[0]/nlines
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

topo_fname = os.path.join(os.path.split(work_dir)[0],'PER3_topo.trn')
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

# Define inversion output domain as foreground area
foreground_nodes = np.hstack([mr2.boundaries['foreground'][ikey][0] for ikey in mr2.boundaries['foreground']['order']])
foreground_nodes = np.hstack([foreground_nodes,foreground_nodes[0]])
foreground_xy = np.array([mr2.mesh_dict['nodes'][inode-1][:2] for inode in foreground_nodes])
mesh_type = 3

inv_options = {'error_mod':2,
               'inverse_type':1}
r2_in_dict = {'electrode_array':mr2.electrode_array,
              'output_domain':None,'mesh_type':mesh_type,
              'reg_elems':[[1,len(mr2.mesh_dict['elements']),1e2]],
              'inv_dict':inv_options,
              'singular_type':0}

meas_temp = meas_data.copy() # keep meas_data in memory, and unchanged                           
cols2d = [0,2,4,6,8,9]# ,10] include last column for error
meas_temp = meas_temp[:,cols2d]
meas_temp[:,5] = -meas_temp[:,5] # correct polarity, or use [4,2,8,6] for p+, p-

# Test effects of 'a_wgt' and 'b_wgt' on resulting tomogram
a_wgts = np.hstack([0,np.logspace(-4,1,11)]) # ohms
b_wgts = [0,
          1e-2,2e-2,3e-2,4e-2,5e-2,6e-2,8e-2,
          1e-1,2e-1,5e-1] # x 100 to give percent error in resistance data

output_data = []
all_er_data = []
icatch_warnings=[]
convg_matrix = []
iter_matrix = []
icount = 0

force_run = False # for restarting anlaysis

for a_wgt in a_wgts:
    b_output_data = []
    inner_convg_matrix = []
    inner_iter_matrix = []
    for b_wgt in b_wgts:
        icount+=1
        print("-------------------------------------------------------")
        print("a_wgt = {0:1.0e}; b_wgt = {1:1.0e}".format(a_wgt,b_wgt))    
        inv_work_dir = os.path.join(work_dir,'inv_a{0:1.0e}_b{1:1.0e}'.format(a_wgt,b_wgt))
        if not os.path.isdir(inv_work_dir):
            os.makedirs(inv_work_dir)
            
        if not os.path.isfile(os.path.join(inv_work_dir,'f001_res.dat')) or force_run:
            if a_wgt==0 and b_wgt==0:
                # Use errors provided by measurements
                meas_temp2 = meas_data.copy()[:,[0,2,4,6,8,9,10]]
                meas_temp2[:,5] = -meas_temp2[:,5] # correct polarity, or use [4,2,8,6] for p+, p-
                meas_temp2[:,6] = meas_temp2[:,6]*meas_temp2[:,5] # estimate of stddev(R)
                protocol_dict = {'meas_data':meas_temp2}
            else:
                protocol_dict = {'meas_data':meas_temp}
                
            # Copy mesh.dat into inv_work_dir
            new_mesh_fname = os.path.join(inv_work_dir,'mesh.dat')
            if not os.path.isfile(new_mesh_fname):
                copyfile(os.path.join(os.path.split(mr2.msh_name)[0],'mesh.dat'),new_mesh_fname)
            
            r2_dict.update({'survey_dir':inv_work_dir})
            iR2 = r2_tools.R2(**r2_dict)
    
            inv_options.update({'a_wgt':a_wgt,'b_wgt':b_wgt})
            r2_in_dict.update({'inv_dict':inv_options})
    
            run_r2_dict = {'protocol_dict':protocol_dict, 'r2_in_dict':r2_in_dict,
                            'run_bool':True}
            
            iR2.run_all(**run_r2_dict)
            
            output = iR2.run_output

        else:
            output = pyres_utils.load_r2out(inv_work_dir)
        
        # Extract data from output
        catch_warnings=[itemp.strip() for itemp in output if 'WARNING' in itemp or 'FATAL' in itemp]
        if len(catch_warnings)>0:
            iwarn = catch_warnings[-1]
        else:
            iwarn = 'None'
            
        icatch_warnings.append(iwarn)
        
        if 'FATAL' in iwarn:
            rms1=np.nan
            n_iter=0
            convg_status=np.nan
        else:
            rms1 = [float(i.split()[-1]) for i in output if 'Final RMS' in i][-1]
            n_iter = [int(i.split()[-1]) for i in output if 'Iteration' in i][-1]
            convg_status = 1
            if iwarn is not None:
                if 'not converged' in iwarn:
                    convg_status = -1
        
        inner_convg_matrix.append(convg_status)
        inner_iter_matrix.append(n_iter)
        if not np.isnan(rms1):
            # Post process runs
            inv_data = pyres_utils.load_inv_output(work_dir=inv_work_dir)
        else:
            inv_data = None
            
        b_output_data.append([icount,a_wgt,b_wgt,convg_status,rms1,n_iter])
        all_er_data.append(inv_data)
        print("-------------------------------------------------------")
    output_data.append(b_output_data)
    convg_matrix.append(inner_convg_matrix)    
    iter_matrix.append(inner_iter_matrix)
#%% Analyze results
# For pdf text saving
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Arial']})
import matplotlib.pyplot as plt
analyze_bool = True
if analyze_bool:
    # Find which models converged
    finished_models = [i for i,data in enumerate(all_er_data) if data is not None]
    
    # Select those models
    all_models = np.array(all_er_data).reshape((len(a_wgts),len(b_wgts)))
    all_models2 = all_models[::-1,:] # organize to match convergence plot
    
    # Find elements inside foreground
    elements = np.array(mr2.mesh_dict['elements'])
    elem_areas,start_inds = np.unique(elements[:,-1], return_index=True)
    end_inds = np.hstack([start_inds[1:],elements.shape[0]])
    foreground_elems = np.arange(start_inds[0],end_inds[0])
    
    errfunc = lambda obs,mod: (obs-mod)/obs
    main_inds = [6,3]
    relative_bool=True
    if relative_bool:
        inv_col=2
        max_val = 2
        main_er_data = all_models2[main_inds[0],main_inds[1]][:,2:].copy()
    else:
        inv_col=3
    main_inv_data = all_models.ravel()[finished_models[0]][foreground_elems,:]
    new_mesh_dict = {'elements':elements[foreground_elems,:],
                     'nodes':mr2.mesh_dict['nodes']}
    plt_dict = {'mesh_dict':new_mesh_dict,'inv_col':inv_col,
                }


    fig,ax = plt.subplots(len(a_wgts),len(b_wgts),sharex=True,sharey=True)

    ax_all = ax.ravel()
    # Plot all data
    icount = 0
    for ia,a_wgt in enumerate(a_wgts[::-1]): # organize to match convergence plot
        for ib,b_wgt in enumerate(b_wgts):
            if relative_bool:
                if all_models2[ia,ib] is not None:
                    temp_data = all_models2[ia,ib].copy()
                    temp_data[:,2:] = errfunc(main_er_data,temp_data[:,2:])
                else:
                    temp_data = None
                plt_dict.update({'plt_dict':{'vmin':-max_val,'vmax':max_val,'cmap':'RdBu'}})
            else:
                temp_data = all_models2[ia,ib]
                plt_dict.update({'plt_dict':{'vmin':5,'vmax':200,'cmap':'rainbow'}})
            if temp_data is not None:          
                if int(temp_data.shape[0])==len(mr2.mesh_dict['elements']):
                    temp_data=temp_data[foreground_elems,:].copy()
                elif int(temp_data.shape[0])>len(new_mesh_dict['elements']):
                    temp_data=temp_data[:len(new_mesh_dict['elements']),:].copy()
                    
                plt_dict.update({'inv_data':temp_data,
                                 'figax':[fig,ax_all[icount]]})    
                fig,ax_temp,tri=plot_utils.plot_tri_mesh(**plt_dict)
                ax_temp.set_ylim([main_inv_data[:,1].min(),main_inv_data[:,1].max()])
                ax_temp.set_xlim([main_inv_data[:,0].min(),main_inv_data[:,0].max()])
            ax_temp.axes.get_yaxis().set_visible(False)
            ax_temp.axes.get_xaxis().set_visible(False)
            ax_temp.axis('off')
            #ax_all[icount].set_title('a={0:1.1e},b={1:1.2f}'.format(a_wgt,b_wgt))
            icount+=1
    
    # Make one colorbar for all plots 
#    ticks = [1.5,5,10,50,100,200.]
    fig2=plt.figure()
    cbar_ax = fig2.add_axes([0.5,0.15,0.05,.7])
    fig.colorbar(tri,cax=cbar_ax)#,ticks=ticks)
    plot_utils.plt.show()    

save_fig=False
if save_fig:
    photo_dir=r'D:\Research\ER\Paper\Figures'
    fig.savefig(os.path.join(photo_dir,'FieldEx_sens_rel_all.png'),dpi=500,
                format='png',)
    
    #%% Plot convergence matrix
    plt.close('all')
    iters = np.array(iter_matrix)
    fig,ax2 = plt.subplots()
    cmap = plot_utils.plt.cm.get_cmap('RdBu',3)
    plt_opts = {'cmap':cmap,
                'vmin':-1,'vmax':1,'edgecolors':'grey','shading':'flat'}
    im = ax2.pcolormesh(np.array(convg_matrix),**plt_opts)
#    im = ax2.imshow(np.array(convg_matrix),**plt_opts)
    xtick_locs = ax2.get_xticks()
    new_xtick_locs = np.arange(xtick_locs[0],xtick_locs[-1],1)+0.5
    ax2.set_xticks(new_xtick_locs)    
    ax2.set_xticklabels(np.hstack([['{0:1.2f}'.format(b) for b in b_wgts]]))
    ytick_locs = ax2.get_yticks()
    new_ytick_locs = np.arange(ytick_locs[0],ytick_locs[-1],1)+0.5
    ax2.set_yticks(new_ytick_locs)     
    ax2.set_yticklabels(np.hstack([['{0:1.1e}'.format(a) for a in a_wgts]]))
    
    for ia,(a_wgt,iy) in enumerate(zip(a_wgts,new_ytick_locs)):
        for ib,(b_wgt,ix) in enumerate(zip(b_wgts,new_xtick_locs)):
            ax2.text(ix,iy,'{0:d}'.format(iters[ia,ib]),color='white')
            
    ax2.set_xlim([0,11])
#    ax2.set_ylim([-1,13])
    ax2.set_ylabel('a_wgt [ohm]')
    ax2.set_xlabel('b_wgt [-]')
    plt.colorbar(im,ax=ax2,fraction=0.046, pad=0.04,ticks=[-1,0,1])
    
#%% ################## Plot inverse model results ##################
#plot_dict = {'work_dir':work_dir,
#             'method':'linear',
#             'topog_xy':mr2.topo_xyz,
#             'plt_opts':{'vmin':1.5,'vmax':200.,'ticks':[1.5,5,10,50,100,200.]},
#             'invert_y':False,
#             'cmap':'rainbow','keep_log':False}
#fig,ax,[X,Y,ER] = plot_utils.plot_res(**plot_dict)
#c1=ax.contour(X,Y,ER,[5,10,50,100], colors='k')
#plot_utils.plt.clabel(c1, inline=0, fontsize=10,fmt='%3.1f')
##ax.set_ylim([-6,2])
##ax.set_xlim([0,68.75])
#
##%% ################## Plot model sensitivity results ##################
#plot_dict = {'fname':os.path.join(work_dir,'f001_sen.dat'),
#             'method':'linear',
#             'topog_xy':mr2.topo_xyz,
#             'plt_opts':{'vmin':1e-6,'vmax':1e1},
#             'invert_y':False,
#             'cmap':'rainbow','keep_log':True}
#fig2,ax2,[X,Y,sens] = plot_utils.plot_res(**plot_dict)
#
#c1=ax2.contour(X,Y,sens,[-5,-4,-3,-2,-1], colors='k')
#plot_utils.plt.clabel(c1, inline=0, fontsize=10,fmt='%1.0f')
##ax.set_ylim([-6,2])
#ax2.set_xlim([0,68.75])
#
##%% ################## Plot ER with sensitivity transparency ##################
#
## Normalize sensitivity matrix to 0-1 alpha (transparency) values
#min_sens = -5
#full_sens = -2
#sens_norm = 1-(full_sens-sens)/(full_sens-min_sens)
#sens_norm[sens_norm<0] = 0
#sens_norm[sens_norm>1] = 1
#sens_norm = np.ma.filled(sens_norm,0.)
#
#ax.cla()
#c2 = ax.pcolormesh(X,Y,np.ma.filled(ER,np.nan),vmin=1.5,vmax=200.,cmap='rainbow')
#fig.canvas.draw()
#
#new_alphas = []
#new_facecolors = []
#for iface,ialpha in zip(c2.get_facecolors(),sens_norm[1:,1:].flatten()):
#    iface[3]=ialpha # set alpha value of RGBA tuple
#fig.canvas.draw()
#plot_utils.plt.show()





