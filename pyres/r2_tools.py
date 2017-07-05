# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 15:24:16 2017

@author: kbefus
"""
from __future__ import print_function
import numpy as np
import time, os
from pyres import pyres_utils

# R2 Inversion defaults
inv_defaults = {'inverse_type':1,
                'target_decrease':0.0,
                'qual_ratio':0,
                'patch_size_xy':[1,1],
                'rho_min':-1e6,'rho_max':1e6,
                'data_type':1, 'reg_mode':0,
                'tolerance':1.,'no_improve':1.0,
                'max_iterations':10,'error_mod':2,'alpha_aniso':1.0,
                'alpha_s':1e-2,'cginv_tol':1e-4, 'cginv_maxits':5e2,
                'a_wgt': 1e-3,'b_wgt':0.03,'num_xy_poly':0.,
                'mesh_scale':1.0}

x_background_growth = [0.5, 1.0, 2., 5., 10., 20., 50., 100.]
y_background_growth = [0.5, 1.0, 2., 5., 10., 20., 50., 100.]
                
class R2(object):
    '''Main R2 object.'''
    def __init__(self,exe_name='R2.exe',exe_dir='',survey_name=None,
                 survey_dir=None, job_type=None, data_fpath=None,
                 mesh_type=None, flux_type=3.0,):
        '''Create R2 object.
        
        Inputs
        --------------------------
        
        job_type: int, np.int
            Job type for R2, where 
                job_type = 0 is a forward model only 
                job_type = 1 includes an inverse model
        '''
        self.exe_fpath = os.path.join(exe_dir,exe_name)
        if survey_name is not None:
            self.survey_name = survey_name
        else:
            self.survey_name = "Unnamed 3D resistivity dataset"
        
        self.survey_dir = survey_dir
        self.job_type = job_type
        
        self.data_fpath = data_fpath   
        
    def r2_in(self,in_fname='R2.in', startingRfile=None, add_dir_bool=False,
              singular_type=1,mesh_type=None,
              electrode_array=None,
              flux_type=3.0, res_matrix=1, fwd_resis=1e2,
              inv_dict=inv_defaults, node_dict=None,
              output_domain=None,reg_elems=None):
        '''Make R2.in file.
        
        Inputs
        --------------------------
        
        res_matrix: int, np.int
            Matrix of resolution/sensistivity to output from R2.
            res_matrix = 0 gives no sensisitvity outputs
            res_matrix = 1 gives diagonal of [JTWTWJ]
            res_matrix = 2 gives true resolution matrix
            res_matrix = 3 gives sensitivity map and outputs Jacobian matrix.
        
        mesh_type: int, np.int
            Define mesh type
            mesh_type=
            mesh_type = 3: Triangular mesh
            mesh_type = 4: Regular quadrilateral mesh
            mesh_type = 5: Generalised quadrilateral mesh
        '''
                # Add directory if needed
        if os.path.dirname(in_fname) is '':
                R2in_output_file = os.path.join(self.survey_dir,in_fname)
                
        # Add directory if needed
#        add_dir_bool=False
#        if startingRfile is not None:
#            if os.path.dirname(startingRfile) is '':
#                add_dir_bool = True
            
        if add_dir_bool:
                startingRfile = os.path.join(self.survey_dir,startingRfile)
                
        # Insert any new parameters into inversion options
        temp_dict = inv_defaults # do not want to overwrite default values
        temp_dict.update(inv_dict)
        self.inv_dict = temp_dict
        
        # Assign volume to output resistivity data
        if output_domain is not None:
            if isinstance(output_domain,list):
                output_domain = np.array(output_domain)
                
            self.inv_dict['num_xy_poly'] = output_domain.shape[0]
            self.inv_dict['xy_poly'] = output_domain.copy()
        
        if reg_elems is None and startingRfile is None:
            num_regions_flag = 1
        elif reg_elems is None and startingRfile is not None:
            num_regions_flag = 0
        else:
            num_regions_flag = len(reg_elems)
            
        self.r2_options = {'job_type':self.job_type, 'mesh_type':mesh_type,
                      'flux_type': flux_type, 'singular_type':singular_type,
                      'res_matrix':res_matrix}
        self.r2_write_dict = {'in_fname':R2in_output_file,'survey_name':self.survey_name,
                      'r2_options':self.r2_options,
                      'num_regions_flag':num_regions_flag,'electrode_array':electrode_array,
                      'startingRfile':startingRfile,'fwd_resis':fwd_resis,'inv_dict':self.inv_dict,
                      'reg_elems':reg_elems,'node_dict':node_dict}
        
        write_r2_in(**self.r2_write_dict)
        
    def r2_protocol(self,protocol_fname='protocol.dat', meas_data=None):
        '''Write protocol file for R2.'''
                # Add directory if needed
        if os.path.dirname(protocol_fname) is '':
                protocol_fname = os.path.join(self.survey_dir,protocol_fname)
        
        if self.r2_write_dict['startingRfile'] is not None and self.r2_write_dict['num_regions_flag']==0:
            diff_inv = True
        else:
            diff_inv = False
        
        protocol_dict = {'protocol_fname':protocol_fname,
                         'job_type':self.job_type,
                         'inv_dict':self.inv_dict,
                         'meas_data':meas_data,
                         'diff_inv':diff_inv}
        write_r2_protocol(**protocol_dict)

    def run_r2(self):
        '''Run R2.exe.'''
        # Need to change working directory to project directory - shouldn't be necessary
        cwd = self.survey_dir
        os.chdir(cwd)
        
        # Construct command
        cmd_list = [self.exe_fpath]
        
        # Then run R3t.exe in that folder
        run_dict = {'cmd_list':cmd_list,'cwd':cwd,
                    'normal_msg':'','report':True}
        self.r2_success,self.run_output = pyres_utils.run_cmd(**run_dict)
        
        pass    
    
    def run_all(self,load_bool=False,load_dict=None,
                protocol_dict=None, r2_in_dict=None,
                run_bool=False):
        '''Run all steps in R2.'''
        if load_bool:
            self.load_data(**load_dict)
        
        self.r2_in(**r2_in_dict)
        self.r2_protocol(**protocol_dict)
        if run_bool:
            self.run_r2()
# --------------------- Helper functions ---------------------

def write_r2_in(in_fname='R2.in', survey_name=None, r2_options=None, num_regions_flag=1,
              electrode_array=None, startingRfile=None, fwd_resis=None, inv_dict=inv_defaults,
              node_dict=None, ncols=10,reg_elems=None):
    '''Write R2.in.'''
    with open(in_fname,'w') as f:
            
        # Header information
        header = "Project: {}, created using pyres on {}\n".format(\
                    survey_name,time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()))
        
        f.write(header)
        f.write("    {}    {}   {}    {}    {}  << job_type, mesh_type, flux_type, singularity_type, res_matrix \n".format(\
                                                    r2_options['job_type'], 
                                                    r2_options['mesh_type'], 
                                                    r2_options['flux_type'],
                                                    r2_options['singular_type'],
                                                    r2_options['res_matrix']))
        f.write("\n") # blank space
        
        # Mesh-dependent inputs
        if int(r2_options['mesh_type']) == 4:
            f.write("{0:10.0f}  {1:10.0f}  << numnp_x, numnp_y \n".format(node_dict['xx'].shape[0],
                                                            node_dict['yy'].shape[0]))
            f.write("\n")
            # Write x coordinates
            for ival,xval in enumerate(node_dict['xx']):
                f.write("{0:10.4f}  ".format(xval))
                if ival==len(node_dict['xx'])-1:
                    f.write(" << xx \n")
                elif np.fmod(ival+1,ncols)==0:
                    f.write("\n")
                    
            f.write("\n") # blank space 
            
            # Write topography
            
            if node_dict['topog'] is None:
                node_dict['topog'] = np.zeros_like(node_dict['xx'])
            
            for ival,tval in enumerate(node_dict['topog']):
                f.write("{0:10.4f}  ".format(tval))
                if ival==len(node_dict['topog'])-1:
                    f.write(" << topog \n")
                elif np.fmod(ival+1,ncols)==0:
                    f.write("\n")
            
            f.write("\n")
            
            # Write y (i.e., depth below topog) coordinates - positive down
            for ival,yval in enumerate(node_dict['yy']):
                f.write("{0:10.4f}  ".format(yval))
                if ival==len(node_dict['yy'])-1:
                    f.write(" << yy \n")
                elif np.fmod(ival+1,ncols)==0:
                    f.write("\n")

        elif int(r2_options['mesh_type']) == 5:
            f.write("{0:10.0f}  {1:10.0f}  << numnp_x, numnp_y \n".format(node_dict['xx'].shape[0],
                                                            node_dict['yy'].shape[0]))
            f.write("\n")
            # Write x coordinates
            for ival,xval in enumerate(node_dict['xx']):
                f.write("{0:10.4f}  ".format(xval))
                if ival==len(node_dict['xx'])-1:
                    f.write(" << xx \n")
                elif np.fmod(ival+1,ncols)==0:
                    f.write("\n")
                    
            f.write("\n") # blank space 
                       
            # Write y (i.e., depth below topog) coordinates - positive down
            for icol in xrange(node_dict['yy'].shape[1]):
                for ival,irow in enumerate(xrange(node_dict['yy'].shape[0])):
                    f.write("{0:10.4f}  ".format(node_dict['yy'][irow,icol]))
                    if ival==node_dict['yy'].shape[0]-1:
                        f.write(" << yy for xx column = {}\n".format(icol))
                    elif np.fmod(ival+1,ncols)==0:
                        f.write("\n")
        
        elif int(r2_options['mesh_type']) == 3:
            f.write("    {}  << mesh scale \n".format(inv_dict['mesh_scale']))
        
        f.write("\n") # blank space 
        # Line 11
        f.write("    {}  << num_regions\n".format(num_regions_flag))
        
        if num_regions_flag == 0 and startingRfile is not None:
            if len(startingRfile) > 15:
                print("Warning: filename must be 15 characters or less!")
                
            f.write("{}\n".format(startingRfile))
            f.write("\n") # blank space
        else:
            f.write("\n") # blank space
            if num_regions_flag==1 and reg_elems is None:
                reg_elems = [[electrode_array[0,1],electrode_array[-1,1],fwd_resis]]
                
            for iregion,reg_row in enumerate(reg_elems):
                f.write("{0:10.0f}  {1:10.0f}  {2:10.2f}  << elem_1, elem_2, res_value\n".format(*reg_row))

        # Method-specific components
        if r2_options['job_type'] == 1 or r2_options['job_type'] in ['inverse','Inverse','inv','I','i']: # Inverse model options
            f.write("\n") # blank space
            
            if int(r2_options['mesh_type']) in [4, 5]:
                f.write("{0:7.0f}  {1:5.0f}  << no. patches in x, no. patches in y\n".format(*inv_dict['patch_size_xy']))
                if inv_dict['patch_size_xy'][0] == 0 and inv_dict['patch_size_xy'][1] == 0:
                    f.write("{0:7.0f}  {1:5.0f}  << num_param_x, num_param_y\n".format(len(inv_dict['npxy'][0]),len(inv_dict['npxy'][1])))
                    
                    # Write mesh parameters in x
                    f.write("{0:7.0f}  ".format(inv_dict['npxystart'][0]))
                    for ival,npxval in enumerate(inv_dict['npxy'][0]):
                        f.write("{}  ".format(npxval))
                        if ival==len(inv_dict['npxy'][0])-1:
                            f.write(" << npxval \n")
                        elif np.fmod(ival+1,ncols)==0 and ival>0:
                            f.write("\n")
                    
                    f.write("\n") # blank space
                    
                    # Write mesh parameters in y
                    f.write("{0:7.0f}  ".format(inv_dict['npxystart'][1]))
                    for ival,npyval in enumerate(inv_dict['npxy'][1]):
                        f.write("{}  ".format(npyval))
                        if ival==len(inv_dict['npxy'][1])-1:
                            f.write(" << npyval \n")
                        elif np.fmod(ival+1,ncols)==0 and ival>0:
                            f.write("\n")
                
            f.write("\n") # blank space
            
            # Line 18
            f.write("    {0:3.0f}  {1:5.2f}  << inverse_type, target_decrease\n".format(inv_dict['inverse_type'],
                                                                           inv_dict['target_decrease']))
            f.write("\n") # blank space
            if inv_dict['inverse_type'] == 3:
                f.write("{} << qual_ratio\n".format(inv_dict['qual_ratio']))
                f.write("\n") # blank space
                
                f.write("{} {} << rho_min, rho_max\n".format(inv_dict['rho_min'],inv_dict['rho_max']))
                f.write("\n") # blank space
                
            else:
                # Inverse model options
                # Regularization, Line 21
                f.write("{0:6.0f}{1:5.0f}  << data_type, reg_mode \n".format(inv_dict['data_type'],inv_dict['reg_mode']))
                f.write("\n") # blank space
                
                if inv_dict['reg_mode'] in [0,2]: # normal regularization
                    inv_type_txt = "     ".join([str(tempval) for tempval in [inv_dict['tolerance'],
                                              inv_dict['max_iterations'],inv_dict['error_mod'],
                                              inv_dict['alpha_aniso']]])
                    
                else: # 1 = background regularisation
                    inv_type_txt = "     ".join([str(tempval) for tempval in [inv_dict['tolerance'],
                                              inv_dict['max_iterations'],inv_dict['error_mod'],
                                              inv_dict['alpha_aniso'],inv_dict['alpha_s']]])
                # Write line 8, inverse type    
                f.write("{}  << tolerance, max_iterations, error_mod, alpha_aniso, (alpha_s)\n".format(inv_type_txt))
                f.write("\n") # blank space
                
           
                # Error variance model parameters, Line 23 (offset error, relative error)
                model_var_params = [inv_dict['a_wgt'],inv_dict['b_wgt'],inv_dict['rho_min'],inv_dict['rho_max']]
                f.write("{}  {}  {}  {}  << a_wgt, b_wgt, rho_min, rho_max\n".format(*model_var_params))
                f.write("\n") # blank space
                
                if 'param_symbol' in inv_dict.keys():
                    for paramy_row in range(inv_dict['param_symbol']):
                        for paramx_val in paramy_row:
                            f.write("{}  ".format(paramx_val))
                        f.write(" \n")
                    f.write(" \n")
        
            
        # Bounding polyline for output region, eventually default to foreground outline, Line 25
        f.write("{0:d}  << number of points in polyline\n".format(int(inv_dict['num_xy_poly'])))
        
        if inv_dict['num_xy_poly'] > 0: # Write poly coordinates if at least one output polyline
            for ipolyxy in inv_dict['xy_poly']:
                f.write("{}  {}  << x_poly, y_poly\n".format(ipolyxy[0],ipolyxy[1]))
        f.write("\n") # blank space
        f.write("\n") # blank space
        
        # Write electrode information
        f.write("{}  << num_electrodes\n".format(electrode_array.shape[0]))
        if r2_options['mesh_type'] == 3:
            for elec_row in electrode_array:
                f.write("{0:8.0f}  {1:8.0f}  << electrode number, mesh node\n".format(*elec_row))
        else:
            for elec_row in electrode_array:
                f.write("{0:8.0f}  {1:8.0f}  {2:8.0f}  << electrode number, column, row\n".format(*elec_row))
            
        f.write("\n") # blank space     

def write_r2_protocol(protocol_fname=None, job_type=None,
                   inv_dict=None, meas_data=None,
                   diff_inv=False):
    '''Write protocol.dat for R2.'''
    if isinstance(meas_data,list):
        meas_data = np.array(meas_data)
        
    with open(protocol_fname,'w') as f:
        
        f.write("{}  << num_ind_meas\n".format(meas_data.shape[0]))
        
        if job_type == 1 or job_type in ['inverse','Inverse','inv','I','i']: # Inverse model options
            
           
            # 0) measurement id, 1) P+ electrode elec_num, 2) P- electrode elec_num,
            # 3) C+ electrode elec_num, 4) C- electrode elec_num,
            # 5) Voltage over Current ratio (v_i_ratio), 6) background v_i_ratio,
            # 7) Data standard deviation (data_sd)
            if inv_dict['a_wgt']==0 and inv_dict['b_wgt']==0 and inv_dict['reg_mode']==2:
                inv_meas_format = "{0:8.0f}  {1:8.0f} {2:8.0f}  {3:8.0f} {4:8.0f}  {5:e} {6:f}  {7:f} \n"
                if meas_data.shape[1] == 7:
                    meas_data_temp = np.hstack([meas_data[:,:-1],np.zeros((meas_data.shape[0],1)),meas_data[:,-1].reshape((-1,1))])
            elif inv_dict['a_wgt']==0 and inv_dict['b_wgt']==0:
                # 0) measurement id, 1) P+ electrode elec_num, 2) P- electrode elec_num,
                # 3) C+ electrode elec_num, 4) C- electrode elec_num,
                # 5) Voltage over Current ratio (v_i_ratio),
                # 6) Data standard deviation (data_sd)
                inv_meas_format = "{0:8.0f}  {1:8.0f} {2:8.0f}  {3:8.0f} {4:8.0f}  {5:10.6f} {6:f} \n"
                meas_data_temp = meas_data.copy()
            elif diff_inv: # Difference inversion
                # 0) measurement id, 1) P+ electrode elec_num, 2) P- electrode elec_num,
                # 3) C+ electrode elec_num, 4) C- electrode elec_num,
                # 5) t0: Voltage over Current ratio (v_i_ratio)
                # 6) t0: Voltage over Current ratio (v_i_ratio)
                inv_meas_format = "{0:8.0f}  {1:8.0f} {2:8.0f}  {3:8.0f} {4:8.0f}  {5:10.6f}  {6:10.6f}  \n"
                meas_data_temp = meas_data[:,:7].copy()
            else:
                # 0) measurement id, 1) P+ electrode elec_num, 2) P- electrode elec_num,
                # 3) C+ electrode elec_num, 4) C- electrode elec_num,
                # 5) Voltage over Current ratio (v_i_ratio)
                inv_meas_format = "{0:8.0f}  {1:8.0f} {2:8.0f}  {3:8.0f} {4:8.0f}  {5:10.6f}  \n"
                meas_data_temp = meas_data[:,:6].copy()
                
        else: # Forward model solution
            # 0) measurement id, 1) P+ electrode elec_num, 2) P- electrode elec_num,
            # 3) C+ electrode elec_num, 4) C- electrode elec_num
            inv_meas_format = "{0:8.0f}  {1:8.0f} {2:8.0f}  {3:8.0f} {4:8.0f}  \n"
            meas_data_temp = meas_data.copy()
            
        for meas_row in meas_data_temp:
            f.write(inv_meas_format.format(*meas_row))    
        
