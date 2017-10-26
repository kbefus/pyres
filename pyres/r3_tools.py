# -*- coding: utf-8 -*-
"""
Created on Sat Apr 01 15:19:23 2017

@author: kbefus
"""
from __future__ import print_function

import numpy as np
import time, os

from pyres import pyres_utils

inv_defaults = {'inverse_type':0,'data_type':1,'tolerance':1.,'no_improve':1.0,
                'max_iterations':10,'error_mod':2,'alpha_aniso':1.0,
                'alpha_s':1.0,'cginv_tol':1e-4, 'cginv_maxits':5e2,
                'alpha_max':1e10,'num_alpha_steps':10,'min_step':1e-3,
                'a_wgt': 0.,'b_wgt':0.03,'num_xy_poly':0.}

class R3(object):
    '''Main R3 object.'''
    def __init__(self,exe_name='R3t.exe',exe_dir='',survey_name=None,
                 survey_dir=None, job_type=None, data_fpath=None,
                 ):
        '''Create R3 object.
        
        Inputs
        --------------------------
        
        job_type: int, np.int
            Job type for R3, where 
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

    
    def load_data(self,data_fmt='AGI',data_fpath=None,
                  inf_electrode_xy = [None,None],
                  string_nums=1):
        '''Load resistivity data file.'''
        
        # File to read is specified either in R3 initiation or as a variable
        # to load_data. The input to this function is given priority.
        if data_fpath is None:
            data_fpath = self.data_fpath
        else:
            self.data_fpath = data_fpath # Save data_fpath to R3 class variable.
        
        load_dict = {'fname':data_fpath,'data_fmt':data_fmt,
                     'inf_electrode_xy':inf_electrode_xy,
                     'string_nums':string_nums}
        self.meas_data,self.electrode_data = pyres_utils.load_er_data(**load_dict)
          
      
        
        
    def r3_in(self, in_fname='R3t.in', singularity_type=0,add_dir_bool=False,
              electrode_array=None, startingRfile=None, fwd_resis=1e2, inv_dict=inv_defaults,
              output_domain=None,reg_elems=None,zminmax=[-50.,0.]):
        '''Make R3t.in file.
        
        
        Inputs:
            
            inv_dict: dictionary
                Dictionary containing all inversion options (lines 6-15) in R3 documentation.
                Default is set to inv_defaults located in R3tools.py

        '''
        # Add directory if needed
        if os.path.dirname(in_fname) is '':
                R3tin_output_file = os.path.join(self.survey_dir,in_fname)
                
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
            self.inv_dict['num_xy_poly'] = output_domain.shape[0]
            self.inv_dict['xy_poly'] = output_domain.copy()
        
        if reg_elems is None and startingRfile is None:
            num_regions_flag = 1
        elif reg_elems is None and startingRfile is not None:
            num_regions_flag = 0
        else:
            num_regions_flag = len(reg_elems)
        
        self.r3_options = {'job_type':self.job_type,'singularity_type':singularity_type}
        self.r3_write_dict = {'in_fname':R3tin_output_file,'survey_name':self.survey_name,
                              'r3_options':self.r3_options,
                              'num_regions_flag':num_regions_flag,'electrode_array':electrode_array,
                              'startingRfile':startingRfile,'fwd_resis':fwd_resis,'inv_dict':self.inv_dict,
                              'zminmax':zminmax,'reg_elems':reg_elems}
        
        write_r3_in(**self.r3_write_dict)
    
    def r3_protocol(self,protocol_fname='protocol.dat', meas_data=None):
        
                # Add directory if needed
        if os.path.dirname(protocol_fname) is '':
                protocol_fname = os.path.join(self.survey_dir,protocol_fname)
                
        protocol_dict = {'protocol_fname':protocol_fname,
                         'job_type':self.job_type,
                         'inv_dict':self.inv_dict,
                         'meas_data':meas_data}
        write_r3_protocol(**protocol_dict)
   
    def run_r3(self):
        
        # Need to change working directory to project directory - shouldn't be necessary
        cwd = self.survey_dir
        os.chdir(cwd)
        
        # Construct command
        cmd_list = [self.exe_fpath]
        
        # Then run R3t.exe in that folder
        run_dict = {'cmd_list':cmd_list,'cwd':cwd,
                    'normal_msg':'','report':True}
        self.r3_success,self.run_output = pyres_utils.run_cmd(**run_dict)
        
        pass    
    
    def run_all(self,load_bool=False,load_dict=None,
                protocol_dict=None, r3_in_dict=None,
                run_bool=False):
        if load_bool:
            self.load_data(**load_dict)
        
        self.r3_in(**r3_in_dict)
        self.r3_protocol(**protocol_dict)
        if run_bool:
            self.run_r3()
                
# --------------------- Helper functions -----------------------                
def write_r3_in(in_fname='R3t.in', survey_name=None ,r3_options=None, num_regions_flag=1,
              electrode_array=None, startingRfile=None, fwd_resis=None, inv_dict=inv_defaults,
              zminmax=[-50.,0.],reg_elems=None):
    
    with open(in_fname,'w') as f:
            
        # Header information
        header = "Project: {}, created using pyres on {}\n".format(\
                    survey_name,time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()))
        
        f.write(header)
        f.write("\n") # blank space
        f.write("{}  {}  << job_type, singularity_type\n".format(r3_options['job_type'],r3_options['singularity_type']))
        f.write("\n") # blank space
        
        f.write("{}  << num_regions\n".format(num_regions_flag))
        if num_regions_flag == 0 and startingRfile is not None:
            f.write("{}\n".format(startingRfile))
            f.write("\n") # blank space
        else:
            f.write("{}  << constant resistivity applied to all elements\n".format(fwd_resis))

        # Method-specific components
        if r3_options['job_type'] == 1 or r3_options['job_type'] in ['inverse','Inverse','inv','I','i']: # Inverse model options
            f.write("\n") # blank space
            f.write("{0:d}  << data type".format(inv_dict['data_type']))
            f.write("\n") # blank space
            f.write("\n") # blank space
        
            f.write("{0:d}  << inverse_type\n".format(inv_dict['inverse_type']))
            f.write("\n") # blank space
            f.write("\n") # blank space
            
            if inv_dict['inverse_type'] == 0: # normal regularization
                inv_type_txt = "  ".join([str(tempval) for tempval in [inv_dict['tolerance'],inv_dict['no_improve'],
                                          inv_dict['max_iterations'],inv_dict['error_mod'],
                                          inv_dict['alpha_aniso']]])
                
            else: # 1 = background regularisation, 2 = difference regularisation
                inv_type_txt = "  ".join([str(tempval) for tempval in [inv_dict['tolerance'],inv_dict['no_improve'],
                                          inv_dict['max_iterations'],inv_dict['error_mod'],
                                          inv_dict['alpha_aniso'],inv_dict['alpha_s']]])
            # Write line 8, inverse type    
            f.write("{}  << tolerance, no_improve, max_iterations, error_mod, alpha_aniso, (alpha_s)\n".format(inv_type_txt))
            f.write("\n") # blank space
            
            # Inverse model options
            # Conjugate gradient options, Line 9
            f.write("{}  {}  << cginv_tolerance, cginv_maxits\n".format(inv_dict['cginv_tol'],inv_dict['cginv_maxits']))
            f.write("\n") # blank space
            
            # Smoothing options, Line 10
            f.write("{0:e}  {1:d}  << alpha_max, num_alpha_steps\n".format(inv_dict['alpha_max'],inv_dict['num_alpha_steps']))
            f.write("\n") # blank space
            
            # Minimum step length, Line 11 (usually 1e-3 to 1e-2)
            f.write("{}  << min_step\n".format(inv_dict['min_step']))
            f.write("\n") # blank space
        
            # Error variance model parameters, Line 12 (offset error, relative error)
            f.write("{}  {}  << a_wgt, b_wgt\n".format(inv_dict['a_wgt'],inv_dict['b_wgt']))
            f.write("\n") # blank space
        
            # Output region definition
            if 'z_min' in inv_dict.keys() and 'z_max' in inv_dict.keys():
                # assume both vertical bounds already defined
                f.write("{}  {}  << z_min, z_max\n".format(inv_dict['z_min'],inv_dict['z_max']))
            else:
                # Define z_min and z_max based on electrode data
                f.write("{}  {}  << z_min, z_max\n".format(zminmax[0],zminmax[1]))
            
            # Bounding polyline for output region, eventually default to foreground outline
            f.write("{0:d}  << number of points in polyline\n".format(int(inv_dict['num_xy_poly'])))
            
            if inv_dict['num_xy_poly'] > 0: # Write poly coordinates if at least one output polyline
                for ipolyxy in inv_dict['xy_poly']:
                    f.write("{}  {}  << x_poly, y_poly\n".format(ipolyxy[0],ipolyxy[1]))
            f.write("\n") # blank space
            f.write("\n") # blank space
        
        # Write electrode information
        f.write("{}  << num_electrodes\n".format(electrode_array.shape[0]))
        icount=0
        for elec_row in electrode_array:
            if icount==0:
                f.write("{0:8.0f} {1:8.0f} {2:8.0f}  << string number, electrode number, mesh node\n".format(*elec_row))
                icount+=1
            else:
                f.write("{0:8.0f} {1:8.0f} {2:8.0f}\n".format(*elec_row))
            
        f.write("\n") # blank space 

def write_r3_protocol(protocol_fname=None, job_type=None,
                   inv_dict=None, meas_data=None):
    
    with open(protocol_fname,'w') as f:
        
        f.write("{}  << num_ind_meas\n".format(meas_data.shape[0]))
        
        if job_type == 1 or job_type in ['inverse','Inverse','inv','I','i']: # Inverse model options
            
            if inv_dict['a_wgt'] == 0 and inv_dict['b_wgt'] == 0:
                # 0) measurement id, 1) P+ electrode "string" and 2) elec_num, 3) P- electrode "string" and 4) elec_num,
                # 5) C+ electrode "string" and 6) elec_num, 7) C- electrode "string" and 8) elec_num,
                # 9) measured resistance, 10) resistance error
                inv_meas_format = "{0:8.0f}  {1:8.0f} {2:8.0f}  {3:8.0f} {4:8.0f}  {5:8.0f} {6:8.0f}  {7:8.0f} {8:8.0f}  {9:f} {10:f}\n"
            else: 
                # 0) measurement id, 1) P+ electrode "string" and 2) elec_num, 3) P- electrode "string" and 4) elec_num,
                # 5) C+ electrode "string" and 6) elec_num, 7) C- electrode "string" and 8) elec_num,
                # 9) measured resistance
                inv_meas_format = "{0:8.0f}  {1:8.0f} {2:8.0f}  {3:8.0f} {4:8.0f}  {5:8.0f} {6:8.0f}  {7:8.0f} {8:8.0f}  {9:f}\n"
            
        else: # Forward model solution
            # 0) measurement id, 1) P+ electrode "string" and 2) elec_num, 3) P- electrode "string" and 4) elec_num,
            # 5) C+ electrode "string" and 6) elec_num, 7) C- electrode "string" and 8) elec_num
            inv_meas_format = "{0:8.0f}  {1:8.0f} {2:8.0f}  {3:8.0f} {4:8.0f}  {5:8.0f} {6:8.0f}  {7:8.0f} {8:8.0f}\n"
            if meas_data.shape[1] > 9:
                meas_data = meas_data[:,:9]
                
        if meas_data[0,0]==0:
            meas_data[:,0] = meas_data[:,0]+1 # start with measurement 1
                
        for meas_row in meas_data:
            f.write(inv_meas_format.format(*meas_row))
        
        
        
        