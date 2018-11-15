# -*- coding: utf-8 -*-
"""
Created on Fri Apr 07 14:29:53 2017

@author: kbefus
"""
from __future__ import print_function
import os, sys
import numpy as np
from scipy.interpolate import griddata

############### Apparent resistivity to Resistance functions ################

# Dipole-dipole
app_res_to_res_dipdip = lambda a,n,app_res: app_res/(np.pi*a*n*(n+1)*(n+2))

############### ER data format variables ###########################
str_fmt = str
float_fmt = float
int_fmt = int

######### Example of formating for field data #########
# A,B are current electrodes
# M,N are potential electrodes
stg_cols_3d = ['id','user','date','time','resis','resis_error',
               'currentma','appres','line_id',
               'Ax','Ay','Az','Bx','By','Bz',
               'Mx','My','Mz','Nx','Ny','Nz',
               'cmdlinetxt','hv','cyk',
               'mtime','gain','mchantxt'] # last 5 columns not always used

stg_fmt_3d = {'id':int_fmt,'user':str_fmt,'date':str_fmt,'time':str_fmt,
              'resis':float_fmt,'resis_error':int_fmt,
              'currentma':int_fmt,'appres':float_fmt,'line_id':str_fmt,
              'Ax':float_fmt,'Ay':float_fmt,'Az':float_fmt,
              'Bx':float_fmt,'By':float_fmt,'Bz':float_fmt,
              'Mx':float_fmt,'My':float_fmt,'Mz':float_fmt,
              'Nx':float_fmt,'Ny':float_fmt,'Nz':float_fmt,
              'cmdlinetxt':str_fmt,'hv':str_fmt,'cyk':str_fmt,
              'mtime':str_fmt,'gain':str_fmt,'mchantxt':str_fmt}
#----------------------------------------------------------------------
# Syscal format assumes Spa.1=Ax, Spa.2=Bx, Spa.3=Mx, Spa.4=Nx
# Note: may not be complete
syscal_cols = ['id','Ax','Bx','Mx','Nx','resis','resis_error','M',
               'Vp','In','Time','	Ay','By','My','Ny','Stack',
               'Check','Vab','Pab','Rab','Latitude','Longitude','	Channel',
               'Bat','Rx-Bat','Temp.','Date','Gapfiller']
syscal_fmt = {'id':str_fmt,'Ax':float_fmt,'Bx':float_fmt,
              'Mx':float_fmt,'Nx':float_fmt,'resis':float_fmt,
              'resis_error':float_fmt,'M':float_fmt,
              'Vp':float_fmt,'In':float_fmt,'Time':float_fmt,
              'Ay':float_fmt,'	By':float_fmt,
              'My':float_fmt,'Ny':float_fmt,'Stack':float_fmt,
              'Check':float_fmt,'Vab':float_fmt,'Pab':float_fmt,
              'Rab':float_fmt,'Latitude':float_fmt,'Longitude':float_fmt,
              'Channel':float_fmt,'Bat':float_fmt,'Rx-Bat':float_fmt,
              'Temp.':float_fmt,'Date':str_fmt,'Gapfiller':float_fmt}

#--------------------------------------------------------
# Lippmann
lipp_cols = ['id','A','B','M','N','I','U',
             'dU','U90','dU90','rho','phi',
             'f','n','nAB','Profile','Spread','PseudoZ',
             'X','Y','Z','date','time','U(Tx)']
lipp_fmt = {'id':int_fmt,'A':int_fmt,'B':int_fmt,'M':int_fmt,'N':int_fmt,'I':float_fmt,'U':float_fmt,
             'dU':float_fmt,'U90':float_fmt,'dU90':float_fmt,'rho':float_fmt,'phi':float_fmt,
             'f':float_fmt,'n':int_fmt,'nAB':int_fmt,'Profile':float_fmt,'Spread':float_fmt,'PseudoZ':float_fmt,
             'X':float_fmt,'Y':float_fmt,'Z':float_fmt,'date':str_fmt,'time':str_fmt,'U(Tx)':float_fmt}
#--------------------------------------------------------
# Geotomo/Res2dinv files

geotomo_types = {1:'wenner',2:'polepole',3:'dipdip',6:'poledip',7:'schlum',8:'eqdipdip',11:'general',}
geotomo_cols={'header_default':{0:'name',1:'e_spacing',2:'array_type',3:'n_data',4:'x_type',5:'ip_flag'},
             'header_general':{0:'name',1:'e_spacing',2:'array_type',3:'sub_array_type',4:'meas_type_label',5:'app_res_flag',
                               6:'n_data',7:'x_type',8:'ip_flag'},
             'dipdip':['x','a','n','rho'],
             'default':['x','a','rho'],
             'general':['nelect','Bx','By','Ax','Ay','Mx','My','Nx','Ny','rho_or_resis']}
geotomo_fmt={'header':{'name':str_fmt,'e_spacing':float_fmt,'array_type':int_fmt,
                        'n_data':int_fmt,'x_type':int_fmt,'ip_flag':int_fmt,
                        'sub_array_type':int_fmt,'app_res_flag':int_fmt,
                        'meas_type_label':str_fmt},
              'data':{'x':float_fmt,'a':float_fmt,'n':int_fmt,'rho':float_fmt,
                      'Ax':float_fmt,'Ay':float_fmt,'Az':float_fmt,
                      'Bx':float_fmt,'By':float_fmt,'Bz':float_fmt,
                      'Mx':float_fmt,'My':float_fmt,'Mz':float_fmt,
                      'Nx':float_fmt,'Ny':float_fmt,'Nz':float_fmt,
                      'resis':float_fmt,'resis_error':float_fmt,
                      'rho_or_resis':float_fmt,'nelect':int_fmt}}
#--------------------------------------------------------
all_fmt_dict = {'agi':{'fmt':stg_fmt_3d,'cols':stg_cols_3d,
                       'delimiter':',','nheader':3},
                'syscal':{'fmt':syscal_fmt,'cols':syscal_cols,
                          'delimiter':',','nheader':1},
                'lipp':{'fmt':lipp_fmt,'cols':lipp_cols,
                        'delimiter':None,'nheader':269},
                'geotomo':{'fmt':geotomo_fmt,'cols':geotomo_cols,
                           'delimiter':[',','\t',' '],'nheader':6}}

################ Post-processing utilities ################
def doi_oldenburg_li_func(res_array1=None,res_array2=None,
             ref_res1=None, ref_res2=None):
    '''Calculate depth of investigation.'''
    doi_array = (res_array1-res_array2)/(ref_res1-ref_res2)
    return doi_array

def doi_analysis(model_fnames=None, model_dirs=None, ref_res_array=None,
             same_discretization=True,doi_func=doi_oldenburg_li_func):
    '''Calculate depth of investigation after Oldenburg and Li (1999).'''

    model_data = []
    if model_fnames is None:
        for model_dir in model_dirs:
            model_data.append(load_inv_output(work_dir=model_dir))
    else:
        for model_fname in model_fnames:
            model_data.append(load_inv_output(fname=model_fname))
            
    if same_discretization: # no interpolation required
        doi_dict = {'res_array1':model_data[0][:,2],
                    'res_array2':model_data[1][:,2],
                    'ref_res1':ref_res_array[0],
                    'ref_res2':ref_res_array[1]}
        doi_array = doi_func(**doi_dict)
    
        xy_array = model_data[0][:,:2]
        
    return xy_array,doi_array
    
    

################# I/O utilities #####################################
def load_fwd_output(fname='R2_forward.dat',work_dir=None,delimiter=' ',
                    nheader=1,ifile=0,app_res_bool=False):
    '''Load R2 or R3 forward model output file.'''
    if fname is None:
        import glob
        fname = glob.glob(os.path.join(work_dir,'*.dat'))[ifile]
    elif os.path.dirname(fname) in ['']:
        fname=os.path.join(work_dir,fname)
    
    with open(fname,'r') as i_file:
        header_info = []
        data_list = []
        irow = 0
        for iline in i_file:
            if irow < nheader:
                header_info.append(iline)
            else:
                data_list.append([float(ipiece.strip('\t').strip('\n')) for ipiece in iline.split(delimiter) if ipiece.strip() not in ['']])
            irow += 1
    data_out = np.array(data_list)
    
    # Remove apparent resistivity column
    if not app_res_bool and 'R2' in fname:
        data_out = data_out[:,:6]
    elif not app_res_bool and 'R3' in fname:
        data_out = data_out[:,:10]
    
    
    return data_out

def load_inv_output(fname='f001_res.dat',work_dir=None,delimiter=' ',
                    nheader=0,ifile=0):
    '''Load R2 inversion output file.'''
    if fname is None:
        import glob
        fname = glob.glob(os.path.join(work_dir,'*.dat'))[ifile]
    elif os.path.dirname(fname) in ['']:
        fname = os.path.join(work_dir,fname)
    
    with open(fname,'r') as r_file:
        header_info = []
        data_list = []
        irow = 0
        for iline in r_file:
            if irow < nheader:
                header_info.append(iline)
            else:
                data_list.append([float(ipiece.strip('\t').strip('\n')) for ipiece in iline.split(delimiter) if ipiece.strip() not in ['']])
            irow += 1
    data_out = np.array(data_list)
    
    return data_out

def read_trn(work_dir=None,fname=None,
             nheader=3,delimiter=',',ifile=0):
    '''Read terrain xyz delimited text file.'''
    if fname is None:
        import glob
        fname = glob.glob(os.path.join(work_dir,'*.trn'))[ifile]
    
    with open(fname,'r') as trn_file:
        header_info = []
        data_list = []
        irow = 0
        for iline in trn_file:
            if irow < nheader:
                header_info.append(iline)
            else:
                data_list.append([float(ipiece.strip(' \t').strip('\n')) for ipiece in iline.split(delimiter)])
            irow += 1
    data_out = np.array(data_list)
    
    return data_out
                           
def load_er_data(fname=None,data_fmt=None, inf_electrode_xy = [None,None],
                 string_nums=1,fmt_dict=None, nan_val=None, skip_vals=['\n'],
                 error_percent=None):
    '''Load resistivity data.'''
    if fmt_dict is None:
        if data_fmt.lower() in ['agi']:
            fmt_dict = all_fmt_dict['agi']
            dtype = 'agi'
        elif data_fmt.lower() in ['geotomo','res2dinv','res2d']:
            fmt_dict = all_fmt_dict['geotomo']
            header_dict = {} # store header information in dictionary
            dtype = 'geotomo'
        elif data_fmt.lower() in ['lipp','lippmann']:
            fmt_dict = all_fmt_dict['lipp']
            electrode_dict={}
            nan_val = '-'
            dtype = 'lipp'         

    # Load data to an array
    with open(fname,'r') as er_file:
        header_info = []
        data_list = []
        topo_data = []
        irow = 0
        for iline in er_file:
            if iline not in skip_vals:
                if irow < fmt_dict['nheader']:
                    header_info.append(iline)
                    if dtype in ['lipp']:
                        # Find electrode positions
                        if 'Electrode [' in iline:
                            temp_Edata = iline.split('=')
                            pos_data = [float(itemp) for itemp in temp_Edata[1].split()]
                            Enum = int(temp_Edata[0].split('[')[-1].split(']')[0])
                            electrode_dict.update({Enum:pos_data})
                    elif dtype in ['geotomo']:
                        if 'array_type' not in list(header_dict,):
                            header_dict.update({fmt_dict['cols']['header_default'][irow]:fmt_dict['fmt']['header'][fmt_dict['cols']['header_default'][irow]](iline.strip('\n'))})
                        else:

                            if header_dict['array_type']==11: # general has more information
                                fmt_dict['nheader'] = len(fmt_dict['cols']['header_general']) # update
                                header_dict.update({fmt_dict['cols']['header_general'][irow]:fmt_dict['fmt']['header'][fmt_dict['cols']['header_general'][irow]](iline.strip('\n'))})
                            else:
                                fmt_dict['nheader'] = len(fmt_dict['cols']['header_default']) # update
                                header_dict.update({fmt_dict['cols']['header_default'][irow]:fmt_dict['fmt']['header'][fmt_dict['cols']['header_default'][irow]](iline.strip('\n'))})
                else:
                    if dtype in ['geotomo']:
                        if irow < len(header_dict.keys())+header_dict['n_data']:
                            if len(fmt_dict['delimiter'])==1:
                                data_list.append([ipiece.strip(' \t') for ipiece in iline.strip('\n').split(fmt_dict['delimiter']) if len(ipiece)>0])
                            else:
                                for idelim in fmt_dict['delimiter'][1:]:
                                    iline = iline.replace(idelim,fmt_dict['delimiter'][0])
                                
                                data_list.append([ipiece.strip(' \t') for ipiece in iline.strip('\n').split(fmt_dict['delimiter'][0]) if len(ipiece)>0])


                        elif int(iline.strip('\n')) != 0:
                            print(iline)
                            # Read in supporting data (e.g., topography)
                            if 'topography' in iline.lower():
                                active_type = 'topo'
                                itopo_row = 0
                                ntopo_header = 2
                                ntopo = 100 # place holder
                            
                            if itopo_row == 1:
                                topo_type = int(iline.strip('\n'))
                            elif itopo_row == 2:
                                ntopo = int(iline.strip('\n'))
                            elif itopo_row < ntopo+2 and itopo_row != 0:
                                if len(fmt_dict['delimiter'])==1:
                                    topo_data.append([ipiece.strip(' \t').strip('\n') for ipiece in iline.split(fmt_dict['delimiter']) if len(ipiece)>0])
                                else:
                                    for idelim in fmt_dict['delimiter'][1:]:
                                        iline = iline.replace(idelim,fmt_dict['delimiter'][0])
                                    topo_data.append([ipiece.strip(' \t').strip('\n') for ipiece in iline.split(fmt_dict['delimiter'][0]) if len(ipiece)>0])
                            
                            
                            
                            if active_type in ['topo']:
                                itopo_row += 1
                            
                    else:
                        data_list.append([ipiece.strip(' \t').strip('\n') for ipiece in iline.split(fmt_dict['delimiter'])])
                irow += 1
    
    # Convert data types to desired formats
    data_array = np.array(data_list)

    if nan_val is not None:
        data_array[data_array==nan_val] = np.nan
    data_dict = {}
    if dtype in ['geotomo']:
        array_name = geotomo_types[header_dict['array_type']]

        if array_name not in list(fmt_dict['cols'].keys()):
            array_name = 'default'
            
        for icol,active_col in enumerate(fmt_dict['cols'][array_name][:len(data_list[0])]):

            data_dict.update({active_col:data_array[:,icol].astype(fmt_dict['fmt']['data'][active_col])})
    else:
        for icol,active_col in enumerate(fmt_dict['cols'][:len(data_list[0])]):
            data_dict.update({active_col:data_array[:,icol].astype(fmt_dict['fmt'][active_col])})
    

    position_keys = ['A','B','M','N']
    
    # For Lippmann data, need to create Ax to Ny
    if dtype in ['lipp']:
        for pos_key in position_keys:
            temp_xyz = np.array([electrode_dict[d1] for d1 in data_dict[pos_key]])
            data_dict['{}x'.format(pos_key)] = temp_xyz[:,0]
            data_dict['{}y'.format(pos_key)] = temp_xyz[:,1]
            data_dict['{}z'.format(pos_key)] = temp_xyz[:,2]
        
        # Need to also define resis (V/I)
        data_dict['resis'] = data_dict['U']/data_dict['I']
    
    # For geotomo data, need to massage data depending on array type of data input
    if dtype in ['geotomo']:
        if array_name in ['dipdip']:
            # Need to assign electrode positions
            data_dict['Ax'] = data_dict['x']
            data_dict['Bx'] = data_dict['Ax'] + data_dict['a']
            data_dict['Mx'] = data_dict['Bx'] + data_dict['a']*data_dict['n']
            data_dict['Nx'] = data_dict['Mx'] + data_dict['a']
            data_dict['resis'] = app_res_to_res_dipdip(data_dict['a'],data_dict['n'],data_dict['rho'])
            if error_percent is not None:
                data_dict['resis_error'] = error_percent*np.ones_like(data_dict['x'])
            else:
                data_dict['resis_error'] = 0*np.ones_like(data_dict['x'])
            # Force all *y locations to 0, need general array to specify
            data_dict['Ay'] = np.zeros_like(data_dict['x'])
            data_dict['By'] = np.zeros_like(data_dict['x'])
            data_dict['My'] = np.zeros_like(data_dict['x'])
            data_dict['Ny'] = np.zeros_like(data_dict['x'])
            
        elif array_name in ['general'] and header_dict['app_res_flag']==1:
            data_dict['a'] = data_dict['Bx']-data_dict['Ax']
            data_dict['n'] = (data_dict['Mx']-data_dict['Bx'])/data_dict['a']
            
            data_dict['resis'] = app_res_to_res_dipdip(data_dict['a'],data_dict['n'],data_dict['rho_or_resis'])
            if error_percent is not None:
                data_dict['resis_error'] = (error_percent/1e2)*np.ones_like(data_dict['Ax'])
            else:
                data_dict['resis_error'] = 0*np.ones_like(data_dict['Ax'])
        
    inf_edict = {}
    if inf_electrode_xy[0] is not None: # find appearances of inf electrode
        all_x = []
        all_y = []
        for pos_key in position_keys:
            match_dict = {'XY':np.column_stack([data_dict['{}x'.format(pos_key)],
                                                data_dict['{}y'.format(pos_key)]]),
                          'xy_to_match':np.array(inf_electrode_xy)}
            
            inf_edict[pos_key] = match_xy(**match_dict)
            
            temp_mask = np.ones_like(data_dict['{}x'.format(pos_key)],dtype=bool)
            temp_mask[inf_edict[pos_key]] = False
            all_x.extend(data_dict['{}x'.format(pos_key)][temp_mask])
            all_y.extend(data_dict['{}y'.format(pos_key)][temp_mask])
                      
    
    else:
        all_x = np.hstack([data_dict['Ax'],data_dict['Bx'],data_dict['Mx'],data_dict['Nx']])
        all_y = np.hstack([data_dict['Ay'],data_dict['By'],data_dict['My'],data_dict['Ny']])
#        all_z = np.hstack([data_dict['Az'],data_dict['Bz'],data_dict['Mz'],data_dict['Nz']])
        
    # Electrodes only have x,y coordinates, no vertical electrodes
    x_elec_positions = np.unique(all_x)
    y_elec_positions = np.unique(all_y)
    X,Y = np.meshgrid(x_elec_positions,y_elec_positions,indexing='xy')
    elec_nums = np.reshape(np.arange(X.size),X.shape)+1.# not specified by some companies, can be any order
    
    # Need to map data to electrode numbers
    electode_position_data = np.hstack([elec_nums.reshape((-1,1)),
                                        X.reshape((-1,1)),
                                        Y.reshape((-1,1)),
                                        np.zeros((Y.size,1))])
    if inf_electrode_xy[0] is not None:
        electode_position_data = np.vstack([electode_position_data,
                                             np.hstack([-9999,inf_electrode_xy,0.]).reshape((1,-1))])
        # Correct infinity electrode node number once mesh is created (np.inf until then)
        
    for ikey in position_keys:
        enum_key = '_'.join([ikey,'enum'])
        data_dict[enum_key] = np.zeros_like(data_dict['Ax'],dtype=int)
        
        # Loop through each electrode position and find matching measurements
        for elect_row in electode_position_data:    
            same_xy = (data_dict['{}x'.format(ikey)]==elect_row[1]) & (data_dict['{}y'.format(ikey)]==elect_row[2])
            data_dict[enum_key][same_xy] = int(elect_row[0]) 
        
    if dtype in ['agi']:
        data_dict['resis_error']=data_dict['resis_error']/10./1e2 # convert from 10ths of percent to ratio
    elif dtype in ['lipp']:
        data_dict['resis_error'] = (data_dict['dU']/1e2)
        
    meas_nums = np.arange(data_dict['resis'].shape[0],dtype=int).reshape((-1,1))+1
    
    if isinstance(string_nums,(int,float)):
        string_num = string_nums*np.ones((meas_nums.shape[0],4))
        
        epos_strnum = string_nums*np.ones((electode_position_data.shape[0],1))
    elif string_nums in ['x','by_x']:
        string_num = np.nan * np.ones((meas_nums.shape[0],4))
        epos_strnum = np.nan * np.ones((electode_position_data.shape[0],1))
        for ix,xval in enumerate(x_elec_positions):
            string_num[data_dict['Mx']==xval,0] = ix+1 # Don't use 0-based indexing
            string_num[data_dict['Nx']==xval,1] = ix+1
            string_num[data_dict['Ax']==xval,2] = ix+1
            string_num[data_dict['Bx']==xval,3] = ix+1
            epos_strnum[electode_position_data[:,1]==xval] = ix+1
    elif string_nums in ['y','by_y']:
        string_num = np.nan * np.ones((meas_nums.shape[0],4))
        epos_strnum = np.nan * np.ones((electode_position_data.shape[0],1))
        for iy,yval in enumerate(y_elec_positions):
            string_num[data_dict['My']==yval,0] = iy+1 # Don't use 0-based indexing
            string_num[data_dict['Ny']==yval,1] = iy+1
            string_num[data_dict['Ay']==yval,2] = iy+1
            string_num[data_dict['By']==yval,3] = iy+1
            epos_strnum[electode_position_data[:,2]==yval] = iy+1
    else:
        string_num = np.array(string_nums)
    
    # Add string numbers to electrode data
    electode_position_data = np.hstack([epos_strnum,
                                        electode_position_data])
        
    # Construct measured data array 
    meas_data = np.hstack([meas_nums.reshape((-1,1)),
                       string_num[:,0].reshape((-1,1)),
                       data_dict['M_enum'].reshape((-1,1)),
                       string_num[:,1].reshape((-1,1)),
                       data_dict['N_enum'].reshape((-1,1)),
                       string_num[:,2].reshape((-1,1)),
                       data_dict['A_enum'].reshape((-1,1)),
                       string_num[:,3].reshape((-1,1)),
                       data_dict['B_enum'].reshape((-1,1)),
                       data_dict['resis'].reshape((-1,1)), 
                       data_dict['resis_error'].reshape((-1,1))])
    return meas_data,electode_position_data

def load_r2out(work_dir=None):
    '''Load R2.out file.
    
    Inputs:
    -------
    
    work_dir: str
            path to directory containing R2.out
    
    Returns:
    -----------
    
    output_lines: list
            list containing all lines of R2.out
    '''
    r2out_fname = os.path.join(work_dir,'R2.out')
    
    with open(r2out_fname,'r') as r2:
        output_lines = []
        for iline in r2:
            output_lines.append(iline.rstrip())
    return output_lines

def write_res_dat(elem_array=None,startingRfile=None,
                          wdir=None):
    if wdir is not None:
        startingRfile = os.path.join(wdir,startingRfile)
    
    with open(startingRfile,'w') as fout:
        for iline in elem_array:
            fout.write('  {0:2.1f}  {1:2.1f}  {2:2.1f}   {3:10.2f}\n'.format(*iline))
            
################ Conversion utilities ################### 

def grid_inv_data(inv_data=None, nxny=[1e2,1e2],
                  maxy=0, miny=None,
                  interp_method='force_grid',inv_col=3):
    '''Make regular grid for plotting surface.'''
    
    if interp_method not in ['force_grid']:
        # Interpolate ER to regular grid
        # X data
        minx,maxx = np.min(inv_data[:,0]),np.max(inv_data[:,0])
        x = np.linspace(minx,maxx,nxny[0])
        dx = x[1]-x[0]
        x = x+dx/2. # to cell centers
        
        # Y data
        if miny is None:
            miny=np.min(inv_data[:,1])
        y = np.linspace(maxy,miny,nxny[1])
        dy = y[1]-y[0]
        y = y+dy/2. # to cell centers
    
        X,Y = np.meshgrid(x,y)
        
        ER = griddata(inv_data[:,:2],
                    inv_data[:,inv_col],
                    (X,Y),
                    method=interp_method)
    else:
        out_dict=force_grid(inv_data=inv_data)
        X = out_dict['x']
        Y = out_dict['y']
        if inv_col==2:
            ER = out_dict['inv']
        else:
            ER = out_dict['log_inv']
        
    return X,Y,ER

def force_grid(inv_data=None,col_names = ['x','y','inv','log_inv']):
    '''Make regular grid directly from inverse output.'''
    
    x = inv_data[:,0]
    dif_x = np.diff(x)
    # Outputs written by column, find each column
    col_ids = np.hstack([0,1+(dif_x>0.).nonzero()[0],len(x)])
    num_yvals = np.diff(col_ids)
    nrows = np.max(num_yvals)
    ncols = len(col_ids)-1
    
    col_names = col_names[:inv_data.shape[1]] # remove unused names
    out_dict = {}
    for iname,col_name in enumerate(col_names):
        new_array = np.nan*np.zeros((nrows,ncols))
        icol = 0
        for icount_col,(col_edge,row_edge) in enumerate(zip(col_ids[1:],num_yvals)):
            new_array[:row_edge,icount_col] = inv_data[icol:col_edge,iname]
            if col_name in ['x','y']:
                new_array[row_edge:,icount_col] = new_array[row_edge-1,icount_col] # fill lower rows with last x or y value
            icol=col_edge
            
        out_dict.update({col_name:new_array.copy()})
    return out_dict
    

def make_synthetic_survey(n_electrodes=None,array_type='dpdp',max_dipole=4,
                          max_separation=8):
    '''Make electrode combinations for survey.'''
    
    
    elect_array = np.arange(n_electrodes)+1
    if array_type in ['dpdp','dipole-dipole']:
        # Dipole-dipole array
        from itertools import combinations
        e_combos = np.array(list(combinations(elect_array,4))).astype(int)
        pair_test = (np.abs(np.diff(e_combos[:,:2],axis=1))<=max_dipole) & (np.abs(np.diff(e_combos[:,2:],axis=1))<=max_dipole)
        e_combos = e_combos[pair_test.squeeze(),:]
        sep_test = np.abs(np.mean(e_combos[:,:2],axis=1)-np.mean(e_combos[:,2:],axis=1)) <= max_separation
        e_combos = e_combos[sep_test.squeeze(),:]
    elif array_type.lower() in ['wenner']:
        e_combos = []
        for isep in range(1,max_separation+1):
            sep_array = np.column_stack([elect_array[isep:-(isep*2)],
                                         elect_array[isep*2:-isep],
                                         elect_array[0:-(isep*3)],
                                         elect_array[isep*3:]])
            e_combos.append(sep_array)
        
        e_combos = np.vstack(e_combos)
        
    # Add measurement identifier
    out_protocol = np.hstack([np.arange(e_combos.shape[0]).reshape((-1,1))+1,
                              e_combos])
    return out_protocol
        
############ Array utilities ###################

def match_xy(XY=None, xy_to_match=None):
    '''Find matching pairs of data by row.'''
    XY = np.ascontiguousarray(XY.copy())
    xy_to_match = np.ascontiguousarray(xy_to_match.copy()).reshape((-1,2))
    
    # Create view of arrays that treats each row as a single value
    pt1 = xy_to_match.view([('', xy_to_match.dtype)]*xy_to_match.shape[1])
    XY1 = XY.view([('',XY.dtype)]*XY.shape[1])
    
    # Find pt1 in XY1
    xy_match_ind = []
    for pt in pt1:
        ifind = (np.array(XY1==pt)).nonzero()[0]
        if len(ifind)>0:
            xy_match_ind.append(ifind)
    
    return np.array(xy_match_ind).astype(int)
    
def extrap(x, xp, yp, method='linear'):
    """np.interp function with linear extrapolation"""
    ind_order = np.argsort(xp)
    xp,yp = xp[ind_order],yp[ind_order]
    y = np.interp(x, xp, yp)
    if method in ['linear']:
        y[x < xp[0]] = yp[0] + (x[x<xp[0]]-xp[0]) * (yp[0]-yp[1]) / (xp[0]-xp[1])
        y[x > xp[-1]]= yp[-1] + (x[x>xp[-1]]-xp[-1]) * (yp[-1]-yp[-2])/(xp[-1]-xp[-2])
    elif method in ['constant']:
        y[x < xp[0]] = yp[0]
        y[x > xp[-1]]= yp[-1]
    return y 
    
    
############ Running external programs ###############

def run_cmd(cmd_list=None,cwd='./',async=False,
            silent=False, pause=False, report=False,
            normal_msg=None,failed_words = ["fail", "error"]):
    '''Run command in DOS.
    
    This function will run the cmd_list using subprocess.Popen.  It
    communicates with the process' stdout asynchronously/syncrhonously and reports
    progress to the screen with timestamps
    
    Parameters
    ----------
    cmd_list : list
        List of [Executable name (with path, if necessary), other args].
    cwd : str
        current working directory, where inputs are stored. (default is the
        current working directory - './')
    silent : boolean
        Echo run information to screen (default is True).
    pause : boolean, optional
        Pause upon completion (default is False).
    report : boolean, optional
        Save stdout lines to a list (buff) which is returned
        by the method . (default is False).
    normal_msg : str
        Normal termination message used to determine if the
        run terminated normally. (default is None)
    failed_words: list
        List of words to search for that indicates problem with running
        command. (default is ["fail","error"])
    async : boolean
        asynchonously read model stdout and report with timestamps.  good for
        models that take long time to run.  not good for models that run
        really fast
    Returns
    -------
    (success, buff)
    success : boolean
    buff : list of lines of executable output (i.e., stdout)
    
    
    Source: after flopy.mbase.run_model'''
    
    # Load libraries
    from datetime import datetime
    import subprocess as sp
    import threading
    if sys.version_info > (3, 0):
        import queue as Queue
    else:
        import Queue
    
    success = False
    buff = []

    # simple function for the thread to target
    def q_output(output, q):
        for line in iter(output.readline, b''):
            q.put(line)
            # time.sleep(1)
            # output.close()

    proc = sp.Popen(cmd_list,
                    stdout=sp.PIPE, stderr=sp.STDOUT, cwd=cwd)

    # Run executable and handle all output at once
    if not async:
        while True:
            line = proc.stdout.readline()
            c = line.decode('utf-8')
            if c != '':
                    if normal_msg in c.lower() and normal_msg is not None:
                        success = True
                    c = c.rstrip('\r\n')
                    if not silent:
                        print('{}'.format(c))
                    if report == True:
                        buff.append(c)
            else:
                break
        return success, buff
    
    # ------- Run exe and collect output while still running -------------
    # some tricks for the async stdout reading
    q = Queue.Queue()
    thread = threading.Thread(target=q_output, args=(proc.stdout, q))
    thread.daemon = True
    thread.start()

    last = datetime.now()
    lastsec = 0.
    while True:
        try:
            line = q.get_nowait()
        except Queue.Empty:
            pass
        else:
            if line == '':
                break
            line = line.decode().lower().strip()
            if line != '':
                now = datetime.now()
                dt = now - last
                tsecs = dt.total_seconds() - lastsec
                line = "(elapsed:{0})-->{1}".format(tsecs, line)
                lastsec = tsecs + lastsec
                buff.append(line)
                if not silent:
                    print(line)
                for fword in failed_words:
                    if fword in line:
                        success = False
                        break
        if proc.poll() is not None:
            break
    proc.wait()
    thread.join(timeout=1)
    buff.extend(proc.stdout.readlines())
    proc.stdout.close()

    # Look for success message in output if one provided
    if normal_msg is not None:
        for line in buff:
            if normal_msg in line:
                print("success")
                success = True
                break

    if pause:
        input('Press Enter to continue...')
    return success, buff 