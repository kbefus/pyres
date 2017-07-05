# -*- coding: utf-8 -*-
"""
Created on Fri May 05 15:42:04 2017

@author: kbefus
"""
from __future__ import print_function
import numpy as np
from pyres import pyres_utils


# ----------------------- Mesh construction functions -----------------------


def make_region(geom_dict=None, start_line_num=None,start_pt_num=None,
                  active_domain=None,boundary_dict=None,extend_to=None,
                  region_xyzpts=None,clen=None):
    '''Extend mesh regions to a boundary.
    
    Inputs
    -----------
    
    geom_dict: dictionary
        Python dictionary containing mesh dictionaries:
                - points
                - boundaries
                - surflines
                - fore_lines
                
    start_line_num: int, np.int
        Starting number for creating additional lines.

    start_pt_num: int, np.int
        Starting number for creating additional points.        
        
    active_dict: str
        Active boundary dictionary: either 'background' or 'foreground'.
    
    '''
    
    points_dict = geom_dict['points'].copy()
    boundaries_dict = geom_dict['boundaries'].copy()
    surflines_dict = geom_dict['surflines'].copy()
    forelines_list = geom_dict['fore_lines']
  
    # Initialize region dictionaries
    region_points = {}
    region_lines = {}
    region_surflines = {}
    region_surfaces = {}

    for x,y,z in region_xyzpts:
        region_points.update({start_pt_num:[x,y,z,clen]})
        if 'order' not in region_points.keys():
            region_points.update({'order':[start_pt_num]})
        else:
            region_points['order'].append(start_pt_num)
        start_pt_num+=1
        
    save_main_points = np.array(region_points['order'])
    
    # Make lines
    for istart,iend in zip(save_main_points[:-1],save_main_points[1:]):
        region_lines.update({start_line_num:[istart,iend]})
        if 'order' not in region_lines.keys():
            region_lines.update({'order':[start_line_num]})
        else:
            region_lines['order'].append(start_line_num)
            
        start_line_num+=1
        
    if extend_to is not None: 
        if extend_to in ['foreground']:
            # Define foreground-specific data
            region_buffer_xy = boundary_dict['fdxdy']
            
        elif extend_to in ['background']:
            # Define background-specific data
            region_buffer_xy = boundary_dict['bdxdy']
        
        # prepare domain points
        all_eline_pts = np.unique([points_dict[ikey]['order'] for ikey in points_dict['order']])
        fdict= {'boundary_dict':boundaries_dict[active_domain],
                'exclude_array':all_eline_pts}
                
        all_domain_pts, domain_lines = find_boundary_pts(**fdict) 
        domain_pts = all_domain_pts[1:-1] # remove first and last
        
        # Make additional points to meet domain boundary
        x1 = region_xyzpts[-1,0]+region_buffer_xy[0]
        x2 = region_xyzpts[0,0]-region_buffer_xy[0]
        
        # First connector point to domain boundary
        region_points.update({start_pt_num:np.hstack([x1,region_xyzpts[-1,1:].flatten(),clen]).tolist()})
        region_points['order'].append(start_pt_num)
        region_edge_pts = [start_pt_num]
        
        start_pt_num+=1
        # Final connector point to domain boundary            
        region_points.update({start_pt_num:np.hstack([x2,region_xyzpts[0,1:].flatten(),clen]).tolist()})
        region_points['order'].append(start_pt_num)
        region_edge_pts.append(start_pt_num)
        start_pt_num+=1
        
        new_line_keys = [start_line_num,start_line_num+1]
        start_line_num+=2
        insert_dict = {'orig_dict':boundaries_dict[active_domain],
                       'in_keys':domain_lines[::2],
                       'new_pts':region_edge_pts,
                       'new_keys':new_line_keys
                       }
        # Update domain lines
        boundaries_dict[active_domain]=insert_key_halfway(**insert_dict)
        

        # Update foreground surface
        surflines_dict.update({surflines_dict['order'][0]:np.roll(boundaries_dict['foreground']['order'],1).tolist()}) # foreground boundary
        
        # Update background surface                
        foreline_order = -1*np.array(boundaries_dict[active_domain]['order'])[np.arange(boundaries_dict['foreground']['order'].index(forelines_list[0]),
                                      boundaries_dict[active_domain]['order'].index(forelines_list[-1])+2,
                                      dtype=int)][::-1]
        background_order = np.hstack([foreline_order,boundaries_dict['background']['order']])
        surflines_dict.update({surflines_dict['order'][1]:background_order.tolist()}) # background line

        # First line connectors to domain
        region_lines.update({start_line_num:[save_main_points[-1],region_edge_pts[0]]})
        region_lines['order'].append(start_line_num)
        start_line_num+=1
        
        region_lines.update({start_line_num:[region_edge_pts[0],domain_pts[0]]})
        region_lines['order'].append(start_line_num)
        start_line_num+=1
        
        # Connect domain points
        for istart,iend in zip(domain_pts[:-1],domain_pts[1:]):
            region_lines.update({start_line_num:[istart,iend]})

            region_lines['order'].append(start_line_num)
            start_line_num+=1
            
        # Final connectors
        region_lines.update({start_line_num:[domain_pts[-1],region_edge_pts[-1]]})
        region_lines['order'].append(start_line_num)
        start_line_num+=1
        
        region_lines.update({start_line_num:[region_edge_pts[-1],save_main_points[0]]})
        region_lines['order'].append(start_line_num)
        start_line_num+=1
    else:
        # Add final connecting line of internal feature
        region_lines.update({start_line_num:[save_main_points[-1],save_main_points[0]]})
        region_lines['order'].append(start_line_num)
            
        start_line_num+=1
    
    # Check region lines with other line dictionaries
    region_lines,ind_changes = check_dict_entries(boundaries_dict[active_domain],region_lines)
        
    # Write surface
    region_surflines.update({start_line_num:region_lines['order']})
    region_surflines.update({'order':[start_line_num]})
    start_line_num+=1
    
    region_surfaces.update({start_line_num:[start_line_num-1]})
    region_surfaces.update({'order':[start_line_num]})
    start_line_num+=1
    
    # Consolidate dictionaries into output dictionary only as needed (if changed)
    output_dict = {}
    output_dict.update({'boundaries':boundaries_dict})
    output_dict.update({'surflines':surflines_dict})
    # Region information
    output_dict.update({'region_dict':{}})
    output_dict['region_dict'].update({'active_domain':active_domain})
    output_dict['region_dict'].update({'region_points':region_points})
    output_dict['region_dict'].update({'region_lines':region_lines})
    output_dict['region_dict'].update({'region_surflines':region_surflines})
    output_dict['region_dict'].update({'region_surfaces':region_surfaces})

    return output_dict,start_line_num,start_pt_num
    
def quad_mesh_region(mesh_xy=[None,None],target_xypos=None, target_res=None,
                     target_firstnrows=None,region_elem_list=None,background_res=None):
    """Define quadrilateral mesh region."""
    xx_mesh,yy_mesh = mesh_xy
    if len(xx_mesh.shape) == 1:
        elem_numbering = np.reshape(np.arange(xx_mesh.shape[0]*yy_mesh.shape[0]),
                                       (yy_mesh.shape[0],xx_mesh.shape[0]),order='F')+1   
        XX,YY = np.meshgrid(xx_mesh,yy_mesh)
    else: # 2d mesh coordinates input
        elem_numbering = np.reshape(np.arange(xx_mesh.shape[0]*xx_mesh.shape[1]),
                                       (xx_mesh.shape[0],xx_mesh.shape[1]),order='F')+1  
        XX,YY = np.array(xx_mesh),np.array(yy_mesh)
    
    fill_cols = None
    
    if target_xypos is not None:
        all_xy_pairs = np.hstack([XX.reshape((-1,1)),YY.reshape((-1,1))])
        matched_inds = pyres_utils.match_xy(XY=all_xy_pairs,xy_to_match=np.array(target_xypos).astype(float)).squeeze()
        unraveled_inds = np.array(np.unravel_index(matched_inds,XX.shape)).T
    
        unique_cols = np.unique(unraveled_inds[:,1])
        fill_cols = np.arange(unique_cols[0],unique_cols[-1]) # can have +1 at end to keep last match.
        start_end_rows = []
        for unq_col in unique_cols:
            rows_found = np.sort(unraveled_inds[unraveled_inds[:,1]==unq_col,0])
            start_end_rows.append([rows_found[0],rows_found[-1]])
        
        start_end_rows = np.array(start_end_rows)
        start_rows_fill = np.round(np.interp(fill_cols,unique_cols,start_end_rows[:,0]),0).astype(int)
        end_rows_fill = np.round(np.interp(fill_cols,unique_cols,start_end_rows[:,1]),0).astype(int)
        
    elif target_firstnrows is not None:
        fill_cols = np.arange(XX.shape[1])
        start_rows_fill = np.zeros(XX.shape[1]-1)
        end_rows_fill = target_firstnrows*np.ones(XX.shape[1]-1)

    # Make background element region first
    if region_elem_list is None:
        # Initialize list with whole domain
        region_elem_list = [[1,int((XX.shape[0]-1)*(XX.shape[1]-1)),background_res]]
        
    # Add other regions if needed
    if fill_cols is not None:
        for icol,irowstart,irowend in zip(fill_cols,start_rows_fill,end_rows_fill):
            elems_found = np.sort(elem_numbering[irowstart:irowend+1,icol])
            if len(elems_found)>1:
                region_elem_list.append([elems_found[0]-icol,elems_found[-1]-icol-1,target_res])
            else:
                region_elem_list.append([elems_found[0]-icol,elems_found[0]-icol,target_res])
        
    return region_elem_list

def make_tri_region_elements(mesh_dict=None, target_zones=None,
                             main_res=None,target_res=None):
    """Define regions for triangular mesh."""
    
    # Make main region
    reg_elems = []
    reg_elems.append([1,len(mesh_dict['elements']),main_res])
    
    # Make target region(s)
    if target_zones is not None:
        if isinstance(target_zones,(int,float)):
            target_zones = [target_zones]
        if isinstance(target_res,(int,float)):
            target_res = [target_res]
            
        for target_zone,temp_res in zip(target_zones,target_res):
            all_elements = np.column_stack([np.arange(len(mesh_dict['elements']))+1,np.array(mesh_dict['elements'])])
            target_elements = np.sort(all_elements[mesh_dict['zones']==target_zone,0])
            reg_elems.append([target_elements[0],target_elements[-1],temp_res]) # if elements are in order
            
    return reg_elems   
        
# ----------------------- Mesh helper functions -------------------------------
def find_boundary_pts(boundary_dict=None,exclude_array=None):
    '''Find mesh points not in exclude_array.'''
    out_points = []
    out_edge_keys = []
    for ikey1 in boundary_dict.keys():
        temp_catch = True
        for ipt in boundary_dict[ikey1]:
            if ipt not in exclude_array and ipt not in out_points:
                out_points.extend([ipt])
            elif ipt in exclude_array:
                temp_catch=False
        if temp_catch:
            out_edge_keys.append(ikey1)
            
    return np.array(out_points),out_edge_keys


def load_region_xyz(fname,nheaders=2,delimiter=',',
                    col_scales = [1,1]):
    '''Load three column delimited text file.'''
    with open(fname,'r') as f:
        header_info = []
        data_list = []
        irow = 0
        for iline in f:
            if irow < nheaders:
                header_info.append(iline)
            else:
                data_list.append([float(ipiece.strip(' \t')) for ipiece in iline.split(delimiter)])
            irow += 1
    
    data_array = np.array(data_list)
    
    for icol,iscale in enumerate(col_scales):
        if icol < data_array.shape[-1]:
            if float(iscale) != 0.0: 
                data_array[:,icol] = data_array[:,icol]*float(iscale)
            
    return data_array,header_info    
    
def insert_key_halfway(orig_dict=None,in_keys=None,
                       new_pts=None,new_keys=None):
    '''Insert key into Python dictionary.'''   
    out_dict = orig_dict.copy()

    # Reassign key entries as ints
    if 'order' in out_dict.keys():
        dict_keys = out_dict['order']
    else:
        dict_keys = np.sort(list(out_dict.keys())).tolist()
    
    for in_key,new_pt,new_key in zip(in_keys,new_pts,new_keys):
        entry_orig = orig_dict[in_key]
        new_entries = dict(zip((in_key,new_key),
                          ([entry_orig[0],new_pt],
                           [new_pt,entry_orig[1]])))
        in_key_ind = dict_keys.index(in_key)
        dict_keys.insert(in_key_ind+1,new_key)
        out_dict.update(new_entries)
        
    out_dict.update({'order':dict_keys})

    return out_dict
    
    
def renumber_dict(in_dict=None,start_num=None):
    '''Renumber Python dictionary keys in a sorted order.'''
    if 'order' in in_dict.keys():
        dict_keys = in_dict['order']
    else:
        dict_keys = np.sort(list(in_dict.keys()))
    final_dict={}
    out_order = []
    icount_keys = start_num
    for ikey in dict_keys:
        final_dict[icount_keys] = in_dict[ikey]
        out_order.append(icount_keys)
        icount_keys+=1
        
    final_dict['order'] = out_order    
    return final_dict,icount_keys
    
def check_dict_entries(main_dict=None,other_dict=None):
    '''Check if dictionary entry exists in main dictionary.'''
    if 'order' in main_dict.keys():
        main_keys = main_dict['order']
    else:
        main_keys = np.sort(list(main_dict.keys())) 
        
    main_entries = np.array([main_dict[main_key] for main_key in main_keys])
    
    if 'order' in other_dict.keys():
        other_keys = other_dict['order']
    else:
        other_keys = np.sort(list(other_dict.keys())) 

    out_key_order =[]
    new_dict = other_dict.copy()
    replaced_keys = []
    for ikey,other_key in enumerate(other_keys):

        other_entry = new_dict[other_key]
        ind_found = pyres_utils.match_xy(main_entries,np.array(other_entry))
        if len(ind_found.squeeze())> 0:
            old_key = new_dict.pop(other_key,None)
            replaced_keys.append([other_key,main_keys[ind_found[0]]])
            out_key_order.append(main_keys[ind_found[0]])
            new_dict[main_keys[ind_found[0]]]=main_dict[main_keys[ind_found[0]]]
        else:
            out_key_order.append(other_key)
            
    new_dict.update({'order':out_key_order})
    return new_dict, replaced_keys
    