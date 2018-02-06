# -*- coding: utf-8 -*-
"""
Created on Fri May 05 15:42:04 2017

@author: kbefus
"""
from __future__ import print_function
import numpy as np
from pyres import pyres_utils


# ----------------------- Mesh construction functions -----------------------


def make_points(mesh_obj=None,clen=None,start_pts=1,ndim=2):
    elec_list = []
    iline_order = []
    for iline in np.arange(mesh_obj.line_strexyz.shape[0]):    
        if iline not in mesh_obj.points.keys():
            mesh_obj.points[iline] = {}
        iline_order.append(iline)
        key_order = []
        for linestr,enum,x,y,z in mesh_obj.line_strexyz[iline]:
            mesh_obj.points[iline].update({start_pts:[x,y,z,clen]})
            key_order.append(start_pts)
            if ndim==2:
                elec_list.append([enum,start_pts]) # save data for protocol.dat
            elif ndim==3:
                elec_list.append([linestr,enum,start_pts])
            start_pts+=1
            
        mesh_obj.points[iline].update({'order':key_order})
        
    mesh_obj.points.update({'order':iline_order})    
    mesh_obj.count_pts = start_pts        
    mesh_obj.electrode_array = np.array(elec_list)

def make_lines(mesh_obj,start_lines=1):
        
        mesh_obj.start_lines=start_lines
        mesh_obj.lines = {}
        mesh_obj.surveyline_startend = {}
        mesh_obj.startpt,mesh_obj.endpt = [],[]
        key_order_se = []
        for ikey in mesh_obj.points['order']:
            # Collect consecutive points into lines in reverse order
            line_keys = mesh_obj.points[ikey]['order']
            
            # record start and end positions of line
            mesh_obj.startpt.append(line_keys[0])
            mesh_obj.endpt.append(line_keys[-1])
            
            # Save first line in section
            mesh_obj.surveyline_startend[ikey] = [mesh_obj.start_lines]
            key_order_se.append(ikey)
            
            key_order = []
            for istart,iend in zip(line_keys[:-1],line_keys[1:]):
                mesh_obj.lines.update({mesh_obj.start_lines:[istart,iend]})
                key_order.append(mesh_obj.start_lines)
                mesh_obj.start_lines+=1
                
            mesh_obj.lines.update({'order':key_order})
            # Save last line in section
            mesh_obj.surveyline_startend[ikey].append(mesh_obj.start_lines-1)
        
        mesh_obj.surveyline_startend.update({'order':key_order_se})
        
        # Connectors for the lines
        if len(mesh_obj.startpt)>1:
            for istart1,istart2 in zip(mesh_obj.startpt[:-1],mesh_obj.startpt[1:]):
                mesh_obj.lines.update({mesh_obj.start_lines:[istart1,istart2]})
                mesh_obj.lines['order'].append(mesh_obj.start_lines)
                mesh_obj.start_lines+=1
            
            for end1,iend2 in zip(mesh_obj.endpt[:-1],mesh_obj.endpt[1:]):
                mesh_obj.lines.update({mesh_obj.start_lines:[end1,iend2]})
                mesh_obj.lines['order'].append(mesh_obj.start_lines)
                mesh_obj.start_lines+=1    

def make_boundaries(mesh_obj,nforeground=4,nbackground=4,ndim=None):
        '''Make boundaries from survey bounds.
        
        
        Can eventually add capability to define non-square boundaries
        '''
        mesh_obj.boundaries = {}
        # First make foreground boundaries by extending some constant beyond
        # the corner points
        fore_corners = np.arange(mesh_obj.count_pts,mesh_obj.count_pts+nforeground)
        mesh_obj.count_pts = mesh_obj.count_pts+nforeground
        
        fore_pairs = np.c_[fore_corners,np.roll(fore_corners,-1)]
        
        mesh_obj.boundaries['foreground'] = {}
        key_order = []
        
        if ndim==2:
            # Include lines connecting electrodes in foreground
            mesh_obj.boundaries['foreground'].update(mesh_obj.lines.copy())
            key_order.extend(mesh_obj.lines['order'])
            
            # Start foreground at end of line
            mesh_obj.boundaries['foreground'].update({mesh_obj.start_lines:[mesh_obj.endpt[0],fore_pairs[-1,1]]})
            key_order.append(mesh_obj.start_lines)
            mesh_obj.start_lines += 1
            nforeground = nforeground-1
            
        mesh_obj.fore_lines=[]
        for ifore in np.arange(nforeground):
            mesh_obj.boundaries['foreground'].update({mesh_obj.start_lines:
                                                    fore_pairs[ifore,:].tolist()})
            key_order.append(mesh_obj.start_lines)
            mesh_obj.fore_lines.append(mesh_obj.start_lines)
            mesh_obj.start_lines += 1
        
        if ndim == 2: 
            # Complete foreground by connecting last foreground line to first electrode
            mesh_obj.boundaries['foreground'].update({mesh_obj.start_lines:[fore_pairs[-1,0],mesh_obj.startpt[0]]})
            key_order.append(mesh_obj.start_lines)
            mesh_obj.start_lines += 1
        
        # Save foreground order to mesh object
        mesh_obj.boundaries['foreground'].update({'order':key_order})
        
        # Background points and lines
        back_corners = np.arange(mesh_obj.count_pts,mesh_obj.count_pts+nbackground)
        mesh_obj.count_pts = mesh_obj.count_pts + nbackground
        
        back_pairs = np.c_[back_corners,np.roll(back_corners,-1)]
        key_order = []
        mesh_obj.boundaries['background'] = {}  
        if ndim==2:
            # Connect foreground to background
            mesh_obj.boundaries['background'].update({mesh_obj.start_lines:[fore_pairs[0,0],back_pairs[0,0]]})
            key_order.append(mesh_obj.start_lines)
            mesh_obj.start_lines += 1
            nbackground=nbackground-1
            
        for iback in np.arange(nbackground):
            mesh_obj.boundaries['background'].update({mesh_obj.start_lines:
                                                    back_pairs[iback,:].tolist()})
            key_order.append(mesh_obj.start_lines)
            mesh_obj.start_lines += 1
        
        if ndim==2:
            # Connect last background line to foreground
            mesh_obj.boundaries['background'].update({mesh_obj.start_lines:[back_pairs[-1,0],fore_pairs[-1,0]]})
            key_order.append(mesh_obj.start_lines)
            mesh_obj.start_lines += 1
        
        # Save background order to mesh object
        mesh_obj.boundaries['background'].update({'order':key_order})
        

def make_surfaces(mesh_obj,ndim=None):
    
    mesh_obj.surfaces = {}   
    seg_dict = {'nelectrodes':mesh_obj.nelectrodes_per_line,
                'nlines':mesh_obj.nlines,
                'topbot':True,
                'iline':np.arange(2*(mesh_obj.nlines-1))}
    mesh_obj.connecting_segs = segment_loc(**seg_dict).reshape((2,-1)).T

    mesh_obj.surflines = {} # will have "Line Loop(#) = {}
    key_order_sl = []
    key_order_s = []
    
    if ndim == 3: # R3 surfaces
        # Make foreground surfaces between individual lines
        for iloop in np.arange(mesh_obj.nlines-1):
            topbot_edges = mesh_obj.connecting_segs[iloop,:]
            l1_startend = mesh_obj.surveyline_startend[iloop+1]
            line1 = np.arange(l1_startend[0],l1_startend[1]+1,1)
            l2_startend = mesh_obj.surveyline_startend[iloop]
            line2 = np.arange(-l2_startend[1],-l2_startend[0]+1,1)
            out_lines = np.hstack([topbot_edges[0],
                                   line1,
                                   -topbot_edges[1],
                                   line2])
            
            mesh_obj.surflines.update({mesh_obj.start_lines:out_lines.tolist()})
            key_order_sl.append(mesh_obj.start_lines)
            mesh_obj.start_lines += 1
            
            mesh_obj.surfaces.update({mesh_obj.start_lines:[mesh_obj.start_lines-1]})
            key_order_s.append(mesh_obj.start_lines)
            mesh_obj.start_lines += 1

        # Make bounding foreground surface from foreground lines and line-bounding lines
        mesh_obj.surflines.update({mesh_obj.start_lines:np.roll(np.sort(list(mesh_obj.boundaries['foreground']['order'])),1).tolist()}) # foreground boundary
        key_order_sl.append(mesh_obj.start_lines)
        iforeground = mesh_obj.start_lines
        mesh_obj.start_lines += 1
        lastsurvey_lines = mesh_obj.surveyline_startend[mesh_obj.nlines-1]
        firstsurvey_lines = mesh_obj.surveyline_startend[0]
        survey_line_bounds = np.hstack([mesh_obj.connecting_segs[:,0].ravel(),
                                        np.arange(lastsurvey_lines[0],lastsurvey_lines[-1]+1),
                                        -mesh_obj.connecting_segs[::-1,1].ravel(),
                                        np.arange(-firstsurvey_lines[-1],-firstsurvey_lines[0]+1)])
        mesh_obj.surflines.update({mesh_obj.start_lines:np.roll(survey_line_bounds,-2).tolist()})# survey line boundary
        key_order_sl.append(mesh_obj.start_lines)
        mesh_obj.start_lines += 1
        mesh_obj.surfaces.update({mesh_obj.start_lines:[mesh_obj.start_lines-1,mesh_obj.start_lines-2]})
        key_order_s.append(mesh_obj.start_lines)
        mesh_obj.start_lines += 1
        
        # Make background surface between foreground and background lines
        mesh_obj.surflines.update({mesh_obj.start_lines:np.roll(np.sort(list(mesh_obj.boundaries['background']['order'])),1).tolist()}) # background line
        key_order_sl.append(mesh_obj.start_lines)
        mesh_obj.start_lines += 1
        mesh_obj.surfaces.update({mesh_obj.start_lines:[iforeground,mesh_obj.start_lines-1]})
        key_order_s.append(mesh_obj.start_lines)
        
    elif ndim==2:
        # Make bounding foreground surface from foreground lines and line-bounding lines
        mesh_obj.surflines.update({mesh_obj.start_lines:np.roll(mesh_obj.boundaries['foreground']['order'],1).tolist()}) # foreground boundary
        key_order_sl.append(mesh_obj.start_lines)
        mesh_obj.start_lines += 1
        mesh_obj.surfaces.update({mesh_obj.start_lines:[mesh_obj.start_lines-1]})
        key_order_s.append(mesh_obj.start_lines)
        mesh_obj.start_lines += 1
        
        # Make background surface between foreground and background lines
        foreline_order = -1*np.array(mesh_obj.fore_lines)[::-1]
        background_order = np.hstack([foreline_order,mesh_obj.boundaries['background']['order']])
        mesh_obj.surflines.update({mesh_obj.start_lines:background_order.tolist()}) # background line
        key_order_sl.append(mesh_obj.start_lines)
        mesh_obj.start_lines += 1
        
        mesh_obj.surfaces.update({mesh_obj.start_lines:[mesh_obj.start_lines-1]})
        key_order_s.append(mesh_obj.start_lines)
        
    # Record order into mesh_obj
    mesh_obj.surflines.update({'order':key_order_sl})
    mesh_obj.surfaces.update({'order':key_order_s})

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

def select_tet_region_nodes(mesh_dict=None,target_zones=None,
                            main_res=None, target_res=None,buffer_size=0.,
                            startingRfile=None):
    
    tet_elem_nodes = element_xyz(mesh_dict)
    tet_elem_centroids = tet_elem_nodes.mean(axis=1)
    
    reg_elems = None
    elem_array = None
    
    if not startingRfile:
        # Make region element list
        
        reg_elems.append([1,len(mesh_dict['elements']),main_res]) # Background
    else:
        # Set all elements to background
        elem_array = np.zeros([tet_elem_nodes.shape[0],4])
        elem_array[:,3] = main_res
        
    for target_zone,temp_res in zip(target_zones,target_res):
        # target zone has following shape for rectuangular prism
        # target_zone = [[x_min,x_max],
        #                [y_min,y_max],
        #                [z_min,z_max]]
        ifound = (tet_elem_centroids[:,0]+buffer_size >= target_zone[0][0]) &\
                 (tet_elem_centroids[:,0]-buffer_size <= target_zone[0][1]) &\
                 (tet_elem_centroids[:,1]+buffer_size >= target_zone[1][0]) &\
                 (tet_elem_centroids[:,1]-buffer_size <= target_zone[1][1]) &\
                 (tet_elem_centroids[:,2]+buffer_size >= target_zone[2][0]) &\
                 (tet_elem_centroids[:,2]-buffer_size <= target_zone[2][1]) 
        
        if not startingRfile:
            # Need to write continuous elements as a region...TBD
            pass
        else:
            elem_array[ifound,3] = temp_res
            
    if elem_array is not None:
        return elem_array
    elif reg_elems is not None:
        return reg_elems
        
# ----------------------- Mesh helper functions -------------------------------
def segment_loc(iline=None,nelectrodes=None,nlines=None,topbot=False):
    if topbot:
        return nlines*(nelectrodes-1) + iline + 1
    else:
        return iline*(nelectrodes-1)

def element_xyz(mesh_dict=None):
    '''Assign xyz coordinates to element nodes.'''
    all_elements = np.array(mesh_dict['elements'])
    all_element_nums = all_elements[:,:mesh_dict['npere']]
    uniq_elements = np.unique(all_element_nums)
    out_element_xyz = {}
    for ielement in uniq_elements:
        out_element_xyz.update({ielement:mesh_dict['nodes'][ielement-1]})
    
    out_xyz=[]
    for element_row in all_element_nums:
        element_xyz_temp = []
        for element_in_row in element_row:
            element_xyz_temp.append(out_element_xyz[element_in_row])
        out_xyz.append(element_xyz_temp)
    return np.array(out_xyz)

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

def find_line_orientation(mesh_obj=None):
    
    line_startpts = mesh_obj.startpt
    line_endpts = mesh_obj.endpt
    
    start_list = []
    end_list = []
    for iline,(stpt,endpt) in enumerate(zip(line_startpts,line_endpts)):
        start_list.append(mesh_obj.points[iline][stpt][:-1]) # don't need clen
        end_list.append(mesh_obj.points[iline][endpt][:-1])
    
    s_array = np.array(start_list)
    e_array = np.array(end_list)
    # Need to find dx dy for first and last lines
    dx_0 = e_array[0,0]-s_array[0,0]
    dy_0 = e_array[1,1]-e_array[0,1]
    dx_end = e_array[-1,0]-s_array[-1,0]
    dy_end = e_array[-1,1]-e_array[-2,1]
    
    out_multipliers = np.array([[-1,-1],[1,-1],[-1,1],[1,1]]) # e.g., 1, 8, 49, 56
    if dx_0 < 0:
        out_multipliers[:2] = out_multipliers[:2][::-1]
    if dy_0 < 0 and dy_end < 0:
        out_multipliers[:,1] = -out_multipliers[:,1]
    if dx_end < 0: 
        out_multipliers[2:3] = out_multipliers[2:3][::-1]

    return out_multipliers
    
    
    

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
    