# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 09:50:23 2017

@author: kbefus
"""
from __future__ import print_function
import numpy as np
import time, os

from pyres import pyres_utils,mesh_utils


class meshR3(object):
    def __init__(self,geom_dict=None,name=None,nelectrodes_per_line=None,
                 nlines=None):
        self.name = name
        self.nelectrodes_per_line = nelectrodes_per_line
        self.nlines = nlines


    def make_points(self, line_strexyz=None, clen=0.1,start_pts=1):
        '''Make point dictionary.'''
        
        if line_strexyz is None:
            # Make line data using default line settings
            self.make_strxyz()
        else:
            self.line_strexyz = line_strexyz
        
        if not hasattr(self,'points'):
            self.points = {}
        
        mesh_utils.make_points(self,clen,start_pts,ndim=3)  
        

    def make_lines(self,start_lines=1):
        
        mesh_utils.make_lines(self,start_lines)
        
    def make_boundaries(self,nforeground=4,nbackground=4):
        '''Make boundaries from survey bounds.
        
        Can eventually add capability to define non-square boundaries
        '''
        
        mesh_utils.make_boundaries(self,nforeground,nbackground)
    
    def make_surfaces(self):
        
        mesh_utils.make_surfaces(self,ndim=3)
        
    def make_regions(self, region_xyzpts=None, extend_to=None,
                    clen=None, boundary_dict=None,region_fname=None,load_dict=None,
                    outside_foreground=False,active_domain='foreground'):
        '''Create geometry components for assigning model region.'''
        
        
    def write_geo(self,out_fname=None,boundary_dict=None,extrude_dict=None,ndigits=2):
       '''Write .geo file.
       
       '''
       gmsh_geo_dict = {'gmsh_obj':self,'out_fname':out_fname,
                        'boundary_dict':boundary_dict,
                        'extrude_dict':extrude_dict,
                        'ndigits':ndigits}
       write_gmsh_geo(**gmsh_geo_dict)
       self.geo_name = out_fname
    
    def run_gmsh(self,silent=True,report=True,
                 gmshdir = r'D:\Research\Software\gmsh-2.16.0-Windows',
                 gmshexe = 'gmsh.exe',optimize=False, num_quad_lloyd_steps=0, dim=3,
                 normal_msg = 'stopped on'):
        '''Run gmsh executable.
        
        Source: after https://github.com/nschloe/pygmsh/blob/master/pygmsh/helper.py#generate_mesh
        '''
        cwd = os.path.dirname(self.geo_name)
        
        cmd_list = [os.path.join(gmshdir,gmshexe), '-{0:d}'.format(dim), self.geo_name, '-o', self.msh_name]
        
        if optimize:
            cmd_list += ['-optimize']
        
        if num_quad_lloyd_steps > 0:
            cmd_list += ['-optimize_lloyd', str(num_quad_lloyd_steps)]
        
        run_dict = {'cmd_list':cmd_list,'cwd':cwd,'silent':silent,
                    'normal_msg':normal_msg,'report':report}
        self.success,self.buff = pyres_utils.run_cmd(**run_dict)
            
    def locate_inf_electrode(self,inf_electrode_xy=[None,None],
                             msh_name=None):
        '''Find mesh node of inf electrode(s).'''
        
        if hasattr(self,'msh_name'):
            # Load .msh file
            msh_dict = read_gmsh(msh_fname=self.msh_name)
        elif msh_name is not None:
            msh_dict = read_gmsh(msh_fname=msh_name)
        else:
            print("msh_name needs to be assigned.")
        # Find closest node to infinity electrode
        
        # 1) Try to find perfect match
        match_inds = pyres_utils.match_xy(XY=msh_dict['nodes'][:,0:2],
                                       xy_to_match = np.array(inf_electrode_xy))
        if len(match_inds)<1.:
            dist_array = np.sqrt((inf_electrode_xy[0]-msh_dict['nodes'][:,0])**2.+\
                                 (inf_electrode_xy[1]-msh_dict['nodes'][:,1])**2.)
            match_ind = np.argmin(dist_array)+1 # account for 1-based indexing
        
        return match_ind,msh_dict['nodes'][match_ind,:]
        
    def msh_to_dat(self,msh_fname=None,out_file='mesh3d.dat', topo_correct=False,
               topo_dict=None, overwrite=True, verbose=False, job_type=1,
               alpha_scale=1.):
        '''Convert .msh to .dat for R3.'''
        if msh_fname is None:
            msh_fname = self.msh_name
            
        gmsh_to_R3_dict= {'msh_fname':msh_fname,'out_file':out_file,
                          'topo_correct':topo_correct,'topo_dict':topo_dict,
                          'overwrite':overwrite,'verbose':verbose,
                          'job_type':job_type,'alpha_scale':alpha_scale}
        self.gmsh_to_R3_msg,self.nregions, self.mesh_dict,self.topo_xyz = gmsh_to_R3(**gmsh_to_R3_dict)
    
    def make_mesh(self,geom_dict=None,write_dict=None,
                  run_gmsh=False, gmsh_dict=None,region_dict=None):
        '''Run individual functions to make mesh.
        '''
        self.make_points(**geom_dict)
        self.make_lines()
        self.make_boundaries()
        self.make_surfaces()
        if region_dict is not None:
            self.make_regions(**region_dict)
            write_dict.update({'region_dict':region_dict})
        self.write_geo(**write_dict)
        
        if run_gmsh:
            self.msh_name = "{}.msh".format(os.path.splitext(self.geo_name)[0])
            if gmsh_dict is None:
                self.run_gmsh()
            else:
                self.run_gmsh(**gmsh_dict)

                
                
class meshR2(object):
    def __init__(self,geom_dict=None,name=None,nelectrodes_per_line=None,
                 nlines=1,electrode_spacing=None):
        self.name = name
        self.nelectrodes_per_line = nelectrodes_per_line
        self.nlines = nlines
        self.electrode_spacing = electrode_spacing


    def make_strxyz(self):
        '''Make 2D line string, electrode number, x, y, z array.
        
        
        Note: Assumes all electrodes at surface (z=0) and are on the 
              y-axis (y=0).
        '''
        self.line_strexyz = np.array([np.hstack([np.ones((self.nelectrodes_per_line,1)),
                                  np.arange(self.nelectrodes_per_line).reshape([-1,1]),
                                  self.electrode_spacing*np.arange(self.nelectrodes_per_line).reshape([-1,1]),
                                  np.zeros((self.nelectrodes_per_line,1)),
                                  np.zeros((self.nelectrodes_per_line,1))])])
        
    def make_points(self, line_strexyz=None, clen=0.1,start_pts=1):
        '''Make point dictionary.'''
        
        if line_strexyz is None:
            # Make line data using default line settings
            self.make_strxyz()
        else:
            self.line_strexyz = line_strexyz
        
        if not hasattr(self,'points'):
            self.points = {}
        
        elec_list = []
        iline_order = []
        for iline in np.arange(self.line_strexyz.shape[0]):    
            if iline not in self.points.keys():
                self.points[iline] = {}
            iline_order.append(iline)
            key_order = []
            for linestr,enum,x,y,z in self.line_strexyz[iline]:
                self.points[iline].update({start_pts:[x,y,z,clen]})
                key_order.append(start_pts)
                elec_list.append([enum,start_pts]) # save data for protocol.dat
                start_pts+=1
                
            self.points[iline].update({'order':key_order})
            
        self.points.update({'order':iline_order})    
        self.count_pts = start_pts        
        self.electrode_array = np.array(elec_list)       

    def make_lines(self,start_lines=1):
        
        self.start_lines=start_lines
        self.lines = {}
        self.surveyline_startend = {}
        self.startpt,self.endpt = [],[]
        key_order_se = []
        for ikey in self.points['order']:
            # Collect consecutive points into lines in reverse order
            line_keys = self.points[ikey]['order']
            
            # record start and end positions of line
            self.startpt.append(line_keys[0])
            self.endpt.append(line_keys[-1])
            
            # Save first line in section
            self.surveyline_startend[ikey] = [self.start_lines]
            key_order_se.append(ikey)
            
            key_order = []
            for istart,iend in zip(line_keys[:-1],line_keys[1:]):
                self.lines.update({self.start_lines:[istart,iend]})
                key_order.append(self.start_lines)
                self.start_lines+=1
                
            self.lines.update({'order':key_order})
            # Save last line in section
            self.surveyline_startend[ikey].append(self.start_lines-1)
        
        self.surveyline_startend.update({'order':key_order_se})
        
        # Connectors for the lines
        if len(self.startpt)>1:
            for istart1,istart2 in zip(self.startpt[:-1],self.startpt[1:]):
                self.lines.update({self.start_lines:[istart1,istart2]})
                self.lines['order'].append(self.start_lines)
                self.start_lines+=1
            
            for end1,iend2 in zip(self.endpt[:-1],self.endpt[1:]):
                self.lines.update({self.start_lines:[end1,iend2]})
                self.lines['order'].append(self.start_lines)
                self.start_lines+=1
        
    def make_boundaries(self,nforeground=4,nbackground=4):
        '''Make boundaries from survey bounds.
        
        
        Can eventually add capability to define non-square boundaries
        '''
        
        self.boundaries = {}
        # First make foreground boundaries by extending some constant beyond
        # the corner points
        fore_corners = np.arange(self.count_pts,self.count_pts+nforeground)
        self.count_pts = self.count_pts+nforeground
        
        fore_pairs = np.c_[fore_corners,np.roll(fore_corners,-1)]
        
        self.boundaries['foreground'] = {}
        self.boundaries['foreground'].update(self.lines.copy())
        key_order = []
        key_order.extend(self.lines['order'])
        self.boundaries['foreground'].update({self.start_lines:[self.endpt[0],fore_pairs[-1,1]]})
        key_order.append(self.start_lines)
        self.start_lines += 1
        self.fore_lines=[]
        for ifore in np.arange(nforeground-1):
            self.boundaries['foreground'].update({self.start_lines:
                                                    fore_pairs[ifore,:].tolist()})
            key_order.append(self.start_lines)
            self.fore_lines.append(self.start_lines)
            self.start_lines += 1
        
        self.boundaries['foreground'].update({self.start_lines:[fore_pairs[-1,0],self.startpt[0]]})
        key_order.append(self.start_lines)
        self.boundaries['foreground'].update({'order':key_order})
        self.start_lines += 1
        
        # Background points and lines
        back_corners = np.arange(self.count_pts,self.count_pts+nbackground)
        self.count_pts = self.count_pts + nbackground
        
        back_pairs = np.c_[back_corners,np.roll(back_corners,-1)]
        key_order = []
        self.boundaries['background'] = {}        
        self.boundaries['background'].update({self.start_lines:[fore_pairs[0,0],back_pairs[0,0]]})
        key_order.append(self.start_lines)
        self.start_lines += 1
            
        for iback in np.arange(nbackground-1):
            self.boundaries['background'].update({self.start_lines:
                                                    back_pairs[iback,:].tolist()})
            key_order.append(self.start_lines)
            self.start_lines += 1
        
        self.boundaries['background'].update({self.start_lines:[back_pairs[-1,0],fore_pairs[-1,0]]})
        key_order.append(self.start_lines)
        self.boundaries['background'].update({'order':key_order})
        self.start_lines += 1
    
    def make_surfaces(self):
        
        self.surfaces = {}   
        seg_dict = {'nelectrodes':self.nelectrodes_per_line,
                    'nlines':self.nlines,
                    'topbot':True,
                    'iline':np.arange(2*(self.nlines-1))}
        self.connecting_segs = mesh_utils.segment_loc(**seg_dict).reshape((2,-1)).T

        self.surflines = {} # will have "Line Loop(#) = {}
        key_order_sl = []
        key_order_s = []
        # Make bounding foreground surface from foreground lines and line-bounding lines
        self.surflines.update({self.start_lines:np.roll(self.boundaries['foreground']['order'],1).tolist()}) # foreground boundary
        key_order_sl.append(self.start_lines)
        self.start_lines += 1
        self.surfaces.update({self.start_lines:[self.start_lines-1]})
        key_order_s.append(self.start_lines)
        self.start_lines += 1
        
        # Make background surface between foreground and background lines
        foreline_order = -1*np.array(self.fore_lines)[::-1]
        background_order = np.hstack([foreline_order,self.boundaries['background']['order']])
        self.surflines.update({self.start_lines:background_order.tolist()}) # background line
        key_order_sl.append(self.start_lines)
        self.surflines.update({'order':key_order_sl})
        self.start_lines += 1
        
        self.surfaces.update({self.start_lines:[self.start_lines-1]})
        key_order_s.append(self.start_lines)
        self.surfaces.update({'order':key_order_s})

    def make_regions(self, region_xyzpts=None, extend_to=None,
                    clen=None, boundary_dict=None,region_fname=None,load_dict=None,
                    outside_foreground=False,active_domain='foreground'):
        '''Create geometry components for assigning model region.'''
        
        
        # ----------------- Load region data ---------------------------
        if region_fname is not None:
            region_xyzpts=[]
            for region_temp in region_fname:
                region_xyz_temp,header_info = load_delim(region_temp,**load_dict)
                if region_xyz_temp.shape[-1]<3:
                    region_xyz_temp = np.hstack([region_xyz_temp,np.zeros((region_xyz_temp.shape[0],1))])
                
                region_xyzpts.append(region_xyz_temp)
        
        if isinstance(region_xyzpts,(list,tuple)):
            try:
                region_xyzpts = np.array(region_xyzpts)
                

            except:
                print("Could not convert region_xyzpts to np.ndarray...Assuming it is a list of region arrays")
        
        if isinstance(region_xyzpts,np.ndarray):
            # Fill elevation, z, column with zeros if not existing   
            if region_xyzpts.shape[-1]<3:
                region_xyzpts = np.hstack([region_xyzpts,np.zeros((region_xyzpts.shape[0],1))])
            
            if len(region_xyzpts.shape)<3:
                # Stack regions into 3rd dimension
                region_xyzpts = [region_xyzpts]

        self.region_xyzpts = region_xyzpts
        # -------------- Make region geometry components ---------------------  
        self.region_dict = {}
        iregion_count = 0        
        for iregion_xyzpts in self.region_xyzpts:
            if active_domain in ['foreground'] and extend_to in ['background']:
                # First create a region in the foreground that will then be 
                # extended into the background as well
                geom_dict = {'points':self.points.copy(),
                             'boundaries':self.boundaries.copy(),
                             'surflines':self.surflines.copy(),
                             'fore_lines':self.fore_lines}
                             
                # Find xyz points inside foreground
                foreground_bounds,foreground_bool = mesh_utils.calc_region_bounds(region='foreground',
                                                                  geom_dict=geom_dict,
                                                                  xyz = iregion_xyzpts)
                
                ipts_foreground = iregion_xyzpts[foreground_bool,:].copy()
                ipts_background = iregion_xyzpts[~foreground_bool,:].copy()
                iregion_xyzpts = ipts_background.copy() # replace with only pts in background region
                
                make_region_dict = {'geom_dict':geom_dict,
                                    'start_line_num':self.start_lines,
                                    'start_pt_num':self.count_pts,
                                    'extend_to':'foreground',
                                    'region_xyzpts':ipts_foreground,
                                    'boundary_dict':boundary_dict,
                                    'clen':clen,
                                    'active_domain':active_domain}
                region_out_dict,self.start_lines,self.count_pts = mesh_utils.make_region(**make_region_dict)
                self.region_dict.update({iregion_count:region_out_dict['region_dict']})
                self.boundaries = region_out_dict['boundaries'].copy()
                self.surflines = region_out_dict['surflines'].copy()
                
                active_domain = 'background' # overwrite to make background region
                
                # Update region order
                if 'order' not in self.region_dict.keys():
                    self.region_dict.update({'order':[iregion_count]})
                else:
                    self.region_dict['order'].append(iregion_count)
            
            # Create region in "active_domain"
            geom_dict = {'points':self.points.copy(),
                             'boundaries':self.boundaries.copy(),
                             'surflines':self.surflines.copy(),
                             'fore_lines':self.fore_lines}
                
            make_region_dict = {'geom_dict':geom_dict,
                                'start_line_num':self.start_lines,
                                'start_pt_num':self.count_pts,
                                'extend_to':extend_to,
                                'region_xyzpts':iregion_xyzpts,
                                'boundary_dict':boundary_dict,
                                'clen':clen,
                                'active_domain':active_domain}
            region_out_dict,self.start_lines,self.count_pts = mesh_utils.make_region(**make_region_dict)
            self.region_dict.update({iregion_count:region_out_dict['region_dict']})
            self.boundaries = region_out_dict['boundaries'].copy()
            self.surflines = region_out_dict['surflines'].copy()
            
            # Update region order
            if 'order' not in self.region_dict.keys():
                self.region_dict.update({'order':[iregion_count]})
            else:
                self.region_dict['order'].append(iregion_count)
                
        
    def write_geo(self,out_fname=None,boundary_dict=None,extrude_dict=None,ndigits=2,
                  region_dict=None):
       '''Write .geo file.
       
       '''
       gmsh_geo_dict = {'gmsh_obj':self,'out_fname':out_fname,
                        'boundary_dict':boundary_dict,
                        'extrude_dict':extrude_dict,
                        'ndigits':ndigits,'mesh_dim':2,
                        'region_dict':region_dict}
       write_gmsh_geo(**gmsh_geo_dict)
       self.geo_name = out_fname
    
    def run_gmsh(self,silent=True,report=True,
                 gmshdir = r'D:\Research\Software\gmsh-2.16.0-Windows',
                 gmshexe = 'gmsh.exe',optimize=False, num_quad_lloyd_steps=0, dim=3,
                 normal_msg = 'stopped on'):
        '''Run gmsh executable.
        
        Source: after https://github.com/nschloe/pygmsh/blob/master/pygmsh/helper.py#generate_mesh
        '''
        cwd = os.path.dirname(self.geo_name)
        
        cmd_list = [os.path.join(gmshdir,gmshexe), '-{0:d}'.format(dim), self.geo_name, '-o', self.msh_name]
        
        if optimize:
            cmd_list += ['-optimize']
        
        if num_quad_lloyd_steps > 0:
            cmd_list += ['-optimize_lloyd', str(num_quad_lloyd_steps)]
        
        run_dict = {'cmd_list':cmd_list,'cwd':cwd,'silent':silent,
                    'normal_msg':normal_msg,'report':report}
        self.success,self.buff = pyres_utils.run_cmd(**run_dict)
            
    def locate_inf_electrode(self,inf_electrode_xy=[None,None],
                             msh_name=None):
        '''Find mesh node of inf electrode(s).'''
        
        if hasattr(self,'msh_name'):
            # Load .msh file
            msh_dict = read_gmsh(msh_fname=self.msh_name)
        elif msh_name is not None:
            msh_dict = read_gmsh(msh_fname=msh_name)
        else:
            print("msh_name needs to be assigned.")
        # Find closest node to infinity electrode
        
        # 1) Try to find perfect match
        match_inds = pyres_utils.match_xy(XY=msh_dict['nodes'][:,0:2],
                                       xy_to_match = np.array(inf_electrode_xy))
        if len(match_inds)<1.:
            dist_array = np.sqrt((inf_electrode_xy[0]-msh_dict['nodes'][:,0])**2.+\
                                 (inf_electrode_xy[1]-msh_dict['nodes'][:,1])**2.)
            match_ind = np.argmin(dist_array)+1 # account for 1-based indexing
        
        return match_ind,msh_dict['nodes'][match_ind,:]
        
    def msh_to_dat(self,msh_fname=None,out_file='mesh.dat', topo_correct=False,
               topo_dict=None, overwrite=True, verbose=False, job_type=1,
               alpha_scale=1.,node_zones=None,param_zones=None):
        '''Convert .msh to .dat for R2.'''
        if msh_fname is None:
            msh_fname = self.msh_name

        gmsh_to_R2_dict= {'msh_fname':msh_fname,'out_file':out_file,
                          'topo_correct':topo_correct,'topo_dict':topo_dict,
                          'overwrite':overwrite,'verbose':verbose,
                          'alpha_scale':alpha_scale,
                          'node_zones':node_zones,'param_zones':param_zones}
        self.gmsh_to_R2_msg,self.nregions, self.mesh_dict,self.topo_xyz = gmsh_to_R2(**gmsh_to_R2_dict)
    
    def make_mesh(self,geom_dict=None,write_dict=None,
                  run_gmsh=False, gmsh_dict=None,
                  region_dict=None):
        '''Run individual functions to make mesh.
        '''
        self.make_points(**geom_dict)
        self.make_lines()
        self.make_boundaries()
        self.make_surfaces()
        if region_dict is not None:
            self.make_regions(**region_dict)
            write_dict.update({'region_dict':region_dict})
            
        self.write_geo(**write_dict)
        
        if run_gmsh:
            self.msh_name = "{}.msh".format(os.path.splitext(self.geo_name)[0])
            if gmsh_dict is None:
                self.run_gmsh()
            else:
                self.run_gmsh(**gmsh_dict)
#%% ------------------ Utilities ------------------            
            
def write_gmsh_geo(gmsh_obj=None,out_fname=None,boundary_dict=None,
                   extrude_dict=None,ndigits=2,mesh_dim=3,region_dict=None):
    '''Write gmsh .geo file.'''
    if extrude_dict is None and mesh_dim==3:
        mesh_dim = 2
    
    # Keep track of written geometries
    all_pts = []
    all_lines = []
    foreground_pt_fmt = "Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }}\n"
    out_fore_mult = mesh_utils.find_line_orientation(gmsh_obj)
    
    with open(out_fname,'w') as f:
        f.write("// Gmsh project created using Python on {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()))) 
        # Write point data by survey row
        for i in gmsh_obj.points['order']:
            f.write("// Points for row {}\n".format(i+1))
            for ipt,lineptkey in enumerate(gmsh_obj.points[i]['order']):
                x,y,z,clen = gmsh_obj.points[i][lineptkey]
                f.write("Point({0}) = {{{1}, {2}, {3}, {4}}};\n".format(lineptkey,x,y,z,clen)) # reverse x,y to match R3 example
                all_pts.append(lineptkey)
            f.write("\n") # blank space
        
        if mesh_dim == 2:
            # Write Translate command for external foreground boundaries
            f.write(foreground_pt_fmt.format(boundary_dict['fdxdy'][0],0.,0,np.max(gmsh_obj.endpt)))
            f.write(foreground_pt_fmt.format(boundary_dict['fdxdy'][0],-boundary_dict['fdxdy'][1],0,np.max(gmsh_obj.endpt)))
            f.write(foreground_pt_fmt.format(-boundary_dict['fdxdy'][0],-boundary_dict['fdxdy'][1],0,np.min(gmsh_obj.startpt)))
            f.write(foreground_pt_fmt.format(-boundary_dict['fdxdy'][0],0.,0,np.min(gmsh_obj.startpt)))
            f.write("\n") # blank space
            
            # For far background
            f.write("p1[]=Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }};\n".format(boundary_dict['bdxdy'][0]+boundary_dict['fdxdy'][0],0.,0,np.max(gmsh_obj.endpt)))
            f.write("p3[]=Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }};\n".format(boundary_dict['bdxdy'][0]+boundary_dict['fdxdy'][0],-boundary_dict['bdxdy'][1],0,np.min(gmsh_obj.endpt)))
            f.write("p4[]=Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }};\n".format(-boundary_dict['bdxdy'][0]-boundary_dict['fdxdy'][0],-boundary_dict['bdxdy'][1],0,np.max(gmsh_obj.startpt)))
            f.write("p2[]=Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }};\n".format(-boundary_dict['bdxdy'][0]-boundary_dict['fdxdy'][0],0.,0,np.min(gmsh_obj.startpt)))
            f.write("\n") # blank space
        else:
            
            # Write Translate command for external boundaries, foreground first
            f.write(foreground_pt_fmt.format(out_fore_mult[2][0]*boundary_dict['fdxdy'][0],out_fore_mult[2][1]*boundary_dict['fdxdy'][1],0,np.max(gmsh_obj.startpt)))
            f.write(foreground_pt_fmt.format(out_fore_mult[0][0]*boundary_dict['fdxdy'][0],out_fore_mult[0][1]*boundary_dict['fdxdy'][1],0,np.min(gmsh_obj.startpt)))
            f.write(foreground_pt_fmt.format(out_fore_mult[1][0]*boundary_dict['fdxdy'][0],out_fore_mult[1][1]*boundary_dict['fdxdy'][1],0,np.min(gmsh_obj.endpt)))
            f.write(foreground_pt_fmt.format(out_fore_mult[3][0]*boundary_dict['fdxdy'][0],out_fore_mult[3][1]*boundary_dict['fdxdy'][1],0,np.max(gmsh_obj.endpt)))
            f.write("\n") # blank space
            
            # For far background
            f.write("p1[]=Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }};\n".format(boundary_dict['bdxdy'][0],-boundary_dict['bdxdy'][1],0,np.max(gmsh_obj.startpt)))
            f.write("p2[]=Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }};\n".format(-boundary_dict['bdxdy'][0],-boundary_dict['bdxdy'][1],0,np.min(gmsh_obj.startpt)))
            f.write("p3[]=Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }};\n".format(-boundary_dict['bdxdy'][0],boundary_dict['bdxdy'][1],0,np.min(gmsh_obj.endpt)))
            f.write("p4[]=Translate {{{0},{1},{2}}} {{ Duplicata {{ Point{{{3}}}; }} }};\n".format(boundary_dict['bdxdy'][0],boundary_dict['bdxdy'][1],0,np.max(gmsh_obj.endpt)))
            f.write("\n") # blank space
        
        # Set characteristic length for background elements
        if isinstance(clen,str):
            f.write("Characteristic Length {{p1[0], p2[0], p3[0], p4[0]}} = {} * {};\n".format(boundary_dict['bclen'],clen)) # assumes clen used above in this function
        else:
            f.write("Characteristic Length {{p1[0], p2[0], p3[0], p4[0]}} = {};\n".format(boundary_dict['bclen']*clen))
        f.write("\n") # blank space
        
        if hasattr(gmsh_obj,'region_dict'):
            if gmsh_obj.region_dict is not None:
                # Write region geometry
                f.write("// Region points \n")
                for i in gmsh_obj.region_dict['order']:
                    # Write region points
                    f.write("// Points for region {}\n".format(i+1))
                    for ipt,lineptkey in enumerate(gmsh_obj.region_dict[i]['region_points']['order']):
                        x,y,z,clen = gmsh_obj.region_dict[i]['region_points'][lineptkey]
                        f.write("Point({0}) = {{{1}, {2}, {3}, {4}}};\n".format(lineptkey,x,y,z,clen)) # must reverse x,y to match R3 example
                        all_pts.append(lineptkey)
                    f.write("\n")
        
        # Write foreground boundary lines
        if mesh_dim == 2:
            # Foreground connects into electrode lines
            foreground_lines = gmsh_obj.boundaries['foreground']['order']
            
            f.write("// Foreground boundary lines\n")
            for i in foreground_lines:
                if i not in all_lines:
                    f.write("Line({0}) = {{{1}}};\n".format(i,', '.join([str(itemp) for itemp in gmsh_obj.boundaries['foreground'][i]])))
                    all_lines.append(i)
        else:
            
            # Write Line data by row
            for iline in gmsh_obj.surveyline_startend['order']:
                f.write("// Line {} of electrodes\n".format(iline+1))
                line_nums = np.arange(gmsh_obj.surveyline_startend[iline][0],gmsh_obj.surveyline_startend[iline][1]+1)
                for line_num in line_nums:
                    if line_num not in all_lines:
                        f.write("Line({0}) = {{{1}}};\n".format(line_num,', '.join([str(itemp) for itemp in gmsh_obj.lines[line_num]])))
                        all_lines.append(line_num)
                f.write("\n") # blank space
            
            # Joining lines between electrode lines
            join_lines = np.sort(gmsh_obj.connecting_segs.ravel())
            f.write("// Lines connecting lines of electrodes\n")
            for join_line in join_lines:
                if join_line not in all_lines:
                    f.write("Line({0}) = {{{1}}};\n".format(join_line,', '.join([str(itemp) for itemp in gmsh_obj.lines[join_line]])))
                    all_lines.append(join_line)
            f.write("\n") # blank space
            
            f.write("// Foreground boundary lines\n")
            for i in gmsh_obj.boundaries['foreground']['order']:
                if i not in all_lines:
                    f.write("Line({0}) = {{{1}}};\n".format(i,', '.join([str(itemp) for itemp in gmsh_obj.boundaries['foreground'][i]])))
                    all_lines.append(i)
        f.write("\n") # blank space
        
        # Write background boundary lines
        f.write("// Background boundary lines\n")
        for i in gmsh_obj.boundaries['background']['order']:
            if i not in all_lines:
                f.write("Line({0}) = {{{1}}};\n".format(i,', '.join([str(itemp) for itemp in gmsh_obj.boundaries['background'][i]])))
                all_lines.append(i)
        f.write("\n") # blank space
        
        
        nregion_surfs = 0
        surf_keys_ordered = gmsh_obj.surfaces['order']
        if hasattr(gmsh_obj,'region_dict'):
            if gmsh_obj.region_dict is not None:
                # Write region geometry
                f.write("// Region data \n")
                region_phys_surfline_list = []
                region_back_phys_surfline_list = []
                region_fore_phys_surfline_list = []
                region_phys_surf_list = []
                for i in gmsh_obj.region_dict['order']:
                    
                    # Write region lines
                    f.write("// Lines for region {}\n".format(i))
                    line_nums = gmsh_obj.region_dict[i]['region_lines']['order']
                    for line_num in line_nums:
                        if line_num not in all_lines:
                            f.write("Line({0}) = {{{1}}};\n".format(line_num,', '.join([str(itemp) for itemp in gmsh_obj.region_dict[i]['region_lines'][line_num]])))
                            all_lines.append(line_num)
                    f.write("\n") # blank space
                    
                    # Write region surface
                    for surfline_key in gmsh_obj.region_dict[i]['region_surflines']['order']:
                        f.write("Line Loop({0}) = {{{1}}};\n".format(surfline_key,', '.join([str(itemp) for itemp in gmsh_obj.region_dict[i]['region_surflines'][surfline_key]])))
                        if gmsh_obj.region_dict[i]['active_domain'] in ['foreground']:
                            region_fore_phys_surfline_list.append(surfline_key)
                        else:
                            region_back_phys_surfline_list.append(surfline_key)
                        
                    for surf_key in gmsh_obj.region_dict[i]['region_surfaces']['order']:
                        f.write("Plane Surface({0}) = {{{1}}};\n".format(surf_key,', '.join([str(itemp) for itemp in gmsh_obj.region_dict[i]['region_surfaces'][surf_key]])))
                        region_phys_surf_list.append(surf_key)
                        
                    f.write("\n") # blank space
                # Assign the physical surface for the region 
                # Assumes only foreground and background assigned earlier
                if len(region_back_phys_surfline_list)>0:
                    gmsh_obj.surfaces[surf_keys_ordered[1]].append(region_back_phys_surfline_list[0]) # add to end of background surface
                
                if len(region_fore_phys_surfline_list)>0:
                    gmsh_obj.surfaces[surf_keys_ordered[0]].append(region_fore_phys_surfline_list[0]) # add to end of foreground surface
                
                for ireg,reg_surf in enumerate(region_phys_surf_list):
                    nregion_surfs +=1
                    f.write("Physical Surface({0}) = {{{1}}};\n".format(nregion_surfs,reg_surf))
        
        
        # Write the lines needed to make surfaces
        f.write("// Foreground and Background surface lines and surfaces\n")
        all_surf_keys = np.sort(np.hstack([gmsh_obj.surflines['order'],gmsh_obj.surfaces['order']]))
        for i in all_surf_keys:
            if i in gmsh_obj.surflines['order']:
                f.write("Line Loop({0}) = {{{1}}};\n".format(i,', '.join([str(itemp) for itemp in gmsh_obj.surflines[i]])))
            elif i in gmsh_obj.surfaces['order']:
                f.write("Plane Surface({0}) = {{{1}}};\n".format(i,', '.join([str(itemp) for itemp in gmsh_obj.surfaces[i]])))
                f.write("\n") # blank space
                
        f.write("\n") # blank space
        
        if mesh_dim == 2:                    
            f.write("\n")
            
            for icount,isurf in enumerate(gmsh_obj.surfaces['order']):
                f.write("Physical Surface({0}) = {{{1}}};\n".format(icount+nregion_surfs+1,isurf))
            
        else:
            # Extrude foreground surfaces
            f.write("// Extrude foreground surfaces\n")
            prep_nlayers = ', '.join([str(round(itemp,ndigits)) for itemp in extrude_dict['nlayer_list']])
            prep_zratios = ', '.join([str(round(itemp,ndigits)) for itemp in extrude_dict['domain_zratios']])
            fg_list = []
            for icount,isurf in enumerate(gmsh_obj.surfaces['order'][:-1]): # all but the background surface
                f.write("fg{0}[] = Extrude {{0,0,{1}}} {{Surface{{{2}}}; Layers{{ {{{3}}}, {{{4}}} }}; Recombine; }};\n".format(\
                        icount,-extrude_dict['max_depth'],
                        isurf,
                        prep_nlayers,
                        prep_zratios))
                fg_list.append("fg{0}[1]".format(icount))
                
            f.write("\n") # blank space
            
            # Extrude background surface
            f.write("// Extrude background surface\n")
            f.write("bg[] = Extrude {{0,0,{1}}} {{Surface{{{2}}}; Layers{{ {{{3}}}, {{{4}}} }}; Recombine; }};\n".format(\
                    icount,-extrude_dict['max_depth'],
                   gmsh_obj.surfaces['order'][-1],
                    prep_nlayers,prep_zratios))
            f.write("\n") # blank space
             
            # Write Physical volumes
            f.write("Physical Volume(1) = {{{0}}};\n".format(", ".join(fg_list)))
            f.write("Physical Volume(2) = {{{0}}};\n".format('bg[1]'))
        

        

def add_topo(msh_fname=None,mesh_dict=None,out_file=None,topo_file=None, 
             surf_elev=0.,
             nxy = 100, method='linear', fill_value=0.,
             write_gmsh_bool=False, gmsh_out_fname=None,
             r2_bool=True,topo_load_dict={'cols2read':[0,1,2],
                                           'nheaders':1}):
    '''Add topography to gmsh .msh file.
    
    Adds topography to a gmsh .msh file. Reads the nodes in 'msh_fname' (the gmsh
    3D msh file) and for each elevation 'plane' it perturbs the node
    elevation according to the values given in z. z holds the topography
    values corresponding to x and y and is interpolated across the gmsh mesh. 
    NOTE: The upper surface of the original gmsh .msh must be flat and set to surf_elev. 
    
    Inputs:
        msh_fname: str
            gmsh 3D .msh file with flat upper surface with an elevation of surf_elev.
        
        topo_file: str
            text file containing x,y,z coordinates in three columns.
        
        method: str
            choice of interpolation for 2D interpolation. Currently does not include
            options available in addTopo.m (gridfit.m).
        
        r2_bool: boolean
            defines whether R2 is being used. If False, assumes 3D mesh for R3.
    
    Output: 
        out_file: str
            the file in which new node elevations and co-ordinates are written 
            in ASCII format. These data should replace the 'nodes' section
            of the original gmsh .msh file.
    
    Source: After Matlab addTopo.m function by James Ramm 2012 in R3 package.'''
    from scipy.interpolate import griddata
    
    # Read topography file
    if topo_file is not None and os.path.isfile(topo_file):
        
        topo_data,header_info = load_delim(topo_file,**topo_load_dict)
        
        # Read in mesh
        if mesh_dict is None: # Load .msh, otherwise use input mesh_dict
            mesh_dict = read_gmsh(msh_fname=msh_fname)
        
            
        # Expand nxy
        if isinstance(nxy,(float,int)):
            nx,ny = nxy,nxy
        else:
            nx,ny = nxy
        if r2_bool:
            # -------- Add topography to 2D mesh -------------
            node_xyz = mesh_dict['nodes']
            if topo_data.shape[-1]<3:
                ydata = topo_data[:,1]
            else:
                ydata = topo_data[:,2]
            # Use z in topofile to adjust node y values
            new_y_diff = surf_elev+pyres_utils.extrap(node_xyz[:,0],topo_data[:,0],
                                              ydata,method=method)
            
            node_xyz[:,1] = node_xyz[:,1]+new_y_diff
            mesh_dict['nodes'] = node_xyz
            
        else:
            # -------- Add topography to 3D mesh -------------
            # Find surface coordinates
            surface_xyz = mesh_dict['nodes'][mesh_dict['nodes'][:,2]==surf_elev,:] # y,x,z
            
    
            
            # Create regularly interpolated grid over full domain
    #        xdomain = np.linspace(surface_xyz[:,0].min(),surface_xyz[:,0].max(),nx)
    #        ydomain = np.linspace(surface_xyz[:,1].min(),surface_xyz[:,1].max(),ny)
    #        X,Y = np.meshgrid(xdomain,ydomain)
            
            # Interpolate topographic data to surface nodes, need to flip x,y in 
            # surface_xyz to match R3 example (surface_xyz[:,:2][:,::-1])
            new_surface_z = griddata(topo_data[:,:2],topo_data[:,2],
                            surface_xyz[:,:2],
                            method=method,fill_value=fill_value)
            
            surface_unique_z,layer_inds = np.unique(mesh_dict['nodes'][:,2],return_inverse=True)
            
            # Calculate node-based amount to vertically shift nodes in each layer
            surface_zdiff = new_surface_z - surface_xyz[:,2]
    
            # Adjust nodes in all layers by a constant amount
            # Note: do not need to loop through individual nodes if all layers are 
            # ordered the same (i.e., built from extruding a single plane)
            output_node_array = mesh_dict['nodes'].copy()
            for ilayer, layerz in enumerate(surface_unique_z):
                if ilayer != 0: # Do not change the elevation of nodes in the lowest layer
                    active_layer = layer_inds == ilayer
                    output_node_array[active_layer,2] = output_node_array[active_layer,2] + surface_zdiff
            
            mesh_dict['nodes'] = output_node_array
        
        
        if write_gmsh_bool:
            if gmsh_out_fname is None:
                gmsh_out_fname = msh_fname
            write_gmsh(msh_fname=gmsh_out_fname, mesh_dict=mesh_dict)
        
        return mesh_dict,topo_data     
        

    else:
        print("topo_file not specified or does not exist: {}\n".format(topo_file))
        
def read_gmsh(msh_fname=None):
    '''Read gmsh .msh file.'''
    # Initialize values of switches
    inFormat, inNodes, inElements = 0, 0, 0
    
    # Read gmsh .msh file
    with open(msh_fname,'r') as fid:
        for line in fid:
            # Check for section headers to determine action
            if line[0] == '$':
                if line.find('Nodes') > 0: inNodes = 1
                if line.find('EndNodes') > 0: inNodes = 0
                if line.find('Elements') > 0: inElements = 1
                if line.find('EndElements') > 0: inElements = 0
                if line.find('MeshFormat') > 0: inFormat = 1
                if line.find('EndMeshFormat') > 0: inFormat = 0
            else: # Not a section header
                if inNodes == 1: # Node section
                    if len(line.split()) == 1: # No spatial information, initiate
                        nodes = []
                    else: # Store spatial information
                        entry = list(map(float, line.split()))[1:] # Convert to float
                        nodes.append(entry)
                elif inElements == 1: # Element section
                    if len(line.split()) == 1: # No spatial information, initiate
                        prisms, tets, zones = [], [], []
                        triags = []
                        elem_info = []
                    else:
                        entry = list(map(int, line.split()))[1:]
                        elem_info.append(entry[:4])
                        # mapping physical region to zone number
                        if entry[0] == 6:
                            prisms.append(entry[-6:] + [entry[2]])
                            zones.append(entry[2])
                        elif entry[0] == 4:
                            tets.append(entry[-4:] + [entry[2]])
                            zones.append(entry[2])
                        elif entry[0] == 2:
                            triags.append(entry[-3:] + [entry[2]])
                            zones.append(entry[2])
                elif inFormat == 1: # Format information (for re-writing .msh)
                    mesh_format = line # save as read
        
    # Check and assign mesh type
    if len(prisms) > len(tets):
        npere = 6
        mesh_type = 'Triangular prisms'
        elements = prisms
    elif len(tets) > len(triags):
        npere = 4
        mesh_type = 'Tetrahedra'
        elements = tets  
    else:
        npere = 2
        mesh_type = 'Triangular'
        elements = triags
    
    mesh_dict = {'npere':npere,'mesh_type':mesh_type,'nodes':np.array(nodes),
                 'elements':np.array(elements),'elem_info':np.array(elem_info),
                 'zones':zones,'mesh_format':mesh_format}
                 
    return mesh_dict
    
def write_gmsh(msh_fname=None,mesh_dict=None):
    '''Write gmsh .msh file.'''
    
    with open(msh_fname,'w') as mfid:
        
        # Write mesh format
        mfid.write('$MeshFormat\n')
        mfid.write("{}".format(mesh_dict['mesh_format']))
        mfid.write('$EndMeshFormat\n')
        
        # Write node data
        mfid.write('$Nodes\n')
        mfid.write('{0:d}\n'.format(mesh_dict['nodes'].shape[0]))
        for inode, node_row in enumerate(mesh_dict['nodes'].tolist()):
            mfid.write('{0:d} {1:f} {2:f} {3:f}\n'.format(inode+1,node_row[0],node_row[1],node_row[2]))
        mfid.write('$EndNodes\n')  
        
        # Write element data
        mfid.write('$Elements\n')
        mfid.write('{0:d}\n'.format(mesh_dict['elements'].shape[0]))
        for ielem,(elem_row,info_row) in enumerate(zip(mesh_dict['elements'].astype(str),mesh_dict['elem_info'].astype(str))):

            # Number of entities, dimension, slave tag, master tag, slave node, master node
            mfid.write('{0} {1} {2}\n'.format(ielem+1," ".join(info_row)," ".join(elem_row[:-1])))

        mfid.write('$EndElements\n')  
        
        
        
def gmsh_to_R3(msh_fname=None,out_file='mesh3d.dat', topo_correct=False,
               topo_dict=None, overwrite=True, verbose=False,alpha_scale=1,
               job_type=1):
    '''Convert .msh to R3 mesh3d.dat.
    
    Source: After Gmsh2R3t.py by Florian Wagner 2012 from R3t package'''
    
    if topo_correct:
        mesh_dict,topo_data = add_topo(msh_fname=msh_fname,**topo_dict)
    else:
        mesh_dict = read_gmsh(msh_fname=msh_fname)
        topo_data = None
    
    # Write new mesh3d.dat file
    if os.path.dirname(out_file) in  ['']:
        # Save to same dirictory as .msh file
        output_file = os.path.join(os.path.dirname(msh_fname),out_file)
    else:
        output_file = out_file # Directory already specified in output
    
    exist_flag=0
    if os.path.isfile(output_file):
        print("Warning: {} is already created. Default is to overwrite.".format(output_file))
        exist_flag = 1
        
    unq_zones = np.unique(mesh_dict['elements'][:,-1])
    
    if exist_flag == 0 or overwrite:
        with open(output_file,'w') as wfid:
            # Define mesh components. Create one dirichlet node, datum = 0
            wfid.write('%d %d 1 0 %d \n' % (len(mesh_dict['elements']), len(mesh_dict['nodes']), mesh_dict['npere']))
            
            if isinstance(mesh_dict['elements'],np.ndarray):
                
                mesh_dict['elements'] = mesh_dict['elements'].tolist()
            
            for i, elem in enumerate(mesh_dict['elements']):
                wfid.write('%d ' % (i + 1))
                elem.insert(-1, i+1) # unique parameter number
                for entry in elem:
                    wfid.write('%d ' % entry)
                wfid.write('\n')
            
            wfid.write('\n')
            for i, node in enumerate(mesh_dict['nodes']):
                node_entry = list(node)
                node_entry.insert(0, i+1)
                wfid.write('%d %10.3f %10.3f %10.3f \n' % tuple(node_entry) )
            wfid.write('1\n')  # dirichlet node
    
            # Define smoothing scale for each mesh zone
            if job_type == 1 or job_type in ['inv','inverse','Inv','Inverse']:
                if len(unq_zones)> 1:
                    if isinstance(alpha_scale,(int,float)):
                        alpha_scale = alpha_scale*np.ones_like(unq_zones)
                    
                    for izone,alpha_scale_temp in enumerate(alpha_scale):
                        wfid.write('{0:d} {1:f}\n'.format(izone+1,alpha_scale_temp))
            
    # Show mesh information
    message = """
    %s has been successfully created.
    
    Mesh Info:
    
    Type of elements: %s
    Nodes: %d
    Elements: %d
    Zones: %d
    """ % (os.path.basename(output_file), mesh_dict['mesh_type'],
           len(mesh_dict['elements']), len(mesh_dict['nodes']),
           len(set(mesh_dict['zones'])))            
    
    if verbose:
        print(message)
        
    return message,len(unq_zones),mesh_dict,topo_data

def gmsh_to_R2(msh_fname=None,out_file='mesh.dat', topo_correct=False,
               topo_dict=None, overwrite=True, verbose=False,alpha_scale=1,
               node_zones=None,param_zones=None,default_zone=1):
    '''Convert .msh to R2 mesh.dat.
    
    Source: After Gmsh2R3t.py by Florian Wagner 2012 from R3t package'''
    
    if topo_correct:
        mesh_dict,topo_data = add_topo(msh_fname=msh_fname,**topo_dict)
    else:
        mesh_dict = read_gmsh(msh_fname=msh_fname)
        topo_data = []
    
    if node_zones is None and len(np.unique(mesh_dict['zones']))==1:
        node_zones = {}
        for icount,elem_num in enumerate(mesh_dict['elements']):
            node_zones.update({icount:1})
    else:
        node_zones = {}
        for icount,(elem_num,elem_zone) in enumerate(zip(mesh_dict['elements'],mesh_dict['zones'])):
            node_zones.update({icount:elem_zone})
        
    # Write new mesh3d.dat file
    if os.path.dirname(out_file) in  ['']:
        # Save to same dirictory as .msh file
        output_file = os.path.join(os.path.dirname(msh_fname),out_file)
    else:
        output_file = out_file # Directory already specified in output
    
    exist_flag=0
    if os.path.isfile(output_file):
        print("Warning: {} is already created. Default is to overwrite.".format(output_file))
        exist_flag = 1
        
    all_zones = []
    
    if exist_flag == 0 or overwrite:
        with open(output_file,'w') as wfid:
            # Define mesh components. Create one dirichlet node, datum = 0
            wfid.write('%d %d \n' % (len(mesh_dict['elements']), len(mesh_dict['nodes'])))
            wfid.write('\n')
            if isinstance(mesh_dict['elements'],np.ndarray):
                
                mesh_dict['elements'] = mesh_dict['elements'].tolist()
            icount=1
            for i, elem in enumerate(mesh_dict['elements']):
                wfid.write('{0:10.0f}'.format(i + 1))
                if param_zones is not None:
                    if node_zones[i] in param_zones.keys():
                        if param_zones[node_zones[i]]['param']: # if it is true
                            icount-=1 # keep icount the same if still in the zone

                elem.insert(-1,icount) # unique parameter number                
                icount+=1
                for entry in elem[:-1]:
                    wfid.write('{0:10.0f}'.format(entry))
                
                temp_node_zone = default_zone
                if node_zones[i] in param_zones.keys():
                    if param_zones[node_zones[i]]['zone']: # if it is true
                        temp_node_zone = node_zones[i]+1 # increase from default if 1
                        
                wfid.write('{0:10.0f}'.format(temp_node_zone))
                    
                all_zones.append(temp_node_zone)
                wfid.write('\n')
            
            wfid.write('\n')
            wfid.write('\n')
            for i, node in enumerate(mesh_dict['nodes']):
                node_entry = list(node)
                node_entry.insert(0, i+1)
                wfid.write('{0:d} {1:10.3f} {2:10.3f} \n'.format(*node_entry[:-1]))
            
    unq_zones = np.unique(all_zones)    
    
    # Show mesh information
    message = """
    %s has been successfully created.
    
    Mesh Info:
    
    Type of elements: %s
    Nodes: %d
    Elements: %d
    Zones: %d
    """ % (os.path.basename(output_file), mesh_dict['mesh_type'],
       len(mesh_dict['elements']), len(mesh_dict['nodes']),
       len(set(mesh_dict['zones']))) 
    if verbose:

        print(message)
        
    return message,len(unq_zones),mesh_dict,topo_data
    
def gridfit(xyz=None,xynodes=None, smoothness=1, interp='triangle',
            regularize_type='gradient', solver='vector', max_iter=None,
            extend='warning', tilesize=np.inf, overlap=0.2, mask=None,
            autoscale=True, xyscale=[1,1]):    
    '''Conform surface to regular grid.
    
    
    Source:  After Matlab gridfit.m function by James Ramm 2012 in R3 package.
    '''
    
    # Unpack xyz
    if isinstance(xyz,list):
        x,y,z = xyz # Unpack xyz list
    elif isinstance(xyz,np.ndarray):
        x,y,z = xyz[:,0], xyz[:,1], xyz[:,2]
    
    print("Error: gridfit function in development and incomplete.")
    return x,y,z
    
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
        
    
def load_delim(fname,nheaders=2,delimiter=' ',
                    col_scales = [1,1,1,1],cols2read=[0,1],
                    drop_na=True):
    '''Load delimited text file.'''
    with open(fname,'r') as f:
        header_info = []
        data_list = []
        irow = 0
        for iline in f:
            if irow < nheaders:
                header_info.append(iline)
            elif iline not in ['\n']:
                data_list.append([float(ipiece.strip().strip(' \t')) for ipiece in np.array(iline.split(delimiter))[cols2read]])
            irow += 1
    
    data_array = np.array(data_list)
    
    for icol,iscale in enumerate(col_scales):
        if icol < data_array.shape[-1]:
            if float(iscale) != 0.0: 
                data_array[:,icol] = data_array[:,icol]*float(iscale)            
    
    if drop_na:
        data_array = data_array[~np.isnan(data_array.sum(axis=1))]            
    
    return data_array,header_info
    
    
def insert_key_halfway(orig_dict=None,in_keys=None,
                       new_pts=None,new_keys=None):
    '''Insert key into unordered dictionary.'''   
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
    '''Renumber unordered dictionary.'''
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
        if len(ind_found)> 0:
            old_key = new_dict.pop(other_key,None)
            replaced_keys.append([other_key,main_keys[ind_found[0]]])
            out_key_order.append(main_keys[ind_found[0]])
            new_dict[main_keys[ind_found[0]]]=main_dict[main_keys[ind_found[0]]]
        else:
            out_key_order.append(other_key)
            
    new_dict.update({'order':out_key_order})
    return new_dict, replaced_keys
            
    
    
    