# -*- coding: utf-8 -*-
"""
Created on Tue May 02 10:17:51 2017

@author: kbefus
"""

from __future__ import print_function
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatter 
from matplotlib.image import NonUniformImage
from scipy.ndimage.filters import uniform_filter1d
from scipy.interpolate import griddata
import os

from pyres.pyres_utils import load_inv_output,grid_inv_data,extrap


def calculate_appres(point_dict=None,meas_data=None,line_num=0,
                     verbose_bool=False,cols2d = [0,2,4,6,8,9,10]):
    '''Calculate apparent resistivity.'''
    
    geom_factor_func = lambda Ax,Bx,Mx,Nx: (2.*np.pi) / \
                                            (((1./((Mx-Ax))) - \
                                              (1./((Nx-Ax))))-\
                                             ((1./((Mx-Bx))) - \
                                              (1./((Nx-Bx)))))

    point_dict = point_dict[line_num]
    point_data=[]
    for point_entry in np.sort(point_dict.keys()):
        point_data.append(point_dict[point_entry])
    
    point_array = np.array(point_data) # x, y, z, clen
    
    
    meas_temp = meas_data.copy()  
    if not meas_temp.shape[1] <= len(cols2d):                          
       meas_temp = meas_temp[:,cols2d] # num, m, n, a, b, V/I, error
    
    abmn_array = meas_temp[:,1:5]
    abmn_inds = abmn_array.astype(int).ravel()-1
    abmn_x = point_array[abmn_inds,0].reshape(abmn_array.shape)
    geom_factors=[]
    for abmn in abmn_x:
        geom_factors.append(geom_factor_func(*abmn))
    
    geom_factors = np.array(geom_factors)
    
    app_res = meas_temp[:,-2]*geom_factors
    
    if verbose_bool:
        from scipy.stats import gmean
        print("Geometric mean apparent resistivity = {0:10.2f} ohm-m".format(gmean(np.abs(app_res))))
                                            
    return app_res, geom_factors
    

def calc_pseudo_locs(meas_data=None):
    '''Calculate pseudo-section data locations.'''
    # Construct x locations for points using dipole centers
    x1 = np.abs(meas_data[:,2]-meas_data[:,1])/2.+np.min(meas_data[:,1:3],axis=1)
    x2 = np.abs(meas_data[:,3]-meas_data[:,4])/2.+np.min(meas_data[:,3:5],axis=1)
    xloc = x2+np.abs(x1-x2)/2.
    
    # Construct y locations based on dipole spacing
    yloc = np.abs(x1-x2)/2.  
    
    yloc_dict = {}
    for iunq,yuniq in enumerate(np.unique(yloc)):
        yloc_dict.update({yuniq:(yloc==yuniq).nonzero()[0]})
    
    return [xloc,yloc],yloc_dict

def plot_raw(meas_data=None, norm_y=True,
             v_i_ratio_col=5, xybuff=2., 
             median_thresh=.1, neighbor_thresh=0.2,
             npt_neighbor=7, plot_bool=True,
             plot_grid=True):
    '''Plot raw resistance data by electrode spacing.'''

    [xloc,yloc],yloc_dict = calc_pseudo_locs(meas_data=meas_data)
    
    if plot_bool:
        fig,ax = plt.subplots()
        if plot_grid:
            for ytemp in yloc:
                ax.plot(xloc,ytemp*np.ones_like(xloc),'-',color='lightgrey',
                        linewidth=0.5)
    else:
        ax = None
                
    bad_inds_dict = {}
    all_bad_inds = []
    for key_val in yloc_dict.keys():
        
        xvals = xloc[yloc_dict[key_val]]
        orig_yvals = meas_data[yloc_dict[key_val],v_i_ratio_col]
        if norm_y:
            median_val = np.median(orig_yvals)
            min_max_dif = np.max(orig_yvals) - \
                            np.min(orig_yvals)
            yvals = key_val + (orig_yvals-median_val)/min_max_dif
        else:
            yvals = key_val + meas_data[yloc_dict[key_val],v_i_ratio_col]
        
        if plot_bool:
            ax.plot(xvals,yvals,'o-',label="edepth = {}".format(key_val))
        
        # Error checking
        if median_thresh is not None and neighbor_thresh is not None:
            median_bad_bool = np.abs(yvals-key_val) > median_thresh
            
            filt_yvals = uniform_filter1d(orig_yvals,size=npt_neighbor)
            filt_ratio = (orig_yvals-filt_yvals)/filt_yvals
            npt_smear = np.ceil(npt_neighbor/2.).astype(int)
            filt_bad_bool = np.abs(filt_ratio)>neighbor_thresh#[npt_smear:-npt_smear]
            bad_inds = (median_bad_bool & filt_bad_bool).nonzero()[0]
            bad_inds_dict.update({key_val:yloc_dict[key_val][bad_inds]})
            all_bad_inds.extend(yloc_dict[key_val][bad_inds])
            if plot_bool:
                ax.plot(xvals[bad_inds],yvals[bad_inds],'kx',markersize=30)            
        if plot_bool:    
            ax.set_xlim([xloc.min()-xybuff,xloc.max()+xybuff])
            ax.set_ylim([yloc.max()+xybuff,yloc.min()-xybuff])
        
    return [xloc,yloc],yloc_dict,bad_inds_dict,all_bad_inds,ax
    
def plot_pseudo(meas_data=None,plt_options=None,
                data_col=5, log_bool=False,
                xybuff=2., app_res_bool=True,
                point_dict=None):
    '''Plot pseudo-section.'''
    [xloc,yloc],_ = calc_pseudo_locs(meas_data=meas_data)
    
    if app_res_bool and point_dict is not None:
        
        app_res,gf = calculate_appres(point_dict=point_dict,meas_data=meas_data)
        plot_data=app_res
        color_label = 'Apparent resistivity [Ohm-m]'
    else:
        plot_data = meas_data[:,data_col]
        color_label = r'V/I [Ohm]'
    
    fig,ax = plt.subplots()
    if log_bool:
        s1 = ax.scatter(xloc,yloc,c=np.log10(np.abs(plot_data)),edgecolors='none',
                    s=40)
    else:
        s1 = ax.scatter(xloc,yloc,c=plot_data,edgecolors='none',
                        s=40)
    cbar = plt.colorbar(s1,ax=ax)
    cbar.ax.set_ylabel(color_label)
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_ylabel('Pseudo depth [n electrode spacings]')
    ax.set_xlabel('Electrode number')
    ax.set_xlim([xloc.min()-xybuff,xloc.max()+xybuff])
    ax.set_ylim([yloc.max()+xybuff,yloc.min()-xybuff])

def plot_quad_mesh(xxyy_mesh=None,invert_y=True,
                   label_nodes=False,label_buffer=.05,
                   nskip_label=2):
    '''Plot quadrilateral mesh.'''
    xx_mesh,yy_mesh=xxyy_mesh
    
    XX,YY = np.meshgrid(xx_mesh,yy_mesh)
    fig,ax = plt.subplots()
    ax.plot(XX,YY,'-',color='gray')
    ax.plot(XX.T,YY.T,'-',color='gray')
    ax.plot(XX.T,YY.T,'.',color='k')
    
    if label_nodes:
        node_numbering = np.reshape(np.arange(xx_mesh.shape[0]*yy_mesh.shape[0]),
                                   (yy_mesh.shape[0],xx_mesh.shape[0]),order='F')+1
        for x,y,label in zip(XX.ravel(),YY.ravel(),node_numbering.ravel())[::nskip_label]:
            ax.annotate(str(int(label)),xy=(x+label_buffer,y-label_buffer))
            
        
    
    if invert_y:
        ax.set_ylim(ax.get_ylim()[::-1])
    ax.set_ylabel('y-axis')
    ax.set_xlabel('x-axis')
    ax.set_title('Quadrilateral mesh')
    return ax

def plot_tri_mesh(mesh_dict=None,plt_dict={'color':'0.8'},inv_data=None,inv_col=None,
                  nticks=8,keep_log=False,unq_domains=False,force_color=False,
                  figax=None):
    '''Plot triangular finite element mesh.'''
    import matplotlib.tri as mtri
    nodes = mesh_dict['nodes']
    elements = np.array(mesh_dict['elements'])[:,:3]-1
    if unq_domains:
        # Use domain element labels to plot differently
        elem_areas,elem_inds = np.unique(np.array(mesh_dict['elements'])[:,-1], return_index=True)
        tri_obj = []
        colors = plt.cm.rainbow(np.linspace(0,1,elem_areas.shape[0]))
        end_inds = np.hstack([elem_inds[1:],elements.shape[0]])
        
        for elem_ind,end_ind in zip(elem_inds,end_inds):
            temp_tri_obj = mtri.Triangulation(nodes[:,0],nodes[:,1],triangles=elements[elem_ind:end_ind,:])
            tri_obj.append(temp_tri_obj)
    else:
        tri_obj = mtri.Triangulation(nodes[:,0],nodes[:,1],triangles=elements)
    
    if figax is None:
        fig,ax = plt.subplots()
    else:
        fig,ax = figax
    
    if inv_data is None:
        if isinstance(tri_obj,list):
            for icolor,tri_obj_temp in zip(colors,tri_obj):
                if not force_color:
                    plt_dict.update({'color':icolor})
                ax.triplot(tri_obj_temp,**plt_dict)
        else:
            ax.triplot(tri_obj,**plt_dict)
    else:
        if 'vmin' in plt_dict.keys() and 'vmax' in plt_dict.keys():
            ticks = np.linspace(plt_dict['vmin'],plt_dict['vmax'],nticks)
        else:
            ticks = np.linspace(inv_data[:,inv_col].min(),inv_data[:,inv_col].max(),nticks)
        
        all_ER = inv_data[:,inv_col]
        
        if not keep_log and inv_col==3:
            all_ER = 10.**(all_ER)
        elif keep_log:
            plt_dict['vmin'] = np.log10(plt_dict['vmin'])
            plt_dict['vmax'] = np.log10(plt_dict['vmax'])
            ticks = np.linspace(plt_dict['vmin'],plt_dict['vmax'],nticks)
        
        tc = ax.tripcolor(tri_obj,all_ER,**plt_dict)
        plt.colorbar(tc,ax=ax,ticks=ticks)
        
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    
    return fig,ax,tc
    
# ----------------- Inversion output plotting -------------
def plot_outnodes(fname=None,work_dir=None,fwd_transfer_bool=False):
    '''Plot model nodes.'''
    if fwd_transfer_bool==True:
        nheader=1
    else:
        nheader=0
    
    if fname is None and work_dir is not None:
        inv_data = load_inv_output(work_dir=work_dir,nheader=nheader)    
    elif os.path.dirname(fname) in ['']:
        inv_data = load_inv_output(fname=os.path.join(work_dir,fname),nheader=nheader)
    else:
        inv_data = load_inv_output(fname=fname,nheader=nheader)
    
    fig,ax = plt.subplots()
#    ax.plot(inv_data[:,0],inv_data[:,1],'.')
    ax.scatter(inv_data[:,0],inv_data[:,1],c=inv_data[:,3],edgecolor='none')
    ax.set_xlabel('Distance')
    ax.set_ylabel('Depth or Elevation')

def plot_tri_res(fname=None,work_dir=None,inv_data=None,
                 plt_opts={},nheader=0,
                 mesh_dict=None,inv_col=3,keep_log=False,
                 xylims=[None,None,None,None]):
    '''Plot inversion results on triangular mesh.'''
    if inv_data is None: # load data
        if fname is None and work_dir is not None:
            inv_data = load_inv_output(work_dir=work_dir,nheader=nheader)    
        elif os.path.dirname(fname) in ['']:
            inv_data = load_inv_output(fname=os.path.join(work_dir,fname),nheader=nheader)
        else:
            inv_data = load_inv_output(fname=fname,nheader=nheader)
    
    print(inv_data.shape)
    tri_dict = {'mesh_dict':mesh_dict,'inv_data':inv_data,
                'inv_col':inv_col,'plt_dict':plt_opts,
                'keep_log':keep_log}
    fig,ax,tri_obj = plot_tri_mesh(**tri_dict)
    if xylims[0] is not None:
        ax.set_xlim(xylims[:2])
        ax.set_ylim(xylims[2:])
        
    return fig,ax
       
def plot_res(fname=None,work_dir=None,
           topog_xy=None,region_xyz=None,
           method='force_grid',nxny=[1e2,1e2],inv_col=3,
           plt_opts={}, miny=None,maxy=0,
           plot_buffer=1.,cmap=plt.cm.inferno_r,save_dict=None,
           fwd_transfer_bool=False,invert_y=False,interp='bicubic',
           aspect=1,keep_log=False,nticks=8):
    
    '''Simple inverse output plot.
    
    
    interp: str
        str defining interpolation method for imshow
        options: https://matplotlib.org/examples/images_contours_and_fields/interpolation_methods.html
    
    '''
    if fwd_transfer_bool==True:
        nheader=1
    else:
        nheader=0
    
    if isinstance(cmap,str):
        cmap = plt.get_cmap(cmap)
        
    if 'cmap' not in plt_opts.keys():
        plt_opts['cmap'] = cmap
    
    if 'aspect' not in plt_opts.keys() and topog_xy is None:
        plt_opts['aspect'] = aspect
        
    if fname is None and work_dir is not None:
        inv_data = load_inv_output(work_dir=work_dir,nheader=nheader)    
    elif os.path.dirname(fname) in ['']:
        inv_data = load_inv_output(fname=os.path.join(work_dir,fname),nheader=nheader)
    else:
        inv_data = load_inv_output(fname=fname,nheader=nheader)
    
    if 'vmin' not in plt_opts.keys():
        plt_opts['vmin'] = np.min(inv_data[:,inv_col])
    
    if 'vmax' not in plt_opts.keys():
        plt_opts['vmax'] = np.max(inv_data[:,inv_col])
        
    # Make regular grid for plotting surface
    grid_dict = {'inv_data':inv_data,'nxny':nxny,
                 'maxy':maxy,'miny':miny,
                 'interp_method':method,'inv_col':inv_col}

    X,Y,all_ER = grid_inv_data(**grid_dict)
    
    minx,maxx = np.nanmin(X),np.nanmax(X)

    
    if topog_xy is not None:
        ydif = extrap(X[0,:],topog_xy[:,0],topog_xy[:,1])-Y[0,:]
    else:
        ydif=0

    if 'ticks' in list(plt_opts.keys()):
        ticks = plt_opts['ticks']
        _=plt_opts.pop('ticks')
    elif inv_col == 3 and not keep_log:
        # Log transform already applied, change to untransformed
        
        ticks = log_steps([plt_opts['vmin'],plt_opts['vmax']])
    elif inv_col == 3 and keep_log:
        plt_opts['vmin'] = np.log10(plt_opts['vmin'])
        plt_opts['vmax'] = np.log10(plt_opts['vmax'])
        ticks = np.linspace(plt_opts['vmin'],plt_opts['vmax'],nticks)
    else:
        ticks = np.linspace(plt_opts['vmin'],plt_opts['vmax'],nticks)
    
    if not keep_log:
        all_ER = 10.**(all_ER)
        
    upper_ER = np.ma.masked_invalid(all_ER)
    if miny is None:
        miny=np.nanmin(Y)

    fig,ax = plt.subplots()
    if topog_xy is not None:
        Y = Y+ydif
        miny=np.nanmin(Y)
        maxy=np.nanmax(Y)
        im = ax.pcolormesh(X,Y,upper_ER,**plt_opts)

    else:
        im = ax.imshow(upper_ER,interpolation=interp,extent=(minx,maxx,miny,maxy),**plt_opts)

    ax.set_ylabel('Elevation')
    ax.set_xlabel('Distance')
    ax.set_xlim(minx,maxx)
    ax.set_ylim(miny,maxy)
    if 'aspect' in plt_opts.keys():
        ax.set_title('Aspect = {}:1'.format(plt_opts['aspect']))
    plt.colorbar(im,ax=ax,ticks=ticks,extend='both',shrink=0.5)
    if invert_y:
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylabel('Depth')
        
    if save_dict is not None:

        fig.savefig(save_dict['fig_fname'],**save_dict['fig_opts'])
        plt.close('all')
        
    else:
        return fig,ax,[X,Y,upper_ER]

def plot_final_res(fname=None,work_dir=None,
                   region_xyz=None, topo_xyz=None,
                   method='constant',nxny=[1e2,1e2],inv_col=3,
                   plt_opts={}, miny=None,plt_opts2={},region_plot=False,
                   plot_buffer=1.,cmap=plt.cm.inferno_r,save_dict=None,
                   fwd_transfer_bool=False,invert_y=False):

    if fwd_transfer_bool==True:
        nheader=1
    else:
        nheader=0
     
    if isinstance(cmap,str):
        cmap = plt.get_cmap(cmap)
        
    if fname is None and work_dir is not None:
        inv_data = load_inv_output(work_dir=work_dir,nheader=nheader)    
    elif os.path.dirname(fname) in ['']:
        inv_data = load_inv_output(fname=os.path.join(work_dir,fname),nheader=nheader)
    else:
        inv_data = load_inv_output(fname=fname,nheader=nheader)
    
    # Make regular grid for plotting surface
    x = np.linspace(np.min(inv_data[:,0]),np.max(inv_data[:,0]),nxny[0])
    if miny is None:
        miny=np.min(inv_data[:,1])
    y = np.linspace(0.,miny,nxny[1])
    X,Y = np.meshgrid(x,y)
    
   
    if 'vmin' in plt_opts.keys():
        vmin=10.**plt_opts['vmin']
        _ = plt_opts.pop('vmin',None)
    else:
        vmin=10.
    
    if 'vmax' in plt_opts.keys():
        vmax=10.**plt_opts['vmax']
        _ = plt_opts.pop('vmax',None)
    else:
        vmax=1e3
    
    if 'cmap' not in plt_opts.keys():
        plt_opts['cmap'] = cmap
    
    formatter = LogFormatter(10, labelOnlyBase=False)
    ticks = log_steps([vmin,vmax])
        
    if region_xyz is not None and topo_xyz is not None and region_plot:
        # Use region data to separate plotting zones
        # 1) Find elevation of divider (originally depth)
        div_y = extrap(region_xyz[:,0],topo_xyz[:,0],topo_xyz[:,2],method)+region_xyz[:,1]
        # 2) Find y locations above region divider
        div_y_all = extrap(X,region_xyz[:,0],div_y,method)
        # 3) Find y shift for topography
        yshift = extrap(X,topo_xyz[:,0],topo_xyz[:,2],method)
        Y = Y+yshift
        lower_mask = (Y > div_y_all) | ((X > region_xyz[:,0].max()) | (X < region_xyz[:,0].min()))
        upper_mask = (Y < div_y_all) | ((X > region_xyz[:,0].max()) | (X < region_xyz[:,0].min()))
        all_ER = griddata(inv_data[:,:2],
                            inv_data[:,inv_col],
                            (X,Y),
                            method='linear')
        if inv_col == 3:
            # Log transform already applied, change to untransformed
            all_ER = 10.**(all_ER)
        
        upper_ER = np.ma.masked_array(all_ER,mask=upper_mask)
        fig,ax = plt.subplots()
            
        s1 = ax.pcolormesh(X,Y,upper_ER,norm=LogNorm(vmin,vmax),**plt_opts)
#        ax.plot(X,Y,'g.')
        ax.plot(topo_xyz[:,0],topo_xyz[:,2],'k-')
        ax.plot(region_xyz[:,0],div_y,'-o',color='lightgrey')
        ax.set_xlim([region_xyz[:,0].min()-plot_buffer,region_xyz[:,0].max()+plot_buffer])
        ax.set_ylim([miny-plot_buffer,np.max(Y)+plot_buffer])
        plt.colorbar(s1,ax=ax,ticks=ticks,extend='both',format=formatter)
        
        ax2 = ax.twinx()
        lower_ER = np.ma.masked_array(all_ER,mask=lower_mask)
        s2 = ax.pcolormesh(X,Y,lower_ER,norm=LogNorm(vmin,vmax),**plt_opts2)
        ax2.set_xlim([region_xyz[:,0].min()-plot_buffer,region_xyz[:,0].max()+plot_buffer])
        ax2.set_ylim([miny-plot_buffer,np.max(Y)+plot_buffer])
        plt.colorbar(s2,ax=ax2,ticks=ticks,extend='both',format=formatter)
    elif topo_xyz is not None and region_xyz is not None:
        # 1) Find elevation of divider (originally depth)
        div_y = extrap(region_xyz[:,0],topo_xyz[:,0],topo_xyz[:,2],method)+region_xyz[:,1]
        # Find y shift for topography
        yshift = extrap(X,topo_xyz[:,0],topo_xyz[:,2],method)
        Y = Y+yshift
        upper_mask =(X > topo_xyz[:,0].max()) | (X < topo_xyz[:,0].min())
        all_ER = griddata(inv_data[:,:2],
                            inv_data[:,inv_col],
                            (X,Y),
                            method='linear')
        if inv_col == 3:
            # Log transform already applied, change to untransformed
            all_ER = 10.**(all_ER)
            
        upper_ER = np.ma.masked_array(all_ER,mask=upper_mask)
        fig,ax = plt.subplots()
        s1 = ax.pcolormesh(X,Y,upper_ER,norm=LogNorm(vmin,vmax),**plt_opts)
#        ax.plot(X,Y,'g.')
        ax.plot(topo_xyz[:,0],topo_xyz[:,2],'k-')
        ax.plot(region_xyz[:,0],div_y,'-o',color='lightgrey')
        ax.set_xlim([topo_xyz[:,0].min()-plot_buffer,topo_xyz[:,0].max()+plot_buffer])
        ax.set_ylim([miny-plot_buffer,np.max(Y)+plot_buffer])
        plt.colorbar(s1,ax=ax,ticks=ticks,extend='both',format=formatter)
    else:
        if topo_xyz is None: # Copy inv data output y coordinates
            topo_xyz = np.zeros_like(inv_data)
            topo_xyz[:,:2] = inv_data[:,:2]
        # Find y shift for topography
        yshift = extrap(X,topo_xyz[:,0],topo_xyz[:,2],method)
        Y = Y+yshift
        all_ER = griddata(inv_data[:,:2],
                            inv_data[:,inv_col],
                            (X,Y),
                            method='linear')
        if inv_col == 3:
            # Log transform already applied, change to untransformed
            all_ER = 10.**(all_ER)
            
        upper_ER = np.ma.masked_invalid(all_ER)
        fig,ax = plt.subplots()
        s1 = ax.pcolormesh(X,Y,upper_ER,norm=LogNorm(vmin,vmax),**plt_opts)
#        ax.plot(X,Y,'g.')
        ax.plot(topo_xyz[:,0],topo_xyz[:,2],'k-')
        ax.set_xlim([topo_xyz[:,0].min()-plot_buffer,topo_xyz[:,0].max()+plot_buffer])
        ax.set_ylim([miny-plot_buffer,np.max(Y)+plot_buffer])
        plt.colorbar(s1,ax=ax,ticks=ticks,extend='both',format=formatter)
    
    if invert_y:
        ax.set_ylim(ax.get_ylim()[::-1])
        
    # put back in dictionary
    plt_opts['vmin'] = np.log10(vmin)
    plt_opts['vmax'] = np.log10(vmax)
    if save_dict is not None:

        fig.savefig(save_dict['fig_fname'],**save_dict['fig_opts'])
        plt.close('all')
    else:
        return fig,ax
        


def log_steps(lims=None,in_log=False,nvals=10):
    if not in_log:
        lims = np.log10(lims)
    
    lims = np.round(lims,0)
    orders = np.arange(np.diff(lims)+1)+lims[0]
    out_values = []
    for order1,order2 in zip(orders[:-1],orders[1:]):
        out_values.extend(np.linspace(10.**order1,10**order2,nvals))
    
    return np.unique(out_values)