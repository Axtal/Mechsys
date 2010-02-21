#!BPY

########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raul D. D. Farfan             #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################


"""
Name: 'MechSys'
Blender: 2.48
Group: 'Themes'
Tooltip: 'Open Library for Mechanics'
"""
__author__  = "Dorival Pedroso"
__version__ = "19.10.09"
__bpydoc__="""\
TODO
"""

# In Mac OS X, the scripts folder is:
#   /Applications/blender.app/Contents/MacOS/.blender/scripts/

# Modules
import Blender
from   Blender import Draw, BGL
from   Blender.Mathutils import Vector
import bpy
import msys_cad  as ca
import msys_mesh as me
import msys_dict as di
import msys_fem  as fem
import msys_dem  as dem
import msys_gui  as gu

#import rpdb2; rpdb2.start_embedded_debugger('msys')

# ================================================================================= Constants

CLOSE_VERT_TOL       = 1.0e-4 # tolerace to decide whether vertices are close to each other or not
EVT_INC              = 10000  # increment used when passing IDs through events
EVT_NONE             = 0      # used for buttons with callbacks
# SETTINGS
EVT_SET_DELPROPS     = 10 # delete all properties
# CAD
EVT_CAD_ADDXYZ       = 20 # type in values to define point(s)
EVT_CAD_FILLET       = 21 # fillet two edges
EVT_CAD_BREAK        = 22 # break edge
EVT_CAD_BREAKM       = 23 # break edge at middle point
EVT_CAD_EINT         = 24 # edge closest distance
EVT_CAD_FPOINT       = 25 # read points from file
EVT_CAD_FSPLINE      = 26 # create a spline from points in a file
EVT_CAD_EXPTS        = 27 # export points to a file
# Mesh
EVT_MESH_SETVTAG     = 30 # set vertex tag
EVT_MESH_SETETAG     = 31 # set edges tag
EVT_MESH_SETFTAG     = 32 # set faces tag
EVT_MESH_DELALLVTAGS = 33 # delete all vertex tags
EVT_MESH_DELALLETAGS = 34 # delete all edge tags
EVT_MESH_DELALLFTAGS = 35 # delete all face tags
EVT_MESH_DELMESH     = 36 # delete mesh object
# Mesh -- structured
EVT_MESH_ADDBLK      = 40 # set block 2D: 4 or 8 edges, 3D: 8 or 20 edges
EVT_MESH_DELALLBLKS  = 41 # set block 2D: 4 or 8 edges, 3D: 8 or 20 edges
EVT_MESH_GENSTRU     = 42 # generate structured mesh using MechSys module
EVT_MESH_GENSTRUS    = 43 # script for structured mesh generation
# Mesh -- unstructured
EVT_MESH_ADDREG      = 50 # add region
EVT_MESH_DELALLREGS  = 51 # delete all regions
EVT_MESH_ADDHOL      = 52 # add hole
EVT_MESH_DELALLHOLS  = 53 # delete all holes
EVT_MESH_GENUNSTRU   = 54 # generate unstructured mesh using MechSys module
EVT_MESH_GENUNSTRUS  = 55 # script for unstructured mesh generation
# Mesh -- extra
EVT_MESH_ADDLINE     = 60 # add linear element
EVT_MESH_DELALLLINES = 61 # delete all linear elements
# Materials
EVT_MAT_ADDMAT       = 70 # add material
EVT_MAT_DELALLMAT    = 71 # delete all materials
# FEM
EVT_FEM_ADDSTAGE     = 80 # add stage
EVT_FEM_DELSTAGE     = 81 # add stage
EVT_FEM_DELALLSTAGES = 82 # add stage
EVT_FEM_ADDNBRY      = 83 # add nodes boundary
EVT_FEM_ADDEBRY      = 84 # add edges boundary
EVT_FEM_ADDFBRY      = 85 # add faces boundary
EVT_FEM_ADDEATT      = 86 # add element attributes
EVT_FEM_DELALLNBRY   = 87 # delete all nodes boundary
EVT_FEM_DELALLEBRY   = 88 # delete all edges boundary
EVT_FEM_DELALLFBRY   = 89 # delete all faces boundary
EVT_FEM_DELALLEATT   = 90 # delete all element attributes
EVT_FEM_RUN          = 91 # run a FE simulation
EVT_FEM_STOP         = 92 # run a FE simulation
EVT_FEM_SCRIPT       = 93 # generate script for FEM 
EVT_FEM_PARAVIEW     = 94 # view in ParaView
EVT_FEM_DPARAVIEW    = 95 # view in ParaView (dynamics)
EVT_FEM_CLEAR        = 96 # clear error message
# DEM
EVT_DEM_GEN_PKG      = 100 # generate packing
EVT_DEM_GEN_TTT      = 101 # generate true triaxial packing
EVT_DEM_GEN_SCR      = 102 # generate script
EVT_DEM_RUN          = 103 # run simulation
EVT_DEM_STOP         = 104 # stop simulation
EVT_DEM_CLEAR        = 105 # clear error message


# ==================================================================================== Events

def hide_all():
    di.set_key('gui_show_cad',  False)
    di.set_key('gui_show_mesh', False)
    di.set_key('gui_show_mat',  False)
    di.set_key('gui_show_fem',  False)
    di.set_key('gui_show_dem',  False)
    di.set_key('gui_show_res',  False)
    di.set_key('gui_inirow',    0)
    Blender.Window.QRedrawAll()

def show_only(key):
    hide_all()
    di.set_key(key, True)
    Blender.Window.QRedrawAll()

# Handle input events
def event(evt, val):
    if   evt==Draw.QKEY and not val: Draw.Exit()
    elif evt==Draw.ESCKEY: Draw.Redraw(1)
    elif evt==Draw.WHEELDOWNMOUSE or evt==Draw.DOWNARROWKEY:
        d = di.load_dict_for_scroll()
        d['gui_inirow'] -= 40
        Blender.Window.QRedrawAll()
    elif evt==Draw.WHEELUPMOUSE or evt==Draw.UPARROWKEY:
        d = di.load_dict_for_scroll()
        if d['gui_inirow']<0: d['gui_inirow'] += 40
        Blender.Window.QRedrawAll()
    elif evt==Draw.PAGEDOWNKEY:
        d = di.load_dict_for_scroll()
        d['gui_inirow'] -= 180
        Blender.Window.QRedrawAll()
    elif evt==Draw.PAGEUPKEY:
        d = di.load_dict_for_scroll()
        if d['gui_inirow']<0: d['gui_inirow'] += 180
        Blender.Window.QRedrawAll()


def try_catch(func):
    def wrapper(*arg):
        try: func(*arg)
        except Exception, inst:
            msg = inst.args[0]
            print '[1;34mMechSys[0m: Error: '+'[1;31m'+msg+'[0m'
            print '>>>>>>>>>>> exception caught by main <<<<<<<<<<<'
            Blender.Draw.PupMenu('ERROR|'+msg)
    return wrapper


def delete_mesh():
    obj = di.get_obj()
    if obj.properties.has_key('msh_name'):
        old_msh = bpy.data.objects[obj.properties['msh_name']]
        scn     = bpy.data.scenes.active
        scn.objects.unlink (old_msh)
        if obj.properties.has_key('msh_type'): obj.properties.pop('msh_type')
        obj.properties.pop('msh_name')

# Handle button events
#@try_catch ##################### comment this out to track errors ######################################################################################################################################
def button_event(evt):

    # ----------------------------------------------------------------------------------- Settings

    if evt==EVT_SET_DELPROPS:
        scn = bpy.data.scenes.active
        obs = scn.objects.selected
        if len(obs)>0:
            message = 'Confirm delete ALL properties?%t|Yes'
            result  = Blender.Draw.PupMenu(message)
            if result>0:
                for o in obs:
                    for k, v in o.properties.iteritems():
                        print '[1;34mMechSys[0m: Object = ', o.name, ' [1;35mDeleted property:[0m ', k
                        if k=='msh_name':
                            old_msh = bpy.data.objects[v]
                            scn.objects.unlink (old_msh)
                        o.properties.pop(k)
                    for p in o.getAllProperties():
                        print '[1;34mMechSys[0m: Object = ', o.name, ' Deleted property: ', p.name, p.data
                        o.removeProperty(p)
                Blender.Window.QRedrawAll()

    # ----------------------------------------------------------------------------------- CAD

    elif evt==EVT_CAD_ADDXYZ:   ca.add_point     (float(di.key('cad_x')), float(di.key('cad_y')), float(di.key('cad_z')))
    elif evt==EVT_CAD_FILLET:   ca.fillet        (float(di.key('cad_rad')), di.key('cad_stp'))
    elif evt==EVT_CAD_BREAK:    ca.break_edge    ()
    elif evt==EVT_CAD_BREAKM:   ca.break_edge    (True)
    elif evt==EVT_CAD_EINT:     ca.edge_intersect()
    elif evt==EVT_CAD_FPOINT:   Blender.Window.FileSelector(ca.add_points_from_file, 'Read X Y Z columns')
    elif evt==EVT_CAD_FSPLINE:  Blender.Window.FileSelector(ca.add_spline_from_file, 'Read X Y Z columns')
    elif evt==EVT_CAD_EXPTS:    Blender.Window.FileSelector(ca.export_points,        'Export Points', 'allpoints.txt')

    # ---------------------------------------------------------------------------------- Mesh

    elif evt==EVT_MESH_SETVTAG:
        tag = di.key('newvtag')
        edm, obj, msh = di.get_msh()
        for vid in msh.verts.selected():
            di.props_set_with_tag ('vtags', str(vid), tag, di.key('newvtag'))
        Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1)

    elif evt==EVT_MESH_SETETAG:
        tag = di.key('newetag')
        edm, obj, msh = di.get_msh()
        for eid in di.get_selected_edges(msh):
            di.props_set_with_tag ('etags', str(eid), tag, di.key('newetag'))
        Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1)

    elif evt==EVT_MESH_SETFTAG:
        tag = di.key('newftag')
        edm, obj, msh = di.get_msh()
        sel           = di.get_selected_edges(msh)
        nedges        = len(sel)
        if nedges==3 or nedges==6 or nedges==4 or nedges==8:
            eds, vds = di.sort_edges_and_verts (msh, sel, msh.edges[sel[0]].v1.index) # will erase sel
            vids = '_'.join([str(vid) for vid in vds])
            di.props_set_with_tag ('ftags', vids, tag, di.key('newftag'))
        else:
            if edm: Blender.Window.EditMode(1)
            raise Exception('To set a face tag (FTAG), 3, 4, 6, or 8 edges must be selected')
        Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1)

    elif evt==EVT_MESH_DELALLVTAGS: di.props_del_all_tags('vtags')
    elif evt==EVT_MESH_DELALLETAGS: di.props_del_all_tags('etags')
    elif evt==EVT_MESH_DELALLFTAGS: di.props_del_all_tags('ftags')
    elif evt==EVT_MESH_DELMESH:
        msg = 'Confirm delete mesh?%t|Yes'
        res = Blender.Draw.PupMenu(msg)
        if res>0:
            delete_mesh()
            Blender.Window.QRedrawAll()

    # -------------------------------------------------------------------- Mesh -- structured

    elif evt==EVT_MESH_ADDBLK:
        props = di.new_blk_props()
        di.props_push_new('blks', props, True, 18, 18+int(props[17]))

    elif evt==EVT_MESH_DELALLBLKS: di.props_del_all('blks')
    elif evt==EVT_MESH_GENSTRU:
        obj  = di.get_obj()
        mesh = me.gen_struct_mesh ()
        me.add_mesh (obj, mesh, 'struct')
    elif evt==EVT_MESH_GENSTRUS:
        obj = di.get_obj()
        txt = me.gen_struct_mesh(True,None,di.key('smsh_cpp'))
        if not di.key('smsh_cpp'):
            txt.write('\nobj = bpy.data.objects["'+obj.name+'"]\n')
            txt.write('me.add_mesh (obj, mesh, "struct")\n')

    # ------------------------------------------------------------------ Mesh -- unstructured

    elif evt==EVT_MESH_ADDREG:     di.props_push_new('regs', di.new_reg_props())
    elif evt==EVT_MESH_DELALLREGS: di.props_del_all ('regs')
    elif evt==EVT_MESH_ADDHOL:     di.props_push_new('hols', di.new_hol_props())
    elif evt==EVT_MESH_DELALLHOLS: di.props_del_all ('hols')
    elif evt==EVT_MESH_GENUNSTRU:  me.gen_unstruct_mesh ()
    elif evt==EVT_MESH_GENUNSTRUS: me.gen_unstruct_mesh (True,None,di.key('umsh_cpp'))

    # ------------------------------------------------------------------------- Mesh -- extra

    elif evt==EVT_MESH_ADDLINE:      di.props_push_new('lins', di.new_lin_props())
    elif evt==EVT_MESH_DELALLLINES:  di.props_del_all ('lins')

    # ---------------------------------------------------------------------------------- Materials

    elif evt==EVT_MAT_ADDMAT:    di.props_push_new_mat()
    elif evt==EVT_MAT_DELALLMAT: di.props_del_all_mats()

    # ----------------------------------------------------------------------------------- FEM 

    elif evt==EVT_FEM_ADDSTAGE: di.props_push_new_stage()
    elif evt==EVT_FEM_DELSTAGE: di.props_del_stage     ()
    elif evt==EVT_FEM_ADDNBRY:  di.props_push_new_fem  ([di.key('fem_stage')], 'nbrys', di.new_nbry_props())
    elif evt==EVT_FEM_ADDEBRY:  di.props_push_new_fem  ([di.key('fem_stage')], 'ebrys', di.new_ebry_props())
    elif evt==EVT_FEM_ADDFBRY:  di.props_push_new_fem  ([di.key('fem_stage')], 'fbrys', di.new_fbry_props())
    elif evt==EVT_FEM_ADDEATT:
        obj  = di.get_obj()
        sids = [k for k in obj.properties['stages']]
        di.props_push_new_fem (sids, 'eatts', di.new_eatt_props())

    elif evt==EVT_FEM_DELALLSTAGES: di.props_del_all_stages()
    elif evt==EVT_FEM_DELALLNBRY:   di.props_del_all_fem   (str(di.key('fem_stage')), 'nbrys')
    elif evt==EVT_FEM_DELALLEBRY:   di.props_del_all_fem   (str(di.key('fem_stage')), 'ebrys')
    elif evt==EVT_FEM_DELALLFBRY:   di.props_del_all_fem   (str(di.key('fem_stage')), 'fbrys')
    elif evt==EVT_FEM_DELALLEATT:
        obj  = di.get_obj()
        sids = [int(k) for k, v in obj.properties['stages'].iteritems()]
        di.props_del_all_fem (sids, 'eatts')

    elif evt==EVT_FEM_SCRIPT:    fem.gen_script ()
    elif evt==EVT_FEM_RUN:       fem.run        ()
    elif evt==EVT_FEM_STOP:      fem.stop       ()
    elif evt==EVT_FEM_PARAVIEW:  fem.paraview   ()
    elif evt==EVT_FEM_DPARAVIEW: fem.paraview   (True)
    elif evt==EVT_FEM_CLEAR:
        di.key('fem_fatal').value = 0
        Blender.Window.QRedrawAll()

    # ----------------------------------------------------------------------------------- DEM 

    elif evt==EVT_DEM_GEN_PKG: dem.gen_pkg    (False)
    elif evt==EVT_DEM_GEN_TTT: dem.gen_pkg    (True)
    elif evt==EVT_DEM_GEN_SCR: dem.gen_script ()
    elif evt==EVT_DEM_RUN:     dem.run        ()
    elif evt==EVT_DEM_STOP:    dem.stop       ()
    elif evt==EVT_DEM_CLEAR:
        di.key('dem_fatal').value = 0
        Blender.Window.QRedrawAll()

# ================================================================================= Callbacks

# ---------------------------------- Show/Hide

@try_catch
def cb_gui_show_cad(evt,val):
    di.set_key ('gui_show_cad', val)
    show_only  ('gui_show_cad')
@try_catch
def cb_gui_show_mesh(evt,val):
    di.set_key ('gui_show_mesh', val)
    show_only  ('gui_show_mesh')
@try_catch
def cb_gui_show_mat(evt,val):
    di.set_key ('gui_show_mat', val)
    show_only  ('gui_show_mat')
@try_catch
def cb_gui_show_fem(evt,val):
    di.set_key ('gui_show_fem', val)
    show_only  ('gui_show_fem')
@try_catch
def cb_gui_show_dem(evt,val):
    di.set_key ('gui_show_dem', val)
    show_only  ('gui_show_dem')

# ---------------------------------- Settings

@try_catch
def cb_show_props(evt,val): di.set_key ('show_props', val)
@try_catch
def cb_show_e_ids(evt,val): di.set_key ('show_e_ids', val)
@try_catch
def cb_show_v_ids(evt,val): di.set_key ('show_v_ids', val)
@try_catch
def cb_show_blks (evt,val): di.set_key ('show_blks',  val)
@try_catch
def cb_show_Bids (evt,val): di.set_key ('show_Bids',  val)
@try_catch
def cb_show_axes (evt,val): di.set_key ('show_axes',  val)
@try_catch
def cb_show_regs (evt,val): di.set_key ('show_regs',  val)
@try_catch
def cb_show_vtags(evt,val): di.set_key ('show_vtags', val)
@try_catch
def cb_show_etags(evt,val): di.set_key ('show_etags', val)
@try_catch
def cb_show_ftags(evt,val): di.set_key ('show_ftags', val)
@try_catch
def cb_over(evt,val):
    edm, obj, msh = di.get_msh()
    if val: # mark
        nedges = len(msh.edges)
        over   = []
        for i in range(nedges):
            v1 = msh.edges[i].v1.co
            v2 = msh.edges[i].v2.co
            uu = v2-v1
            for j in range(i+1,nedges):
                v3 = msh.edges[j].v1.co
                v4 = msh.edges[j].v2.co
                t3 = (v2[2]*v3[2]-v1[2]*v3[2]-v1[2]*v2[2]+v1[2]**2.0+v2[1]*v3[1]-v1[1]*v3[1]-v1[1]*v2[1]+v1[1]**2.0+v2[0]*v3[0]-v1[0]*v3[0]-v1[0]*v2[0]+v1[0]**2.0)/(v2[2]**2.0-2.0*v1[2]*v2[2]+v1[2]**2.0+v2[1]**2.0-2.0*v1[1]*v2[1]+v1[1]**2.0+v2[0]**2.0-2.0*v1[0]*v2[0]+v1[0]**2.0)
                t4 = (v2[2]*v4[2]-v1[2]*v4[2]-v1[2]*v2[2]+v1[2]**2.0+v2[1]*v4[1]-v1[1]*v4[1]-v1[1]*v2[1]+v1[1]**2.0+v2[0]*v4[0]-v1[0]*v4[0]-v1[0]*v2[0]+v1[0]**2.0)/(v2[2]**2.0-2.0*v1[2]*v2[2]+v1[2]**2.0+v2[1]**2.0-2.0*v1[1]*v2[1]+v1[1]**2.0+v2[0]**2.0-2.0*v1[0]*v2[0]+v1[0]**2.0)
                w3 = v1+t3*uu
                w4 = v1+t4*uu
                if ((v3-w3).length<CLOSE_VERT_TOL) and ((v4-w4).length<CLOSE_VERT_TOL): # parallel at the same supporting line
                    if (t3>0.0 and t3<1.0) or (t4>0.0 and t4<1.0): # overlapping
                        if not i in over: over.append(i)
                        if not j in over: over.append(j)
        if len(over)>0:
            obj.properties['over'] = over
            Blender.Draw.PupMenu('Found %d overlapping edges'%len(over))
        else: Blender.Draw.PupMenu('Could not find any overlapping edges')
    else:
        if obj.properties.has_key('over'): obj.properties.pop('over')
    Blender.Window.QRedrawAll()
@try_catch
def cb_hide_mesh(evt,val):
    obj = di.get_obj()
    if obj.properties.has_key('msh_name'):
        msh_obj = bpy.data.objects[obj.properties['msh_name']]
        if val: msh_obj.layers = []
        else:   msh_obj.layers = [1]
        di.set_key('hide_mesh', val)
        Blender.Window.QRedrawAll()

# ---------------------------------- CAD

def cb_set_x     (evt,val): di.set_key('cad_x',   val)
@try_catch
def cb_set_y     (evt,val): di.set_key('cad_y',   val)
@try_catch
def cb_set_z     (evt,val): di.set_key('cad_z',   val)
@try_catch
def cb_fillet_rad(evt,val): di.set_key('cad_rad', val)
@try_catch
def cb_fillet_stp(evt,val): di.set_key('cad_stp', val)

# ---------------------------------- Mesh

@try_catch
def cb_is3d(evt,val): di.props_set_val('is3d', val)
@try_catch
def cb_iso2(evt,val):
    msg = 'This action will delete current mesh. Confirm?%t|Yes'
    res = Blender.Draw.PupMenu(msg)
    if res>0:
        delete_mesh()
        di.props_set_val('iso2', val)
@try_catch
def cb_frame (evt,val):
    if val==1: di.props_set_val('msh_type', 'frame')
    else:
        obj = di.get_obj()
        if obj.properties.has_key('msh_type'): obj.properties.pop('msh_type')
        Blender.Window.QRedrawAll()
@try_catch
def cb_vtag  (evt,val): di.set_key      ('newvtag', val)
@try_catch
def cb_etag  (evt,val): di.set_key      ('newetag', val)
@try_catch
def cb_ftag  (evt,val): di.set_key      ('newftag', val)
@try_catch
def cb_smsh_cpp (evt,val): di.set_key('smsh_cpp', val)
@try_catch
def cb_umsh_cpp (evt,val): di.set_key('umsh_cpp', val)

# ---------------------------------- Mesh -- structured

@try_catch
def cb_blk_tag(evt,val): di.props_set_item ('blks', evt-EVT_INC,  0, val)
@try_catch
def cb_blk_xax(evt,val): di.blk_set_axis   (        evt-EVT_INC,  1)
@try_catch
def cb_blk_yax(evt,val): di.blk_set_axis   (        evt-EVT_INC,  2)
@try_catch
def cb_blk_zax(evt,val): di.blk_set_axis   (        evt-EVT_INC,  3)
@try_catch
def cb_blk_nx (evt,val): di.props_set_item ('blks', evt-EVT_INC,  8, val)
@try_catch
def cb_blk_ny (evt,val): di.props_set_item ('blks', evt-EVT_INC,  9, val)
@try_catch
def cb_blk_nz (evt,val): di.props_set_item ('blks', evt-EVT_INC, 10, val)
@try_catch
def cb_blk_ax (evt,val): di.props_set_item ('blks', evt-EVT_INC, 11, float(val))
@try_catch
def cb_blk_ay (evt,val): di.props_set_item ('blks', evt-EVT_INC, 12, float(val))
@try_catch
def cb_blk_az (evt,val): di.props_set_item ('blks', evt-EVT_INC, 13, float(val))
@try_catch
def cb_blk_nlx(evt,val): di.props_set_item ('blks', evt-EVT_INC, 14, val)
@try_catch
def cb_blk_nly(evt,val): di.props_set_item ('blks', evt-EVT_INC, 15, val)
@try_catch
def cb_blk_nlz(evt,val): di.props_set_item ('blks', evt-EVT_INC, 16, val)
@try_catch
def cb_blk_del(evt,val): di.props_del      ('blks', evt-EVT_INC)

# ---------------------------------- Mesh -- unstructured

@try_catch
def cb_minang (evt,val): di.props_set_val('minang', float(val))
@try_catch
def cb_maxarea(evt,val): di.props_set_val('maxarea',float(val))

# ---------------------------------- Mesh -- unstructured -- regions

@try_catch
def cb_regs_tag    (evt,val): di.props_set_item ('regs', evt-EVT_INC, 0, val)
@try_catch
def cb_regs_maxarea(evt,val): di.props_set_item ('regs', evt-EVT_INC, 1, float(val))
@try_catch
def cb_regs_setx   (evt,val): di.props_set_item ('regs', evt-EVT_INC, 2, float(val))
@try_catch
def cb_regs_sety   (evt,val): di.props_set_item ('regs', evt-EVT_INC, 3, float(val))
@try_catch
def cb_regs_setz   (evt,val): di.props_set_item ('regs', evt-EVT_INC, 4, float(val))
@try_catch
def cb_regs_del    (evt,val): di.props_del      ('regs', evt-EVT_INC)

# ---------------------------------- Mesh -- unstructured -- holes

@try_catch
def cb_hols_setx(evt,val): di.props_set_item ('hols', evt-EVT_INC, 0, float(val))
@try_catch
def cb_hols_sety(evt,val): di.props_set_item ('hols', evt-EVT_INC, 1, float(val))
@try_catch
def cb_hols_setz(evt,val): di.props_set_item ('hols', evt-EVT_INC, 2, float(val))
@try_catch
def cb_hols_del (evt,val): di.props_del      ('hols', evt-EVT_INC)

# ---------------------------------- Mesh -- extra -- linear elements

@try_catch
def cb_show_lins(evt,val): di.set_key ('show_lins', val)
@try_catch
def cb_mesh_letag (evt,val): di.props_set_item ('lins', evt-EVT_INC, 0, int(val))
@try_catch
def cb_mesh_getN0  (evt,val):
    msh_obj = di.get_msh_obj (di.get_obj())
    di.props_set_item ('lins', evt-EVT_INC, 1, di.find_node(msh_obj))
@try_catch
def cb_mesh_leN0  (evt,val): di.props_set_item ('lins', evt-EVT_INC, 1, float(val))
@try_catch
def cb_mesh_getN1  (evt,val):
    msh_obj = di.get_msh_obj (di.get_obj())
    di.props_set_item ('lins', evt-EVT_INC, 2, di.find_node(msh_obj))
@try_catch
def cb_mesh_leN1  (evt,val): di.props_set_item ('lins', evt-EVT_INC, 2, float(val))
@try_catch
def cb_mesh_ledel (evt,val): di.props_del      ('lins', evt-EVT_INC)

# ---------------------------------- Materials -- mats

@try_catch
def cb_mat_setname  (evt,val): di.props_set_text ('texts', evt-EVT_INC, val)
@try_catch
def cb_mat_setmodel (evt,val): di.props_set_item ('mats', evt-EVT_INC, 0, val-1)
@try_catch
def cb_mat_del      (evt,val): di.props_del      ('mats', evt-EVT_INC)

@try_catch
def cb_mat_setE  (evt,val): di.props_set_item ('mats', evt-EVT_INC, 2, float(val))
@try_catch
def cb_mat_setA  (evt,val): di.props_set_item ('mats', evt-EVT_INC, 3, float(val))
@try_catch
def cb_mat_setIzz(evt,val): di.props_set_item ('mats', evt-EVT_INC, 4, float(val))
@try_catch
def cb_mat_setnu (evt,val): di.props_set_item ('mats', evt-EVT_INC, 5, float(val))
@try_catch
def cb_mat_setfc (evt,val): di.props_set_item ('mats', evt-EVT_INC, 6, val-1)
@try_catch
def cb_mat_setsY (evt,val): di.props_set_item ('mats', evt-EVT_INC, 7, float(val))
@try_catch
def cb_mat_setcu (evt,val): di.props_set_item ('mats', evt-EVT_INC, 8, float(val))
@try_catch
def cb_mat_setlam(evt,val): di.props_set_item ('mats', evt-EVT_INC, 9, float(val))
@try_catch
def cb_mat_setkap(evt,val): di.props_set_item ('mats', evt-EVT_INC, 10, float(val))
@try_catch
def cb_mat_setphi(evt,val): di.props_set_item ('mats', evt-EVT_INC, 11, float(val))
@try_catch
def cb_mat_setk  (evt,val): di.props_set_item ('mats', evt-EVT_INC, 12, float(val))

# ---------------------------------- FEM -- nbrys

@try_catch
def cb_nbry_tag(evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC, 0, int(val))
@try_catch
def cb_nbry_key(evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC, 1, val-1)
@try_catch
def cb_nbry_val(evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC, 2, float(val))
@try_catch
def cb_nbry_del(evt,val): di.props_del_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC)

# ---------------------------------- FEM -- ebrys

@try_catch
def cb_ebry_tag(evt,val): di.props_set_fem ([di.key('fem_stage')], 'ebrys', evt-EVT_INC, 0, int(val))
@try_catch
def cb_ebry_key(evt,val): di.props_set_fem ([di.key('fem_stage')], 'ebrys', evt-EVT_INC, 1, val-1)
@try_catch
def cb_ebry_val(evt,val): di.props_set_fem ([di.key('fem_stage')], 'ebrys', evt-EVT_INC, 2, float(val))
@try_catch
def cb_ebry_del(evt,val): di.props_del_fem ([di.key('fem_stage')], 'ebrys', evt-EVT_INC)

# ---------------------------------- FEM -- fbrys

@try_catch
def cb_fbry_tag(evt,val): di.props_set_fem ([di.key('fem_stage')], 'fbrys', evt-EVT_INC, 0, int(val))
@try_catch
def cb_fbry_key(evt,val): di.props_set_fem ([di.key('fem_stage')], 'fbrys', evt-EVT_INC, 1, val-1)
@try_catch
def cb_fbry_val(evt,val): di.props_set_fem ([di.key('fem_stage')], 'fbrys', evt-EVT_INC, 2, float(val))
@try_catch
def cb_fbry_del(evt,val): di.props_del_fem ([di.key('fem_stage')], 'fbrys', evt-EVT_INC)

# ---------------------------------- FEM -- eatts

@try_catch
def cb_eatt_tag  (evt,val): di.props_set_fem_all_stg ('eatts', evt-EVT_INC, 0, int(val))
@try_catch
def cb_eatt_gebp (evt,val): di.props_set_fem_all_stg ('eatts', evt-EVT_INC, 2, int(val)-1)
@try_catch
def cb_eatt_gty  (evt,val): di.props_set_fem_all_stg ('eatts', evt-EVT_INC, 3, int(val)-1)
@try_catch
def cb_eatt_mat  (evt,val): di.props_set_fem_all_stg ('eatts', evt-EVT_INC, 4, int(val)-1)
@try_catch
def cb_eatt_h    (evt,val): di.props_set_fem_all_stg ('eatts', evt-EVT_INC, 5, float(val))
@try_catch
def cb_eatt_isact(evt,val): di.props_set_fem ([di.key('fem_stage')], 'eatts', evt-EVT_INC, 6, int(val))
@try_catch
def cb_eatt_act  (evt,val): di.props_set_fem ([di.key('fem_stage')], 'eatts', evt-EVT_INC, 7, int(val))
@try_catch
def cb_eatt_deact(evt,val): di.props_set_fem ([di.key('fem_stage')], 'eatts', evt-EVT_INC, 8, int(val))
@try_catch
def cb_eatt_del  (evt,val): di.props_del_fem_all_stg ('eatts', evt-EVT_INC)

# ---------------------------------- FEM

@try_catch
def cb_fem_cpp   (evt,val): di.set_key('fem_cpp', val)

# ---------------------------------- FEM -- stages

@try_catch
def cb_fem_prob  (evt,val): di.props_set_val ('pty', val-1)
@try_catch
def cb_fem_stage (evt,val):
    obj = di.get_obj()
    for k, v in obj.properties['stages'].iteritems():
        if val==v[0]: di.set_key('fem_stage', int(k))
    Blender.Window.QRedrawAll()
@try_catch
def cb_fem_stage_desc       (evt,val): di.props_set_text ('texts', evt-EVT_INC, val)
@try_catch
def cb_fem_stage_apply_bf   (evt,val): di.props_set_item ('stages', evt-EVT_INC, 2, int(val))
@try_catch
def cb_fem_stage_clear_disp (evt,val): di.props_set_item ('stages', evt-EVT_INC, 3, int(val))
@try_catch
def cb_fem_stage_ndiv       (evt,val): di.props_set_item ('stages', evt-EVT_INC, 4, int(val))
@try_catch
def cb_fem_stage_act        (evt,val): di.props_set_item ('stages', evt-EVT_INC, 5, int(val))
@try_catch
def cb_fem_stage_type       (evt,val): di.props_set_item ('stages', evt-EVT_INC, 6, val-1)
@try_catch
def cb_fem_stage_tf         (evt,val): di.props_set_item ('stages', evt-EVT_INC, 7, float(val))
@try_catch
def cb_fem_stage_dt         (evt,val): di.props_set_item ('stages', evt-EVT_INC, 8, float(val))
@try_catch
def cb_fem_stage_dtOut      (evt,val): di.props_set_item ('stages', evt-EVT_INC, 9, float(val))
@try_catch
def cb_fem_stage_output     (evt,val): di.props_set_item ('stages', evt-EVT_INC, 31, int(val))
@try_catch
def cb_fem_stage_tolR       (evt,val): di.props_set_item ('stages', evt-EVT_INC, 10, float(val))
@try_catch
def cb_fem_stage_calcWork   (evt,val): di.props_set_item ('stages', evt-EVT_INC, 11, int(val))
@try_catch
def cb_fem_stage_scheme     (evt,val): di.props_set_item ('stages', evt-EVT_INC, 12, val-1)
@try_catch
def cb_fem_stage_dscheme    (evt,val): di.props_set_item ('stages', evt-EVT_INC, 24, val-1)
@try_catch
def cb_fem_stage_DampTy     (evt,val): di.props_set_item ('stages', evt-EVT_INC, 25, val-1)
@try_catch
def cb_fem_stage_DampAm     (evt,val): di.props_set_item ('stages', evt-EVT_INC, 26, float(val))
@try_catch
def cb_fem_stage_DampAk     (evt,val): di.props_set_item ('stages', evt-EVT_INC, 27, float(val))
@try_catch
def cb_fem_stage_DynTh1     (evt,val): di.props_set_item ('stages', evt-EVT_INC, 28, float(val))
@try_catch
def cb_fem_stage_DynTh2     (evt,val): di.props_set_item ('stages', evt-EVT_INC, 29, float(val))

# ---------------------------------- DEM

@try_catch
def cb_dem_Lx(evt,val): di.set_key ('dem_Lx', float(val))
@try_catch
def cb_dem_Ly(evt,val): di.set_key ('dem_Ly', float(val))
@try_catch
def cb_dem_Lz(evt,val): di.set_key ('dem_Lz', float(val))
@try_catch
def cb_dem_Nx(evt,val): di.set_key ('dem_Nx', val)
@try_catch
def cb_dem_Ny(evt,val): di.set_key ('dem_Ny', val)
@try_catch
def cb_dem_Nz(evt,val): di.set_key ('dem_Nz', val)
@try_catch
def cb_dem_R(evt,val): di.set_key ('dem_R', float(val))
@try_catch
def cb_dem_rho(evt,val): di.set_key ('dem_rho', float(val))
@try_catch
def cb_dem_seed(evt,val): di.set_key ('dem_seed', val)
@try_catch
def cb_dem_prob(evt,val): di.set_key ('dem_prob', float(val))
@try_catch
def cb_dem_pkg(evt,val): di.set_key ('dem_pkg', val-1)
@try_catch
def cb_dem_res(evt,val): di.set_key ('dem_res', val)
@try_catch
def cb_dem_draw_verts(evt,val): di.set_key ('dem_draw_verts', val)
@try_catch
def cb_dem_draw_edges(evt,val): di.set_key ('dem_draw_edges', val)

@try_catch
def cb_dem_Kn (evt,val): di.set_key ('dem_Kn', float(val))
@try_catch
def cb_dem_mu (evt,val): di.set_key ('dem_mu', float(val))
@try_catch
def cb_dem_Kt (evt,val): di.set_key ('dem_Kt', float(val))
@try_catch
def cb_dem_beta (evt,val): di.set_key ('dem_beta', float(val))
@try_catch
def cb_dem_Gn (evt,val): di.set_key ('dem_Gn', float(val))
@try_catch
def cb_dem_eta (evt,val): di.set_key ('dem_eta', float(val))
@try_catch
def cb_dem_Gt (evt,val): di.set_key ('dem_Gt', float(val))

@try_catch
def cb_dem_iso_pf (evt,val): di.set_key ('dem_iso_pf', float(val))
@try_catch
def cb_dem_iso_timef (evt,val): di.set_key ('dem_iso_timef', float(val))
@try_catch
def cb_dem_iso_dt (evt,val): di.set_key ('dem_iso_dt', float(val))
@try_catch
def cb_dem_iso_dtout (evt,val): di.set_key ('dem_iso_dtout', float(val))
@try_catch
def cb_dem_iso_render (evt,val): di.set_key ('dem_iso_render', val)

@try_catch
def cb_dem_ttt_pf (evt,val): di.set_key ('dem_ttt_pf', float(val))
@try_catch
def cb_dem_ttt_qf (evt,val): di.set_key ('dem_ttt_qf', float(val))
@try_catch
def cb_dem_ttt_thf (evt,val): di.set_key ('dem_ttt_thf', float(val))
@try_catch
def cb_dem_ttt_pex (evt,val): di.set_key ('dem_ttt_pex', val)
@try_catch
def cb_dem_ttt_pey (evt,val): di.set_key ('dem_ttt_pey', val)
@try_catch
def cb_dem_ttt_pez (evt,val): di.set_key ('dem_ttt_pez', val)
@try_catch
def cb_dem_ttt_exf (evt,val): di.set_key ('dem_ttt_exf', float(val))
@try_catch
def cb_dem_ttt_eyf (evt,val): di.set_key ('dem_ttt_eyf', float(val))
@try_catch
def cb_dem_ttt_ezf (evt,val): di.set_key ('dem_ttt_ezf', float(val))
@try_catch
def cb_dem_ttt_timef (evt,val): di.set_key ('dem_ttt_timef', float(val))
@try_catch
def cb_dem_ttt_dt (evt,val): di.set_key ('dem_ttt_dt', float(val))
@try_catch
def cb_dem_ttt_dtout (evt,val): di.set_key ('dem_ttt_dtout', float(val))
@try_catch
def cb_dem_ttt_render (evt,val): di.set_key ('dem_ttt_render', val)

@try_catch
def cb_dem_cpp_script (evt,val): di.set_key ('dem_cpp_script', val)

# ======================================================================================= GUI

# Draw GUI
def gui():
    # load dictionary
    d = di.load_dict()

    # Get current selected mesh
    try:                    edm, obj, msh = di.get_msh (False)
    except Exception, inst: edm, obj, msh = 0, None, None

    # Data from current selected object
    markover = False
    is3d     = False
    iso2     = False
    isframe  = False
    blks     = {}
    minang   = -1.0
    maxarea  = -1.0
    regs     = {}
    hols     = {}
    lins     = {}
    mats     = {}
    texts    = {}
    stages   = {}
    nstages  = 0
    nbrys    = {}
    ebrys    = {}
    fbrys    = {}
    eatts    = {}
    pty      = 0
    if obj!=None:
        if obj.properties.has_key('over'):  markover   = True
        if obj.properties.has_key('is3d'):  is3d       = obj.properties['is3d']
        if obj.properties.has_key('iso2'):  iso2       = obj.properties['iso2']
        if obj.properties.has_key('msh_type'):
            if obj.properties['msh_type']=='frame': isframe = True
        if obj.properties.has_key('blks'):    blks     = obj.properties['blks']
        if obj.properties.has_key('minang'):  minang   = obj.properties['minang']
        if obj.properties.has_key('maxarea'): maxarea  = obj.properties['maxarea']
        if obj.properties.has_key('regs'):    regs     = obj.properties['regs']
        if obj.properties.has_key('hols'):    hols     = obj.properties['hols']
        if obj.properties.has_key('lins'):    lins     = obj.properties['lins']
        if obj.properties.has_key('mats'):    mats     = obj.properties['mats']
        if obj.properties.has_key('texts'):   texts    = obj.properties['texts']
        if obj.properties.has_key('stages'):  stages   = obj.properties['stages']
        nstages = len(stages)
        if nstages>0:
            stg = 'stg_'+str(d['fem_stage'])
            if obj.properties[stg].has_key('nbrys'): nbrys = obj.properties[stg]['nbrys']
            if obj.properties[stg].has_key('ebrys'): ebrys = obj.properties[stg]['ebrys']
            if obj.properties[stg].has_key('fbrys'): fbrys = obj.properties[stg]['fbrys']
            if obj.properties[stg].has_key('eatts'): eatts = obj.properties[stg]['eatts']
        if obj.properties.has_key('pty'):  pty = obj.properties['pty']

    # materials menu
    matmnu   = 'Materials %t'
    matnames = {}
    for k, v in mats.iteritems():
        tid              = int(v[1])       # text_id
        desc             = texts[str(tid)] # description
        matmnu          += '|'+desc+' %x'+str(int(k)+1)
        matnames[int(k)] = desc

    # restore EditMode
    if edm: Blender.Window.EditMode(1)

    # constants
    rg  = 12    # row gap
    cg  = 10    # column gap
    rh  = 20    # row height
    srg =  6    # small row gap

    # geometry
    W,H = Blender.Window.GetAreaSize()
    r   = H-rh-rg-d['gui_inirow']
    c   = cg
    w   = W-2*c
    r0  = r # initial r

    # materials: extra rows
    mat_extra_rows = 0
    for k, v in mats.iteritems():
        mdl = int(v[0])
        mat_extra_rows += d['mdl2nrow'][mdl]

    # stages: extra rows
    stg_extra_rows = 0
    if len(stages)>0:
        sid = str(d['fem_stage'])  # stage_id
        if int(stages[sid][6])==0: # equilib/steady
            stg_extra_rows = 2
        elif int(stages[sid][6])==1: # transient
            stg_extra_rows = 3
        elif int(stages[sid][6])==2: # dynamics
            stg_extra_rows = 6

    # height of boxes
    h_set           = 6*rh+2*srg+2*rg
    h_cad           = 4*rh+2*rg
    h_msh_stru_blks = rh+srg+2*rh*len(blks) if len(blks)>0 else 0
    h_msh_stru      = 2*rh+srg+h_msh_stru_blks+2*rg
    h_msh_unst_regs = rh+srg+rh*len(regs) if len(regs)>0 else 0
    h_msh_unst_hols = rh+srg+rh*len(hols) if len(hols)>0 else 0
    h_msh_unst      = 4*rh+3*srg+h_msh_unst_regs+h_msh_unst_hols+2*rg
    h_msh_ext_lins  = rh+srg+rh*len(lins) if len(lins)>0 else 0
    h_msh_ext       = rh+2*rg+h_msh_ext_lins
    h_msh           = 11*rh+srg+h_msh_stru+h_msh_unst+2*rg+h_msh_ext
    h_mat_mats      = rg+(srg+rh)*len(mats)+rh*mat_extra_rows if len(mats)>0 else 0
    h_mat           = rh+h_mat_mats+2*rg
    h_fem_nbrys     = rh+srg+rh*len(nbrys)                  if len(nbrys)>0 else 0
    h_fem_ebrys     = rh+srg+rh*len(ebrys)                  if len(ebrys)>0 else 0
    h_fem_fbrys     = rh+srg+rh*len(fbrys)                  if len(fbrys)>0 else 0
    h_fem_eatts     = rg+srg+rh*len(eatts)*2+srg*len(eatts) if len(eatts)>0 else 0
    h_fem_stage     = 13*rh+srg+h_fem_nbrys+h_fem_ebrys+h_fem_fbrys+h_fem_eatts+rh*stg_extra_rows
    h_fem           = 5*rh+h_fem_stage+3*rg + (rh if (d['fem_running'].value or d['fem_fatal'].value) else 0)
    h_dem_pkg       = 8*rh+3*rg
    h_dem_cte       = 5*rh+srg
    h_dem_ttt_iso   = 3*rh+2*srg
    h_dem_ttt_she   = 5*rh+2*srg
    h_dem_ttt       = 6*rh+2*srg+h_dem_ttt_iso+h_dem_ttt_she + (rh if (d['dem_running'].value or d['dem_fatal'].value) else 0)
    h_dem           = h_dem_pkg+h_dem_cte+srg+6*rh+h_dem_ttt

    # clear background
    gu.background()

    # ======================================================== Settings

    gu.caption1_(c,r0,w,rh,'SETTINGS')
    r = r0
    r, c, w = gu.box1_in(W,cg,rh,rg, c,r,w,h_set)
    Draw.Toggle ('Show/Hide All', EVT_NONE, c, r, 320, rh, d['show_props'], 'Show mesh properties'   , cb_show_props)
    r -= srg
    r -= rh
    gu.text(c,r,"Geometry:")
    Draw.Toggle ('Vertices',   EVT_NONE, c+ 80, r, 80, rh, d['show_v_ids'], 'Show Vertex IDs'        , cb_show_v_ids)
    Draw.Toggle ('Edges',      EVT_NONE, c+160, r, 80, rh, d['show_e_ids'], 'Show Edges IDs'         , cb_show_e_ids)
    r -= rh
    gu.text(c,r,"Mesh:")
    Draw.Toggle ('Local Axes', EVT_NONE, c+ 80, r, 80, rh, d['show_axes'],  'Show local system axes'  , cb_show_axes)
    Draw.Toggle ('Blocks',     EVT_NONE, c+160, r, 80, rh, d['show_blks'],  'Show blocks tags'        , cb_show_blks)
    Draw.Toggle ('B Loc IDs',  EVT_NONE, c+240, r, 80, rh, d['show_Bids'],  'Show local IDs of vertices in blocks', cb_show_Bids)
    r -= rh
    Draw.Toggle ('Regs/Holes', EVT_NONE, c+ 80, r, 80, rh, d['show_regs'],  'Show regions and holes' , cb_show_regs)
    Draw.Toggle ('Lines',      EVT_NONE, c+160, r, 80, rh, d['show_lins'],  'Show linear cells'      , cb_show_lins)
    r -= rh
    gu.text(c,r,"Tags:")
    Draw.Toggle ('V Tags',     EVT_NONE, c+ 80, r, 80, rh, d['show_vtags'], 'Show vertex tags'       , cb_show_vtags)
    Draw.Toggle ('E Tags',     EVT_NONE, c+160, r, 80, rh, d['show_etags'], 'Show edge tags'         , cb_show_etags)
    Draw.Toggle ('F Tags',     EVT_NONE, c+240, r, 80, rh, d['show_ftags'], 'Show face tags'         , cb_show_ftags)
    r -= rh
    r -= srg
    Draw.PushButton ('Delete all properties',  EVT_SET_DELPROPS, c,     r, 160, rh,           'Delete all properties')
    Draw.Toggle     ('Mark overlapping edges', EVT_NONE,         c+160, r, 160, rh, markover, 'Mark overlapping edges', cb_over)
    r, c, w = gu.box1_out(W,cg,rh,rg, c,r)
    r -= rh
    r0 = r

    # ======================================================== CAD

    gu.caption1(c,r0,w,rh,'CAD',c,d['gui_show_cad'],cb_gui_show_cad)
    if d['gui_show_cad']:
        r = r0-rh
        r, c, w = gu.box1_in(W,cg,rh,rg, c,r,w,h_cad)
        Draw.String     ('X=',           EVT_NONE,        c,     r, 80, rh, d['cad_x'],128,     'Set the X value of the new point to be added',cb_set_x)
        Draw.String     ('Y=',           EVT_NONE,        c+ 80, r, 80, rh, d['cad_y'],128,     'Set the Y value of the new point to be added',cb_set_y)
        Draw.String     ('Z=',           EVT_NONE,        c+160, r, 80, rh, d['cad_z'],128,     'Set the Z value of the new point to be added',cb_set_z)
        Draw.PushButton ('Add point',    EVT_CAD_ADDXYZ,  c+240, r, 80, rh,                     'Add point from coordinates')
        r-=rh
        Draw.String     ('Rad=',         EVT_NONE,        c,     r, 80, rh, d['cad_rad'],  128, 'Set radius for fillet operation',        cb_fillet_rad)
        Draw.Number     ('Stp=',         EVT_NONE,        c+ 80, r, 80, rh, d['cad_stp'],1,100, 'Set steps for fillet operation' ,        cb_fillet_stp)
        Draw.PushButton ('Fillet',       EVT_CAD_FILLET,  c+160, r, 80, rh,                     'Create a fillet between two edges')
        Draw.PushButton ('Read points',  EVT_CAD_FPOINT,  c+240, r, 80, rh,                     'Add points by reading a list from file')
        r-=rh
        Draw.PushButton ('Break',        EVT_CAD_BREAK ,  c,     r, 80, rh,                     'Break an edge at a previously selected point')
        Draw.PushButton ('Break at mid', EVT_CAD_BREAKM,  c+ 80, r, 80, rh,                     'Break an edge at its middle point')
        Draw.PushButton ('Edge int',     EVT_CAD_EINT,    c+160, r, 80, rh,                     'Find the intersection (smaller distance) between two edges')
        Draw.PushButton ('Read spline',  EVT_CAD_FSPLINE, c+240, r, 80, rh,                     'Add a spline by reading its points from file')
        r-=rh
        Draw.PushButton ('Export points',EVT_CAD_EXPTS,   c,     r, 80, rh,                     'Export points to a file')
        r, c, w = gu.box1_out(W,cg,rh,rg, c,r)
    r -= rh
    r -= rg

    # ======================================================== Mesh

    gu.caption1(c,r0,w,rh,'MESH',c+80,d['gui_show_mesh'],cb_gui_show_mesh)
    if d['gui_show_mesh']:
        r = r0-rh
        r, c, w = gu.box1_in(W,cg,rh,rg, c,r,w,h_msh)
        gu.text (c,r,'Vertex tags:')
        Draw.Number      ('',          EVT_NONE,             c+ 80, r, 80, rh, d['newvtag'],   -1000, 0, 'New vertex tag',                      cb_vtag)
        Draw.PushButton  ('Set V tag', EVT_MESH_SETVTAG,     c+160, r, 80, rh,                           'Set vertices tag (0 => remove tag)')
        Draw.PushButton  ('Del VTags', EVT_MESH_DELALLVTAGS, c+240, r, 80, rh, 'Delete all vertex tags')
        r -= rh
        gu.text (c,r,'Edge tags:')
        Draw.Number      ('',          EVT_NONE,             c+ 80, r, 80, rh, d['newetag'],   -1000, 0, 'New edge tag',                      cb_etag)
        Draw.PushButton  ('Set E tag', EVT_MESH_SETETAG,     c+160, r, 80, rh,                           'Set edges tag (0 => remove tag)')
        Draw.PushButton  ('Del ETags', EVT_MESH_DELALLETAGS, c+240, r, 80, rh, 'Delete all edge tags')
        r -= rh
        gu.text (c,r,'Face tags:')
        Draw.Number      ('',          EVT_NONE,             c+ 80, r, 80, rh, d['newftag'],   -1000, 0, 'New face tag',                      cb_ftag)
        Draw.PushButton  ('Set F tag', EVT_MESH_SETFTAG,     c+160, r, 80, rh,                           'Set faces tag (0 => remove tag)')
        Draw.PushButton  ('Del FTags', EVT_MESH_DELALLFTAGS, c+240, r, 80, rh, 'Delete all face tags')
        r -= rh
        gu.text (c,r,'Settings:')
        Draw.Toggle      ('3D mesh',     EVT_NONE, c+80,  r, 80, rh, is3d,                     'Set 3D mesh',                       cb_is3d)
        Draw.Toggle      ('Frame Mesh',  EVT_NONE, c+160, r, 80, rh, isframe,         'Set frame (truss/beams only) mesh', cb_frame)
        Draw.Toggle      ('O2 Elements', EVT_NONE, c+240, r, 80, rh, iso2,            'Generate quadratic (o2) elements' , cb_iso2)
        r -= srg
        r -= rh
        Draw.Toggle      ('Hide mesh',   EVT_NONE,         c,     r, 160, rh, d['hide_mesh'], 'Hide mesh', cb_hide_mesh)
        Draw.PushButton  ('Delete mesh', EVT_MESH_DELMESH, c+160, r, 160, rh, 'Delete current mesh')
        r -= rh
        r -= rh

        # ----------------------- Mesh -- structured

        gu.caption2(c,r,w,rh,'Structured mesh')
        r, c, w = gu.box2_in(W,cg,rh,rg, c,r,w,h_msh_stru)

        # ----------------------- Mesh -- structured -- blocks

        gu.caption3(c,r,w,rh, 'Blocks', EVT_MESH_ADDBLK,EVT_MESH_DELALLBLKS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_stru_blks)
        if len(blks)>0: gu.text(c,r,'   ID:Tag    Local axes      nX           nY          nZ')
        else: r += (rh+srg)
        for k, v in blks.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     (str(i)+':', EVT_INC+i, c     , r-rh, 60, 2*rh, int(v[0]), -1000, -1,'Block tag',                       cb_blk_tag)
            Draw.PushButton ('X',        EVT_INC+i, c+ 60 , r,    20,   rh,                      'Set X-axis',                      cb_blk_xax)
            Draw.PushButton ('Y',        EVT_INC+i, c+ 80 , r,    20,   rh,                      'Set Y-axis',                      cb_blk_yax)
            Draw.PushButton ('Z',        EVT_INC+i, c+100 , r,    20,   rh,                      'Set Z-axis',                      cb_blk_zax)
            Draw.Number     ('',         EVT_INC+i, c+120 , r,    60,   rh, int(v[ 8]), 1, 1000, 'Number of divisions along X',     cb_blk_nx)
            Draw.Number     ('',         EVT_INC+i, c+180 , r,    60,   rh, int(v[ 9]), 1, 1000, 'Number of divisions along Y',     cb_blk_ny)
            Draw.Number     ('',         EVT_INC+i, c+240 , r,    60,   rh, int(v[10]), 1, 1000, 'Number of divisions along Z',     cb_blk_nz)
            r -= rh
            Draw.Toggle     ('X',        EVT_INC+i, c+ 60 , r,    20,   rh, int(v[14]),          'Set nonlinear divisions along X', cb_blk_nlx)
            Draw.Toggle     ('Y',        EVT_INC+i, c+ 80 , r,    20,   rh, int(v[15]),          'Set nonlinear divisions along Y', cb_blk_nly)
            Draw.Toggle     ('Z',        EVT_INC+i, c+100 , r,    20,   rh, int(v[16]),          'Set nonlinear divisions along Z', cb_blk_nlz)
            Draw.String     ('aX=',      EVT_INC+i, c+120 , r,    60,   rh, str(v[11]),      32, 'Set division coefficient aX',     cb_blk_ax)
            Draw.String     ('aY=',      EVT_INC+i, c+180 , r,    60,   rh, str(v[12]),      32, 'Set division coefficient aY',     cb_blk_ay)
            Draw.String     ('aZ=',      EVT_INC+i, c+240 , r,    60,   rh, str(v[13]),      32, 'Set division coefficient aZ',     cb_blk_az)
            Draw.PushButton ('Del',      EVT_INC+i, c+300 , r,    40, 2*rh,                      'Delete this row',                 cb_blk_del)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- structured -- END

        r -= srg
        Draw.Toggle     ('C++',          EVT_NONE,          c,     r,  40, rh, d['smsh_cpp'], 'Generate C++ script', cb_smsh_cpp)
        Draw.PushButton ('Write Script', EVT_MESH_GENSTRUS, c+40,  r, 100, rh, 'Create script for structured mesh generation')
        Draw.PushButton ('Generate',     EVT_MESH_GENSTRU,  c+140, r, 120, rh, 'Generated structured mesh')
        r, c, w = gu.box2_out(W,cg,rh,rg, c,r)

        # ----------------------- Mesh -- unstructured

        r -= rh
        gu.caption2(c,r,w,rh,'Unstructured mesh')
        r, c, w = gu.box2_in(W,cg,rh,rg, c,r,w,h_msh_unst)
        gu.text     (c,r,'Quality:')
        Draw.String ('q=', EVT_NONE, c+ 50, r, 80, rh, '%g'%minang,  32, 'Set the minimum angle between edges/faces (-1 => use default)',                        cb_minang)
        Draw.String ('a=', EVT_NONE, c+130, r, 80, rh, '%g'%maxarea, 32, 'Set the maximum area/volume (uniform) for triangles/tetrahedrons (-1 => use default)', cb_maxarea)
        r -= rh
        r -= srg

        # ----------------------- Mesh -- unstructured -- regions

        gu.caption3(c,r,w,rh, 'Regions', EVT_MESH_ADDREG,EVT_MESH_DELALLREGS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_unst_regs)
        if len(regs)>0: gu.text(c,r,'  ID:Tag     Amax        X           Y          Z')
        else: r += (rh+srg)
        for k, v in regs.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     (str(i)+':', EVT_INC+i, c,     r, 60, rh, int(v[0]),-1000,-1,'Region tag',                  cb_regs_tag)
            Draw.String     ('',         EVT_INC+i, c+ 60, r, 50, rh, '%g'%v[1],   32,   'Max area (-1 => use default)',cb_regs_maxarea)
            Draw.String     ('',         EVT_INC+i, c+110, r, 50, rh, '%g'%v[2],   32,   'X of the region',             cb_regs_setx)
            Draw.String     ('',         EVT_INC+i, c+160, r, 50, rh, '%g'%v[3],   32,   'Y of the region',             cb_regs_sety)
            Draw.String     ('',         EVT_INC+i, c+210, r, 50, rh, '%g'%v[4],   32,   'Z of the region',             cb_regs_setz)
            Draw.PushButton ('Del',      EVT_INC+i, c+260, r, 30, rh,                    'Delete this row',             cb_regs_del)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- unstructured -- holes

        r -= srg
        gu.caption3(c,r,w,rh, 'Holes', EVT_MESH_ADDHOL,EVT_MESH_DELALLHOLS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_unst_hols)
        if len(hols)>0: gu.text(c,r,'   ID         X             Y             Z')
        else: r += (rh+srg)
        for k, v in hols.iteritems():
            r -= rh
            i  = int(k)
            gu.label        (str(i),            c    , r, 40, rh)
            Draw.String     ('',     EVT_INC+i, c+ 40, r, 60, rh, '%g'%v[0], 32, 'X of the hole',   cb_hols_setx)
            Draw.String     ('',     EVT_INC+i, c+100, r, 60, rh, '%g'%v[1], 32, 'Y of the hole',   cb_hols_sety)
            Draw.String     ('',     EVT_INC+i, c+160, r, 60, rh, '%g'%v[2], 32, 'Z of the hole',   cb_hols_setz)
            Draw.PushButton ('Del',  EVT_INC+i, c+220, r, 40, rh,                'Delete this row', cb_hols_del)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- unstructured -- END

        r -= srg
        Draw.Toggle     ('C++',          EVT_NONE,            c,     r,  40, rh, d['umsh_cpp'], 'Generate C++ script', cb_umsh_cpp)
        Draw.PushButton ('Write Script', EVT_MESH_GENUNSTRUS, c+40,  r, 100, rh, 'Create script for unstructured mesh generation')
        Draw.PushButton ('Generate',     EVT_MESH_GENUNSTRU , c+140, r, 120, rh, 'Generated unstructured mesh')
        r, c, w = gu.box2_out(W,cg,rh,rg, c,r+rh)

        # ----------------------- Mesh -- extra elements

        r -= rh
        r -= rh
        gu.caption2(c,r,w,rh,'Extra elements')
        r, c, w = gu.box2_in(W,cg,rh,rg, c,r,w,h_msh_ext)

        # ----------------------- Mesh -- extra elements -- linear elements

        gu.caption3(c,r,w,rh, 'Linear elements', EVT_MESH_ADDLINE,EVT_MESH_DELALLLINES)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_ext_lins)
        if len(lins)>0: gu.text(c,r,'     Tag                        N0                          N1')
        else: r += (rh+srg)
        for k, v in lins.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     ('',       EVT_INC+i, c,     r, 60, rh, v[0],-1000, 0, 'New linear element tag',   cb_mesh_letag)
            Draw.PushButton ('Get N0', EVT_INC+i, c+ 60, r, 60, rh,                'Get N0 of linear element', cb_mesh_getN0)
            Draw.Number     ('',       EVT_INC+i, c+120, r, 60, rh, v[1], 0, 1e8,  'N0: start point',          cb_mesh_leN0)
            Draw.PushButton ('Get N1', EVT_INC+i, c+180, r, 60, rh,                'Get N1 of linear element', cb_mesh_getN1)
            Draw.Number     ('',       EVT_INC+i, c+240, r, 60, rh, v[2], 0, 1e8,  'N1: end point',            cb_mesh_leN1)
            Draw.PushButton ('Del',    EVT_INC+i, c+300, r, 40, rh,                'Delete this row',          cb_mesh_ledel)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- extra elements -- linear elements -- END

        r, c, w = gu.box2_out(W,cg,rh,rg, c,r+rh)
        r, c, w = gu.box1_out(W,cg,rh,rg, c,r)

    r -= rh
    r -= rg

    # ======================================================== Materials

    gu.caption1(c,r0,w,rh,'MATERIALS',c+160,d['gui_show_mat'],cb_gui_show_mat)
    if d['gui_show_mat']:
        r = r0-rh
        r, c, w = gu.box1_in(W,cg,rh,rg, c,r,w,h_mat)

        # ----------------------- Mat -- parameters

        gu.caption2_(c,r,w,rh,'Materials', EVT_MAT_ADDMAT,EVT_MAT_DELALLMAT)
        r, c, w = gu.box2__in(W,cg,rh,rg, c,r,w,h_mat_mats)
        for k, v in mats.iteritems():
            r  -= rh
            i   = int(k)
            tid = int(v[1])       # text_id
            des = texts[str(tid)] # description
            Draw.Menu       (d['mdlmnu'], EVT_INC+i,   c,     r, 120, rh, int(v[0])+1, 'Constitutive model: ex.: LinElastic', cb_mat_setmodel)
            Draw.String     ('',          EVT_INC+tid, c+120, r, 120, rh, des, 32,     'Material name/description',           cb_mat_setname)
            Draw.PushButton ('Delete',    EVT_INC+i,   c+240, r,  60, rh,              'Delete this row',                     cb_mat_del)
            r -= rh
            if d['mdl'][int(v[0])]=='Rod': # Rod
                Draw.String ('E = ', EVT_INC+i, c,     r, 100, rh, '%g'%v[2], 32, 'E: Young modulus',        cb_mat_setE)
                Draw.String ('A = ', EVT_INC+i, c+100, r, 100, rh, '%g'%v[3], 32, 'A: Cross-sectional area', cb_mat_setA)
            elif d['mdl'][int(v[0])]=='Beam': # Beam
                Draw.String ('E = ',   EVT_INC+i, c,     r, 100, rh, '%g'%v[2], 32, 'E: Young modulus',             cb_mat_setE)
                Draw.String ('A = ',   EVT_INC+i, c+100, r, 100, rh, '%g'%v[3], 32, 'A: Cross-sectional area',      cb_mat_setA)
                Draw.String ('Izz = ', EVT_INC+i, c+200, r, 100, rh, '%g'%v[4], 32, 'Izz: Cross-sectional inertia', cb_mat_setIzz)
            elif d['mdl'][int(v[0])]=='LinElastic': # LinElastic
                Draw.String ('E = ',  EVT_INC+i, c,     r, 100, rh, '%g'%v[2], 32, 'E: Young modulus',  cb_mat_setE)
                Draw.String ('nu = ', EVT_INC+i, c+100, r, 100, rh, '%g'%v[5], 32, 'nu: Poisson ratio', cb_mat_setnu)
            elif d['mdl'][int(v[0])]=='ElastoPlastic': # ElastoPlastic
                Draw.String ('E = ',     EVT_INC+i, c,     r, 100, rh, '%g'%v[2], 32, 'E: Young modulus',      cb_mat_setE)
                Draw.String ('nu = ',    EVT_INC+i, c+100, r, 100, rh, '%g'%v[5], 32, 'nu: Poisson ratio',     cb_mat_setnu)
                Draw.Menu   (d['fcmnu'], EVT_INC+i, c+200, r, 100, rh, int(v[6])+1,   'fc: failure criterion', cb_mat_setfc)
                r -= rh
                Draw.String ('sY = ', EVT_INC+i, c,     r, 100, rh, '%d'%v[7], 32, 'sY: Yield stress',       cb_mat_setsY)
                Draw.String ('cu = ', EVT_INC+i, c+100, r, 100, rh, '%d'%v[8], 32, 'cu: Undrained cohesion', cb_mat_setcu)
            elif d['mdl'][int(v[0])]=='CamClay': # CamClay
                Draw.String ('lam = ', EVT_INC+i, c,     r, 100, rh, '%g'%v[ 9], 32, 'lam: Lambda',       cb_mat_setlam)
                Draw.String ('kap = ', EVT_INC+i, c+100, r, 100, rh, '%g'%v[10], 32, 'kap: Kappa',        cb_mat_setkap)
                Draw.String ('nu = ',  EVT_INC+i, c+200, r, 100, rh, '%g'%v[ 5], 32, 'nu: Poisson ratio', cb_mat_setnu)
                r -= rh
                Draw.String ('phi = ', EVT_INC+i, c, r, 100, rh, '%g'%v[11], 32, 'phi: Shear angle at critical state (compression/degrees)', cb_mat_setphi)
            if d['mdl'][int(v[0])]=='LinFlow': # LinFlow
                Draw.String ('k = ', EVT_INC+i, c, r, 100, rh, '%g'%v[12], 32, 'k: Diffusion coefficient', cb_mat_setk)
            r -= srg
        r -= srg
        r -= rh
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- Mat -- END
        r, c, w = gu.box1_out(W,cg,rh,rg, c,r+rh)
    else: r -= rh
    r -= rg

    # ======================================================== FEM

    gu.caption1(c,r0-rh,w,rh,'FEM',c,d['gui_show_fem'],cb_gui_show_fem)
    if d['gui_show_fem']:
        r = r0-rh
        r, c, w = gu.box1_in(W,cg,rh,rg, c,r,w,h_fem)
        gu.text (c,r,"Problem:")
        Draw.Menu (d['ptymnu'], EVT_NONE, c+80, r, 200, rh, pty+1, 'Problem type', cb_fem_prob)
        r -= rh
        r -= rg

        # ----------------------- FEM -- stages

        gu.caption2__(c,r,w,rh,'Stage #                  / %d'%(nstages),EVT_FEM_ADDSTAGE,EVT_FEM_DELSTAGE,EVT_FEM_DELALLSTAGES)
        if len(stages)>0:
            sid = str(d['fem_stage'])  # stage_id
            num = int(stages[sid][0])  # num
            Draw.Number ('', EVT_NONE, c+55, r+2, 60, rh-4, num, 0,99,'Show stage', cb_fem_stage)
            r, c, w = gu.box2_in(W,cg,rh,rg, c,r,w,h_fem_stage)

            # ----------------------- FEM -- eatts

            r -= srg
            if num==1: gu.caption3 (c,r,w,rh,'Elements attributes', EVT_FEM_ADDEATT,EVT_FEM_DELALLEATT)
            else:      gu.caption3_(c,r,w,rh,'Elements attributes')
            r, c, w = gu.box3__in(W,cg,rh, c,r,w,h_fem_eatts)
            if len(eatts)==0: r += (rh+srg)
            else:             r += srg
            if num>1:
                for k, v in obj.properties['stages'].iteritems():
                    if int(v[0])==1: # first stage
                        fstg = 'stg_'+k
                        break
                featts = obj.properties[fstg]['eatts'] if obj.properties[fstg].has_key('eatts') else {} # first stage eatts
            etags = [(int(v[0]),k) for k, v in eatts.iteritems()]
            #etags.sort(reverse=True)
            for tag, k in etags:
                v  = eatts[k]
                r -= rh
                i  = int(k)
                if num==1:
                    Draw.Number     ('',                  EVT_INC+i,   c,     r-rh, 60, 2*rh, int(v[0]),-1000,0, 'Set tag',                 cb_eatt_tag)
                    Draw.Menu       (d['pty2gpmnu'][pty], EVT_INC+i,   c+ 60, r,    70,   rh, int(v[2])+1,       'Type: Beam, Quad4, Hex8', cb_eatt_gebp)
                    Draw.Menu       (d['pty2gmnu'][pty],  EVT_INC+i,   c+130, r,    60,   rh, int(v[3])+1,       'Geometry type: psa, pse', cb_eatt_gty)
                    Draw.String     ('',                  EVT_INC+i,   c+190, r,    70,   rh, 'h=%g'%v[5],32,    'h: Thickness',            cb_eatt_h)
                    Draw.PushButton ('Del',               EVT_INC+i,   c+260, r-rh, 30, 2*rh,                    'Delete this row',         cb_eatt_del)
                    r -= rh
                    Draw.Menu       (matmnu,       EVT_INC+i, c+60,  r, 100, rh, int(v[4])+1,  'Choose material', cb_eatt_mat)
                    Draw.Toggle     ('Is Active',  EVT_INC+i, c+160, r, 100, rh, int(v[6]),    'Active',          cb_eatt_isact)
                else:
                    tag   =                str(int(featts[k][0]))
                    geom  = d['pty2geom'][pty][int(featts[k][2])]
                    gty   = d['pty2gty'] [pty][int(featts[k][3])]
                    mat   = matnames          [int(featts[k][4])] if matnames.has_key(int(featts[k][4])) else ''
                    thick =                    str(featts[k][5])
                    gu.label_   (tag,                     c,  r-rh,  60, 2*rh)
                    gu.label    (geom,                    c+ 60, r,  70,   rh)
                    gu.label    (gty,                     c+130, r,  60,   rh)
                    gu.label    (thick,                   c+190, r,  70,   rh)
                    r -= rh
                    gu.label    (mat,                     c+ 60, r,   70,  rh)
                    Draw.Toggle ('Activate',   EVT_INC+i, c+130, r,   60,  rh, int(v[7]), 'Activate this element at this stage?',   cb_eatt_act)
                    Draw.Toggle ('Deactivate', EVT_INC+i, c+190, r,   70,  rh, int(v[8]), 'Deactivate this element at this stage?', cb_eatt_deact)
                r -= srg
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- nbrys

            r -= rh
            gu.caption3(c,r,w,rh,'Nodes BCs', EVT_FEM_ADDNBRY,EVT_FEM_DELALLNBRY)
            r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_fem_nbrys)
            if len(nbrys)>0: gu.text(c,r,'       Tag         Key         Value')
            else: r += (rh+srg)
            for k, v in nbrys.iteritems():
                r -= rh
                i  = int(k)
                Draw.Number     ('',                 EVT_INC+i, c    , r, 80, rh, int(v[0]),-1000,-1,'Set tag',                     cb_nbry_tag)
                Draw.Menu       (d['pty2Nmnu'][pty], EVT_INC+i, c+ 80, r, 60, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz',  cb_nbry_key)
                Draw.String     ('',                 EVT_INC+i, c+140, r, 80, rh, str(v[2]), 128,    'Value of boundary condition', cb_nbry_val)
                Draw.PushButton ('Del',              EVT_INC+i, c+220, r, 40, rh,                    'Delete this row',             cb_nbry_del)
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- ebrys

            r -= rh
            gu.caption3(c,r,w,rh,'Edges BCs', EVT_FEM_ADDEBRY,EVT_FEM_DELALLEBRY)
            r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_fem_ebrys)
            if len(ebrys)>0: gu.text(c,r,'       Tag         Key         Value')
            else: r += (rh+srg)
            etags = [(int(v[0]),k) for k, v in ebrys.iteritems()]
            #etags.sort(reverse=True)
            for tag, k in etags:
                v  = ebrys[k]
                r -= rh
                i  = int(k)
                Draw.Number     ('',                 EVT_INC+i, c,     r, 80, rh, int(v[0]),-1000,-1,'Set tag',                     cb_ebry_tag)
                Draw.Menu       (d['pty2Fmnu'][pty], EVT_INC+i, c+ 80, r, 60, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz',  cb_ebry_key)
                Draw.String     ('',                 EVT_INC+i, c+140, r, 80, rh, str(v[2]), 128,    'Value of boundary condition', cb_ebry_val)
                Draw.PushButton ('Del',              EVT_INC+i, c+220, r, 40, rh,                    'Delete this row',             cb_ebry_del)
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- fbrys

            r -= rh
            gu.caption3(c,r,w,rh,'Faces BCs', EVT_FEM_ADDFBRY,EVT_FEM_DELALLFBRY)
            r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_fem_fbrys)
            if len(fbrys)>0: gu.text(c,r,'       Tag         Key         Value')
            else: r += (rh+srg)
            ftags = [(int(v[0]),k) for k, v in fbrys.iteritems()]
            #ftags.sort(reverse=True)
            for tag, k in ftags:
                v   = fbrys[k]
                r  -= rh
                i   = int(k)
                Draw.Number      ('',                 EVT_INC+i, c,     r, 80, rh, int(v[0]),-1000,-1,'Set tag',                    cb_fbry_tag)
                Draw.Menu        (d['pty2Fmnu'][pty], EVT_INC+i, c+ 80, r, 60, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz', cb_fbry_key)
                Draw.String      ('',                 EVT_INC+i, c+140, r, 80, rh, str(v[2]), 128,    'Value boundary condition',   cb_fbry_val)
                Draw.PushButton  ('Del',              EVT_INC+i, c+220, r, 40, rh,                    'Delete this row',            cb_fbry_del)
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- stages -- END

            i   = d['fem_stage']       # stage_id
            tid = int(stages[sid][1])  # text_id
            r  -= rh
            Draw.Toggle ('Active', EVT_INC+i,   c,    r,  80, rh,           int(stages[sid][5]),        'Set stage Active/Inactive during simulations?', cb_fem_stage_act)
            Draw.String ('Desc: ', EVT_INC+tid, c+80, r, 220, rh, texts[str(int(stages[sid][1]))], 256, 'Description of this stage',                     cb_fem_stage_desc)
            r -= rh
            Draw.Toggle ('Apply body forces',   EVT_INC+i, c,     r, 150, rh, int(stages[sid][2]),      'Apply body forces ?',                cb_fem_stage_apply_bf)
            Draw.Toggle ('Clear displacements', EVT_INC+i, c+150, r, 150, rh, int(stages[sid][3]),      'Clear displacements (and strains)?', cb_fem_stage_clear_disp)
            r -= rh
            Draw.Menu (d['stagetypemnu'], EVT_INC+i, c, r, 300, rh, int(stages[sid][6])+1, 'Stage type: equilib/steady, transient, dynamics',  cb_fem_stage_type)
            r -= rh
            r -= srg
            if int(stages[sid][6])==0: # equilib/steady
                gu.text(c,r,'Scheme:'); Draw.Menu (d['schememnu'], EVT_INC+i, c+60, r, 150, rh, int(stages[sid][12])+1, 'Global solution scheme',  cb_fem_stage_scheme)
                r -= rh
                gu.text(c,r,'NDiv:'); Draw.Number ('', EVT_INC+i, c+60, r, 60, rh, int(stages[sid][4]), 1,10000, 'Number of divisions', cb_fem_stage_ndiv)
            elif int(stages[sid][6])==1: # transient
                gu.text(c,    r,'tf:');    Draw.String ('', EVT_INC+i, c+60, r, 60, rh, '%g'%stages[sid][7], 128, 'Final time',                cb_fem_stage_tf);    r -= rh
                gu.text(c,    r,'dt:');    Draw.String ('', EVT_INC+i, c+60, r, 60, rh, '%g'%stages[sid][8], 128, 'Time increment',            cb_fem_stage_dt);    r -= rh
                gu.text(c,    r,'dtOut:'); Draw.String ('', EVT_INC+i, c+60, r, 60, rh, '%g'%stages[sid][9], 128, 'Time increment for output', cb_fem_stage_dtOut)
            elif int(stages[sid][6])==2: # dynamics
                gu.text(c, r,'Scheme:');  Draw.Menu (d['dschememnu'], EVT_INC+i, c+60, r, 120, rh, int(stages[sid][24])+1, 'Dynamics scheme', cb_fem_stage_dscheme)
                r -= rh
                gu.text(c, r,'Damping:'); Draw.Menu (d['damptymnu'],  EVT_INC+i, c+60, r, 120, rh, int(stages[sid][25])+1, 'Damping type',    cb_fem_stage_DampTy);
                r -= rh
                gu.text(c,     r,'DampAm:'); Draw.String ('', EVT_INC+i, c+60,  r, 60, rh, '%g'%stages[sid][26], 128, 'Damping coefficient Am (multiplies M matrix)', cb_fem_stage_DampAm);
                gu.text(c+130, r,'DynTh1:'); Draw.String ('', EVT_INC+i, c+190, r, 60, rh, '%g'%stages[sid][28], 128, 'Dynamic coefficient Theta1',                   cb_fem_stage_DynTh1);
                r -= rh
                gu.text(c,     r,'DampAk:'); Draw.String ('', EVT_INC+i, c+60,  r, 60, rh, '%g'%stages[sid][27], 128, 'Damping coefficient Ak (multiplies K matrix)', cb_fem_stage_DampAk);
                gu.text(c+130, r,'DynTh2:'); Draw.String ('', EVT_INC+i, c+190, r, 60, rh, '%g'%stages[sid][29], 128, 'Dynamic coefficient Theta2',                   cb_fem_stage_DynTh2);
                r -= rh
                gu.text(c,    r,'timef');
                gu.text(c+60, r,'dt');
                gu.text(c+120,r,'dtout');
                r -= rh
                Draw.String('',       EVT_INC+i, c,     r, 60, rh, '%g'%stages[sid][7], 128, 'Final time',                  cb_fem_stage_tf)
                Draw.String('',       EVT_INC+i, c+60,  r, 60, rh, '%g'%stages[sid][8], 128, 'Time increment',              cb_fem_stage_dt)
                Draw.String('',       EVT_INC+i, c+120, r, 60, rh, '%g'%stages[sid][9], 128, 'Time increment for output',   cb_fem_stage_dtOut)
                Draw.Toggle('Output', EVT_INC+i, c+180, r, 60, rh, int(stages[sid][31]),     'Output VTU file each dtOut?', cb_fem_stage_output)
            r -= rh
            r -= srg

            r, c, w = gu.box2_out(W,cg,rh,rg, c,r+rh)

        # ----------------------- FEM -- END
        r -= rh
        Draw.Toggle     ('C++',             EVT_NONE,       c,     r, 60,  rh, d['fem_cpp'], 'Generate C++ script instead of Python ?', cb_fem_cpp)
        Draw.PushButton ('Generate Script', EVT_FEM_SCRIPT, c+60,  r, 120, rh, 'Generate Script for FE simulation')
        Draw.PushButton ('Run',             EVT_FEM_RUN,    c+180, r, 60,  rh, 'Run simulation')
        Draw.PushButton ('Stop',            EVT_FEM_STOP,   c+240, r, 60,  rh, 'Stop simulation')
        if d['fem_running'].value:
            r -= rh
            gu.text (c, r, 'Simulation running ....................................................')
        if d['fem_fatal'].value:
            r -= rh
            gu.text (c, r, 'Simulation FAILED (check terminal)')
            Draw.PushButton ('clr', EVT_FEM_CLEAR, c+240, r, 60,  rh, 'Clear error message')
        r -= rh
        Draw.PushButton ('ParaView: stages',     EVT_FEM_PARAVIEW,  c,     r, 150, rh, 'View results in ParaView for the end of stages')
        Draw.PushButton ('ParaView: transient',  EVT_FEM_DPARAVIEW, c+150, r, 150, rh, 'View transient results in ParaView')

        r, c, w = gu.box1_out(W,cg,rh,rg, c,r)
    r -= rg

    # ======================================================== DEM

    gu.caption1(c,r0-rh,w,rh,'DEM',c+80,d['gui_show_dem'],cb_gui_show_dem)
    if d['gui_show_dem']:
        r = r0-rh
        r, c, w = gu.box1_in(W,cg,rh,rg, c,r,w,h_dem)

        # ----------------------- DEM -- particle generation

        gu.caption2(c,r,w,rh,'Particle generation')
        r, c, w = gu.box2_in(W,cg,rh,rg, c,r,w,h_dem_pkg)

        gu.text(c,    r,'Lx:'); Draw.String('', EVT_NONE, c+40,  r, 60, rh, str(d['dem_Lx']), 128, 'Length of box', cb_dem_Lx)
        gu.text(c+120,r,'R:');  Draw.String('', EVT_NONE, c+160, r, 60, rh, str(d['dem_R']),  128, 'Radius',        cb_dem_R)
        r -= rh
        gu.text(c,    r,'Ly:');  Draw.String('', EVT_NONE, c+40,  r, 60, rh, str(d['dem_Ly']),  128, 'Length of box', cb_dem_Ly)
        gu.text(c+120,r,'rho:'); Draw.String('', EVT_NONE, c+160, r, 60, rh, str(d['dem_rho']), 128, 'Density',       cb_dem_rho)
        r -= rh
        gu.text(c,    r,'Lz:');   Draw.String('', EVT_NONE, c+40,  r, 60, rh, str(d['dem_Lz']),  128,   'Length of box', cb_dem_Lz)
        gu.text(c+120,r,'Seed:'); Draw.Number('', EVT_NONE, c+160, r, 60, rh, int(d['dem_seed']),1,1000,'Seed',          cb_dem_seed)
        r -= rh
        gu.text(c,    r,'Nx:');   Draw.Number('', EVT_NONE, c+40,  r, 60, rh, int(d['dem_Nx']),1,1000,'Number of x divisions', cb_dem_Nx)
        gu.text(c+120,r,'Prob:'); Draw.String('', EVT_NONE, c+160, r, 60, rh, str(d['dem_prob']), 128,'Probability/fraction',  cb_dem_prob)
        r -= rh
        gu.text(c,r,'Ny:'); Draw.Number('', EVT_NONE,  c+40,  r, 60,  rh, int(d['dem_Ny']),1,1000,'Number of y divisions', cb_dem_Ny)
        Draw.Menu (d['dem_pkgs_mnu'],       EVT_NONE,  c+120, r, 100, rh, int(d['dem_pkg'])+1,    'Packing type',          cb_dem_pkg)
        r -= rh
        gu.text(c,r,'Nz:'); Draw.Number('', EVT_NONE,  c+40, r, 60, rh, int(d['dem_Nz']),1,1000,'Number of z divisions', cb_dem_Nz)
        r -= rh+rg
        gu.text(c,r,'Resolution:'); Draw.Number('', EVT_NONE, c+80,  r, 60, rh, int(d['dem_res']),3,128,'Resolution',    cb_dem_res)
        Draw.Toggle     ('Draw Verts',              EVT_NONE, c+140, r, 80, rh, d['dem_draw_verts'],    'Draw Vertices', cb_dem_draw_verts)
        Draw.Toggle     ('Draw Edges',              EVT_NONE, c+220, r, 80, rh, d['dem_draw_edges'],    'Draw Edges',    cb_dem_draw_edges)
        r -= rh
        Draw.PushButton ('Generate Particles',   EVT_DEM_GEN_PKG, c,     r, 150, rh, 'Generate particles')
        Draw.PushButton ('Generate TTT Packing', EVT_DEM_GEN_TTT, c+150, r, 150, rh, 'Generate true triaxial (TTT) packing')

        r, c, w = gu.box2_out(W,cg,rh,rg, c,r)

        # ----------------------- DEM -- constants

        r -= rh
        gu.caption2(c,r,w,rh,'Internal constants')
        r, c, w = gu.box2_in(W,cg,rh,rg, c,r,w,h_dem_cte)

        gu.text(c,     r,'Kn:'  ); Draw.String('', EVT_NONE, c+40,  r, 60, rh, str(d['dem_Kn']),   128, 'Normal stiffness between particles',                   cb_dem_Kn)
        gu.text(c+120, r,'mu:'  ); Draw.String('', EVT_NONE, c+160, r, 60, rh, str(d['dem_mu']),   128, 'Microscopic friction coefficient',                     cb_dem_mu)
        r -= rh
        gu.text(c,     r,'Kt:'  ); Draw.String('', EVT_NONE, c+40,  r, 60, rh, str(d['dem_Kt']),   128, 'Tangential stiffness between particles',               cb_dem_Kt)
        gu.text(c+120, r,'beta:'); Draw.String('', EVT_NONE, c+160, r, 60, rh, str(d['dem_beta']), 128, 'Spheres only: rolling stiffness coefficient',          cb_dem_beta)
        r -= rh
        gu.text(c,     r,'Gn:'  ); Draw.String('', EVT_NONE, c+40,  r, 60, rh, str(d['dem_Gn']),   128, 'Normal dissipative coefficient between particles',     cb_dem_Gn)
        gu.text(c+120, r,'eta:' ); Draw.String('', EVT_NONE, c+160, r, 60, rh, str(d['dem_eta']),  128, 'Spheres only: plastic moment coefficient (0 indicates that the rolling resistance will not be used)', cb_dem_eta)
        r -= rh
        gu.text(c,     r,'Gt:'  ); Draw.String('', EVT_NONE, c+40, r, 60, rh, str(d['dem_Gt']),   128, 'Tangential dissipative coefficient between particles', cb_dem_Gt)
        r -= rh

        r, c, w = gu.box2_out(W,cg,rh,rg, c,r)

        # ----------------------- DEM -- True Triaxial Test

        gu.caption2(c,r,w,rh,'True Triaxial Test (TTT)')
        r, c, w = gu.box2_in(W,cg,rh,rg, c,r,w,h_dem_ttt)

        # ----------------------- DEM -- True Triaxial Test --- Stage 1

        gu.caption3_(c,r,w,rh, 'Stage 1: Isotropic Compression')
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_dem_ttt_iso)
        r -= srg

        gu.text(c, r,'pf:'); Draw.String('', EVT_NONE, c+30, r, 60, rh, str(d['dem_iso_pf']), 128, 'p at the and of isotropic compression', cb_dem_iso_pf)
        r -= rh
        gu.text(c,    r,'timef');
        gu.text(c+60, r,'dt');
        gu.text(c+120,r,'dtout');
        r -= rh
        Draw.String('', EVT_NONE, c,    r, 60, rh,     str(d['dem_iso_timef']), 128, 'final time',                cb_dem_iso_timef)
        Draw.String('', EVT_NONE, c+60, r, 60, rh,     str(d['dem_iso_dt']),    128, 'time increment',            cb_dem_iso_dt)
        Draw.String('', EVT_NONE, c+120,r, 60, rh,     str(d['dem_iso_dtout']), 128, 'time increment for output', cb_dem_iso_dtout)
        Draw.Toggle ('Render', EVT_NONE, c+180, r, 60, rh, d['dem_iso_render'],      'Render video ?',            cb_dem_iso_render)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- DEM -- True Triaxial Test --- Stage 2

        r -= rh
        gu.caption3_(c,r,w,rh, 'Stage 2: Shear Loading')
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_dem_ttt_she)
        r -= srg

        gu.text(c,     r,'pf:');  Draw.String('', EVT_NONE, c+30,  r, 60, rh, str(d['dem_ttt_pf']),  128, 'p at the and of shearing',  cb_dem_ttt_pf)
        gu.text(c+110, r,'Exf:'); Draw.String('', EVT_NONE, c+150, r, 60, rh, str(d['dem_ttt_exf']), 128, 'Final Ex (strain)',         cb_dem_ttt_exf)
        Draw.Toggle     ('pEx:',                  EVT_NONE, c+210, r, 60, rh,     d['dem_ttt_pex'],       'Prescribed Ex (strain)',    cb_dem_ttt_pex)
        r -= rh
        gu.text(c,     r,'qf:');  Draw.String('', EVT_NONE, c+30,  r, 60, rh, str(d['dem_ttt_qf']),  128, 'q at the and of shearing',  cb_dem_ttt_qf)
        gu.text(c+110, r,'Eyf:'); Draw.String('', EVT_NONE, c+150, r, 60, rh, str(d['dem_ttt_eyf']), 128, 'Final Ey (strain)',         cb_dem_ttt_eyf)
        Draw.Toggle     ('pEy:',                  EVT_NONE, c+210, r, 60, rh,     d['dem_ttt_pey'],       'Prescribed Ey (strain)',    cb_dem_ttt_pey)
        r -= rh
        gu.text(c,     r,'thf:'); Draw.String('', EVT_NONE, c+30,  r, 60, rh, str(d['dem_ttt_thf']), 128, 'th at the and of shearing', cb_dem_ttt_thf)
        gu.text(c+110, r,'Ezf:'); Draw.String('', EVT_NONE, c+150, r, 60, rh, str(d['dem_ttt_ezf']), 128, 'Final Ez (strain)',         cb_dem_ttt_ezf)
        Draw.Toggle     ('pEz:',                  EVT_NONE, c+210, r, 60, rh,     d['dem_ttt_pez'],       'Prescribed Ez (strain)',    cb_dem_ttt_pez)
        r -= rh
        gu.text(c,    r,'timef');
        gu.text(c+60, r,'dt');
        gu.text(c+120,r,'dtout');
        r -= rh
        Draw.String('', EVT_NONE, c,    r, 60, rh,     str(d['dem_ttt_timef']), 128, 'final time',                cb_dem_ttt_timef)
        Draw.String('', EVT_NONE, c+60, r, 60, rh,     str(d['dem_ttt_dt']),    128, 'time increment',            cb_dem_ttt_dt)
        Draw.String('', EVT_NONE, c+120,r, 60, rh,     str(d['dem_ttt_dtout']), 128, 'time increment for output', cb_dem_ttt_dtout)
        Draw.Toggle ('Render', EVT_NONE, c+180, r, 60, rh, d['dem_ttt_render'],      'Render video ?',            cb_dem_ttt_render)
        r -= rh

        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- DEM -- True Triaxial Test --- Stage 2 -- END

        Draw.Toggle     ('C++',             EVT_NONE,        c,     r, 60,  rh, d['dem_cpp_script'], 'Generate C++ script instead of Python ?', cb_dem_cpp_script)
        Draw.PushButton ('Generate Script', EVT_DEM_GEN_SCR, c+60,  r, 120, rh, 'Generate Script')
        Draw.PushButton ('Run',             EVT_DEM_RUN,     c+180, r, 60,  rh, 'Run simulation')
        Draw.PushButton ('Stop',            EVT_DEM_STOP,    c+240, r, 60,  rh, 'Stop simulation')
        if d['dem_running'].value:
            r -= rh
            gu.text (c, r, 'Simulation running ....................................................')
        if d['dem_fatal'].value:
            r -= rh
            gu.text (c, r, 'Simulation FAILED (check terminal)')
            Draw.PushButton ('clr', EVT_DEM_CLEAR, c+240, r, 60,  rh, 'Clear error message')

        r, c, w = gu.box2_out(W,cg,rh,rg, c,r)
        r, c, w = gu.box1_out(W,cg,rh,rg, c,r)
    r -= rh
    r -= rg


# Register GUI
Draw.Register (gui, event, button_event)


# ================================================================================ ScriptLink

# Load script for View3D drawing
msdict     = di.load_dict()
script_lnk = 'msys_3dlink.py'
script_dir = Blender.Get('scriptsdir')+Blender.sys.sep
if script_lnk in [t.name for t in Blender.Text.Get()]:
    Blender.Text.unlink(Blender.Text.Get(script_lnk))
print '[1;34mMechSys[0m: loading <'+script_dir+script_lnk+'>'
Blender.Text.Load(script_dir+script_lnk)


# Connect View3D drawing script
scn = bpy.data.scenes.active
scn.clearScriptLinks()
if scn.getScriptLinks('Redraw')==None:
    scn.addScriptLink(script_lnk,'Redraw')
elif script_lnk not in scn.getScriptLinks('Redraw'):
    scn.addScriptLink(script_lnk,'Redraw')


# ============================================================================== SpaceHandler

# Load script for View3D space handler
script_shandler = 'msys_shandler.py'
if script_shandler in [t.name for t in Blender.Text.Get()]:
    Blender.Text.unlink(Blender.Text.Get(script_shandler))
print '[1;34mMechSys[0m: loading <'+script_dir+script_shandler+'>'
Blender.Text.Load(script_dir+script_shandler)

# TODO: Connect Space Handler to View3D
