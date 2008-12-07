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
Blender: 2.46
Group: 'Themes'
Tooltip: 'CAD, Mesh generator, and FEM in Blender'
"""
__author__  = "Dorival Pedroso"
__version__ = "06.05.08"
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
import msys_fem  as fe
import msys_gui  as gu
import msys_res  as re

#import rpdb2; rpdb2.start_embedded_debugger('msys')

# ================================================================================= Constants

EVT_INC              = 10000 # increment used when passing IDs through events
EVT_NONE             =  0 # used for buttons with callbacks
EVT_REFRESH          =  1 # refresh all windows
# SETTINGS
EVT_SET_SHOWHIDE     = 10 # show/hide SET box
EVT_SET_SHOWONLY     = 11 # show only this box
EVT_SET_DELPROPS     = 12 # delete all properties
EVT_SET_HIDEALL      = 13 # hide all boxes
# CAD
EVT_CAD_SHOWHIDE     = 20 # show/hide CAD box
EVT_CAD_SHOWONLY     = 21 # show/hide CAD box
EVT_CAD_ADDXYZ       = 22 # type in values to define point(s)
EVT_CAD_FILLET       = 23 # fillet two edges
EVT_CAD_BREAK        = 24 # break edge
EVT_CAD_BREAKM       = 25 # break edge at middle point
EVT_CAD_EINT         = 26 # edge closest distance
EVT_CAD_FPOINT       = 27 # read points from file
EVT_CAD_FSPLINE      = 28 # create a spline from points in a file
# Mesh
EVT_MESH_SHOWHIDE    = 40 # show/hide MESH box
EVT_MESH_SHOWONLY    = 41 # show only this box
EVT_MESH_SETETAG     = 42 # set edges tag
EVT_MESH_SETFTAG     = 43 # set faces tag
EVT_MESH_SETRTAG     = 44 # set reinforcement
EVT_MESH_DELALLETAGS = 45 # delete all edge tags
EVT_MESH_DELALLFTAGS = 46 # delete all face tags
EVT_MESH_DELALLRTAGS = 47 # delete all reinforcement tags
EVT_MESH_DELMESH     = 48 # delete mesh object
# Mesh -- structured
EVT_MESH_ADDBLK      = 60 # set block 2D: 4 or 8 edges, 3D: 8 or 20 edges
EVT_MESH_DELALLBLKS  = 61 # set block 2D: 4 or 8 edges, 3D: 8 or 20 edges
EVT_MESH_GENSTRU     = 62 # generate structured mesh using MechSys module
EVT_MESH_GENSTRUS    = 63 # script for structured mesh generation
# Mesh -- unstructured
EVT_MESH_ADDREG      = 70 # add region
EVT_MESH_DELALLREGS  = 71 # delete all regions
EVT_MESH_ADDHOL      = 72 # add hole
EVT_MESH_DELALLHOLS  = 73 # delete all holes
EVT_MESH_GENUNSTRU   = 74 # generate unstructured mesh using MechSys module
EVT_MESH_GENUNSTRUS  = 75 # script for unstructured mesh generation
# Materials
EVT_MAT_SHOWHIDE     = 80 # show/hide CAD box
EVT_MAT_SHOWONLY     = 81 # show/hide CAD box
EVT_MAT_ADDMAT       = 82 # add material
EVT_MAT_DELALLMAT    = 83 # delete all materials
# FEM
EVT_FEM_SHOWHIDE     =  90 # show/hide FEM box
EVT_FEM_SHOWONLY     =  91 # show only this box
EVT_FEM_ADDSTAGE     =  92 # add stage
EVT_FEM_DELSTAGE     =  93 # add stage
EVT_FEM_DELALLSTAGES =  94 # add stage
EVT_FEM_ADDNBRY      =  95 # add nodes boundary (given coordinates)
EVT_FEM_ADDNBID      =  96 # add nodes boundary (given nodes IDs)
EVT_FEM_ADDEBRY      = 100 # add edges boundary
EVT_FEM_ADDFBRY      = 101 # add faces boundary
EVT_FEM_ADDEATT      = 102 # add element attributes
EVT_FEM_DELALLNBRY   = 103 # delete all nodes boundary
EVT_FEM_DELALLNBID   = 104 # delete all nodes boundary (IDs)
EVT_FEM_DELALLEBRY   = 105 # delete all edges boundary
EVT_FEM_DELALLFBRY   = 106 # delete all faces boundary
EVT_FEM_DELALLEATT   = 200 # delete all element attributes
EVT_FEM_RUN          = 201 # run a FE simulation
EVT_FEM_SCRIPT       = 202 # generate script for FEM 
EVT_FEM_PARAVIEW     = 203 # view in ParaView
EVT_FEM_SAVESTAGES   = 204 # save stage info to a new object
EVT_FEM_READSTAGES   = 205 # read stage info from another object
# Results
EVT_RES_SHOWHIDE     = 300 # show/hide results box
EVT_RES_SHOWONLY     = 301 # show only this box
EVT_RES_STATS        = 302 # show statistics
EVT_RES_REPORT       = 303 # show statistics


# ==================================================================================== Events

def hide_all():
    di.set_key('gui_show_set',  False)
    di.set_key('gui_show_cad',  False)
    di.set_key('gui_show_mesh', False)
    di.set_key('gui_show_mat',  False)
    di.set_key('gui_show_fem',  False)
    di.set_key('gui_show_res',  False)
    di.set_key('gui_inirow',    0)
    Blender.Window.QRedrawAll()

def show_only(key):
    di.set_key('gui_show_set',  False)
    di.set_key('gui_show_cad',  False)
    di.set_key('gui_show_mesh', False)
    di.set_key('gui_show_mat',  False)
    di.set_key('gui_show_fem',  False)
    di.set_key('gui_show_res',  False)
    di.set_key(key,             True)
    di.set_key('gui_inirow',    0)
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
            Blender.Draw.PupMenu('ERROR|'+msg)
    return wrapper

# Handle button events
@try_catch
def button_event(evt):
    if evt==EVT_REFRESH: Blender.Window.QRedrawAll()

    # ----------------------------------------------------------------------------------- Settings

    if evt==EVT_SET_SHOWHIDE: di.toggle_key('gui_show_set')
    if evt==EVT_SET_SHOWONLY: show_only    ('gui_show_set')

    elif evt==EVT_SET_DELPROPS:
        scn = bpy.data.scenes.active
        obs = scn.objects.selected
        if len(obs)>0:
            message = 'Confirm delete ALL properties?%t|Yes'
            result  = Blender.Draw.PupMenu(message)
            if result>0:
                for o in obs:
                    for k, v in o.properties.iteritems():
                        print '[1;34mMechSys[0m: Object = ', o.name, ' [1;35mDeleted property:[0m ', k, str(v)
                        o.properties.pop(k)
                    for p in o.getAllProperties():
                        print '[1;34mMechSys[0m: Object = ', o.name, ' Deleted property: ', p.name, p.data
                        o.removeProperty(p)
                Blender.Window.QRedrawAll()

    elif evt==EVT_SET_HIDEALL: hide_all()

    # ----------------------------------------------------------------------------------- CAD

    elif evt==EVT_CAD_SHOWHIDE: di.toggle_key('gui_show_cad')
    elif evt==EVT_CAD_SHOWONLY: show_only    ('gui_show_cad')
    elif evt==EVT_CAD_ADDXYZ:   ca.add_point     (float(di.key('cad_x')), float(di.key('cad_y')), float(di.key('cad_z')))
    elif evt==EVT_CAD_FILLET:   ca.fillet        (float(di.key('cad_rad')), di.key('cad_stp'))
    elif evt==EVT_CAD_BREAK:    ca.break_edge    ()
    elif evt==EVT_CAD_BREAKM:   ca.break_edge    (True)
    elif evt==EVT_CAD_EINT:     ca.edge_intersect()
    elif evt==EVT_CAD_FPOINT:   Blender.Window.FileSelector(ca.add_points_from_file, 'Read X Y Z columns')
    elif evt==EVT_CAD_FSPLINE:  Blender.Window.FileSelector(ca.add_spline_from_file, 'Read X Y Z columns')

    # ---------------------------------------------------------------------------------- Mesh

    elif evt==EVT_MESH_SHOWHIDE: di.toggle_key('gui_show_mesh')
    elif evt==EVT_MESH_SHOWONLY: show_only    ('gui_show_mesh')

    elif evt==EVT_MESH_SETETAG:
        tag = di.key('newetag')[0]
        edm, obj, msh = di.get_msh()
        for eid in di.get_selected_edges(msh):
            di.props_set_with_tag ('etags', str(eid), tag, di.key('newetag'))
        Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1)

    elif evt==EVT_MESH_SETFTAG:
        tag = di.key('newftag')[0]
        edm, obj, msh = di.get_msh()
        sel           = di.get_selected_edges(msh)
        nedges        = len(sel)
        if nedges==3 or nedges==6 or nedges==4 or nedges==8:
            eids = '_'.join([str(eid) for eid in sel])
            di.props_set_with_tag ('ftags', eids, tag, di.key('newftag'))
        else: raise Exception('To set a face tag (FTAG), 3, 4, 6, or 8 edges must be selected')
        Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1)

    elif evt==EVT_MESH_SETRTAG:
        tag = di.key('newrtag')[0]
        edm, obj, msh = di.get_msh()
        for eid in di.get_selected_edges(msh):
            di.props_set_with_tag ('rtags', str(eid), tag, di.key('newrtag'))
        Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1)

    elif evt==EVT_MESH_DELALLETAGS: di.props_del_all_tags('etags')
    elif evt==EVT_MESH_DELALLFTAGS: di.props_del_all_tags('ftags')
    elif evt==EVT_MESH_DELALLRTAGS: di.props_del_all_tags('rtags')
    elif evt==EVT_MESH_DELMESH:
        obj = di.get_obj()
        if obj.properties.has_key('msh_name'):
            msg = 'Confirm delete mesh?%t|Yes'
            res = Blender.Draw.PupMenu(msg)
            if res>0:
                old_msh = bpy.data.objects[obj.properties['msh_name']]
                scn     = bpy.data.scenes.active
                scn.objects.unlink (old_msh)
                obj.properties.pop('mesh_type')
                obj.properties.pop('msh_name')
                obj.properties.pop('elems')
                if obj.properties.has_key('res'): obj.properties.pop('res')
                Blender.Window.QRedrawAll()


    # ---------------------------------------------------------------------------------- Materials

    elif evt==EVT_MAT_SHOWHIDE:  di.toggle_key        ('gui_show_mat')
    elif evt==EVT_MAT_SHOWONLY:  show_only            ('gui_show_mat')
    elif evt==EVT_MAT_ADDMAT:    di.props_push_new_mat()
    elif evt==EVT_MAT_DELALLMAT: di.props_del_all_mats()

    # -------------------------------------------------------------------- Mesh -- structured

    elif evt==EVT_MESH_ADDBLK:
        props = di.new_blk_props()
        di.props_push_new('blks', props, True, 18, 18+int(props[17]))

    elif evt==EVT_MESH_DELALLBLKS: di.props_del_all('blks')
    elif evt==EVT_MESH_GENSTRU:
        obj  = di.get_obj()
        mesh = me.gen_struct_mesh()
        me.add_mesh (obj, mesh, 'struct')
    elif evt==EVT_MESH_GENSTRUS:
        obj = di.get_obj()
        txt = me.gen_struct_mesh(True)
        txt.write('obj = bpy.data.objects["'+obj.name+'"]\n')
        txt.write('me.add_mesh     (obj, mesh, "struct")\n')

    # -------------------------------------------------------------------------- Unstructured

    elif evt==EVT_MESH_ADDREG:     di.props_push_new('regs', di.new_reg_props())
    elif evt==EVT_MESH_DELALLREGS: di.props_del_all ('regs')
    elif evt==EVT_MESH_ADDHOL:     di.props_push_new('hols', di.new_hol_props())
    elif evt==EVT_MESH_DELALLHOLS: di.props_del_all ('hols')
    elif evt==EVT_MESH_GENUNSTRU:
        obj  = di.get_obj()
        mesh = me.gen_unstruct_mesh()
        me.add_mesh (obj, mesh, 'unstruct')
    elif evt==EVT_MESH_GENUNSTRUS:
        obj = di.get_obj()
        txt = me.gen_unstruct_mesh(True)
        txt.write('obj = bpy.data.objects["'+obj.name+'"]\n')
        txt.write('me.add_mesh           (obj, mesh, "unstruct")\n')

    # ----------------------------------------------------------------------------------- FEM 

    elif evt==EVT_FEM_SHOWHIDE: di.toggle_key('gui_show_fem')
    elif evt==EVT_FEM_SHOWONLY: show_only    ('gui_show_fem')

    elif evt==EVT_FEM_ADDSTAGE: di.props_push_new_stage()
    elif evt==EVT_FEM_DELSTAGE: di.props_del_stage     ()
    elif evt==EVT_FEM_ADDNBRY:  di.props_push_new_fem  ([di.key('fem_stage')], 'nbrys', di.new_nbry_props())
    elif evt==EVT_FEM_ADDNBID:  di.props_push_new_fem  ([di.key('fem_stage')], 'nbsID', di.new_nbID_props())
    elif evt==EVT_FEM_ADDEBRY:  di.props_push_new_fem  ([di.key('fem_stage')], 'ebrys', di.new_ebry_props())
    elif evt==EVT_FEM_ADDFBRY:  di.props_push_new_fem  ([di.key('fem_stage')], 'fbrys', di.new_fbry_props())
    elif evt==EVT_FEM_ADDEATT:
        obj  = di.get_obj()
        sids = [k for k in obj.properties['stages']]
        di.props_push_new_fem (sids, 'eatts', di.new_eatt_props())

    elif evt==EVT_FEM_DELALLSTAGES: di.props_del_all_stages()
    elif evt==EVT_FEM_DELALLNBRY:   di.props_del_all_fem   (str(di.key('fem_stage')), 'nbrys')
    elif evt==EVT_FEM_DELALLNBID:   di.props_del_all_fem   (str(di.key('fem_stage')), 'nbsID')
    elif evt==EVT_FEM_DELALLEBRY:   di.props_del_all_fem   (str(di.key('fem_stage')), 'ebrys')
    elif evt==EVT_FEM_DELALLFBRY:   di.props_del_all_fem   (str(di.key('fem_stage')), 'fbrys')
    elif evt==EVT_FEM_DELALLEATT:
        obj  = di.get_obj()
        sids = [int(k) for k, v in obj.properties['stages'].iteritems()]
        di.props_del_all_fem (sids, 'eatts')

    elif evt==EVT_FEM_RUN:        fe.run_analysis         ()
    elif evt==EVT_FEM_SCRIPT:     fe.run_analysis         (True)
    elif evt==EVT_FEM_PARAVIEW:   fe.paraview             ()
    elif evt==EVT_FEM_SAVESTAGES: fe.save_stage_mats_info ()
    elif evt==EVT_FEM_READSTAGES: fe.read_stage_mats_info ()

    # ----------------------------------------------------------------------------------- RES 

    elif evt==EVT_RES_SHOWHIDE: di.toggle_key  ('gui_show_res')
    elif evt==EVT_RES_SHOWONLY: show_only      ('gui_show_res')
    elif evt==EVT_RES_STATS:    re.stage_stats (di.key('res_stage'))
    elif evt==EVT_RES_REPORT:   re.report      ()


# ================================================================================= Callbacks

# ---------------------------------- Settings

@try_catch
def cb_show_props(evt,val): di.set_key ('show_props', val)
@try_catch
def cb_show_e_ids(evt,val): di.set_key ('show_e_ids', val)
@try_catch
def cb_show_v_ids(evt,val): di.set_key ('show_v_ids', val)
@try_catch
def cb_show_n_ids(evt,val): di.set_key ('show_n_ids', val)
@try_catch
def cb_show_blks (evt,val): di.set_key ('show_blks',  val)
@try_catch
def cb_show_axes (evt,val): di.set_key ('show_axes',  val)
@try_catch
def cb_show_regs (evt,val): di.set_key ('show_regs',  val)
@try_catch
def cb_show_etags(evt,val): di.set_key ('show_etags', val)
@try_catch
def cb_show_ftags(evt,val): di.set_key ('show_ftags', val)
@try_catch
def cb_show_elems(evt,val): di.set_key ('show_elems', val)
@try_catch
def cb_show_opac (evt,val): di.set_key ('show_opac',  val)
@try_catch
def cb_show_rtags(evt,val): di.set_key ('show_rtags', val)

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
def cb_iso2(evt,val): di.props_set_val('iso2', val)
@try_catch
def cb_frame (evt,val):
    if val==1: di.props_set_val('mesh_type', 'frame')
    else:
        obj = di.get_obj()
        if obj.properties.has_key('mesh_type'): obj.properties.pop('mesh_type')
        Blender.Window.QRedrawAll()
@try_catch
def cb_etag  (evt,val): di.set_key      ('newetag', [val, di.key('newetag')[1]])
@try_catch
def cb_ftag  (evt,val): di.set_key      ('newftag', [val, di.key('newftag')[1]])
@try_catch
def cb_fclr  (evt,val): di.set_key      ('newftag', [di.key('newftag')[0], di.rgb2hex(val)])
@try_catch
def cb_rtag  (evt,val): di.set_key      ('newrtag', [val, di.key('newrtag')[1]])

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

# ---------------------------------- Materials -- mats

@try_catch
def cb_mat_setname  (evt,val): di.props_set_text ('texts', evt-EVT_INC, val)
@try_catch
def cb_mat_setmodel (evt,val): di.props_set_item ('mats', evt-EVT_INC, 0, val-1)
@try_catch
def cb_mat_setE     (evt,val): di.props_set_item ('mats', evt-EVT_INC, 1, float(val))
@try_catch
def cb_mat_setnu    (evt,val): di.props_set_item ('mats', evt-EVT_INC, 2, float(val))
@try_catch
def cb_mat_setk     (evt,val): di.props_set_item ('mats', evt-EVT_INC, 3, float(val))
@try_catch
def cb_mat_setlam   (evt,val): di.props_set_item ('mats', evt-EVT_INC, 4, float(val))
@try_catch
def cb_mat_setkap   (evt,val): di.props_set_item ('mats', evt-EVT_INC, 5, float(val))
@try_catch
def cb_mat_setphics (evt,val): di.props_set_item ('mats', evt-EVT_INC, 6, float(val))
@try_catch
def cb_mat_setG     (evt,val): di.props_set_item ('mats', evt-EVT_INC, 7, float(val))
@try_catch
def cb_mat_setv     (evt,val): di.props_set_item ('mats', evt-EVT_INC, 8, float(val))
@try_catch
def cb_mat_setA     (evt,val): di.props_set_item ('mats', evt-EVT_INC, 9, float(val))
@try_catch
def cb_mat_setIzz   (evt,val): di.props_set_item ('mats', evt-EVT_INC, 10, float(val))
@try_catch
def cb_mat_setgw    (evt,val): di.props_set_item ('mats', evt-EVT_INC, 12, float(val))
@try_catch
def cb_mat_setAr    (evt,val): di.props_set_item ('mats', evt-EVT_INC, 13, float(val))
@try_catch
def cb_mat_setAt    (evt,val): di.props_set_item ('mats', evt-EVT_INC, 14, float(val))
@try_catch
def cb_mat_setks    (evt,val): di.props_set_item ('mats', evt-EVT_INC, 15, float(val))
@try_catch
def cb_mat_setc     (evt,val): di.props_set_item ('mats', evt-EVT_INC, 16, float(val))
@try_catch
def cb_mat_setphi   (evt,val): di.props_set_item ('mats', evt-EVT_INC, 17, float(val))
@try_catch
def cb_mat_del      (evt,val): di.props_del      ('mats', evt-EVT_INC)

# ---------------------------------- FEM -- nbrys

@try_catch
def cb_nbry_setx  (evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC, 0, float(val))
@try_catch
def cb_nbry_sety  (evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC, 1, float(val))
@try_catch
def cb_nbry_setz  (evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC, 2, float(val))
@try_catch
def cb_nbry_setkey(evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC, 3, val-1)
@try_catch
def cb_nbry_setval(evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC, 4, float(val))
@try_catch
def cb_nbry_del   (evt,val): di.props_del_fem ([di.key('fem_stage')], 'nbrys', evt-EVT_INC)

# ---------------------------------- FEM -- nbsID

@try_catch
def cb_nbID_setID (evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbsID', evt-EVT_INC, 0, int(val))
@try_catch
def cb_nbID_setkey(evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbsID', evt-EVT_INC, 1, val-1)
@try_catch
def cb_nbID_setval(evt,val): di.props_set_fem ([di.key('fem_stage')], 'nbsID', evt-EVT_INC, 2, float(val))
@try_catch
def cb_nbID_del   (evt,val): di.props_del_fem ([di.key('fem_stage')], 'nbsID', evt-EVT_INC)

# ---------------------------------- FEM -- ebrys

@try_catch
def cb_ebry_settag(evt,val): di.props_set_fem ([di.key('fem_stage')], 'ebrys', evt-EVT_INC, 0, int(val))
@try_catch
def cb_ebry_setkey(evt,val): di.props_set_fem ([di.key('fem_stage')], 'ebrys', evt-EVT_INC, 1, val-1)
@try_catch
def cb_ebry_setval(evt,val): di.props_set_fem ([di.key('fem_stage')], 'ebrys', evt-EVT_INC, 2, float(val))
@try_catch
def cb_ebry_del   (evt,val): di.props_del_fem ([di.key('fem_stage')], 'ebrys', evt-EVT_INC)

# ---------------------------------- FEM -- fbrys

@try_catch
def cb_fbry_settag(evt,val): di.props_set_fem ([di.key('fem_stage')], 'fbrys', evt-EVT_INC, 0, int(val))
@try_catch
def cb_fbry_setkey(evt,val): di.props_set_fem ([di.key('fem_stage')], 'fbrys', evt-EVT_INC, 1, val-1)
@try_catch
def cb_fbry_setval(evt,val): di.props_set_fem ([di.key('fem_stage')], 'fbrys', evt-EVT_INC, 2, float(val))
@try_catch
def cb_fbry_del   (evt,val): di.props_del_fem ([di.key('fem_stage')], 'fbrys', evt-EVT_INC)
@try_catch
def cb_fbry_setclr(evt,val):
    obj = di.get_obj()
    id  = str(evt-EVT_INC)
    tag = obj.properties['fbrys'][id][0]
    clr = di.rgb2hex(val)
    if obj.properties.has_key('ftags'):
        for k, v in obj.properties['ftags'].iteritems():
            if int(v[0])==tag: obj.properties['ftags'][k][1] = clr
    obj.properties['fbrys'][id][3] = clr
    Blender.Window.QRedrawAll()

# ---------------------------------- FEM -- eatts

@try_catch
def cb_eatt_settag  (evt,val):
    obj  = di.get_obj()
    sids = [int(k) for k, v in obj.properties['stages'].iteritems()]
    di.props_set_fem (sids, 'eatts', evt-EVT_INC, 0, int(val))
@try_catch
def cb_eatt_settype(evt,val):
    obj  = di.get_obj()
    sids = [int(k) for k, v in obj.properties['stages'].iteritems()]
    di.props_set_fem (sids, 'eatts', evt-EVT_INC, 1, val-1)
@try_catch
def cb_eatt_setmat  (evt,val):
    obj  = di.get_obj()
    sids = [int(k) for k, v in obj.properties['stages'].iteritems()]
    di.props_set_fem (sids, 'eatts', evt-EVT_INC, 2, val-1)
@try_catch
def cb_eatt_setprops(evt,val): di.props_set_text(                       'texts', evt-EVT_INC, val)
@try_catch
def cb_eatt_isact   (evt,val): di.props_set_fem ([di.key('fem_stage')], 'eatts', evt-EVT_INC, 4, int(val))
@try_catch
def cb_eatt_act     (evt,val): di.props_set_fem ([di.key('fem_stage')], 'eatts', evt-EVT_INC, 5, int(val))
@try_catch
def cb_eatt_deact   (evt,val): di.props_set_fem ([di.key('fem_stage')], 'eatts', evt-EVT_INC, 6, int(val))
@try_catch
def cb_eatt_del     (evt,val):
    obj  = di.get_obj()
    sids = [int(k) for k, v in obj.properties['stages'].iteritems()]
    di.props_del_fem (sids, 'eatts', evt-EVT_INC)

# ---------------------------------- FEM

@try_catch
def cb_fem_fullsc(evt,val): di.set_key('fullsc', val)
@try_catch
def cb_fem_stage (evt,val):
    obj = di.get_obj()
    for k, v in obj.properties['stages'].iteritems():
        if val==v[0]: di.set_key('fem_stage', int(k))
    Blender.Window.QRedrawAll()
@try_catch
def cb_fem_stage_act  (evt,val): di.props_set_item ('stages', evt-EVT_INC, 6, int(val))
@try_catch
def cb_fem_stage_desc (evt,val): di.props_set_text ('texts', evt-EVT_INC, val)
@try_catch
def cb_fem_apply_bf   (evt,val): di.props_set_item ('stages', evt-EVT_INC, 2, int(val))
@try_catch
def cb_fem_clear_disp (evt,val): di.props_set_item ('stages', evt-EVT_INC, 3, int(val))
@try_catch
def cb_fem_ndiv       (evt,val): di.props_set_item ('stages', evt-EVT_INC, 4, int(val))
@try_catch
def cb_fem_dtime      (evt,val): di.props_set_item ('stages', evt-EVT_INC, 5, float(val))

# ---------------------------------- Results

@try_catch
def cb_res_show       (evt,val): di.set_key ('show_res',        val)
@try_catch
def cb_res_stage      (evt,val):
    obj         = di.get_obj()
    res_nstages = len(obj.properties['res'])
    if   res_nstages==0:   di.set_key ('res_stage', 1)
    elif val<=res_nstages: di.set_key ('res_stage', val)
    else: Blender.Window.QRedrawAll()
@try_catch
def cb_res_lbl        (evt,val): di.set_key ('res_lbl',         val-1)
@try_catch
def cb_res_show_scalar(evt,val): di.set_key ('res_show_scalar', val)
@try_catch
def cb_res_warp_scale (evt,val): di.set_key ('res_warp_scale',  val)
@try_catch
def cb_res_show_warp  (evt,val): di.set_key ('res_show_warp',   val)

@try_catch
def cb_res_ext       (evt,val): di.set_key ('res_ext',         val-1)
@try_catch
def cb_res_show_ext  (evt,val): di.set_key ('res_show_extra',  val)
@try_catch
def cb_res_ext_scale (evt,val): di.set_key ('res_ext_scale',   val)
@try_catch
def cb_res_ext_txt   (evt,val): di.set_key ('res_ext_txt',     val)
@try_catch
def cb_res_nodes     (evt,val):
    obj = di.get_obj()
    obj.properties['res_nodes'] = val
    Blender.Window.QRedrawAll()


# ======================================================================================= GUI

# Draw GUI
def gui():
    # load dictionary
    d = di.load_dict()

    # Get current selected mesh
    try:                    edm, obj, msh = di.get_msh (False)
    except Exception, inst: edm, obj, msh = 0, None, None

    # Data from current selected object
    is3d        = False
    iso2        = False
    isframe     = False
    blks        = {}
    minang      = -1.0
    maxarea     = -1.0
    regs        = {}
    hols        = {}
    mats        = {}
    texts       = {}
    stages      = {}
    nstages     = 0
    nbrys       = {}
    nbsID       = {}
    ebrys       = {}
    fbrys       = {}
    eatts       = {}
    res_nstages = 0
    res_nodes   = ''
    lblmnu      = 'Labels %t'
    if obj!=None:
        if obj.properties.has_key('is3d'):  is3d       = obj.properties['is3d']
        if obj.properties.has_key('iso2'):  iso2       = obj.properties['iso2']
        if obj.properties.has_key('mesh_type'):
            if obj.properties['mesh_type']=='frame': isframe = True
        if obj.properties.has_key('blks'):    blks     = obj.properties['blks']
        if obj.properties.has_key('minang'):  minang   = obj.properties['minang']
        if obj.properties.has_key('maxarea'): maxarea  = obj.properties['maxarea']
        if obj.properties.has_key('regs'):    regs     = obj.properties['regs']
        if obj.properties.has_key('hols'):    hols     = obj.properties['hols']
        if obj.properties.has_key('mats'):    mats     = obj.properties['mats']
        if obj.properties.has_key('texts'):   texts    = obj.properties['texts']
        if obj.properties.has_key('stages'):  stages   = obj.properties['stages']
        nstages = len(stages)
        if nstages>0:
            stg = 'stg_'+str(d['fem_stage'])
            if obj.properties[stg].has_key('nbrys'): nbrys = obj.properties[stg]['nbrys']
            if obj.properties[stg].has_key('nbsID'): nbsID = obj.properties[stg]['nbsID']
            if obj.properties[stg].has_key('ebrys'): ebrys = obj.properties[stg]['ebrys']
            if obj.properties[stg].has_key('fbrys'): fbrys = obj.properties[stg]['fbrys']
            if obj.properties[stg].has_key('eatts'): eatts = obj.properties[stg]['eatts']
        if obj.properties.has_key('res'):
            res_nstages = len(obj.properties['res'])
            stage_num   = str(d['res_stage'])
            if obj.properties['res'].has_key(stage_num):
                if obj.properties['res'][stage_num].has_key('lblmnu'):
                    lblmnu = obj.properties['res'][stage_num]['lblmnu']
        if obj.properties.has_key('res_nodes'): res_nodes = obj.properties['res_nodes']

    # materials menu
    matmnu   = 'Materials %t'
    matnames = {}
    for k, v in mats.iteritems():
        tid              = int(v[11])      # text_id
        desc             = texts[str(tid)] # description
        matmnu          += '|'+desc+' %x'+str(int(k)+1)
        matnames[int(k)] = desc

    # restore EditMode
    if edm: Blender.Window.EditMode(1)

    # constants
    rg  = 20    # row gap
    cg  = 14    # column gap
    rh  = 20    # row height
    srg =  6    # small row gap

    # geometry
    W,H = Blender.Window.GetAreaSize()
    r   = H-rh-rg-d['gui_inirow']
    c   = cg
    w   = W-2*c

    # number of extra rows
    mat_extra_rows = 0
    for k, v in mats.iteritems():
        mdl = int(v[0])
        if mdl==5: mat_extra_rows += 1 # Reinforcement

    # height of boxes
    h_set           = 5*rh+srg
    h_cad           = 5*rh
    h_msh_stru_blks = rh+srg+2*rh*len(blks)
    h_msh_stru      = 4*rh+srg+h_msh_stru_blks
    h_msh_unst_regs = rh+srg+rh*len(regs)
    h_msh_unst_hols = rh+srg+rh*len(hols)
    h_msh_unst      = 6*rh+3*srg+h_msh_unst_regs+h_msh_unst_hols
    h_msh           = 9*rh+srg+h_msh_stru+h_msh_unst
    h_mat_mats      = rh+srg+2*rh*(len(mats))+rh*mat_extra_rows+srg*len(mats)
    h_mat           = 3*rh+h_mat_mats
    h_fem_nbrys     = rh+srg+rh*len(nbrys)
    h_fem_nbsID     = rh+srg+rh*len(nbsID)
    h_fem_ebrys     = rh+srg+rh*len(ebrys)
    h_fem_fbrys     = rh+srg+rh*len(fbrys)
    h_fem_eatts     = rh+srg+rh*len(eatts)*2+srg*len(eatts)
    h_fem_stage     = 9*rh+5*srg+h_fem_nbrys+h_fem_nbsID+h_fem_ebrys+h_fem_fbrys+h_fem_eatts
    h_fem           = 6*rh+srg+h_fem_stage if nstages>0 else 6*rh+srg
    h_res_stage     = 4*rh
    h_res           = 4*rh+srg+h_res_stage if res_nstages>0 else 4*rh+srg

    # clear background
    gu.background()

    # ======================================================== Settings

    gu.caption1(c,r,w,rh,'SETTINGS',EVT_REFRESH,EVT_SET_HIDEALL,EVT_SET_SHOWONLY,EVT_SET_SHOWHIDE)
    if d['gui_show_set']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_set)
        Draw.Toggle ('ON/OFF',     EVT_NONE, c    , r-rh, 60, 2*rh, d['show_props'], 'Show mesh properties'   , cb_show_props)
        Draw.Toggle ('V IDs',      EVT_NONE, c+ 60, r,    60,   rh, d['show_v_ids'], 'Show Vertex IDs'        , cb_show_v_ids)
        Draw.Toggle ('N IDs',      EVT_NONE, c+120, r,    60,   rh, d['show_n_ids'], 'Show mesh Node IDs'     , cb_show_n_ids)
        Draw.Toggle ('Blocks',     EVT_NONE, c+180, r,    60,   rh, d['show_blks'],  'Show blocks tags'       , cb_show_blks )
        Draw.Toggle ('Local Axes', EVT_NONE, c+240, r,    80,   rh, d['show_axes'],  'Show local system axes' , cb_show_axes )
        Draw.Toggle ('Regions',    EVT_NONE, c+320, r,    60,   rh, d['show_regs'],  'Show regions and holes' , cb_show_regs )
        r -= rh
        Draw.Toggle ('E Tags',     EVT_NONE, c+ 60, r,    60,   rh, d['show_etags'], 'Show edge tags'         , cb_show_etags)
        Draw.Toggle ('F Tags',     EVT_NONE, c+120, r,    60,   rh, d['show_ftags'], 'Show face tags'         , cb_show_ftags)
        Draw.Toggle ('Elems',      EVT_NONE, c+180, r,    60,   rh, d['show_elems'], 'Show elements tags'     , cb_show_elems)
        Draw.Slider ('',           EVT_NONE, c+240, r,    80,   rh, d['show_opac'],0.0,1.0,0,'Set opacitity to paint faces with tags', cb_show_opac)
        Draw.Toggle ('Reinf',      EVT_NONE, c+320, r,    60,   rh, d['show_rtags'], 'Show reinforcement tags', cb_show_rtags)
        r -= rh
        r -= srg
        Draw.PushButton ('Delete all properties', EVT_SET_DELPROPS, c, r, 200, rh, 'Delete all properties')
        r, c, w = gu.box1_out(W,cg,rh, c,r)
    r -= rh
    r -= rg

    # ======================================================== CAD

    gu.caption1(c,r,w,rh,'CAD',EVT_REFRESH,EVT_SET_HIDEALL,EVT_CAD_SHOWONLY,EVT_CAD_SHOWHIDE)
    if d['gui_show_cad']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_cad)
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
        r, c, w = gu.box1_out(W,cg,rh, c,r)
    r -= rh
    r -= rg

    # ======================================================== Mesh

    gu.caption1(c,r,w,rh,'MESH',EVT_REFRESH,EVT_SET_HIDEALL,EVT_MESH_SHOWONLY,EVT_MESH_SHOWHIDE)
    if d['gui_show_mesh']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_msh)
        Draw.Toggle      ('3D mesh',    EVT_NONE,          c,     r, 80, rh, is3d,                        'Set 3D mesh',                       cb_is3d)
        Draw.Number      ('',           EVT_NONE,          c+ 80, r, 60, rh, d['newetag'][0],    -100, 0, 'New edge tag',                      cb_etag)
        Draw.PushButton  ('Edge',       EVT_MESH_SETETAG,  c+140, r, 60, rh,                              'Set edges tag (0 => remove tag)')
        Draw.Number      ('',           EVT_NONE,          c+200, r, 60, rh, d['newftag'][0],   -1000, 0, 'New face tag',                      cb_ftag)
        Draw.ColorPicker (              EVT_NONE,          c+260, r, 60, rh, di.hex2rgb(d['newftag'][1]), 'Select color to paint tagged face', cb_fclr)
        Draw.PushButton  ('Face',       EVT_MESH_SETFTAG,  c+320, r, 60, rh,                              'Set faces tag (0 => remove tag)')
        r -= rh                         
        Draw.Toggle      ('Frame Mesh',         EVT_NONE,          c,     r,  80, rh, isframe,                     'Set frame (truss/beams only) mesh', cb_frame)
        Draw.Toggle      ('Quadratic Elements', EVT_NONE,          c+ 80, r, 120, rh, iso2,                        'Generate quadratic (o2) elements' , cb_iso2)
        Draw.Number      ('',                   EVT_NONE,          c+200, r,  60, rh, d['newrtag'][0],    -100, 0, 'New reinforcement tag',             cb_rtag)
        Draw.PushButton  ('Reinforcement',      EVT_MESH_SETRTAG,  c+260, r, 120, rh,                              'Set reinforcements tag (0 => remove tag)')
        r -= rh
        r -= srg
        Draw.PushButton  ('Del all ETags', EVT_MESH_DELALLETAGS, c,     r, 90, rh, 'Delete all edge tags')
        Draw.PushButton  ('Del all FTags', EVT_MESH_DELALLFTAGS, c+ 90, r, 90, rh, 'Delete all face tags')
        Draw.PushButton  ('Del all Reinf', EVT_MESH_DELALLRTAGS, c+180, r, 90, rh, 'Delete all reinforcement tags')
        Draw.PushButton  ('Delete mesh',   EVT_MESH_DELMESH,     c+270, r, 80, rh, 'Delete current mesh')
        r -= rh
        r -= rh

        # ----------------------- Mesh -- structured

        gu.caption2(c,r,w,rh,'Structured mesh')
        r, c, w = gu.box2_in(W,cg,rh, c,r,w,h_msh_stru)

        # ----------------------- Mesh -- structured -- blocks

        gu.caption3(c,r,w,rh, 'Blocks', EVT_MESH_ADDBLK,EVT_MESH_DELALLBLKS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_stru_blks)
        gu.text(c,r,'   ID:Tag    Local axes      nX           nY          nZ')
        for k, v in blks.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     (str(i)+':', EVT_INC+i, c     , r-rh, 60, 2*rh, int(v[0]), -100, -1, 'Block tag',                       cb_blk_tag)
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
        Draw.PushButton ('Generate (quadrilaterals/hexahedrons)', EVT_MESH_GENSTRU,  c,     r, 240, rh, 'Generated structured mesh')
        Draw.PushButton ('Write Script',                          EVT_MESH_GENSTRUS, c+240, r, 100, rh, 'Create script for structured mesh generation')
        r, c, w = gu.box2_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- unstructured

        r -= rh
        gu.caption2(c,r,w,rh,'Unstructured mesh')
        r, c, w = gu.box2_in(W,cg,rh, c,r,w,h_msh_unst)
        gu.text     (c,r,'Quality:')
        Draw.String ('q=', EVT_NONE, c+ 50, r, 80, rh, str(minang),  32, 'Set the minimum angle between edges/faces (-1 => use default)',                        cb_minang)
        Draw.String ('a=', EVT_NONE, c+130, r, 80, rh, str(maxarea), 32, 'Set the maximum area/volume (uniform) for triangles/tetrahedrons (-1 => use default)', cb_maxarea)
        r -= rh
        r -= srg

        # ----------------------- Mesh -- unstructured -- regions

        gu.caption3(c,r,w,rh, 'Regions', EVT_MESH_ADDREG,EVT_MESH_DELALLREGS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_unst_regs)
        gu.text(c,r,'  ID:Tag    max area        X             Y            Z')
        for k, v in regs.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     (str(i)+':', EVT_INC+i, c,     r, 60, rh, int(v[0]), -100,-1,'Region tag',                  cb_regs_tag)
            Draw.String     ('',         EVT_INC+i, c+ 60, r, 60, rh, '%g'%v[1],   32,   'Max area (-1 => use default)',cb_regs_maxarea)
            Draw.String     ('',         EVT_INC+i, c+120, r, 60, rh, '%g'%v[2],   32,   'X of the region',             cb_regs_setx)
            Draw.String     ('',         EVT_INC+i, c+180, r, 60, rh, '%g'%v[3],   32,   'Y of the region',             cb_regs_sety)
            Draw.String     ('',         EVT_INC+i, c+240, r, 60, rh, '%g'%v[4],   32,   'Z of the region',             cb_regs_setz)
            Draw.PushButton ('Del',      EVT_INC+i, c+300, r, 40, rh,                    'Delete this row',             cb_regs_del)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- unstructured -- holes

        r -= srg
        gu.caption3(c,r,w,rh, 'Holes', EVT_MESH_ADDHOL,EVT_MESH_DELALLHOLS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_unst_hols)
        gu.text(c,r,'   ID         X             Y             Z')
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
        Draw.PushButton ('Generate (triangles/tetrahedrons)', EVT_MESH_GENUNSTRU , c,     r, 240, rh, 'Generated unstructured mesh')
        Draw.PushButton ('Write Script',                      EVT_MESH_GENUNSTRUS, c+240, r, 100, rh, 'Create script for unstructured mesh generation')
        r, c, w = gu.box2_out(W,cg,rh, c,r+rh)
        r, c, w = gu.box1_out(W,cg,rh, c,r)
    r -= rh
    r -= rg

    # ======================================================== Materials

    gu.caption1(c,r,w,rh,'Materials',EVT_REFRESH,EVT_SET_HIDEALL,EVT_MAT_SHOWONLY,EVT_MAT_SHOWHIDE)
    if d['gui_show_mat']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_mat)

        # ----------------------- Mat -- parameters

        gu.caption2_(c,r,w,rh,'Materials', EVT_MAT_ADDMAT,EVT_MAT_DELALLMAT)
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_mat_mats)
        gu.text(c,r,'          Model                   Name/Description')
        for k, v in mats.iteritems():
            r  -= rh
            i   = int(k)
            tid = int(v[11])      # text_id
            des = texts[str(tid)] # description
            Draw.Menu       (d['mdlmnu'], EVT_INC+i,   c    , r, 140, rh, int(v[0])+1, 'Constitutive model: ex.: LinElastic', cb_mat_setmodel)
            Draw.String     ('',          EVT_INC+tid, c+140, r, 220, rh, des, 32,     'Material name/description',           cb_mat_setname)
            Draw.PushButton ('Del',       EVT_INC+i,   c+360, r,  40, rh,              'Delete this row',                     cb_mat_del)
            r -= rh
            if int(v[0])==0: # LinElastic
                Draw.String ('E=',     EVT_INC+i, c    , r, 80, rh, '%g'%v[1],   32, 'E: Young modulus',                                      cb_mat_setE)
                Draw.String ('nu=',    EVT_INC+i, c+ 80, r, 80, rh, '%g'%v[2],   32, 'nu: Poisson ratio',                                     cb_mat_setnu)
            elif int(v[0])==1: # LinDiffusion
                Draw.String ('k=',     EVT_INC+i, c    , r, 80, rh, '%g'%v[3],   32, 'k: Diffusion coefficient',                              cb_mat_setk)
            elif int(v[0])==2: # CamClay
                Draw.String ('lam=',   EVT_INC+i, c    , r, 80, rh, '%.4f'%v[4], 32, 'lam: Lambda',                                           cb_mat_setlam)
                Draw.String ('kap=',   EVT_INC+i, c+ 80, r, 80, rh, '%.4f'%v[5], 32, 'kap: Kappa',                                            cb_mat_setkap)
                Draw.String ('phics=', EVT_INC+i, c+160, r, 80, rh, '%.2f'%v[6], 32, 'phics: Shear angle at critical state',                  cb_mat_setphics)
                Draw.String ('G=',     EVT_INC+i, c+240, r, 80, rh, '%g'  %v[7], 32, 'G: Shear modulus',                                      cb_mat_setG)
                Draw.String ('v=',     EVT_INC+i, c+320, r, 80, rh, '%.4f'%v[8], 32, 'v: Specific volume',                                    cb_mat_setv)
            elif int(v[0])==3: # BeamElastic
                Draw.String ('E=',     EVT_INC+i, c    , r, 80, rh, '%g'%v[ 1],  32, 'E: Young modulus',                                       cb_mat_setE)
                Draw.String ('A=',     EVT_INC+i, c+ 80, r, 80, rh, '%g'%v[ 9],  32, 'A: Cross-sectional area',                                cb_mat_setA)
                Draw.String ('Izz=',   EVT_INC+i, c+160, r, 80, rh, '%g'%v[10],  32, 'Izz: Cross-sectional inertia',                           cb_mat_setIzz)
            elif int(v[0])==4: # BiotElastic                                                                                                         
                Draw.String ('E=',     EVT_INC+i, c    , r, 80, rh, '%g'%v[ 1],  32, 'E: Young modulus',                                       cb_mat_setE)
                Draw.String ('nu=',    EVT_INC+i, c+ 80, r, 80, rh, '%g'%v[ 2],  32, 'nu: Poisson ratio',                                      cb_mat_setnu)
                Draw.String ('k=',     EVT_INC+i, c+160, r, 80, rh, '%g'%v[ 3],  32, 'k: Isotropic permeability',                              cb_mat_setk)
                Draw.String ('gw=',    EVT_INC+i, c+240, r, 80, rh, '%g'%v[12],  32, 'gw: Water specific weight',                              cb_mat_setgw)
            elif int(v[0])==5: # Reinforcement
                Draw.String ('E=',     EVT_INC+i, c    , r, 80, rh, '%g'%v[ 1],  32, 'E: Young modulus',                                       cb_mat_setE)
                Draw.String ('Ar=',    EVT_INC+i, c+ 80, r, 80, rh, '%g'%v[13],  32, 'Ar: Area of reinforcement (steel cross sectional area)', cb_mat_setAr)
                Draw.String ('At=',    EVT_INC+i, c+160, r, 80, rh, '%g'%v[14],  32, 'At: Total area of reinforcement (steel + covering)',     cb_mat_setAt)
                Draw.String ('ks=',    EVT_INC+i, c+240, r, 80, rh, '%g'%v[15],  32, 'ks: Interfacial spring stiffness',                       cb_mat_setks)
                r -= rh
                Draw.String ('c=',     EVT_INC+i, c    , r, 80, rh, '%g'%v[16],  32, 'c: Cohesion',                                            cb_mat_setc)
                Draw.String ('phi=',   EVT_INC+i, c+ 80, r, 80, rh, '%g'%v[17],  32, 'phi: friction angle',                                    cb_mat_setphi)
            r -= srg
        r -= srg
        r -= rh
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- Mat -- END
        r, c, w = gu.box1_out(W,cg,rh, c,r+rh)
    else: r -= rh
    r -= rg

    # ======================================================== FEM

    gu.caption1(c,r,w,rh,'FEM',EVT_REFRESH,EVT_SET_HIDEALL,EVT_FEM_SHOWONLY,EVT_FEM_SHOWHIDE)
    if d['gui_show_fem']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_fem)

        Draw.PushButton ('Save stage/mats info', EVT_FEM_SAVESTAGES, c,     r, 140, rh, 'Save stage and materials information to a new object')
        Draw.PushButton ('Read stage/mats info', EVT_FEM_READSTAGES, c+140, r, 140, rh, 'Read stage and materials information from another object')
        r -= rh
        r -= rh

        # ----------------------- FEM -- stages

        gu.caption2__(c,r,w,rh,'Stage #                  / %d'%(nstages),EVT_FEM_ADDSTAGE,EVT_FEM_DELSTAGE,EVT_FEM_DELALLSTAGES)
        if (len(stages)>0):
            i   = d['fem_stage']
            sid = str(i)               # stage_id
            num = int(stages[sid][0])  # num
            tid = int(stages[sid][1])  # text_id
            abf = int(stages[sid][2])  # apply body forces
            cdi = int(stages[sid][3])  # clear displacements
            ndi = int(stages[sid][4])  # ndiv
            dti = '%g'%stages[sid][5]  # dtime
            act = int(stages[sid][6])  # active?
            des = texts[str(tid)]      # description
            Draw.Number ('', EVT_NONE, c+55, r+2, 60, rh-4, num, 0,99,'Show stage', cb_fem_stage)
            r, c, w = gu.box2_in(W,cg,rh, c,r,w,h_fem_stage)
            Draw.Toggle ('Active', EVT_INC+i,   c,    r,  80, rh, act,      'Set stage Active/Inactive during simulations?', cb_fem_stage_act)
            Draw.String ('Desc: ', EVT_INC+tid, c+80, r, 280, rh, des, 256, 'Description of this stage',                     cb_fem_stage_desc)
            r -= rh
            Draw.Toggle ('Apply body forces',   EVT_INC+i, c,     r, 120, rh, abf,          'Apply body forces ?',                cb_fem_apply_bf)
            Draw.Toggle ('Clear displacements', EVT_INC+i, c+120, r, 120, rh, cdi,          'Clear displacements (and strains)?', cb_fem_clear_disp)
            Draw.Number ('',                    EVT_INC+i, c+240, r,  60, rh, ndi, 1,10000, 'Number of divisions',                cb_fem_ndiv)
            Draw.String ('',                    EVT_INC+i, c+300, r,  60, rh, dti, 23,      'Delta time',                         cb_fem_dtime)
            r -= rh
            r -= srg

            # ----------------------- FEM -- nbrys

            gu.caption3(c,r,w,rh,'Nodes boundary conditions (X-Y-Y)', EVT_FEM_ADDNBRY,EVT_FEM_DELALLNBRY)
            r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_fem_nbrys)
            gu.text(c,r,'     X             Y             Z        Key      Value')
            for k, v in nbrys.iteritems():
                r -= rh
                i  = int(k)
                Draw.String     ('',          EVT_INC+i, c    , r, 60, rh, str(v[0]), 32, 'X of the node with boundary condition',                                      cb_nbry_setx)
                Draw.String     ('',          EVT_INC+i, c+ 60, r, 60, rh, str(v[1]), 32, 'Y of the node with boundary condition',                                      cb_nbry_sety)
                Draw.String     ('',          EVT_INC+i, c+120, r, 60, rh, str(v[2]), 32, 'Z of the node with boundary condition',                                      cb_nbry_setz)
                Draw.Menu       (d['dfvmnu'], EVT_INC+i, c+180, r, 60, rh, int(v[3])+1,   'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_nbry_setkey)
                Draw.String     ('',          EVT_INC+i, c+240, r, 60, rh, str(v[4]), 32, 'Value of essential/natural boundary condition',                              cb_nbry_setval)
                Draw.PushButton ('Del',       EVT_INC+i, c+300, r, 40, rh,                'Delete this row',                                                            cb_nbry_del)
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- nbsID

            r -= srg
            gu.caption3(c,r,w,rh,'Nodes boundary conditions (given IDs)', EVT_FEM_ADDNBID,EVT_FEM_DELALLNBID)
            r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_fem_nbsID)
            gu.text(c,r,'     ID        Key     Value')
            for k, v in nbsID.iteritems():
                r -= rh
                i  = int(k)
                Draw.Number     ('',          EVT_INC+i, c    , r, 60, rh, int(v[0]),0,10000, 'Set the node ID',                                                            cb_nbID_setID)
                Draw.Menu       (d['dfvmnu'], EVT_INC+i, c+ 60, r, 60, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_nbID_setkey)
                Draw.String     ('',          EVT_INC+i, c+120, r, 60, rh, str(v[2]), 32,     'Value of essential/natural boundary condition',                              cb_nbID_setval)
                Draw.PushButton ('Del',       EVT_INC+i, c+180, r, 40, rh,                    'Delete this row',                                                            cb_nbID_del)
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- ebrys

            r -= srg
            gu.caption3(c,r,w,rh,'Edges boundary conditions', EVT_FEM_ADDEBRY,EVT_FEM_DELALLEBRY)
            r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_fem_ebrys)
            gu.text(c,r,'    Tag       Key     Value')
            for k, v in ebrys.iteritems():
                r -= rh
                i  = int(k)
                Draw.Number     ('',          EVT_INC+i, c,     r, 60, rh, int(v[0]),-1000,-1,'Set tag',                                                                    cb_ebry_settag)
                Draw.Menu       (d['dfvmnu'], EVT_INC+i, c+ 60, r, 60, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_ebry_setkey)
                Draw.String     ('',          EVT_INC+i, c+120, r, 60, rh, str(v[2]), 128,    'Value of essential/natural boundary condition',                              cb_ebry_setval)
                Draw.PushButton ('Del',       EVT_INC+i, c+180, r, 40, rh,                    'Delete this row',                                                            cb_ebry_del)
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- fbrys

            r -= srg
            gu.caption3(c,r,w,rh,'Faces boundary conditions', EVT_FEM_ADDFBRY,EVT_FEM_DELALLFBRY)
            r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_fem_fbrys)
            gu.text(c,r,'    Tag        Colour     Key     Value')
            for k, v in fbrys.iteritems():
                r  -= rh
                i   = int(k)
                clr = di.hex2rgb(v[3])
                Draw.Number      ('',          EVT_INC+i, c,     r, 60, rh, int(v[0]),-1000,-1,'Set tag',                                                                    cb_fbry_settag)
                Draw.ColorPicker (             EVT_INC+i, c+ 60, r, 60, rh, clr,               'Select color to paint tagged face',                                          cb_fbry_setclr)
                Draw.Menu        (d['dfvmnu'], EVT_INC+i, c+120, r, 60, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_fbry_setkey)
                Draw.String      ('',          EVT_INC+i, c+180, r, 60, rh, str(v[2]), 128,    'Value of essential/natural boundary condition',                              cb_fbry_setval)
                Draw.PushButton  ('Del',       EVT_INC+i, c+240, r, 40, rh,                    'Delete this row',                                                            cb_fbry_del)
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- eatts

            r -= srg
            if num==1: gu.caption3 (c,r,w,rh,'Elements attributes', EVT_FEM_ADDEATT,EVT_FEM_DELALLEATT)
            else:      gu.caption3_(c,r,w,rh,'Elements attributes')
            r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_fem_eatts)
            gu.text(c,r,'    Tag               Type                        Material')
            if num>1:
                for k, v in obj.properties['stages'].iteritems():
                    if int(v[0])==1: # first stage
                        fstg = 'stg_'+k
                        break
                featts = obj.properties[fstg]['eatts'] if obj.properties[fstg].has_key('eatts') else {} # first stage eatts
            for k, v in eatts.iteritems():
                r    -= rh
                i     = int(k)
                tid   = int(v[3])       # text id
                props = texts[str(tid)] # properties
                if num==1:
                    Draw.Number     ('',           EVT_INC+i,   c,     r-rh,  60, 2*rh, int(v[0]),-1000,0, 'Set tag',                         cb_eatt_settag)
                    Draw.Menu       (d['etymnu'],  EVT_INC+i,   c+ 60, r,    120,   rh, int(v[1])+1,       'Element type: ex.: Quad4PStrain', cb_eatt_settype)
                    Draw.Menu       (matmnu,       EVT_INC+i,   c+180, r,    140,   rh, int(v[2])+1,       'Choose material',                 cb_eatt_setmat)
                    Draw.PushButton ('Del',        EVT_INC+i,   c+320, r-rh,  40, 2*rh,                    'Delete this row',                 cb_eatt_del)
                    r -= rh                         
                    Draw.String     ('',           EVT_INC+tid, c+ 60, r,    120,   rh, props,  128,       'Additional properties (gam=specific weight, cq=correct moment due to distributed load in beams...)', cb_eatt_setprops)
                    Draw.Toggle     ('Is Active',  EVT_INC+i,   c+180, r,    140,   rh, int(v[4]),         'Is active ?',                            cb_eatt_isact)
                else:
                    etag =          str(featts[k][0])
                    etyp = d['ety'][int(featts[k][1])]
                    emat = matnames[int(featts[k][2])] if matnames.has_key(int(featts[k][2])) else ''
                    prps = texts[str(str(int(featts[k][3])))] # properties
                    gu.label_ (etag,  c,     r-rh,  60, 2*rh)
                    gu.label  (etyp,  c+ 60, r,    120,   rh)
                    gu.label  (emat,  c+180, r,    140,   rh)
                    r -= rh                         
                    gu.label    (prps,                    c+ 60, r, 120, rh)
                    Draw.Toggle ('Activate',   EVT_INC+i, c+180, r,  70, rh, int(v[5]),   'Activate this element at this stage?',   cb_eatt_act)
                    Draw.Toggle ('Deactivate', EVT_INC+i, c+250, r,  70, rh, int(v[6]),   'Deactivate this element at this stage?', cb_eatt_deact)
                r -= srg
            r -= srg
            r, c, w = gu.box3_out(W,cg,rh, c,r)

            # ----------------------- FEM -- stages -- END
            r, c, w = gu.box2_out(W,cg,rh, c,r+rh)
        else: r -= rh

        # ----------------------- FEM -- END
        r -= srg
        Draw.PushButton ('Run analysis',     EVT_FEM_RUN,      c    , r, 100, rh, 'Run a FE analysis directly (without script)')
        Draw.Toggle     ('full',             EVT_NONE,         c+100, r,  40, rh, d['fullsc'], 'Generate full script (including mesh setting up)', cb_fem_fullsc)
        Draw.PushButton ('Write script',     EVT_FEM_SCRIPT,   c+140, r,  80, rh, 'Generate script for FEM')
        Draw.PushButton ('View in ParaView', EVT_FEM_PARAVIEW, c+220, r, 120, rh, 'View results in ParaView')
        r, c, w = gu.box1_out(W,cg,rh, c,r)
    r -= rh
    r -= rg

    # ======================================================== Results

    gu.caption1(c,r,w,rh,'RESULTS',EVT_REFRESH,EVT_SET_HIDEALL,EVT_RES_SHOWONLY,EVT_RES_SHOWHIDE)
    if d['gui_show_res']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_res)
        gu.caption2(c,r,w,rh,'Stage #                  / %d'%(res_nstages))
        if (res_nstages>0):
            Draw.Number ('', EVT_NONE, c+55, r+2, 60, rh-4, d['res_stage'], 1,100,'Show results of specific stage', cb_res_stage)
            r, c, w = gu.box2_in(W,cg,rh, c,r,w,h_res_stage)
            Draw.Toggle     ('ON/OFF',    EVT_NONE,      c    , r-rh, 60, 2*rh, d['show_res'],             'Show results'               , cb_res_show)
            Draw.Menu       (lblmnu,      EVT_NONE,      c+ 60, r,    60,   rh, d['res_lbl']+1,            'Key such as ux, uy, fx, fz' , cb_res_lbl)
            Draw.Toggle     ('Scalar',    EVT_NONE,      c+120, r,    60,   rh, d['res_show_scalar'] ,     'Show scalar values'         , cb_res_show_scalar)
            Draw.String     ('sf=' ,      EVT_NONE,      c+180, r,    60,   rh, d['res_warp_scale']  , 32, 'Set warp (deformed) scale'  , cb_res_warp_scale)
            Draw.Toggle     ('Warp',      EVT_NONE,      c+240, r,    60,   rh, d['res_show_warp']   ,     'Show warped (deformed) mesh', cb_res_show_warp)
            Draw.PushButton ('Stats',     EVT_RES_STATS, c+300, r,    60,   rh,                             'Show statistics')
            r -= rh
            Draw.Menu       (d['extmnu'], EVT_NONE,      c+ 60, r,    60,   rh, d['res_ext']+1,            'Key such as N, M, V'     , cb_res_ext)
            Draw.String     ('sf=' ,      EVT_NONE,      c+120, r,    60,   rh, d['res_ext_scale']  , 32,  'Set extra drawing scale' , cb_res_ext_scale)
            Draw.Toggle     ('Extra',     EVT_NONE,      c+180, r,    60,   rh, d['res_show_extra'] ,      'Show extra output'       , cb_res_show_ext)
            Draw.Toggle     ('Values',    EVT_NONE,      c+240, r,    60,   rh, d['res_ext_txt'] ,         'Show extra values'       , cb_res_ext_txt)
            r, c, w = gu.box1_out(W,cg,rh, c,r)
        r -= rh
        r -= srg
        Draw.String     ('',                 EVT_NONE,       c,     r, 200, rh, res_nodes,  256, 'List of node (separated by commas) to generate full output in report', cb_res_nodes)
        Draw.PushButton ('Generate report',  EVT_RES_REPORT, c+200, r, 120, rh,                  'Generate report')
    r -= rg


# Register GUI
Draw.Register (gui, event, button_event)


# ================================================================================ ScriptLink

# Load script for View3D drawing
dict       = di.load_dict()
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
