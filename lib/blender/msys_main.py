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


# ================================================================================= Constants

EVT_INC              = 10000 # increment used when passing IDs through events
EVT_NONE             =  0 # used for buttons with callbacks
EVT_REFRESH          =  1 # refresh all windows
# SETTINGS
EVT_SET_SHOWHIDE     =  2 # show/hide SET box
EVT_SET_DELPROPS     =  3 # delete all properties
# CAD
EVT_CAD_SHOWHIDE     =  4 # show/hide CAD box
EVT_CAD_ADDXYZ       =  5 # type in values to define point(s)
EVT_CAD_FILLET       =  6 # fillet two edges
EVT_CAD_BREAK        =  7 # break edge
EVT_CAD_BREAKM       =  8 # break edge at middle point
EVT_CAD_EINT         =  9 # edge closest distance
EVT_CAD_FPOINT       = 10 # read points from file
EVT_CAD_FSPLINE      = 11 # create a spline from points in a file
# Mesh
EVT_MESH_SHOWHIDE    = 12 # show/hide MESH box
EVT_MESH_SETETAG     = 13 # set edges tag
EVT_MESH_SETFTAG     = 14 # set faces tag
# Mesh -- structured
EVT_MESH_ADDBLK      = 15 # set block 2D: 4 or 8 edges, 3D: 8 or 20 edges
EVT_MESH_DELALLBLKS  = 16 # set block 2D: 4 or 8 edges, 3D: 8 or 20 edges
EVT_MESH_GENSTRU     = 17 # generate structured mesh using MechSys module
# Mesh -- unstructured
EVT_MESH_ADDREG      = 18 # add region
EVT_MESH_DELALLREGS  = 19 # delete all regions
EVT_MESH_ADDHOL      = 20 # add hole
EVT_MESH_DELALLHOLS  = 21 # delete all holes
EVT_MESH_GENUNSTRU   = 22 # generate structured mesh using MechSys module
# FEM
EVT_FEM_SHOWHIDE     = 23 # show/hide FEM box
EVT_FEM_ADDNBRY      = 24 # add nodes boundary (given coordinates)
EVT_FEM_DELALLNBRY   = 25 # delete all nodes boundary
EVT_FEM_ADDNBRYID    = 26 # add nodes boundary (given nodes IDs)
EVT_FEM_DELALLNBRYID = 27 # delete all nodes boundary
EVT_FEM_RUN          = 28 # run a FE simulation
EVT_FEM_SCRIPT       = 29 # generate script for FEM 
EVT_FEM_PARAVIEW     = 30 # view in ParaView
# Results
EVT_RES_SHOWHIDE     = 31 # show/hide results box


# ==================================================================================== Events

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

def props_del_all(key):
    obj = di.get_obj()
    msg = 'Confirm delete ALL?%t|Yes'
    res  = Blender.Draw.PupMenu(msg)
    if res>0:
        obj.properties.pop(key)
        Blender.Window.QRedrawAll()

# Handle button events
def button_event(evt):
    if evt==EVT_REFRESH: Blender.Window.QRedrawAll()

    if True:#try:
        # ----------------------------------------------------------------------------------- Settings

        if evt==EVT_SET_SHOWHIDE: di.toggle_key_and_redraw('gui_show_set')

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

        # ----------------------------------------------------------------------------------- CAD

        elif evt==EVT_CAD_SHOWHIDE: di.toggle_key_and_redraw('gui_show_cad')
        elif evt==EVT_CAD_ADDXYZ:   ca.add_point     (float(di.key('cad_x')), float(di.key('cad_y')), float(di.key('cad_z')))
        elif evt==EVT_CAD_FILLET:   ca.fillet        (float(di.key('cad_rad')), di.key('cad_stp'))
        elif evt==EVT_CAD_BREAK:    ca.break_edge    ()
        elif evt==EVT_CAD_BREAKM:   ca.break_edge    (True)
        elif evt==EVT_CAD_EINT:     ca.edge_intersect()
        elif evt==EVT_CAD_FPOINT:   Blender.Window.FileSelector(ca.add_points_from_file, 'Read X Y Z columns')
        elif evt==EVT_CAD_FSPLINE:  Blender.Window.FileSelector(ca.add_spline_from_file, 'Read X Y Z columns')

        # ---------------------------------------------------------------------------------- Mesh

        elif evt==EVT_MESH_SHOWHIDE: di.toggle_key_and_redraw('gui_show_mesh')

        elif evt==EVT_MESH_SETETAG:
            edm, obj, msh = di.get_msh()
            for eid in di.get_selected_edges(msh):
                if not obj.properties.has_key('etags'): obj.properties['etags'] = {}
                if di.key('newetag')[0]==0: obj.properties['etags'].pop     (str(eid))
                else:                       obj.properties['etags'].update ({str(eid):di.key('newetag')})
                if len(obj.properties['etags'])==0: obj.properties.pop('etags')
            Blender.Window.QRedrawAll()
            if edm: Blender.Window.EditMode(1)

        elif evt==EVT_MESH_SETFTAG:
            edm, obj, msh = di.get_msh()
            sel           = di.get_selected_edges(msh)
            nedges        = len(sel)
            if nedges==3 or nedges==6 or nedges==4 or nedges==8:
                eids = '_'.join([str(id) for id in sel])
                if not obj.properties.has_key('ftags'): obj.properties['ftags'] = {}
                if di.key('newftag')[0]==0: obj.properties['ftags'].pop     (eids)
                else:                       obj.properties['ftags'].update ({eids:di.key('newftag')})
                if len(obj.properties['ftags'])==0: obj.properties.pop('ftags')
            else: raise Exception('To set a face tag (FTAG), 3, 4, 6, or 8 edges must be selected')
            Blender.Window.QRedrawAll()
            if edm: Blender.Window.EditMode(1)

        # -------------------------------------------------------------------- Mesh -- structured

        elif evt==EVT_MESH_ADDBLK:
            edm, obj, msh = di.get_msh()
            sel           = di.get_selected_edges(msh)
            nedges        = len(sel)
            if di.key('newblk_3d'):
                if not (nedges==8 or nedges==20):
                    raise Exception('To set a 3D block, 8 or 20 edges must be selected')
            else:
                if not (nedges==4 or nedges==8):
                    raise Exception('To set a 2D block, 4 or 8 edges must be selected')
            if not obj.properties.has_key('blks'): obj.properties['blks'] = {}
            eids = '_'.join([str(id) for id in sel])
            if obj.properties['blks'].has_key(eids): raise Exception('Block with edges '+eids.replace('_',' ')+' was already added')
            nblks = len(obj.properties['blks'])
            #                                   0   1   2  3  4  5 6 7   8   9  10  11 12 13 14 15 16 17     0   1      2    3     4     5  6  7   8  9 10    11   12   13     14     15    16      17
            obj.properties['blks'][eids] = [nblks, -1, -1,-1,-1, 2,2,2, 0.0,0.0,0.0,  0,0,0, -1,-1,-1,-1] # ID, tag, x_eid,y_eid,z_eid, nx,ny,nz, ax,ay,az, linx,liny,linz, origin,x_plus,y_plus,z_plus
            Blender.Window.QRedrawAll()
            if edm: Blender.Window.EditMode(1)

        elif evt==EVT_MESH_DELALLBLKS: props_del_all('blks')
        elif evt==EVT_MESH_GENSTRU:    me.gen_struct_mesh()

        # -------------------------------------------------------------------------- Unstructured

        elif evt==EVT_MESH_ADDREG:
            obj     = di.get_obj()
            x, y, z = Blender.Window.GetCursorPos()
            if not obj.properties.has_key('regs'): obj.properties['regs'] = {}
            id = 0
            while obj.properties['regs'].has_key(str(id)): id += 1
            obj.properties['regs'][str(id)] = [-1, x,y,z, -1.0] # tag, x,y,z, max_area
            Blender.Window.QRedrawAll()

        elif evt==EVT_MESH_DELALLREGS: props_del_all('regs')

        elif evt==EVT_MESH_ADDHOL:
            obj     = di.get_obj()
            x, y, z = Blender.Window.GetCursorPos()
            if not obj.properties.has_key('hols'): obj.properties['hols'] = {}
            id = 0
            while obj.properties['hols'].has_key(str(id)): id += 1
            obj.properties['hols'][str(id)] = [x,y,z] # x,y,z
            Blender.Window.QRedrawAll()

        elif evt==EVT_MESH_DELALLHOLS: props_del_all('hols')
        elif evt==EVT_MESH_GENUNSTRU:   me.gen_unstruct_mesh()

        # ----------------------------------------------------------------------------------- FEM 

        elif evt==EVT_FEM_SHOWHIDE: di.toggle_key_and_redraw('gui_show_fem')

        elif evt==EVT_FEM_ADDNBRY:
            obj = di.get_obj()
            if not obj.properties.has_key('nbrys'): obj.properties['nbrys'] = {}
            id = 0
            while obj.properties['nbrys'].has_key(str(id)): id += 1
            obj.properties['nbrys'][str(id)] = [0.0,0.0,0.0, 0, 0.0] # x,y,z, DOFVar=uz..., Val
            Blender.Window.QRedrawAll()

        elif evt==EVT_FEM_DELALLNBRY: props_del_all('nbrys')

        elif evt==EVT_FEM_ADDNBRYID:
            obj = di.get_obj()
            if not obj.properties.has_key('nbrysID'): obj.properties['nbrysID'] = {}
            id = 0
            while obj.properties['nbrysID'].has_key(str(id)): id += 1
            obj.properties['nbrysID'][str(id)] = [0, 0, 0.0] # ID, DOFVar=uz..., Val
            Blender.Window.QRedrawAll()

        elif evt==EVT_FEM_DELALLNBRYID: props_del_all  ('nbrysID')
        elif evt==EVT_FEM_RUN:          fe.run_analysis()
        elif evt==EVT_FEM_SCRIPT:       fe.gen_script  ()
        elif evt==EVT_FEM_PARAVIEW:     fe.paraview    ()

        # ----------------------------------------------------------------------------------- RES 

        elif evt==EVT_RES_SHOWHIDE: di.toggle_key_and_redraw('gui_show_res')


    #except Exception, inst:
        #msg = inst.args[0]
        #print '[1;34mMechSys[0m: Error: '+'[1;31m'+msg+'[0m'
        #Blender.Draw.PupMenu('ERROR|'+msg)


# ================================================================================= Callbacks

# ---------------------------------- Settings

def cb_show_props(evt,val): di.set_key_and_redraw('show_props', val)
def cb_show_e_ids(evt,val): di.set_key_and_redraw('show_e_ids', val)
def cb_show_v_ids(evt,val): di.set_key_and_redraw('show_v_ids', val)
def cb_show_blks (evt,val): di.set_key_and_redraw('show_blks',  val)
def cb_show_axes (evt,val): di.set_key_and_redraw('show_axes',  val)
def cb_show_etags(evt,val): di.set_key_and_redraw('show_etags', val)
def cb_show_ftags(evt,val): di.set_key_and_redraw('show_ftags', val)
def cb_show_elems(evt,val): di.set_key_and_redraw('show_elems', val)
def cb_show_opac (evt,val): di.set_key_and_redraw('show_opac',  val)

# ---------------------------------- CAD

def cb_set_x     (evt,val): di.set_key('cad_x',   val)
def cb_set_y     (evt,val): di.set_key('cad_y',   val)
def cb_set_z     (evt,val): di.set_key('cad_z',   val)
def cb_fillet_rad(evt,val): di.set_key('cad_rad', val)
def cb_fillet_stp(evt,val): di.set_key('cad_stp', val)

# ---------------------------------- Mesh

def cb_etag(evt,val): di.set_key_and_redraw('newetag', [val, di.key('newetag')[1]])
def cb_ftag(evt,val): di.set_key_and_redraw('newftag', [val, di.key('newftag')[1]])
def cb_fclr(evt,val): di.set_key_and_redraw('newftag', [di.key('newftag')[0], di.rgb2hex(val)])

# ---------------------------------- Mesh -- structured

def set_local_system(obj,msh,blk):
    # set local system: origin and vertices on positive directions
    iex = int(obj.properties['blks'][blk][2])
    iey = int(obj.properties['blks'][blk][3])
    if iex>=0 and iey>=0:
        ex     = msh.edges[iex]
        ey     = msh.edges[iey]
        origin = -1
        x_plus = -1
        y_plus = -1
        if ex.v1.index==ey.v1.index:
            origin = ex.v1.index
            x_plus = ex.v2.index
            y_plus = ey.v2.index
        elif ex.v1.index==ey.v2.index:
            origin = ex.v1.index
            x_plus = ex.v2.index
            y_plus = ey.v1.index
        elif ex.v2.index==ey.v1.index:
            origin = ex.v2.index
            x_plus = ex.v1.index
            y_plus = ey.v2.index
        elif ex.v2.index==ey.v2.index:
            origin = ex.v2.index
            x_plus = ex.v1.index
            y_plus = ey.v1.index
        obj.properties['blks'][blk][14] = origin
        obj.properties['blks'][blk][15] = x_plus
        obj.properties['blks'][blk][16] = y_plus
        iez = int(obj.properties['blks'][blk][4])
        if iez>=0:
            z_plus = -1
            ez = msh.edges[iez]
            if ez.v1.index==origin:
                z_plus = ez.v2.index
            elif ez.v2.index==origin:
                z_plus = ez.v1.index
            obj.properties['blks'][blk][17] = z_plus

def set_blk_prop(id,item,val):
    obj = di.get_obj()
    for k, v in obj.properties['blks'].iteritems():
        if int(v[0])==id:
            obj.properties['blks'][k][item] = val
            Blender.Window.QRedrawAll()
            break

def set_blk_axis(id,item):
    edm, obj, msh = di.get_msh()
    sel           = di.get_selected_edges(msh)
    if len(sel)==1:
        eid = sel[0]
        for k, v in obj.properties['blks'].iteritems():
            if int(v[0])==id:
                eds = [int(i) for i in k.split('_')]
                if eid in eds:
                    obj.properties['blks'][k][item] = eid
                    set_local_system(obj,msh,k)
                else: raise Exception('This edge # '+str(eid)+' is not part of the block defined by edges '+k.replace('_',' '))
                Blender.Window.QRedrawAll()
                break
    else: raise Exception('Please, select (only) one edge to define a local axis')
    if edm: Blender.Window.EditMode(1)

def cb_blk_tag(evt,val): set_blk_prop (evt-EVT_INC,  1, val)
def cb_blk_xax(evt,val): set_blk_axis (evt-EVT_INC,  2)
def cb_blk_yax(evt,val): set_blk_axis (evt-EVT_INC,  3)
def cb_blk_zax(evt,val): set_blk_axis (evt-EVT_INC,  4)
def cb_blk_nx (evt,val): set_blk_prop (evt-EVT_INC,  5, val)
def cb_blk_ny (evt,val): set_blk_prop (evt-EVT_INC,  6, val)
def cb_blk_nz (evt,val): set_blk_prop (evt-EVT_INC,  7, val)
def cb_blk_ax (evt,val): set_blk_prop (evt-EVT_INC,  8, float(val))
def cb_blk_ay (evt,val): set_blk_prop (evt-EVT_INC,  9, float(val))
def cb_blk_az (evt,val): set_blk_prop (evt-EVT_INC, 10, float(val))
def cb_blk_nlx(evt,val): set_blk_prop (evt-EVT_INC, 11, val)
def cb_blk_nly(evt,val): set_blk_prop (evt-EVT_INC, 12, val)
def cb_blk_nlz(evt,val): set_blk_prop (evt-EVT_INC, 13, val)

def cb_blk_del(evt,val):
    msg = 'Confirm delete this block?%t|Yes'
    res = Blender.Draw.PupMenu(msg)
    if res>0:
        obj = di.get_obj()
        id  = evt-EVT_INC
        for k, v in obj.properties['blks'].iteritems():
            if int(v[0])==id:
                obj.properties['blks'].pop(k)
                break
        for k, v in obj.properties['blks'].iteritems():
            if int(v[0])>id:
                obj.properties['blks'][k][0] = int(v[0])-1
        Blender.Window.QRedrawAll()

# ---------------------------------- Mesh -- unstructured

def set_obj_prop_and_redraw(key,val):
    obj = di.get_obj()
    obj.properties[key] = val
    Blender.Window.QRedrawAll()

def cb_minang (evt,val): set_obj_prop_and_redraw('minang', val)
def cb_maxarea(evt,val): set_obj_prop_and_redraw('maxarea',val)

# ---------------------------------- Mesh -- unstructured -- regions

def set_obj_subprop_and_redraw(key,subkey,item,val):
    obj = di.get_obj()
    obj.properties[key][subkey][item] = val
    Blender.Window.QRedrawAll()

def del_obj_subprop_and_redraw(key,subkey):
    msg = 'Confirm delete this?%t|Yes'
    res = Blender.Draw.PupMenu(msg)
    if res>0:
        obj = di.get_obj()
        obj.properties[key].pop(subkey)
        Blender.Window.QRedrawAll()

def cb_regs_tag    (evt,val): set_obj_subprop_and_redraw('regs', str(evt-EVT_INC), 0, val)
def cb_regs_setx   (evt,val): set_obj_subprop_and_redraw('regs', str(evt-EVT_INC), 1, float(val))
def cb_regs_sety   (evt,val): set_obj_subprop_and_redraw('regs', str(evt-EVT_INC), 2, float(val))
def cb_regs_setz   (evt,val): set_obj_subprop_and_redraw('regs', str(evt-EVT_INC), 3, float(val))
def cb_regs_maxarea(evt,val): set_obj_subprop_and_redraw('regs', str(evt-EVT_INC), 4, float(val))
def cb_regs_del    (evt,val): del_obj_subprop_and_redraw('regs', str(evt-EVT_INC))

# ---------------------------------- Mesh -- unstructured -- holes

def cb_hols_setx(evt,val): set_obj_subprop_and_redraw('hols', str(evt-EVT_INC), 0, float(val))
def cb_hols_sety(evt,val): set_obj_subprop_and_redraw('hols', str(evt-EVT_INC), 1, float(val))
def cb_hols_setz(evt,val): set_obj_subprop_and_redraw('hols', str(evt-EVT_INC), 2, float(val))
def cb_hols_del (evt,val): del_obj_subprop_and_redraw('hols', str(evt-EVT_INC))

# ---------------------------------- FEM -- nbrys

def cb_nbry_setx  (evt,val): set_obj_subprop_and_redraw('nbrys', str(evt-EVT_INC), 0, float(val))
def cb_nbry_sety  (evt,val): set_obj_subprop_and_redraw('nbrys', str(evt-EVT_INC), 1, float(val))
def cb_nbry_setz  (evt,val): set_obj_subprop_and_redraw('nbrys', str(evt-EVT_INC), 2, float(val))
def cb_nbry_setkey(evt,val): set_obj_subprop_and_redraw('nbrys', str(evt-EVT_INC), 3, val-1)
def cb_nbry_setval(evt,val): set_obj_subprop_and_redraw('nbrys', str(evt-EVT_INC), 4, float(val))
def cb_nbry_del   (evt,val): del_obj_subprop_and_redraw('nbrys', str(evt-EVT_INC))

# ---------------------------------- FEM -- nbrysID

def cb_nbryID_setID (evt,val): set_obj_subprop_and_redraw('nbrysID', str(evt-EVT_INC), 0, int(val))
def cb_nbryID_setkey(evt,val): set_obj_subprop_and_redraw('nbrysID', str(evt-EVT_INC), 1, val-1)
def cb_nbryID_setval(evt,val): set_obj_subprop_and_redraw('nbrysID', str(evt-EVT_INC), 2, float(val))
def cb_nbryID_del   (evt,val): del_obj_subprop_and_redraw('nbrysID', str(evt-EVT_INC))

# ---------------------------------- FEM -- ebrys

def cb_ebry_setkey(evt,val): set_obj_subprop_and_redraw('ebrys', str(evt-EVT_INC), 0, val-1)
def cb_ebry_setval(evt,val): set_obj_subprop_and_redraw('ebrys', str(evt-EVT_INC), 1, float(val))

# ---------------------------------- FEM -- fbrys

def cb_fbry_setkey(evt,val): set_obj_subprop_and_redraw('fbrys', str(evt-EVT_INC), 0, val-1)
def cb_fbry_setval(evt,val): set_obj_subprop_and_redraw('fbrys', str(evt-EVT_INC), 1, float(val))

# ---------------------------------- FEM -- eatts

def set_eatt_and_redraw(tag,item,val):
    obj = di.get_obj()
    tag = str(evt-EVT_INC)
    res = di.sarray_set_val(obj.properties['eatts'][tag], item, val)
    obj.properties['eatts'].pop(tag)
    obj.properties['eatts'][tag] = res
    Blender.Window.QRedrawAll()

def cb_eatt_settype (evt,val): set_eatt_and_redraw(str(evt-EVT_INC), 0, str(val-1))
def cb_eatt_setmodel(evt,val): set_eatt_and_redraw(str(evt-EVT_INC), 1, str(val-1))
def cb_eatt_setprms (evt,val): set_eatt_and_redraw(str(evt-EVT_INC), 2, '_'.join(val.split()))
def cb_eatt_setinis (evt,val): set_eatt_and_redraw(str(evt-EVT_INC), 3, '_'.join(val.split()))

# ---------------------------------- Results

def cb_res_show       (evt,val): di.set_key_and_redraw('res_show',        val)
def cb_res_scalar     (evt,val): di.set_key_and_redraw('res_scalar',      val)
def cb_res_show_scalar(evt,val): di.set_key_and_redraw('res_show_scalar', val)
def cb_res_warp_scale (evt,val): di.set_key_and_redraw('res_warp_scale',  val)
def cb_res_show_warp  (evt,val): di.set_key_and_redraw('res_show_warp',   val)


# ======================================================================================= GUI

def draw_label(x,y,w,h,txt):
    BGL.glColor3f     (0.663, 0.663, 0.663)
    BGL.glRecti       (x, y, x+w, y+h)
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (x+5, y+5)
    Draw.Text         (txt)

# Draw GUI
def gui():
    # load dictionary
    d = di.load_dict()

    # Get current selected mesh
    try:                    edm, obj, msh = di.get_msh (False)
    except Exception, inst: edm, obj, msh = 0, None, None

    # Data from current selected object
    blks    = {}
    minang  = -1.0
    maxarea = -1.0
    regs    = {}
    hols    = {}
    nbrys   = {}
    nbrysID = {}
    ebrys   = {}
    fbrys   = {}
    eatts   = {}
    if obj!=None:
        blks     = obj.properties['blks']    if obj.properties.has_key('blks')    else {}
        minang   = obj.properties['minang']  if obj.properties.has_key('minang')  else -1.0
        maxarea  = obj.properties['maxarea'] if obj.properties.has_key('maxarea') else -1.0
        regs     = obj.properties['regs']    if obj.properties.has_key('regs')    else {}
        hols     = obj.properties['hols']    if obj.properties.has_key('hols')    else {}
        nbrys    = obj.properties['nbrys']   if obj.properties.has_key('nbrys')   else {}
        nbrysID  = obj.properties['nbrysID'] if obj.properties.has_key('nbrysID') else {}
        ebrys    = obj.properties['ebrys']   if obj.properties.has_key('ebrys')   else {}
        fbrys    = obj.properties['fbrys']   if obj.properties.has_key('fbrys')   else {}
        eatts    = obj.properties['eatts']   if obj.properties.has_key('eatts')   else {}

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

    # height of boxes
    h_set           = 5*rh+srg
    h_cad           = 5*rh
    h_msh_stru_blks = rh+srg+2*rh*len(blks)
    h_msh_stru      = 4*rh+srg+h_msh_stru_blks
    h_msh_unst_regs = rh+srg+rh*len(regs)
    h_msh_unst_hols = rh+srg+rh*len(hols)
    h_msh_unst      = 6*rh+3*srg+h_msh_unst_regs+h_msh_unst_hols
    h_msh           = 7*rh+h_msh_stru+h_msh_unst
    h_fem_nbrys     = rh+srg+rh*len(nbrys)
    h_fem_nbrysID   = rh+srg+rh*len(nbrysID)
    h_fem_ebrys     = rh+srg+rh*len(ebrys)
    h_fem_fbrys     = rh+srg+rh*len(fbrys)
    h_fem_eatts     = rh+srg+rh*len(eatts)
    h_fem           = 8*rh+5*srg+h_fem_nbrys+h_fem_nbrysID+h_fem_ebrys+h_fem_fbrys+h_fem_eatts
    h_res           = 3*rh

    # clear background
    gu.background()

    # ======================================================== Settings

    gu.caption1(c,r,w,rh,'SETTINGS',EVT_REFRESH,EVT_SET_SHOWHIDE)
    if d['gui_show_set']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_set)
        Draw.Toggle ('ON/OFF',     EVT_NONE, c    , r-rh, 60, 2*rh, d['show_props'], 'Show mesh properties'   , cb_show_props)
        Draw.Toggle ('E IDs',      EVT_NONE, c+ 60, r,    60,   rh, d['show_e_ids'], 'Show edge IDs'          , cb_show_e_ids)
        Draw.Toggle ('V IDs',      EVT_NONE, c+120, r,    60,   rh, d['show_v_ids'], 'Show vertex IDs'        , cb_show_v_ids)
        Draw.Toggle ('Blocks',     EVT_NONE, c+180, r,    60,   rh, d['show_blks'],  'Show blocks tags'       , cb_show_blks )
        Draw.Toggle ('Local Axes', EVT_NONE, c+240, r,    80,   rh, d['show_axes'],  'Show local system axes' , cb_show_axes )
        r -= rh
        Draw.Toggle ('E Tags',     EVT_NONE, c+ 60, r,    60,   rh, d['show_etags'], 'Show edge tags'         , cb_show_etags)
        Draw.Toggle ('F Tags',     EVT_NONE, c+120, r,    60,   rh, d['show_ftags'], 'Show face tags'         , cb_show_ftags)
        Draw.Toggle ('Elems',      EVT_NONE, c+180, r,    60,   rh, d['show_elems'], 'Show elements tags'     , cb_show_elems)
        Draw.Slider ('',           EVT_NONE, c+240, r,    80,   rh, d['show_opac'],0.0,1.0,0,'Set opacitity to paint faces with tags', cb_show_opac)
        r -= rh
        r -= srg
        Draw.PushButton ('Delete all properties', EVT_SET_DELPROPS, c, r, 200, rh, 'Delete all properties')
        r, c, w = gu.box1_out(W,cg,rh, c,r)
    r -= rh
    r -= rg

    # ======================================================== CAD

    gu.caption1(c,r,w,rh,'CAD',EVT_REFRESH,EVT_CAD_SHOWHIDE)
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

    gu.caption1(c,r,w,rh,'MESH',EVT_REFRESH,EVT_MESH_SHOWHIDE)
    if d['gui_show_mesh']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_msh)
        Draw.Number      ('',     EVT_NONE,          c,     r, 60, rh, d['newetag'][0],    -100, 0, 'New edge tag',                      cb_etag)
        Draw.PushButton  ('Edge', EVT_MESH_SETETAG,  c+ 60, r, 60, rh,                              'Set edges tag (0 => remove tag)')
        Draw.Number      ('',     EVT_NONE,          c+120, r, 80, rh, d['newftag'][0],   -1000, 0, 'New face tag',                      cb_ftag)
        Draw.ColorPicker (        EVT_NONE,          c+200, r, 60, rh, di.hex2rgb(d['newftag'][1]), 'Select color to paint tagged face', cb_fclr)
        Draw.PushButton  ('Face', EVT_MESH_SETFTAG,  c+260, r, 60, rh,                              'Set faces tag (0 => remove tag)')
        r -= rh
        r -= rh

        # ----------------------- Mesh -- structured

        gu.caption2(c,r,w,rh,'Structured mesh')
        r, c, w = gu.box2_in(W,cg,rh, c,r,w,h_msh_stru)

        # ----------------------- Mesh -- structured -- blocks

        gu.caption3(c,r,w,rh, 'Blocks', EVT_MESH_ADDBLK,EVT_MESH_DELALLBLKS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_stru_blks)
        gu.text(c,r,'     Tag     Local axes    nX         nY        nZ')
        for k, v in blks.iteritems():
            r    -= rh
            i     = int(v[0])
            coefs = 'aX=%d aY=%d aZ=%d' % (v[8],v[9],v[10])
            Draw.Number     ('',    EVT_INC+i, c     , r-rh, 60, 2*rh, int(v[1]), -100, 0,  'Block tag',                       cb_blk_tag)
            Draw.PushButton ('X',   EVT_INC+i, c+ 60 , r,    20,   rh,                      'Set X-axis',                      cb_blk_xax)
            Draw.PushButton ('Y',   EVT_INC+i, c+ 80 , r,    20,   rh,                      'Set Y-axis',                      cb_blk_yax)
            Draw.PushButton ('Z',   EVT_INC+i, c+100 , r,    20,   rh,                      'Set Z-axis',                      cb_blk_zax)
            Draw.Number     ('',    EVT_INC+i, c+120 , r,    60,   rh, int(v[ 5]), 1, 1000, 'Number of divisions along X',     cb_blk_nx)
            Draw.Number     ('',    EVT_INC+i, c+180 , r,    60,   rh, int(v[ 6]), 1, 1000, 'Number of divisions along X',     cb_blk_ny)
            Draw.Number     ('',    EVT_INC+i, c+240 , r,    60,   rh, int(v[ 7]), 1, 1000, 'Number of divisions along X',     cb_blk_nz)
            r -= rh                                          
            Draw.Toggle     ('X',   EVT_INC+i, c+ 60 , r,    20,   rh, int(v[11]),          'Set nonlinear divisions along X', cb_blk_nlx)
            Draw.Toggle     ('Y',   EVT_INC+i, c+ 80 , r,    20,   rh, int(v[12]),          'Set nonlinear divisions along Y', cb_blk_nly)
            Draw.Toggle     ('Z',   EVT_INC+i, c+100 , r,    20,   rh, int(v[13]),          'Set nonlinear divisions along Z', cb_blk_nlz)
            Draw.String     ('aX=', EVT_INC+i, c+120 , r,    60,   rh, str(v[ 8]),     128, 'Set division coefficient aX',     cb_blk_ax)
            Draw.String     ('aY=', EVT_INC+i, c+180 , r,    60,   rh, str(v[ 9]),     128, 'Set division coefficient aY',     cb_blk_ay)
            Draw.String     ('aZ=', EVT_INC+i, c+240 , r,    60,   rh, str(v[10]),     128, 'Set division coefficient aZ',     cb_blk_az)
            Draw.PushButton ('Del', EVT_INC+i, c+300 , r,    40, 2*rh,                      'Delete this row',                 cb_blk_del)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- structured -- END

        r -= srg
        Draw.PushButton ('Generate (quadrilaterals/hexahedrons)', EVT_MESH_GENSTRU, c, r, 300, rh, 'Generated structured mesh')
        r, c, w = gu.box2_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- unstructured

        r -= rh
        gu.caption2(c,r,w,rh,'Unstructured mesh')
        r, c, w = gu.box2_in(W,cg,rh, c,r,w,h_msh_unst)
        gu.text     (c,r,'Quality:')
        Draw.String ('q=', EVT_NONE, c+ 50, r, 80, rh, str(minang),  128, 'Set the minimum angle between edges/faces (-1 => use default)',                        cb_minang)
        Draw.String ('a=', EVT_NONE, c+130, r, 80, rh, str(maxarea), 128, 'Set the maximum area/volume (uniform) for triangles/tetrahedrons (-1 => use default)', cb_maxarea)
        r -= rh
        r -= srg

        # ----------------------- Mesh -- unstructured -- regions

        gu.caption3(c,r,w,rh, 'Regions', EVT_MESH_ADDREG,EVT_MESH_DELALLREGS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_unst_regs)
        gu.text(c,r,'    Tag           X             Y            Z        Max area')
        for k, v in regs.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     ('',    EVT_INC+i, c,     r, 60, rh, int(v[0]), -100, 0,'Region tag',                  cb_regs_tag)
            Draw.String     ('',    EVT_INC+i, c+ 60, r, 60, rh, str(v[1]),   32,   'X of the region',             cb_regs_setx)
            Draw.String     ('',    EVT_INC+i, c+120, r, 60, rh, str(v[2]),   32,   'Y of the region',             cb_regs_sety)
            Draw.String     ('',    EVT_INC+i, c+180, r, 60, rh, str(v[3]),   32,   'Z of the region',             cb_regs_setz)
            Draw.String     ('',    EVT_INC+i, c+240, r, 60, rh, str(v[4]),   32,   'Max area (-1 => use default)',cb_regs_maxarea)
            Draw.PushButton ('Del', EVT_INC+i, c+300, r, 40, rh,                    'Delete this row',             cb_regs_del)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- unstructured -- holes

        r -= srg
        gu.caption3(c,r,w,rh, 'Holes', EVT_MESH_ADDHOL,EVT_MESH_DELALLHOLS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_unst_hols)
        gu.text(c,r,'      X            Y             Z')
        for k, v in hols.iteritems():
            r -= rh
            i  = int(k)
            Draw.String     ('',    EVT_INC+i, c    , r, 60, rh, str(v[0]), 32, 'X of the hole',   cb_hols_setx)
            Draw.String     ('',    EVT_INC+i, c+ 60, r, 60, rh, str(v[1]), 32, 'Y of the hole',   cb_hols_sety)
            Draw.String     ('',    EVT_INC+i, c+120, r, 60, rh, str(v[2]), 32, 'Z of the hole',   cb_hols_setz)
            Draw.PushButton ('Del', EVT_INC+i, c+180, r, 40, rh,                'Delete this row', cb_hols_del)
        r -= srg
        r, c, w = gu.box3_out(W,cg,rh, c,r)

        # ----------------------- Mesh -- unstructured -- END

        r -= srg
        Draw.PushButton ('Generate (triangles/tetrahedrons)', EVT_MESH_GENUNSTRU ,  c, r, 300, rh, 'Generated unstructured mesh')
        r, c, w = gu.box2_out(W,cg,rh, c,r+rh)
        r, c, w = gu.box1_out(W,cg,rh, c,r)
    r -= rh
    r -= rg

    # ======================================================== FEM

    gu.caption1(c,r,w,rh,'FEM',EVT_REFRESH,EVT_FEM_SHOWHIDE)
    if d['gui_show_fem']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_fem)

        # ----------------------- FEM -- nbrys

        gu.caption2_(c,r,w,rh,'Nodes boundary conditions (given X-Y-Y)', EVT_FEM_ADDNBRY,EVT_FEM_DELALLNBRY)
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_nbrys)
        gu.text(c,r,'     X             Y             Z        Key      Value')
        for k, v in nbrys.iteritems():
            r -= rh
            i  = int(k)
            Draw.String     ('',                EVT_INC+i, c    , r, 60, rh, str(v[0]), 32, 'X of the node with boundary condition',                                      cb_nbry_setx)
            Draw.String     ('',                EVT_INC+i, c+ 60, r, 60, rh, str(v[1]), 32, 'Y of the node with boundary condition',                                      cb_nbry_sety)
            Draw.String     ('',                EVT_INC+i, c+120, r, 60, rh, str(v[2]), 32, 'Z of the node with boundary condition',                                      cb_nbry_setz)
            Draw.Menu       (d['dofvars_menu'], EVT_INC+i, c+180, r, 40, rh, int(v[3])+1,   'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_nbry_setkey)
            Draw.String     ('',                EVT_INC+i, c+220, r, 60, rh, str(v[4]), 32, 'Value of essential/natural boundary condition',                              cb_nbry_setval)
            Draw.PushButton ('Del',             EVT_INC+i, c+280, r, 40, rh,                'Delete this row',                                                            cb_nbry_del)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- nbrysID

        r -= srg
        gu.caption2_(c,r,w,rh,'Nodes boundary conditions (given IDs)', EVT_FEM_ADDNBRYID,EVT_FEM_DELALLNBRYID)
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_nbrysID)
        gu.text(c,r,'         ID         Key        Value')
        for k, v in nbrysID.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     ('',                EVT_INC+i, c    , r, 60, rh, int(v[0]),0,10000, 'Set the node ID',                                                            cb_nbryID_setID)
            Draw.Menu       (d['dofvars_menu'], EVT_INC+i, c+ 60, r, 40, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_nbryID_setkey)
            Draw.String     ('',                EVT_INC+i, c+100, r, 60, rh, str(v[2]), 32,     'Value of essential/natural boundary condition',                              cb_nbryID_setval)
            Draw.PushButton ('Del',             EVT_INC+i, c+160, r, 40, rh,                    'Delete this row',                                                            cb_nbryID_del)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- ebrys

        r -= srg
        gu.caption2(c,r,w,rh,'Edges boundary conditions')
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_ebrys)
        gu.text(c,r,'        Tag         Key        Value')
        for k, v in ebrys.iteritems():
            r  -= rh
            tag = int(k)
            gu.label    (c,r,40,rh, k)
            Draw.Menu   (d['dofvars_menu'], EVT_INC+tag, c+40, r, 40, rh, int(v[0])+1,    'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_ebry_setkey)
            Draw.String ('',                EVT_INC+tag, c+80, r, 80, rh, str(v[1]), 128, 'Value of essential/natural boundary condition',                              cb_ebry_setval)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- fbrys

        r -= srg
        gu.caption2(c,r,w,rh,'Faces boundary conditions')
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_fbrys)
        gu.text(c,r,'        Tag         Key        Value')
        for k, v in fbrys.iteritems():
            r  -= rh
            tag = int(k)
            clr = di.hex2rgb(v[2])
            gu.label      (c,r,40,rh, k)
            BGL.glColor3f (clr[0], clr[1], clr[2])
            BGL.glRecti   (c+40, r, c+80, r+rh)
            Draw.Menu     (d['dofvars_menu'], EVT_INC+tag, c+ 80, r, 40, rh, int(v[0])+1,    'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_fbry_setkey)
            Draw.String   ('',                EVT_INC+tag, c+120, r, 80, rh, str(v[1]), 128, 'Value of essential/natural boundary condition',                              cb_fbry_setval)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- eatts

        r -= srg
        gu.caption2(c,r,w,rh,'Elements attributes')
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_eatts)
        gu.text(c,r,'      Tag           Type                Model           Parameters        Initial Vals')
        for k, v in eatts.iteritems():
            r  -= rh
            tag = int(k)
            m   = v.split()
            gu.label    (c,r,40,rh, k)
            Draw.Menu   (d['etypes_menu'], EVT_INC+tag, c+ 40, r, 120, rh, int(m[0])+1,               'Element type: ex.: Quad4PStrain',           cb_eatt_settype)
            Draw.Menu   (d['models_menu'], EVT_INC+tag, c+160, r, 100, rh, int(m[1])+1,               'Constitutive model: ex.: LinElastic',       cb_eatt_setmodel)
            Draw.String ('',               EVT_INC+tag, c+260, r, 100, rh, m[2].replace('_',' '),128, 'Parameters: ex.: E=200 nu=0.25',            cb_eatt_setprms)
            Draw.String ('',               EVT_INC+tag, c+360, r,  80, rh, m[3].replace('_',' '),128, 'Initial values: ex.: Sx=0 Sy=0 Sz=0 Sxy=0', cb_eatt_setinis)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- END
        r -= srg
        Draw.PushButton ('Run analysis',     EVT_FEM_RUN,      c,     r, 150, rh, 'Run a FE analysis directly (without script)')
        Draw.PushButton ('Generate script',  EVT_FEM_SCRIPT,   c+150, r, 150, rh, 'Generate script for FEM')
        Draw.PushButton ('View in ParaView', EVT_FEM_PARAVIEW, c+300, r, 150, rh, 'View results in ParaView')
        r, c, w = gu.box1_out(W,cg,rh, c,r)
    r -= rh
    r -= rg

    # ======================================================== Results

    gu.caption1(c,r,w,rh,'RESULTS',EVT_REFRESH,EVT_RES_SHOWHIDE)
    if d['gui_show_res']:
        r, c, w = gu.box1_in(W,cg,rh, c,r,w,h_res)
        Draw.Toggle ('ON/OFF',          EVT_NONE, c    , r, 60, rh, d['res_show'],             'Show results'               , cb_res_show)
        Draw.Menu   (d['dofvars_menu'], EVT_NONE, c+ 60, r, 40, rh, d['res_scalar']+1,         'Key such as ux, uy, fx, fz' , cb_res_scalar)
        Draw.Toggle ('Scalar',          EVT_NONE, c+100, r, 60, rh, d['res_show_scalar'] ,     'Show scalar values'         , cb_res_show_scalar)
        Draw.String ('M=' ,             EVT_NONE, c+160, r, 60, rh, d['res_warp_scale']  , 32, 'Set warp (deformed) scale'  , cb_res_warp_scale)
        Draw.Toggle ('Warp',            EVT_NONE, c+220, r, 60, rh, d['res_show_warp']   ,     'Show warped (deformed) mesh', cb_res_show_warp)
        r, c, w = gu.box1_out(W,cg,rh, c,r)
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
