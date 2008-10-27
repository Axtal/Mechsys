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
import timeit
import Blender
from   Blender import Draw, BGL
from   Blender.Mathutils import Vector
import bpy
import msys_cad  as ca
import msys_mesh as me
import msys_dict as di
import msys_fem  as fe


# ================================================================================= Constants

EVT_INC              = 10000 # increment used when passing IDs through events
EVT_NONE             =  0 # used for buttons with callbacks
EVT_REFRESH          =  1 # refresh all windows
# SETTINGS
EVT_SHOWHIDE_SET     =  2 # show/hide SET box
EVT_SET_DELPROPS     =  3 # delete all properties
# CAD
EVT_SHOWHIDE_CAD     =  4 # show/hide CAD box
EVT_CAD_ADDXYZ       =  5 # type in values to define point(s)
EVT_CAD_FILLET       =  6 # fillet two edges
EVT_CAD_EDGEBREAK    =  7 # break edge
EVT_CAD_EDGEBREAKM   =  8 # break edge at middle point
EVT_CAD_EDGEINTERS   =  9 # edge closest distance
EVT_CAD_FPOINT       = 10 # read points from file
EVT_CAD_FSPLINE      = 11 # create a spline from points in a file
# Mesh
EVT_SHOWHIDE_MESH    = 12 # show/hide MESH box
EVT_MESH_SETETAG     = 13 # set edges tag
EVT_MESH_SETFTAG     = 14 # set faces tag
# Mesh -- structured
EVT_MESH_ADDBLK      = 15 # set block 2D: 4 or 8 edges, 3D: 8 or 20 edges
EVT_MESH_DELALLBLKS  = 16 # set block 2D: 4 or 8 edges, 3D: 8 or 20 edges
EVT_MESH_GENSTRU     = 17 # generate structured mesh using MechSys module
# Mesh -- unstructured
EVT_MESH_ADDREG      = 18 # add region
EVT_MESH_DELALLREGS  = 19 # delete all regions
EVT_MESH_ADDHOLE     = 20 # add hole
EVT_MESH_DELALLHOLES = 21 # delete all holes
EVT_MESH_GENUNSTRU   = 22 # generate structured mesh using MechSys module
# FEM
EVT_SHOWHIDE_FEM     = 23 # show/hide FEM box
EVT_FEM_ADDNBRY      = 24 # add nodes boundary (given coordinates)
EVT_FEM_DELALLNBRY   = 25 # delete all nodes boundary
EVT_FEM_ADDNBRYID    = 26 # add nodes boundary (given nodes IDs)
EVT_FEM_DELALLNBRYID = 27 # delete all nodes boundary
EVT_FEM_RUN          = 28 # run a FE simulation
EVT_FEM_SCRIPT       = 29 # generate script for FEM 
EVT_FEM_PARAVIEW     = 30 # view in ParaView
# Results
EVT_SHOWHIDE_RES     = 31 # show/hide results box


# ==================================================================================== Events

# Handle input events
def event(evt, val):
    if   evt==Draw.QKEY and not val: Draw.Exit()
    elif evt==Draw.ESCKEY: Draw.Redraw(1)
    elif evt==Draw.WHEELDOWNMOUSE or evt==Draw.DOWNARROWKEY:
        d = di.load_dict()
        d['gui_inirow'] -= 40
        Blender.Window.QRedrawAll()
    elif evt==Draw.WHEELUPMOUSE or evt==Draw.UPARROWKEY:
        d = di.load_dict()
        if d['gui_inirow']<0:
            d['gui_inirow'] += 40
            Blender.Window.QRedrawAll()
    elif evt==Draw.PAGEDOWNKEY:
        d = di.load_dict()
        d['gui_inirow'] -= 100
        Blender.Window.QRedrawAll()
    elif evt==Draw.PAGEUPKEY:
        d = di.load_dict()
        if d['gui_inirow']<0:
            d['gui_inirow'] += 100
            Blender.Window.QRedrawAll()


# Handle button events
def button_event(evt):
    d = di.load_dict()
    if evt==EVT_REFRESH: Blender.Window.QRedrawAll()

    if True:#try:
        # ----------------------------------------------------------------------------------- SETTINGS

        # show/hide CAD box
        if evt==EVT_SHOWHIDE_SET: di.toggle_key_and_redraw('gui_show_set')

        # delete all properties
        elif evt==EVT_SET_DELPROPS:
            scn = bpy.data.scenes.active
            obs = scn.objects.selected
            if len(obs)>0:
                message = 'Confirm delete ALL properties?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    for o in obs:
                        for k, v in o.properties.iteritems():
                            print '[1;34mMechSys[0m: >> Object = ', o.name, ' [1;35mDeleted property:[0m ', k, str(v)
                            o.properties.pop(k)
                        for p in o.getAllProperties():
                            print '[1;34mMechSys[0m: Object = ', o.name, ' Deleted property: ', p.name, p.data
                            o.removeProperty(p)
                Blender.Window.QRedrawAll()

        # ----------------------------------------------------------------------------------- CAD

        # show/hide CAD box
        elif evt==EVT_SHOWHIDE_CAD: di.toggle_key_and_redraw('gui_show_cad')

        # add vertices
        elif evt==EVT_CAD_ADDXYZ:
            x = d['newpoint_x']
            y = d['newpoint_y']
            z = d['newpoint_z']
            ca.add_point (float(x), float(y), float(z))

        elif evt==EVT_CAD_FILLET:     ca.fillet(float(d['fillet_radius']),d['fillet_steps'])
        elif evt==EVT_CAD_EDGEBREAK:  ca.break_edge()
        elif evt==EVT_CAD_EDGEBREAKM: ca.break_edge(True)
        elif evt==EVT_CAD_EDGEINTERS: ca.edge_intersect()
        elif evt==EVT_CAD_FPOINT:     Blender.Window.FileSelector(ca.add_points_from_file, "Read X Y Z cols")
        elif evt==EVT_CAD_FSPLINE:    Blender.Window.FileSelector(ca.add_spline_from_file, "Read X Y Z cols")

        # ---------------------------------------------------------------------------------- Mesh

        # show/hide MESH box
        elif evt==EVT_SHOWHIDE_MESH: di.toggle_key_and_redraw('gui_show_mesh')

        # set edges tag
        elif evt==EVT_MESH_SETETAG:
            edm, obj, msh = di.get_msh()
            for eid in di.get_selected_edges(msh):
                if not obj.properties.has_key('etags'): obj.properties['etags'] = {}
                if d['newetag'][0]==0: obj.properties['etags'].pop(str(eid))
                else:                  obj.properties['etags'].update({str(eid):d['newetag']})
                if len(obj.properties['etags'])==0: obj.properties.pop('etags')
            Blender.Window.QRedrawAll()
            if edm: Blender.Window.EditMode(1) # return to EditMode

        # set faces tag
        elif evt==EVT_MESH_SETFTAG:
            edm, obj, msh = di.get_msh()
            sel    = di.get_selected_edges(msh)
            nedges = len(sel)
            if nedges==3 or nedges==6 or nedges==4 or nedges==8:
                eids = '_'.join([str(id) for id in sel])
                if not obj.properties.has_key('ftags'): obj.properties['ftags'] = {}
                if d['newftag'][0]==0: obj.properties['ftags'].pop(eids)
                else:                  obj.properties['ftags'].update({eids:d['newftag']})
                if len(obj.properties['ftags'])==0: obj.properties.pop('ftags')
            else: raise Exception('To set a face tag (FTAG), 3, 4, 6, or 8 edges must be selected')
            Blender.Window.QRedrawAll()
            if edm: Blender.Window.EditMode(1) # return to EditMode

        # ---------------------------------------------------------------------------- Structured

        # add block
        elif evt==EVT_MESH_ADDBLK:
            edm, obj, msh = di.get_msh()
            sel    = di.get_selected_edges(msh)
            nedges = len(sel)
            if d['newblk_3d']:
                if nedges==8 or nedges==20: pass
                else: raise Exception('To set a 3D block, 8 or 20 edges must be selected')
            else:
                if nedges==4 or nedges==8: pass
                else: raise Exception('To set a 2D block, 4 or 8 edges must be selected')
            if not obj.properties.has_key('blks'): obj.properties['blks'] = {}
            eids  = '_'.join([str(id) for id in sel])
            if obj.properties['blks'].has_key(eids): raise Exception('Block with edges '+eids.replace('_',' ')+' was already added')
            nblks = len(obj.properties['blks'])
            #                                   0   1   2  3  4  5 6 7   8   9  10  11 12 13 14 15 16 17     0   1      2    3     4     5  6  7   8  9 10    11   12   13     14     15    16      17
            obj.properties['blks'][eids] = [nblks, -1, -1,-1,-1, 2,2,2, 0.0,0.0,0.0,  0,0,0, -1,-1,-1,-1] # ID, tag, x_eid,y_eid,z_eid, nx,ny,nz, ax,ay,az, linx,liny,linz, origin,x_plus,y_plus,z_plus
            Blender.Window.QRedrawAll()
            if edm: Blender.Window.EditMode(1) # return to EditMode

        # delete all blocks
        elif evt==EVT_MESH_DELALLBLKS:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    obj.properties.pop('blks')
                    Blender.Window.QRedrawAll()

        # generate structured mesh via MechSys
        elif evt==EVT_MESH_GENSTRU:
            tt = timeit.Timer('msys_mesh.gen_struct_mesh()', 'import msys_mesh')
            t  = tt.timeit(number=1)
            print '[1;34mMechSys[0m: time spent on generation and drawing = [1;31m',t,'[0m [1;32mseconds[0m'

        # -------------------------------------------------------------------------- Unstructured

        # add region
        elif evt==EVT_MESH_ADDREG:
            obj = di.get_obj()
            if obj!=None:
                x, y, z = Blender.Window.GetCursorPos()
                regs = di.get_regs (obj)
                di.set_reg (obj, len(regs), -1, '-1', str(x),str(y),str(z))
                Blender.Window.QRedrawAll()

        # delete all regions
        elif evt==EVT_MESH_DELALLREGS:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_regs(obj)
                    Blender.Window.QRedrawAll()

        # add hole
        elif evt==EVT_MESH_ADDHOLE:
            obj = di.get_obj()
            if obj!=None:
                x, y, z = Blender.Window.GetCursorPos()
                hols = di.get_hols (obj)
                di.set_hol (obj, len(hols), str(x),str(y),str(z))
                Blender.Window.QRedrawAll()

        # delete all holes
        elif evt==EVT_MESH_DELALLHOLES:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_hols(obj)
                    Blender.Window.QRedrawAll()

        # generate unstructured mesh via MechSys
        elif evt==EVT_MESH_GENUNSTRU:
            tt = timeit.Timer('msys_mesh.gen_unstruct_mesh()', 'import msys_mesh')
            t  = tt.timeit(number=1)
            print '[1;34mMechSys[0m: time spent on generation and drawing = [1;31m',t,'[0m [1;32mseconds[0m'

        # ----------------------------------------------------------------------------------- FEM 

        # show/hide FEM box
        elif evt==EVT_SHOWHIDE_FEM: di.toggle_key_and_redraw('gui_show_fem')

        # add nodes boundary
        elif evt==EVT_FEM_ADDNBRY:
            obj = di.get_obj()
            if obj!=None:
                nbrys = di.get_nbrys (obj)
                di.set_nbry (obj, len(nbrys), '0.0','0.0','0.0', 'uy', '0.0')
                Blender.Window.QRedrawAll()

        # delete all nodes boundary
        elif evt==EVT_FEM_DELALLNBRY:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_nbrys(obj)
                    Blender.Window.QRedrawAll()

        # add nodes boundary (given tags)
        elif evt==EVT_FEM_ADDNBRYID:
            obj = di.get_obj()
            if obj!=None:
                nbryids = di.get_nbryids (obj)
                di.set_nbryid (obj, len(nbryids), '0', 'ux', '0.0')
                Blender.Window.QRedrawAll()

        # delete all nodes boundary (given tags)
        elif evt==EVT_FEM_DELALLNBRYID:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_nbryids(obj)
                    Blender.Window.QRedrawAll()

        # run a FE simulation
        elif evt==EVT_FEM_RUN:
            obj = di.get_obj ()
            fe.run_analysis (obj)

        # generate FE script
        elif evt==EVT_FEM_SCRIPT:
            obj = di.get_obj ()
            fe.gen_script   (obj)
            Blender.Window.Redraw(Blender.Window.Types.TEXT)

        # view results in ParaView
        elif evt==EVT_FEM_PARAVIEW:
            obj = di.get_obj ()
            fe.paraview (obj)

        # ----------------------------------------------------------------------------------- RES 

        # show/hide RES box
        elif evt==EVT_SHOWHIDE_RES: di.toggle_key_and_redraw('gui_show_res')


    #except Exception, inst:
        #msg = inst.args[0]
        #print '[1;34mMechSys[0m: Error: '+'[1;31m'+msg+'[0m'
        #Blender.Draw.PupMenu('ERROR|'+msg)


# ================================================================================= Callbacks

# ---------------------------------- CAD

def setx_callback         (evt,val): di.set_key('newpoint_x',    val)
def sety_callback         (evt,val): di.set_key('newpoint_y',    val)
def setz_callback         (evt,val): di.set_key('newpoint_z',    val)
def fillet_radius_callback(evt,val): di.set_key('fillet_radius', val)
def fillet_steps_callback (evt,val): di.set_key('fillet_steps',  val)


# ---------------------------------- Mesh

def show_props_callback   (evt,val): di.set_key_and_redraw('show_props',    val)
def show_axes_callback    (evt,val): di.set_key_and_redraw('show_axes',     val)
def show_blks_callback    (evt,val): di.set_key_and_redraw('show_blk_ids',  val)
def show_v_ids_callback   (evt,val): di.set_key_and_redraw('show_v_ids',    val)
def show_e_ids_callback   (evt,val): di.set_key_and_redraw('show_e_ids',    val)
def show_f_ids_callback   (evt,val): di.set_key_and_redraw('show_f_ids',    val)
def ftags_opac_callback   (evt,val): di.set_key_and_redraw('ftags_opac',    val)
def show_etags_callback   (evt,val): di.set_key_and_redraw('show_etags',    val)
def show_ftags_callback   (evt,val): di.set_key_and_redraw('show_ftags',    val)
def show_tags_txt_callback(evt,val): di.set_key_and_redraw('show_tags_txt', val)
def show_elems_callback   (evt,val): di.set_key_and_redraw('show_elems',    val)


# ---------------------------------- etag, ftag, btag

def etag_callback(evt,val): di.set_key_and_redraw('newetag', [val, d['newetag'][1]])
def ftag_callback(evt,val): di.set_key_and_redraw('newftag', [val, d['newftag'][1]])
def fclr_callback(evt,val): di.set_key_and_redraw('newftag', [d['newftag'][0], di.rgb2hex(val)])


# ---------------------------------- blocks

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
    if obj!=None:
        for k, v in obj.properties['blks'].iteritems():
            if int(v[0])==id:
                obj.properties['blks'][k][item] = val
                Blender.Window.QRedrawAll()
                break

def set_blk_axis(id,item):
    edm, obj, msh = di.get_msh()
    sel = di.get_selected_edges(msh)
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
    if edm: Blender.Window.EditMode(1) # return to EditMode

def blk_tag_callback(evt,val): set_blk_prop (evt-EVT_INC,  1, val)
def blk_xax_callback(evt,val): set_blk_axis (evt-EVT_INC,  2)
def blk_yax_callback(evt,val): set_blk_axis (evt-EVT_INC,  3)
def blk_zax_callback(evt,val): set_blk_axis (evt-EVT_INC,  4)
def blk_nx_callback (evt,val): set_blk_prop (evt-EVT_INC,  5, val)
def blk_ny_callback (evt,val): set_blk_prop (evt-EVT_INC,  6, val)
def blk_nz_callback (evt,val): set_blk_prop (evt-EVT_INC,  7, val)
def blk_ax_callback (evt,val): set_blk_prop (evt-EVT_INC,  8, float(val))
def blk_ay_callback (evt,val): set_blk_prop (evt-EVT_INC,  9, float(val))
def blk_az_callback (evt,val): set_blk_prop (evt-EVT_INC, 10, float(val))
def blk_nlx_callback(evt,val): set_blk_prop (evt-EVT_INC, 11, val)
def blk_nly_callback(evt,val): set_blk_prop (evt-EVT_INC, 12, val)
def blk_nlz_callback(evt,val): set_blk_prop (evt-EVT_INC, 13, val)

def blk_del_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        id = evt-EVT_INC
        for k, v in obj.properties['blks'].iteritems():
            if int(v[0])==id:
                obj.properties['blks'].pop(k)
                break
        for k, v in obj.properties['blks'].iteritems():
            if int(v[0])>id:
                obj.properties['blks'][k][0] = int(v[0])-1
    Blender.Window.QRedrawAll()


# ---------------------------------- Unstructured

def minangle_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_minangle (obj, val)
        Blender.Window.QRedrawAll()

def maxarea_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_maxarea (obj, val)
        Blender.Window.QRedrawAll()


# ---------------------------------- Regions

def regs_tag_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_tag (obj, evt-EVT_INC, val)

def regs_maxarea_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_maxarea (obj, evt-EVT_INC, val)

def regs_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_x (obj, evt-EVT_INC, val)
        Blender.Window.QRedrawAll()

def regs_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_y (obj, evt-EVT_INC, val)
        Blender.Window.QRedrawAll()

def regs_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_z (obj, evt-EVT_INC, val)
        Blender.Window.QRedrawAll()

def regs_delonereg_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_reg (obj, evt-EVT_INC)
        Blender.Window.QRedrawAll()


# ---------------------------------- Holes

def holes_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_x (obj, evt-EVT_INC, val)
        Blender.Window.QRedrawAll()

def holes_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_y (obj, evt-EVT_INC, val)
        Blender.Window.QRedrawAll()

def holes_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_z (obj, evt-EVT_INC, val)
        Blender.Window.QRedrawAll()

def holes_delonehole_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_hol (obj, evt-EVT_INC)
        Blender.Window.QRedrawAll()


# ---------------------------------- nbrys

def nodebry_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_x (obj, evt-EVT_INC, val)

def nodebry_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_y (obj, evt-EVT_INC, val)

def nodebry_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_z (obj, evt-EVT_INC, val)

def nodebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_key (obj, evt-EVT_INC, val)

def nodebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_val (obj, evt-EVT_INC, val)

def nodebry_delonenbry_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_nbry (obj,evt-EVT_INC)
        Blender.Window.QRedrawAll()


# ---------------------------------- nbrys (given nodes IDs)

def nodebryid_settag_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbryid_nid (obj, evt-EVT_INC, str(val))

def nodebryid_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbryid_key (obj, evt-EVT_INC, val)

def nodebryid_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbryid_val (obj, evt-EVT_INC, val)

def nodebryid_deloneebry_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_nbryid (obj,evt-EVT_INC)
        Blender.Window.QRedrawAll()


# ---------------------------------- ebrys

def edgebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_INC)
        obj.properties['ebrys'][tag][0] = val-1 # DOFVar==uz,fx,...

def edgebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_INC)
        obj.properties['ebrys'][tag][1] = float(val) # Val


# ---------------------------------- fbrys

def facebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_INC)
        obj.properties['fbrys'][tag][0] = val-1 # DOFVar==uz,fx,...

def facebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_INC)
        obj.properties['fbrys'][tag][1] = float(val) # Val


# ---------------------------------- eatts

def elematt_settype_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_INC)
        res = di.sarray_set_val(obj.properties['eatts'][tag], 0, str(val-1))
        obj.properties['eatts'].pop(tag)
        obj.properties['eatts'][tag] = res

def elematt_setmodel_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_INC)
        res = di.sarray_set_val(obj.properties['eatts'][tag], 1, str(val-1))
        obj.properties['eatts'].pop(tag)
        obj.properties['eatts'][tag] = res

def elematt_setprms_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_INC)
        val = '_'.join(val.split())
        res = di.sarray_set_val(obj.properties['eatts'][tag], 2, val)
        obj.properties['eatts'].pop(tag)
        obj.properties['eatts'][tag] = res

def elematt_setinis_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_INC)
        val = '_'.join(val.split())
        res = di.sarray_set_val(obj.properties['eatts'][tag], 3, val)
        obj.properties['eatts'].pop(tag)
        obj.properties['eatts'][tag] = res


# ---------------------------------- Results

def show_results_callback(evt,val): di.set_key_and_redraw('show_results', val)
def set_scalar_callback  (evt,val): di.set_key_and_redraw('scalar_key',   val)
def show_scalar_callback (evt,val): di.set_key_and_redraw('show_scalar',  val)
def set_warp_callback    (evt,val): di.set_key_and_redraw('warp_scale',   val)
def show_warp_callback   (evt,val): di.set_key_and_redraw('show_warp',    val)


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
    try: edm, obj, msh = di.get_msh (False)
    except Exception, inst:
        edm = 0
        obj = None
        msh = None

    # Data from current selected object
    blks     = {}
    minangle = '-1.0'
    maxarea  = '-1.0'
    regs     = []
    hols     = []
    nbrys    = []
    nbryids  = []
    ebrys    = {}
    fbrys    = {}
    eatts    = {}
    # Set default values
    if obj!=None:
        blks     = obj.properties['blks'] if obj.properties.has_key('blks') else {}
        minangle = di.get_minangle (obj)
        maxarea  = di.get_maxarea  (obj)
        regs     = di.get_regs     (obj)
        hols     = di.get_hols     (obj)
        nbrys    = di.get_nbrys    (obj)
        nbryids  = di.get_nbryids  (obj)
        ebrys    = obj.properties['ebrys'] if obj.properties.has_key('ebrys') else {}
        fbrys    = obj.properties['fbrys'] if obj.properties.has_key('fbrys') else {}
        eatts    = obj.properties['eatts'] if obj.properties.has_key('eatts') else {}

    # restore EditMode
    if edm: Blender.Window.EditMode(1)

    # Width and Height
    wid, hei = Blender.Window.GetAreaSize()

    # Col and Row to add entities
    gx  = 15   # x gap
    gy  = 20   # y gap
    sgy = 10   # small gy
    ggx = 2*gx 
    ggy = 2*gy
    cw  = 100  # column width
    rh  = 20   # row height
    ch  = rh+2 # caption height

    # Heights of windows
    h_set         = 3*rh+sgy+ggy
    h_cad         = 4*rh+ggy
    h_struct_blks = (1+len(blks)*2)*rh+sgy
    h_struct      = 3*rh+sgy + rh+h_struct_blks
    h_unst_regs   = (1+len(regs))*rh+sgy
    h_unst_hols   = (1+len(hols))*rh+sgy
    h_unstruct    = ggy+3*sgy+2*rh + rh+h_unst_regs + rh+h_unst_hols
    h_mesh        = 5*rh + rh+h_struct + rh+h_unstruct
    h_fea_nbrys   = (1+len(nbrys))*rh+sgy
    h_fea_nbryids = (1+len(nbryids))*rh+sgy
    h_fea_ebrys   = (1+len(ebrys))*rh+sgy
    h_fea_fbrys   = (1+len(fbrys))*rh+sgy
    h_fea_eatts   = (1+len(eatts))*rh+sgy
    h_fea         = sgy+ggy+3*rh + rh+h_fea_nbrys + rh+h_fea_nbryids + rh+h_fea_ebrys + rh+h_fea_fbrys + rh+h_fea_eatts
    h_res         = rh+ggy

    # Background color
    BGL.glClearColor (0.531, 0.543, 0.614, 0.0)
    BGL.glClear      (Blender.BGL.GL_COLOR_BUFFER_BIT)
    

    # SETTINGS ========================================================================================================================================

    dx  = 2*gx+50;
    row = hei-gy-d['gui_inirow']
    h   = h_set
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('SETTINGS')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,      wid-ggx-80-70, row+3, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_SET, wid-ggx-80,    row+3, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_set']:
        BGL.glColor3f     (0.85, 0.85, 0.85)
        BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Show:')
        Draw.Toggle       ('ON/OFF',    EVT_NONE, dx    , row-rh, 60, 2*rh, d['show_props'],   'Show mesh properties'      , show_props_callback)
        Draw.Toggle       ('E IDs',     EVT_NONE, dx+ 60, row,    60,   rh, d['show_e_ids'],   'Show edge IDs'             , show_e_ids_callback)
        Draw.Toggle       ('V IDs',     EVT_NONE, dx+120, row,    60,   rh, d['show_v_ids'],   'Show vertex IDs'           , show_v_ids_callback)
        Draw.Toggle       ('Blocks',    EVT_NONE, dx+180, row,    60,   rh, d['show_blk_ids'], 'Show blocks tags'          , show_blks_callback )
        Draw.Toggle       ('Local Axes',EVT_NONE, dx+240, row,    80,   rh, d['show_axes'],    'Show local system axes'    , show_axes_callback ); row -= rh
        Draw.Toggle       ('E Tags',    EVT_NONE, dx+ 60, row,    60,   rh, d['show_etags'],   'Show edge tags'            , show_etags_callback)
        Draw.Toggle       ('F Tags',    EVT_NONE, dx+120, row,    60,   rh, d['show_ftags'],   'Show face tags'            , show_ftags_callback)
        Draw.Toggle       ('Elems',     EVT_NONE, dx+180, row,    60,   rh, d['show_elems'],   'Show elements tags'        , show_elems_callback)
        Draw.Slider       ('',          EVT_NONE, dx+240, row,    80,   rh, d['ftags_opac'],0.0,1.0,0,'Set opacitity to paint faces with tags', ftags_opac_callback); row -= rh+sgy
        Draw.PushButton   ('Delete all properties', EVT_SET_DELPROPS, ggx, row, 200, rh , 'Delete all properties')
        row -= ggy

    else: row -= gy


    # CAD ==============================================================================================================================================

    h = h_cad
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('CAD')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,      wid-ggx-80-70, row+2, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_CAD, wid-ggx-80,    row+2, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_cad']:
        BGL.glColor3f     (0.85, 0.85, 0.85)
        BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Points:')
        Draw.String       ('X=',              EVT_NONE      , dx    , row,  80, rh, d['newpoint_x'],128,'Set the X value of the new point to be added',setx_callback)
        Draw.String       ('Y=',              EVT_NONE      , dx+ 80, row,  80, rh, d['newpoint_y'],128,'Set the Y value of the new point to be added',sety_callback)
        Draw.String       ('Z=',              EVT_NONE      , dx+160, row,  80, rh, d['newpoint_z'],128,'Set the Z value of the new point to be added',setz_callback)
        Draw.PushButton   ('Add x y z point', EVT_CAD_ADDXYZ, dx+240, row, 120, rh, 'Add point from coordinates'); row -= rh
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Edges:');
        Draw.String       ('Radius=',       EVT_NONE      , dx     , row , 120 , rh , d['fillet_radius'],  128,'Set radius for fillet operation',fillet_radius_callback)
        Draw.Number       ('Steps=',        EVT_NONE      , dx+120 , row , 120 , rh , d['fillet_steps' ],1,100,'Set steps for fillet operation' ,fillet_steps_callback)
        Draw.PushButton   ('Fillet edges',  EVT_CAD_FILLET, dx+240 , row , 120 , rh , 'Create a fillet between two edges'); row -= rh
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Edges:');
        Draw.PushButton   ('Break edge',        EVT_CAD_EDGEBREAK , dx,      row ,  120, rh , 'Break an edge at a previously selected point')
        Draw.PushButton   ('Break at middle',   EVT_CAD_EDGEBREAKM, dx+120 , row ,  120, rh , 'Break an edge at its middle point')
        Draw.PushButton   ('Edge intersection', EVT_CAD_EDGEINTERS, dx+240 , row ,  120, rh , 'Find the intersection (smaller distance) between two edges'); row -= rh
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('File:')
        Draw.PushButton   ('Read Points', EVT_CAD_FPOINT,  dx,     row, 120, rh, 'Add points by reading a list from file')
        Draw.PushButton   ('Read Spline', EVT_CAD_FSPLINE, dx+120, row, 120, rh, 'Add a spline by reading its points from file')
        dx   = 2*gx+55;
        row -= ggy

    else: row -= gy


    # Mesh ==============================================================================================================================================

    h = h_mesh
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('MESH')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,       wid-ggx-80-70, row+2, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_MESH, wid-ggx-80,    row+2, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_mesh']:
        netag = d['newetag'][0]
        nftag = d['newftag'][0]
        nfclr = d['newftag'][1]
        BGL.glColor3f     (0.85, 0.85, 0.85)
        BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Set tags:')
        Draw.Number       ('',     EVT_NONE,          dx,     row, 80, rh, netag,    -100, 0, 'New edge tag', etag_callback)
        Draw.PushButton   ('Edge', EVT_MESH_SETETAG,  dx+80,  row, 80, rh,                    'Set edges tag (0 => remove tag)')
        Draw.Number       ('',     EVT_NONE,          dx+160, row, 80, rh, nftag,   -1000, 0, 'New face tag', ftag_callback)
        Draw.ColorPicker  (        EVT_NONE,          dx+240, row, 80, rh, di.hex2rgb(nfclr), 'Select color to paint tagged face',   fclr_callback)
        Draw.PushButton   ('Face', EVT_MESH_SETFTAG,  dx+320, row, 80, rh,                    'Set faces tag (0 => remove tag)'); row -= rh

        # Mesh -- structured
        dx   = 2*gx+70
        ggx += gx
        dx  += gx
        h    = h_struct
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Structured')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h); row -= rh

        # Mesh -- structured -- blocks
        ggx += gx
        h    = h_struct_blks
        BGL.glColor3f     (0.53, 0.54, 0.6)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Blocks')
        BGL.glColor3f     (0.82, 0.82, 0.9)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_MESH_ADDBLK,     wid-ggx-70-80, row+2, 60, rh-4, 'Add block')
        Draw.PushButton   ('Delete all', EVT_MESH_DELALLBLKS, wid-ggx-80,    row+2, 80, rh-4, 'Delete all blocks'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('     Tag     Local axes    nX         nY        nZ'); row -= rh
        for k, v in blks.iteritems():
            i     = int(v[0])
            coefs = 'aX=%d aY=%d aZ=%d' % (v[8],v[9],v[10])
            Draw.Number     ('',    EVT_INC+i, ggx     , row, 60, rh, int(v[1]), -100, 0,  'Block tag',                       blk_tag_callback)
            Draw.PushButton ('X',   EVT_INC+i, ggx+ 60 , row, 20, rh,                      'Set X-axis',                      blk_xax_callback)
            Draw.PushButton ('Y',   EVT_INC+i, ggx+ 80 , row, 20, rh,                      'Set Y-axis',                      blk_yax_callback)
            Draw.PushButton ('Z',   EVT_INC+i, ggx+100 , row, 20, rh,                      'Set Z-axis',                      blk_zax_callback)
            Draw.Number     ('',    EVT_INC+i, ggx+120 , row, 60, rh, int(v[ 5]), 1, 1000, 'Number of divisions along X',     blk_nx_callback)
            Draw.Number     ('',    EVT_INC+i, ggx+180 , row, 60, rh, int(v[ 6]), 1, 1000, 'Number of divisions along X',     blk_ny_callback)
            Draw.Number     ('',    EVT_INC+i, ggx+240 , row, 60, rh, int(v[ 7]), 1, 1000, 'Number of divisions along X',     blk_nz_callback); row -= rh
            Draw.String     ('aX=', EVT_INC+i, ggx+120 , row, 60, rh, str(v[ 8]),     128, 'Set division coefficient aX',     blk_ax_callback)
            Draw.String     ('aY=', EVT_INC+i, ggx+180 , row, 60, rh, str(v[ 9]),     128, 'Set division coefficient aY',     blk_ay_callback)
            Draw.String     ('aZ=', EVT_INC+i, ggx+240 , row, 60, rh, str(v[10]),     128, 'Set division coefficient aZ',     blk_az_callback)
            Draw.Toggle     ('X',   EVT_INC+i, ggx+300 , row, 20, rh, int(v[11]),          'Set nonlinear divisions along X', blk_nlx_callback)
            Draw.Toggle     ('Y',   EVT_INC+i, ggx+320 , row, 20, rh, int(v[12]),          'Set nonlinear divisions along Y', blk_nly_callback)
            Draw.Toggle     ('Z',   EVT_INC+i, ggx+340 , row, 20, rh, int(v[13]),          'Set nonlinear divisions along Z', blk_nlz_callback)
            Draw.PushButton ('Del', EVT_INC+i, ggx+360 , row, 40, rh,                      'Delete this row',                 blk_del_callback); row -= rh

        # Mesh -- structured -- main
        row -= 2*sgy
        ggx -= gx
        Draw.PushButton   ('Generate (quadrilaterals/hexahedrons)', EVT_MESH_GENSTRU,  ggx, row, 300, rh, 'Generated structured mesh')

        # Mesh -- unstructured
        row -= ggy
        dx   = 2*gx+70;
        dx  += gx
        h    = h_unstruct
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Unstructured')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h); row -= rh+gy
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Quality:')
        Draw.String       ('q=',  EVT_NONE, dx,     row, 100, rh, minangle, 128, 'Set the minimum angle between edges/faces (-1 => use default)',                       minangle_callback)
        Draw.String       ('a=',  EVT_NONE, dx+100, row, 100, rh, maxarea,  128, 'Set the maximum area/volume (uniform) for triangles/tetrahedrons (-1 => use default)',maxarea_callback)

        # Mesh -- unstructured -- regs
        row -= sgy
        ggx += gx
        h    = h_unst_regs
        BGL.glColor3f     (0.53, 0.54, 0.6)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Regions')
        BGL.glColor3f     (0.82, 0.82, 0.9)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_MESH_ADDREG,     wid-ggx-70-80, row+2, 60, rh-4, 'Add region')
        Draw.PushButton   ('Delete all', EVT_MESH_DELALLREGS, wid-ggx-80,    row+2, 80, rh-4, 'Delete all regions'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('     Tag        Max area            X                  Y                  Z'); row -= rh
        for i, reg in enumerate(regs):
            Draw.Number     ('',    EVT_INC+i, ggx    , row, 60, rh, int(reg[0]), -100, 0,'Region tag',                  regs_tag_callback)
            Draw.String     ('',    EVT_INC+i, ggx+ 60, row, 80, rh,     reg[1],   128,   'Max area (-1 => use default)',regs_maxarea_callback)
            Draw.String     ('',    EVT_INC+i, ggx+140, row, 80, rh,     reg[2],   128,   'X of the rege',               regs_setx_callback)
            Draw.String     ('',    EVT_INC+i, ggx+220, row, 80, rh,     reg[3],   128,   'Y of the rege',               regs_sety_callback)
            Draw.String     ('',    EVT_INC+i, ggx+300, row, 80, rh,     reg[4],   128,   'Z of the rege',               regs_setz_callback)
            Draw.PushButton ('Del', EVT_INC+i, ggx+380, row, 40, rh,                      'Delete this row',             regs_delonereg_callback); row -= rh

        # Mesh -- unstructured -- holes
        h = h_unst_hols
        BGL.glColor3f     (0.53, 0.54, 0.6)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Holes')
        BGL.glColor3f     (0.82, 0.82, 0.9)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_MESH_ADDHOLE,     wid-ggx-70-80, row+2, 60, rh-4, 'Add hole')
        Draw.PushButton   ('Delete all', EVT_MESH_DELALLHOLES, wid-ggx-80,    row+2, 80, rh-4, 'Delete all holes'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('         X                  Y                  Z'); row -= rh
        for i, hol in enumerate(hols):
            Draw.String     ('',    EVT_INC+i, ggx    , row, 100, rh, hol[0],128,'X of the hole',holes_setx_callback)
            Draw.String     ('',    EVT_INC+i, ggx+100, row, 100, rh, hol[1],128,'Y of the hole',holes_sety_callback)
            Draw.String     ('',    EVT_INC+i, ggx+200, row, 100, rh, hol[2],128,'Z of the hole',holes_setz_callback)
            Draw.PushButton ('Del', EVT_INC+i, ggx+300, row,  40, rh, 'Delete this row', holes_delonehole_callback); row -= rh

        # Mesh -- unstructured -- main
        row -= 2*sgy
        ggx -= gx
        Draw.PushButton ('Generate (triangles/tetrahedrons)', EVT_MESH_GENUNSTRU ,  ggx, row, 300, rh, 'Generated unstructured mesh')

        # gaps/increments
        dx   = 2*gx+55;
        ggx  = 2*gx
        row -= ggy+gy

    else: row -= gy


    # FEM ==============================================================================================================================================

    h = h_fea
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('FEM')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,      wid-ggx-80-70, row+2, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_FEM, wid-ggx-80,    row+2, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_fem']:
        BGL.glColor3f     (0.8, 0.8, 0.8)
        BGL.glRecti       (gx, row, wid-gx, row-h)

        # FEM -- nodes bry (given coordinates)
        ggx += gx
        row -= gy
        h    = h_fea_nbrys
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Nodes boundaries (given X-Y-Z coordinates)')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_FEM_ADDNBRY,    wid-ggx-70-80, row+2, 60, rh-4, 'Add nodes boundary (given X-Y-Z coordinates)')
        Draw.PushButton   ('Delete all', EVT_FEM_DELALLNBRY, wid-ggx-80,    row+2, 80, rh-4, 'Delete all nodes boundary (X-Y-Z)'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('         X                  Y                  Z           Key        Value'); row -= rh
        for i, nb in enumerate(nbrys):
            Draw.String     ('',    EVT_INC+i, ggx    , row, 80, rh, nb[0],128,'X of the node with boundary condition',                                      nodebry_setx_callback)
            Draw.String     ('',    EVT_INC+i, ggx+ 80, row, 80, rh, nb[1],128,'Y of the node with boundary condition',                                      nodebry_sety_callback)
            Draw.String     ('',    EVT_INC+i, ggx+160, row, 80, rh, nb[2],128,'Z of the node with boundary condition',                                      nodebry_setz_callback)
            Draw.String     ('',    EVT_INC+i, ggx+240, row, 40, rh, nb[3],  2,'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', nodebry_setkey_callback)
            Draw.String     ('',    EVT_INC+i, ggx+280, row, 80, rh, nb[4],128,'Value of essential/natural boundary condition',                              nodebry_setval_callback)
            Draw.PushButton ('Del', EVT_INC+i, ggx+360, row, 40, rh, 'Delete this row',                                                                      nodebry_delonenbry_callback); row -= rh

        # FEM -- nodes bry (given nodes ID)
        h = h_fea_nbryids
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Nodes boundaries (given nodes ID)')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_FEM_ADDNBRYID,    wid-ggx-70-80, row+2, 60, rh-4, 'Add nodes boundary (given node ID)')
        Draw.PushButton   ('Delete all', EVT_FEM_DELALLNBRYID, wid-ggx-80,    row+2, 80, rh-4, 'Delete all nodes boundary (node ID)'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('         ID         Key        Value'); row -= rh
        for i, nbid in enumerate(nbryids):
            Draw.Number     ('',    EVT_INC+i, ggx    , row, 80, rh, int(nbid[0]),0,10000,'Set the node ID',                                                           nodebryid_settag_callback)
            Draw.String     ('',    EVT_INC+i, ggx+ 80, row, 40, rh, nbid[1],  2,         'Key such as ux, uy, fx, fz corresponding to the essential/natural variable',nodebryid_setkey_callback)
            Draw.String     ('',    EVT_INC+i, ggx+120, row, 80, rh, nbid[2],128,         'Value of essential/natural boundary condition',                             nodebryid_setval_callback)
            Draw.PushButton ('Del', EVT_INC+i, ggx+200, row, 40, rh, 'Delete this row',                                                                               nodebryid_deloneebry_callback); row -= rh

        # FEM -- edges bry
        h = h_fea_ebrys
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Edges boundaries')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('        Tag         Key        Value'); row -= rh
        for k, v in ebrys.iteritems():
            tag = int(k)
            draw_label  (ggx,row,40,rh, k)
            Draw.Menu   (d['dofvars_menu'], EVT_INC+tag, ggx+40, row, 40, rh, int(v[0])+1,    'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', edgebry_setkey_callback)
            Draw.String ('',                EVT_INC+tag, ggx+80, row, 80, rh, str(v[1]), 128, 'Value of essential/natural boundary condition',                              edgebry_setval_callback); row -= rh

        # FEM -- faces bry
        h = h_fea_fbrys
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Faces boundaries')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('        Tag         Key        Value'); row -= rh
        for k, v in fbrys.iteritems():
            tag = int(k)
            clr = di.hex2rgb(v[2])
            draw_label        (ggx,row,40,rh, k)
            BGL.glColor3f     (clr[0], clr[1], clr[2])
            BGL.glRecti       (ggx+40, row, ggx+80, row+rh)
            Draw.Menu   (d['dofvars_menu'], EVT_INC+tag, ggx+ 80, row, 40, rh, int(v[0])+1,    'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', facebry_setkey_callback)
            Draw.String ('',                EVT_INC+tag, ggx+120, row, 80, rh, str(v[1]), 128, 'Value of essential/natural boundary condition',                              facebry_setval_callback); row -= rh

        # FEM -- elems atts
        h = h_fea_eatts
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Elements attributes')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('      Tag           Type                Model           Parameters        Initial Vals'); row -= rh
        for k, v in eatts.iteritems():
            tag = int(k)
            r   = v.split()
            draw_label      (ggx,row,40,rh, k)
            Draw.Menu       (d['etypes_menu'], EVT_INC+tag, ggx+ 40, row, 120, rh, int(r[0])+1,               'Element type: ex.: Quad4PStrain',           elematt_settype_callback)
            Draw.Menu       (d['models_menu'], EVT_INC+tag, ggx+160, row, 100, rh, int(r[1])+1,               'Constitutive model: ex.: LinElastic',       elematt_setmodel_callback)
            Draw.String     ('',               EVT_INC+tag, ggx+260, row, 100, rh, r[2].replace('_',' '),128, 'Parameters: ex.: E=200 nu=0.25',            elematt_setprms_callback)
            Draw.String     ('',               EVT_INC+tag, ggx+360, row,  80, rh, r[3].replace('_',' '),128, 'Initial values: ex.: Sx=0 Sy=0 Sz=0 Sxy=0', elematt_setinis_callback); row -= rh

        # FEM -- main
        row -= 2*sgy
        ggx -= gx
        Draw.PushButton ('Run analysis',     EVT_FEM_RUN,      ggx,     row, 150, rh, 'Run a FE analysis directly (without script)')
        Draw.PushButton ('Generate script',  EVT_FEM_SCRIPT,   ggx+150, row, 150, rh, 'Generate script for FEM')
        Draw.PushButton ('View in ParaView', EVT_FEM_PARAVIEW, ggx+300, row, 150, rh, 'View results in ParaView'); row -= rh

        # increments
        row = row + rh - ggy

    else: row -= gy


    # Results ===========================================================================================================================================

    h = h_res
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('RESULTS')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,      wid-ggx-80-70, row+2, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_RES, wid-ggx-80,    row+2, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_res']:
        BGL.glColor3f     (0.85, 0.85, 0.85)
        BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Show:')
        Draw.Toggle       ('Results',EVT_NONE, dx    , row, 60, rh, d['show_results'],     'Show results'               , show_results_callback)
        Draw.String       ('Val='   ,EVT_NONE, dx+60 , row, 60, rh, d['scalar_key']  , 32, 'Set value to be shown'      , set_scalar_callback)
        Draw.Toggle       ('Scalar' ,EVT_NONE, dx+120, row, 60, rh, d['show_scalar'] ,     'Show scalar values'         , show_scalar_callback)
        Draw.String       ('M='     ,EVT_NONE, dx+180, row, 60, rh, d['warp_scale']  , 32, 'Set warp (deformed) scale'  , set_warp_callback)
        Draw.Toggle       ('Warp'   ,EVT_NONE, dx+240, row, 60, rh, d['show_warp']   ,     'Show warped (deformed) mesh', show_warp_callback)


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
