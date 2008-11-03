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

#import rpdb2; rpdb2.start_embedded_debugger('msys')

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
EVT_MESH_GENSTRUS    = 18 # script for structured mesh generation
# Mesh -- unstructured
EVT_MESH_ADDREG      = 19 # add region
EVT_MESH_DELALLREGS  = 20 # delete all regions
EVT_MESH_ADDHOL      = 21 # add hole
EVT_MESH_DELALLHOLS  = 22 # delete all holes
EVT_MESH_GENUNSTRU   = 23 # generate unstructured mesh using MechSys module
EVT_MESH_GENUNSTRUS  = 24 # script for unstructured mesh generation
# FEM
EVT_FEM_SHOWHIDE     = 25 # show/hide FEM box
EVT_FEM_ADDNBRY      = 26 # add nodes boundary (given coordinates)
EVT_FEM_ADDNBID      = 27 # add nodes boundary (given nodes IDs)
EVT_FEM_ADDEBRY      = 28 # add edges boundary
EVT_FEM_ADDFBRY      = 29 # add faces boundary
EVT_FEM_ADDEATT      = 30 # add element attributes
EVT_FEM_DELALLNBRY   = 31 # delete all nodes boundary
EVT_FEM_DELALLNBID   = 32 # delete all nodes boundary (IDs)
EVT_FEM_DELALLEBRY   = 33 # delete all edges boundary
EVT_FEM_DELALLFBRY   = 34 # delete all faces boundary
EVT_FEM_DELALLEATT   = 35 # delete all element attributes
EVT_FEM_RUN          = 36 # run a FE simulation
EVT_FEM_SCRIPT       = 37 # generate script for FEM 
EVT_FEM_PARAVIEW     = 38 # view in ParaView
# Results
EVT_RES_SHOWHIDE     = 39 # show/hide results box


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

# Handle button events
def button_event(evt):
    if evt==EVT_REFRESH: Blender.Window.QRedrawAll()

    if True:#try:
        # ----------------------------------------------------------------------------------- Settings

        if evt==EVT_SET_SHOWHIDE: di.toggle_key('gui_show_set')

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

        elif evt==EVT_CAD_SHOWHIDE: di.toggle_key('gui_show_cad')
        elif evt==EVT_CAD_ADDXYZ:   ca.add_point     (float(di.key('cad_x')), float(di.key('cad_y')), float(di.key('cad_z')))
        elif evt==EVT_CAD_FILLET:   ca.fillet        (float(di.key('cad_rad')), di.key('cad_stp'))
        elif evt==EVT_CAD_BREAK:    ca.break_edge    ()
        elif evt==EVT_CAD_BREAKM:   ca.break_edge    (True)
        elif evt==EVT_CAD_EINT:     ca.edge_intersect()
        elif evt==EVT_CAD_FPOINT:   Blender.Window.FileSelector(ca.add_points_from_file, 'Read X Y Z columns')
        elif evt==EVT_CAD_FSPLINE:  Blender.Window.FileSelector(ca.add_spline_from_file, 'Read X Y Z columns')

        # ---------------------------------------------------------------------------------- Mesh

        elif evt==EVT_MESH_SHOWHIDE: di.toggle_key_('gui_show_mesh')

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

        # -------------------------------------------------------------------- Mesh -- structured

        elif evt==EVT_MESH_ADDBLK:
            props = di.new_blk_props()
            di.props_push_new('blks', props, True, 18, 18+int(props[17]))

        elif evt==EVT_MESH_DELALLBLKS: di.props_del_all('blks')
        elif evt==EVT_MESH_GENSTRU:    me.gen_struct_mesh()
        elif evt==EVT_MESH_GENSTRUS:   me.gen_struct_mesh(True)

        # -------------------------------------------------------------------------- Unstructured

        elif evt==EVT_MESH_ADDREG:     di.props_push_new('regs', di.new_reg_props())
        elif evt==EVT_MESH_DELALLREGS: di.props_del_all ('regs')
        elif evt==EVT_MESH_ADDHOL:     di.props_push_new('hols', di.new_hol_props())
        elif evt==EVT_MESH_DELALLHOLS: di.props_del_all ('hols')
        elif evt==EVT_MESH_GENUNSTRU:  me.gen_unstruct_mesh()
        elif evt==EVT_MESH_GENUNSTRUS: me.gen_unstruct_mesh(True)

        # ----------------------------------------------------------------------------------- FEM 

        elif evt==EVT_FEM_SHOWHIDE: di.toggle_key('gui_show_fem')

        elif evt==EVT_FEM_ADDNBRY: di.props_push_new('nbrys', di.new_nbry_props())
        elif evt==EVT_FEM_ADDNBID: di.props_push_new('nbsID', di.new_nbID_props())
        elif evt==EVT_FEM_ADDEBRY: di.props_push_new('ebrys', di.new_ebry_props())
        elif evt==EVT_FEM_ADDFBRY: di.props_push_new('fbrys', di.new_fbry_props())
        elif evt==EVT_FEM_ADDEATT: di.props_push_new('eatts', di.new_eatt_props())

        elif evt==EVT_FEM_DELALLNBRY: di.props_del_all('nbrys')
        elif evt==EVT_FEM_DELALLNBID: di.props_del_all('nbsID')
        elif evt==EVT_FEM_DELALLEBRY: di.props_del_all('ebrys')
        elif evt==EVT_FEM_DELALLFBRY: di.props_del_all('fbrys')
        elif evt==EVT_FEM_DELALLEATT: di.props_del_all('eatts')

        elif evt==EVT_FEM_RUN:      fe.run_analysis()
        elif evt==EVT_FEM_SCRIPT:   fe.run_analysis(True)
        elif evt==EVT_FEM_PARAVIEW: fe.paraview    ()

        # ----------------------------------------------------------------------------------- RES 

        elif evt==EVT_RES_SHOWHIDE: di.toggle_key('gui_show_res')


    #except Exception, inst:
        #msg = inst.args[0]
        #print '[1;34mMechSys[0m: Error: '+'[1;31m'+msg+'[0m'
        #Blender.Draw.PupMenu('ERROR|'+msg)


# ================================================================================= Callbacks

# ---------------------------------- Settings

def cb_show_props(evt,val): di.set_key ('show_props', val)
def cb_show_e_ids(evt,val): di.set_key ('show_e_ids', val)
def cb_show_v_ids(evt,val): di.set_key ('show_v_ids', val)
def cb_show_blks (evt,val): di.set_key ('show_blks',  val)
def cb_show_axes (evt,val): di.set_key ('show_axes',  val)
def cb_show_etags(evt,val): di.set_key ('show_etags', val)
def cb_show_ftags(evt,val): di.set_key ('show_ftags', val)
def cb_show_elems(evt,val): di.set_key ('show_elems', val)
def cb_show_opac (evt,val): di.set_key ('show_opac',  val)

# ---------------------------------- CAD

def cb_set_x     (evt,val): di.set_key('cad_x',   val)
def cb_set_y     (evt,val): di.set_key('cad_y',   val)
def cb_set_z     (evt,val): di.set_key('cad_z',   val)
def cb_fillet_rad(evt,val): di.set_key('cad_rad', val)
def cb_fillet_stp(evt,val): di.set_key('cad_stp', val)

# ---------------------------------- Mesh

def cb_3dmesh(evt,val): di.props_set_val('3dmesh', val)
def cb_etag  (evt,val): di.set_key      ('newetag', [val, di.key('newetag')[1]])
def cb_ftag  (evt,val): di.set_key      ('newftag', [val, di.key('newftag')[1]])
def cb_fclr  (evt,val): di.set_key      ('newftag', [di.key('newftag')[0], di.rgb2hex(val)])

# ---------------------------------- Mesh -- structured

def cb_blk_tag(evt,val): di.props_set_item ('blks', evt-EVT_INC,  0, val)
def cb_blk_xax(evt,val): di.blk_set_axis   (        evt-EVT_INC,  1)
def cb_blk_yax(evt,val): di.blk_set_axis   (        evt-EVT_INC,  2)
def cb_blk_zax(evt,val): di.blk_set_axis   (        evt-EVT_INC,  3)
def cb_blk_nx (evt,val): di.props_set_item ('blks', evt-EVT_INC,  8, val)
def cb_blk_ny (evt,val): di.props_set_item ('blks', evt-EVT_INC,  9, val)
def cb_blk_nz (evt,val): di.props_set_item ('blks', evt-EVT_INC, 10, val)
def cb_blk_ax (evt,val): di.props_set_item ('blks', evt-EVT_INC, 11, float(val))
def cb_blk_ay (evt,val): di.props_set_item ('blks', evt-EVT_INC, 12, float(val))
def cb_blk_az (evt,val): di.props_set_item ('blks', evt-EVT_INC, 13, float(val))
def cb_blk_nlx(evt,val): di.props_set_item ('blks', evt-EVT_INC, 14, val)
def cb_blk_nly(evt,val): di.props_set_item ('blks', evt-EVT_INC, 15, val)
def cb_blk_nlz(evt,val): di.props_set_item ('blks', evt-EVT_INC, 16, val)
def cb_blk_del(evt,val): di.props_del      ('blks', evt-EVT_INC)

# ---------------------------------- Mesh -- unstructured

def cb_minang (evt,val): di.props_set_val('minang', float(val))
def cb_maxarea(evt,val): di.props_set_val('maxarea',float(val))

# ---------------------------------- Mesh -- unstructured -- regions

def cb_regs_tag    (evt,val): di.props_set_item ('regs', str(evt-EVT_INC), 0, val)
def cb_regs_maxarea(evt,val): di.props_set_item ('regs', str(evt-EVT_INC), 1, float(val))
def cb_regs_setx   (evt,val): di.props_set_item ('regs', str(evt-EVT_INC), 2, float(val))
def cb_regs_sety   (evt,val): di.props_set_item ('regs', str(evt-EVT_INC), 3, float(val))
def cb_regs_setz   (evt,val): di.props_set_item ('regs', str(evt-EVT_INC), 4, float(val))
def cb_regs_del    (evt,val): di.props_del      ('regs', str(evt-EVT_INC))

# ---------------------------------- Mesh -- unstructured -- holes

def cb_hols_setx(evt,val): di.props_set_item ('hols', str(evt-EVT_INC), 0, float(val))
def cb_hols_sety(evt,val): di.props_set_item ('hols', str(evt-EVT_INC), 1, float(val))
def cb_hols_setz(evt,val): di.props_set_item ('hols', str(evt-EVT_INC), 2, float(val))
def cb_hols_del (evt,val): di.props_del      ('hols', str(evt-EVT_INC))

# ---------------------------------- FEM -- nbrys

def cb_nbry_setx  (evt,val): di.props_set_item ('nbrys', str(evt-EVT_INC), 0, float(val))
def cb_nbry_sety  (evt,val): di.props_set_item ('nbrys', str(evt-EVT_INC), 1, float(val))
def cb_nbry_setz  (evt,val): di.props_set_item ('nbrys', str(evt-EVT_INC), 2, float(val))
def cb_nbry_setkey(evt,val): di.props_set_item ('nbrys', str(evt-EVT_INC), 3, val-1)
def cb_nbry_setval(evt,val): di.props_set_item ('nbrys', str(evt-EVT_INC), 4, float(val))
def cb_nbry_del   (evt,val): di.props_del      ('nbrys', str(evt-EVT_INC))

# ---------------------------------- FEM -- nbsID

def cb_nbID_setID (evt,val): di.props_set_item ('nbsID', str(evt-EVT_INC), 0, int(val))
def cb_nbID_setkey(evt,val): di.props_set_item ('nbsID', str(evt-EVT_INC), 1, val-1)
def cb_nbID_setval(evt,val): di.props_set_item ('nbsID', str(evt-EVT_INC), 2, float(val))
def cb_nbID_del   (evt,val): di.props_del      ('nbsID', str(evt-EVT_INC))

# ---------------------------------- FEM -- ebrys

def cb_ebry_settag(evt,val): di.props_set_item ('ebrys', str(evt-EVT_INC), 0, int(val))
def cb_ebry_setkey(evt,val): di.props_set_item ('ebrys', str(evt-EVT_INC), 1, val-1)
def cb_ebry_setval(evt,val): di.props_set_item ('ebrys', str(evt-EVT_INC), 2, float(val))
def cb_ebry_del   (evt,val): di.props_del      ('ebrys', str(evt-EVT_INC))

# ---------------------------------- FEM -- fbrys

def cb_fbry_settag(evt,val): di.props_set_item ('fbrys', str(evt-EVT_INC), 0, int(val))
def cb_fbry_setkey(evt,val): di.props_set_item ('fbrys', str(evt-EVT_INC), 1, val-1)
def cb_fbry_setval(evt,val): di.props_set_item ('fbrys', str(evt-EVT_INC), 2, float(val))
def cb_fbry_setclr(evt,val): di.props_set_item ('fbrys', str(evt-EVT_INC), 3, di.rgb2hex(val))
def cb_fbry_del   (evt,val): di.props_del      ('fbrys', str(evt-EVT_INC))

# ---------------------------------- FEM -- eatts

def cb_eatt_settag  (evt,val): di.props_set_item ('eatts', str(evt-EVT_INC), 0, int(val))
def cb_eatt_settype (evt,val): di.props_set_item ('eatts', str(evt-EVT_INC), 1, val-1)
def cb_eatt_setmodel(evt,val): di.props_set_item ('eatts', str(evt-EVT_INC), 2, val-1)
def cb_eatt_setmatID(evt,val): di.props_set_item ('eatts', str(evt-EVT_INC), 3, val)
def cb_eatt_setiniID(evt,val): di.props_set_item ('eatts', str(evt-EVT_INC), 4, val)
def cb_eatt_del     (evt,val): di.props_del      ('eatts', str(evt-EVT_INC))

# ---------------------------------- FEM

def cb_fem_fullsc(evt,val): di.set_key('fem_fullsc', val)

def cb_fem_struct(evt,val):
    di.set_key('fem_struct', val)
    di.set_key('fem_unstru', not val)
    Blender.Window.QRedrawAll()

def cb_fem_unstru(evt,val):
    di.set_key('fem_unstru', val)
    di.set_key('fem_struct', not val)
    Blender.Window.QRedrawAll()

# ---------------------------------- Results

def cb_res_show       (evt,val): di.set_key ('res_show',        val)
def cb_res_scalar     (evt,val): di.set_key ('res_scalar',      val)
def cb_res_show_scalar(evt,val): di.set_key ('res_show_scalar', val)
def cb_res_warp_scale (evt,val): di.set_key ('res_warp_scale',  val)
def cb_res_show_warp  (evt,val): di.set_key ('res_show_warp',   val)


# ======================================================================================= GUI

# Draw GUI
def gui():
    # load dictionary
    d = di.load_dict()

    # Get current selected mesh
    try:                    edm, obj, msh = di.get_msh (False)
    except Exception, inst: edm, obj, msh = 0, None, None

    # Data from current selected object
    is3d    = False
    blks    = {}
    minang  = -1.0
    maxarea = -1.0
    regs    = {}
    hols    = {}
    nbrys   = {}
    nbsID   = {}
    ebrys   = {}
    fbrys   = {}
    eatts   = {}
    if obj!=None:
        if obj.properties.has_key('3dmesh'): is3d = obj.properties['3dmesh']
        else:
            obj.properties['3dmesh'] = False
            is3d                     = False
        blks     = obj.properties['blks']    if obj.properties.has_key('blks')    else {}
        minang   = obj.properties['minang']  if obj.properties.has_key('minang')  else -1.0
        maxarea  = obj.properties['maxarea'] if obj.properties.has_key('maxarea') else -1.0
        regs     = obj.properties['regs']    if obj.properties.has_key('regs')    else {}
        hols     = obj.properties['hols']    if obj.properties.has_key('hols')    else {}
        nbrys    = obj.properties['nbrys']   if obj.properties.has_key('nbrys')   else {}
        nbsID    = obj.properties['nbsID']   if obj.properties.has_key('nbsID')   else {}
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
    h_fem_nbsID     = rh+srg+rh*len(nbsID)
    h_fem_ebrys     = rh+srg+rh*len(ebrys)
    h_fem_fbrys     = rh+srg+rh*len(fbrys)
    h_fem_eatts     = rh+srg+rh*len(eatts)
    h_fem           = 8*rh+5*srg+h_fem_nbrys+h_fem_nbsID+h_fem_ebrys+h_fem_fbrys+h_fem_eatts
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
        Draw.Toggle      ('3D mesh', EVT_NONE,          c,     r, 60, rh, is3d,                        'Set 3D mesh',                       cb_3dmesh)
        Draw.Number      ('',        EVT_NONE,          c+ 60, r, 60, rh, d['newetag'][0],    -100, 0, 'New edge tag',                      cb_etag)
        Draw.PushButton  ('Edge',    EVT_MESH_SETETAG,  c+120, r, 60, rh,                              'Set edges tag (0 => remove tag)')
        Draw.Number      ('',        EVT_NONE,          c+180, r, 60, rh, d['newftag'][0],   -1000, 0, 'New face tag',                      cb_ftag)
        Draw.ColorPicker (           EVT_NONE,          c+240, r, 60, rh, di.hex2rgb(d['newftag'][1]), 'Select color to paint tagged face', cb_fclr)
        Draw.PushButton  ('Face',    EVT_MESH_SETFTAG,  c+300, r, 60, rh,                              'Set faces tag (0 => remove tag)')
        r -= rh
        r -= rh

        # ----------------------- Mesh -- structured

        gu.caption2(c,r,w,rh,'Structured mesh')
        r, c, w = gu.box2_in(W,cg,rh, c,r,w,h_msh_stru)

        # ----------------------- Mesh -- structured -- blocks

        gu.caption3(c,r,w,rh, 'Blocks', EVT_MESH_ADDBLK,EVT_MESH_DELALLBLKS)
        r, c, w = gu.box3_in(W,cg,rh, c,r,w,h_msh_stru_blks)
        gu.text(c,r,'     ID:Tag     Local axes      nX           nY          nZ')
        for k, v in blks.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     (str(i)+':', EVT_INC+i, c     , r-rh, 80, 2*rh, int(v[0]), -100, -1, 'Block tag',                       cb_blk_tag)
            Draw.PushButton ('X',        EVT_INC+i, c+ 80 , r,    20,   rh,                      'Set X-axis',                      cb_blk_xax)
            Draw.PushButton ('Y',        EVT_INC+i, c+100 , r,    20,   rh,                      'Set Y-axis',                      cb_blk_yax)
            Draw.PushButton ('Z',        EVT_INC+i, c+120 , r,    20,   rh,                      'Set Z-axis',                      cb_blk_zax)
            Draw.Number     ('',         EVT_INC+i, c+140 , r,    60,   rh, int(v[ 8]), 1, 1000, 'Number of divisions along X',     cb_blk_nx)
            Draw.Number     ('',         EVT_INC+i, c+200 , r,    60,   rh, int(v[ 9]), 1, 1000, 'Number of divisions along Y',     cb_blk_ny)
            Draw.Number     ('',         EVT_INC+i, c+260 , r,    60,   rh, int(v[10]), 1, 1000, 'Number of divisions along Z',     cb_blk_nz)
            r -= rh                                          
            Draw.Toggle     ('X',        EVT_INC+i, c+ 80 , r,    20,   rh, int(v[14]),          'Set nonlinear divisions along X', cb_blk_nlx)
            Draw.Toggle     ('Y',        EVT_INC+i, c+100 , r,    20,   rh, int(v[15]),          'Set nonlinear divisions along Y', cb_blk_nly)
            Draw.Toggle     ('Z',        EVT_INC+i, c+120 , r,    20,   rh, int(v[16]),          'Set nonlinear divisions along Z', cb_blk_nlz)
            Draw.String     ('aX=',      EVT_INC+i, c+140 , r,    60,   rh, str(v[11]),      32, 'Set division coefficient aX',     cb_blk_ax)
            Draw.String     ('aY=',      EVT_INC+i, c+200 , r,    60,   rh, str(v[12]),      32, 'Set division coefficient aY',     cb_blk_ay)
            Draw.String     ('aZ=',      EVT_INC+i, c+260 , r,    60,   rh, str(v[13]),      32, 'Set division coefficient aZ',     cb_blk_az)
            Draw.PushButton ('Del',      EVT_INC+i, c+320 , r,    40, 2*rh,                      'Delete this row',                 cb_blk_del)
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
        gu.text(c,r,'     ID:Tag      max area        X             Y            Z')
        for k, v in regs.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     (str(i)+':', EVT_INC+i, c,     r, 80, rh, int(v[0]), -100,-1,'Region tag',                  cb_regs_tag)
            Draw.String     ('',         EVT_INC+i, c+ 80, r, 60, rh, str(v[1]),   32,   'Max area (-1 => use default)',cb_regs_maxarea)
            Draw.String     ('',         EVT_INC+i, c+140, r, 60, rh, str(v[2]),   32,   'X of the region',             cb_regs_setx)
            Draw.String     ('',         EVT_INC+i, c+200, r, 60, rh, str(v[3]),   32,   'Y of the region',             cb_regs_sety)
            Draw.String     ('',         EVT_INC+i, c+260, r, 60, rh, str(v[4]),   32,   'Z of the region',             cb_regs_setz)
            Draw.PushButton ('Del',      EVT_INC+i, c+320, r, 40, rh,                    'Delete this row',             cb_regs_del)
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
            Draw.String     ('',     EVT_INC+i, c+ 40, r, 60, rh, str(v[0]), 32, 'X of the hole',   cb_hols_setx)
            Draw.String     ('',     EVT_INC+i, c+100, r, 60, rh, str(v[1]), 32, 'Y of the hole',   cb_hols_sety)
            Draw.String     ('',     EVT_INC+i, c+160, r, 60, rh, str(v[2]), 32, 'Z of the hole',   cb_hols_setz)
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

        # ----------------------- FEM -- nbsID

        r -= srg
        gu.caption2_(c,r,w,rh,'Nodes boundary conditions (given IDs)', EVT_FEM_ADDNBID,EVT_FEM_DELALLNBID)
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_nbsID)
        gu.text(c,r,'     ID        Key     Value')
        for k, v in nbsID.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     ('',                EVT_INC+i, c    , r, 60, rh, int(v[0]),0,10000, 'Set the node ID',                                                            cb_nbID_setID)
            Draw.Menu       (d['dofvars_menu'], EVT_INC+i, c+ 60, r, 40, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_nbID_setkey)
            Draw.String     ('',                EVT_INC+i, c+100, r, 60, rh, str(v[2]), 32,     'Value of essential/natural boundary condition',                              cb_nbID_setval)
            Draw.PushButton ('Del',             EVT_INC+i, c+160, r, 40, rh,                    'Delete this row',                                                            cb_nbID_del)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- ebrys

        r -= srg
        gu.caption2_(c,r,w,rh,'Edges boundary conditions', EVT_FEM_ADDEBRY,EVT_FEM_DELALLEBRY)
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_ebrys)
        gu.text(c,r,'    Tag       Key     Value')
        for k, v in ebrys.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     ('',                EVT_INC+i, c,     r, 60, rh, int(v[0]),-1000,-1,'Set tag',                                                                    cb_ebry_settag)
            Draw.Menu       (d['dofvars_menu'], EVT_INC+i, c+ 60, r, 40, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_ebry_setkey)
            Draw.String     ('',                EVT_INC+i, c+100, r, 60, rh, str(v[2]), 128,    'Value of essential/natural boundary condition',                              cb_ebry_setval)
            Draw.PushButton ('Del',             EVT_INC+i, c+160, r, 40, rh,                    'Delete this row',                                                            cb_ebry_del)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- fbrys

        r -= srg
        gu.caption2_(c,r,w,rh,'Faces boundary conditions', EVT_FEM_ADDFBRY,EVT_FEM_DELALLFBRY)
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_fbrys)
        gu.text(c,r,'    Tag        Colour     Key     Value')
        for k, v in fbrys.iteritems():
            r  -= rh
            i   = int(k)
            clr = di.hex2rgb(v[3])
            Draw.Number      ('',                EVT_INC+i, c,     r, 60, rh, int(v[0]),-1000,-1,'Set tag',                                                                    cb_fbry_settag)
            Draw.ColorPicker (                   EVT_INC+i, c+ 60, r, 60, rh, clr,               'Select color to paint tagged face',                                          cb_fbry_setclr)
            Draw.Menu        (d['dofvars_menu'], EVT_INC+i, c+120, r, 40, rh, int(v[1])+1,       'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', cb_fbry_setkey)
            Draw.String      ('',                EVT_INC+i, c+160, r, 60, rh, str(v[2]), 128,    'Value of essential/natural boundary condition',                              cb_fbry_setval)
            Draw.PushButton  ('Del',             EVT_INC+i, c+220, r, 40, rh,                    'Delete this row',                                                            cb_fbry_del)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- eatts

        r -= srg
        gu.caption2_(c,r,w,rh,'Elements attributes', EVT_FEM_ADDEATT,EVT_FEM_DELALLEATT)
        r, c, w = gu.box2__in(W,cg,rh, c,r,w,h_fem_eatts)
        gu.text(c,r,'    Tag               Type                    Model')
        for k, v in eatts.iteritems():
            r -= rh
            i  = int(k)
            Draw.Number     ('',               EVT_INC+i, c,     r,  60, rh, int(v[0]),-1000,0,         'Set tag',                                   cb_eatt_settag)
            Draw.Menu       (d['etypes_menu'], EVT_INC+i, c+ 60, r, 120, rh, int(v[1])+1,               'Element type: ex.: Quad4PStrain',           cb_eatt_settype)
            Draw.Menu       (d['models_menu'], EVT_INC+i, c+180, r, 100, rh, int(v[2])+1,               'Constitutive model: ex.: LinElastic',       cb_eatt_setmodel)
            Draw.PushButton ('Del',            EVT_INC+i, c+280, r,  40, rh,                            'Delete this row',                           cb_eatt_del)
            #Draw.String     ('',               EVT_INC+i, c+280, r, 100, rh, m[2].replace('_',' '),128, 'Parameters: ex.: E=200 nu=0.25',            cb_eatt_setprms)
            #Draw.String     ('',               EVT_INC+i, c+380, r,  80, rh, m[3].replace('_',' '),128, 'Initial values: ex.: Sx=0 Sy=0 Sz=0 Sxy=0', cb_eatt_setinis)
        r -= srg
        r, c, w = gu.box2__out(W,cg,rh, c,r)

        # ----------------------- FEM -- END
        r -= srg
        Draw.PushButton ('Run analysis',     EVT_FEM_RUN,      c,     r, 100, rh, 'Run a FE analysis directly (without script)')
        Draw.Toggle     ('Full script',      EVT_NONE,         c+100, r,  60, rh, d['fem_fullsc'], 'Show results', cb_fem_fullsc)
        Draw.Toggle     ('Struct',           EVT_NONE,         c+160, r,  40, rh, d['fem_struct'], 'Show results', cb_fem_struct)
        Draw.Toggle     ('Unstruct',         EVT_NONE,         c+200, r,  40, rh, d['fem_unstru'], 'Show results', cb_fem_unstru)
        Draw.PushButton ('Write script',     EVT_FEM_SCRIPT,   c+240, r,  80, rh, 'Generate script for FEM')
        Draw.PushButton ('View in ParaView', EVT_FEM_PARAVIEW, c+320, r, 120, rh, 'View results in ParaView')
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
