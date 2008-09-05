#!BPY

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
import msys_draw as dr
import msys_dict as di
import msys_fem  as fem


# ================================================================================= Constants

EVT_NONE             =  0 # used for buttons with callbacks
EVT_REFRESH          =  1 # refresh all windows
# CAD
EVT_CAD_ADDXYZ       =  2 # type in values to define point(s)
EVT_CAD_FILLET       =  3 # 3D fillet
EVT_CAD_EDGEBREAK    =  4 # break edge
EVT_CAD_EDGEBREAKM   =  5 # break edge at middle point
EVT_CAD_EDGEINTERS   = 55 # edge closest distance
EVT_CAD_FPOINT       =  6 # read points from file
EVT_CAD_FSPLINE      =  7 # create a spline from points in a file
# Mesh
EVT_MESH_SETETAG     =  8 # set edges tag
EVT_MESH_SETFTAG     =  9 # set faces tag
EVT_MESH_DELPROP     = 10 # delete all properties
# Mesh -- structured
EVT_MESH_SETX        = 11 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_SETY        = 12 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_SETZ        = 13 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_GENSTRU     = 14 # generate structured mesh using MechSys module
# Mesh -- unstructured
EVT_MESH_GENUNSTRU   = 15 # generate structured mesh using MechSys module
EVT_MESH_ADDREG      = 16 # add region
EVT_MESH_DELALLREGS  = 17 # delete all regions
EVT_MESH_ADDHOLE     = 18 # add hole
EVT_MESH_DELALLHOLES = 19 # delete all holes
# FEA
EVT_FEA_ADDNBRY      = 20 # add nodes boundary
EVT_FEA_DELALLNBRY   = 21 # delete all nodes boundary
EVT_FEA_ADDEBRY      = 22 # add edges boundary
EVT_FEA_DELALLEBRY   = 23 # delete all edges boundary
EVT_FEA_ADDFBRY      = 24 # add faces boundary
EVT_FEA_DELALLFBRY   = 25 # delete all faces boundary
EVT_FEA_ADDEATT      = 26 # add element attributes
EVT_FEA_DELALLEATT   = 27 # delete all elements attributes
EVT_FEA_RUN          = 28 # run a FE simulation
EVT_FEA_SCRIPT       = 29 # generate script for FEA 


# ==================================================================================== Events

# Handle input events
def event(evt, val):
    if   evt==Draw.QKEY and not val: Draw.Exit()
    elif evt==Draw.ESCKEY: Draw.Redraw(1)
    elif evt==Draw.WHEELDOWNMOUSE:
        d = di.load_dict()
        d['inirow'] -= 40
        Blender.Window.QRedrawAll()
    elif evt==Draw.WHEELUPMOUSE:
        d = di.load_dict()
        if d['inirow']<0:
            d['inirow'] += 40
            Blender.Window.QRedrawAll()


# Handle button events
@dr.print_timing
def button_event(evt):

    # load dictionary
    dict = di.load_dict()

    # refresh all windows
    if evt==EVT_REFRESH: Blender.Window.QRedrawAll()

    # ----------------------------------------------------------------------------------- CAD

    # add vertices
    if evt==EVT_CAD_ADDXYZ:
        x = dict['newpoint_x']
        y = dict['newpoint_y']
        z = dict['newpoint_z']
        dr.add_point (float(x), float(y), float(z))

    # fillet two edges
    elif evt==EVT_CAD_FILLET: dr.fillet(float(dict['fillet_radius']),dict['fillet_steps'])

    # break an edge
    elif evt==EVT_CAD_EDGEBREAK: dr.break_edge()

    # break an edge at middle point
    elif evt==EVT_CAD_EDGEBREAKM: dr.break_edge(True)

    # edge closest distance
    elif evt==EVT_CAD_EDGEINTERS: dr.edge_intersect()

    # read vertices from file
    elif evt==EVT_CAD_FPOINT: Blender.Window.FileSelector(dr.add_points_from_file, "Read X Y Z cols")

    # read spline from file
    elif evt==EVT_CAD_FSPLINE: Blender.Window.FileSelector(dr.add_spline_from_file, "Read X Y Z cols")

    # ---------------------------------------------------------------------------------- Mesh

    # set edges tag
    elif evt==EVT_MESH_SETETAG:
        edm, obj, msh = di.get_msh()
        for id in msh.edges.selected():
            di.set_etag (obj, id, dict['newetag'])
        Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set faces tag
    elif evt==EVT_MESH_SETFTAG:
        edm, obj, msh = di.get_msh()
        nedges = len(msh.edges.selected())
        if nedges==3 or nedges==6 or nedges==4 or nedges==8:
            di.set_ftag (obj, msh.edges.selected(), dict['newftag'])
        Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # delete all properties
    elif evt==EVT_MESH_DELPROP:
        scn = bpy.data.scenes.active
        obs = scn.objects.selected
        for o in obs:
            for p in o.getAllProperties():
                print '[1;34mMechSys[0m: Object = ', o.name, ' Deleted property: ', p.name, p.data
                o.removeProperty(p)
            Blender.Window.QRedrawAll()

    # ---------------------------------------------------------------------------- Structured

    # set local x-axis
    elif evt==EVT_MESH_SETX:
        edm, obj, msh = di.get_msh()
        if len(msh.edges.selected())==1:
            di.set_local_axis (obj, 'x', msh.edges.selected()[0])
            Blender.Window.QRedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select only one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set local y-axis
    elif evt==EVT_MESH_SETY:
        edm, obj, msh = di.get_msh()
        if len(msh.edges.selected())==1:
            di.set_local_axis (obj, 'y', msh.edges.selected()[0])
            Blender.Window.QRedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select only one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set local z-axis
    elif evt==EVT_MESH_SETZ:
        edm, obj, msh = di.get_msh()
        if len(msh.edges.selected())==1:
            di.set_local_axis (obj, 'z', msh.edges.selected()[0])
            Blender.Window.QRedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select onlz one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # generate structured mesh via MechSys
    elif evt==EVT_MESH_GENSTRU:
        tt = timeit.Timer('msys_draw.gen_struct_mesh()', 'import msys_draw')
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
            di.del_all_hols(obj)
            Blender.Window.QRedrawAll()

    # generate unstructured mesh via MechSys
    elif evt==EVT_MESH_GENUNSTRU:
        tt = timeit.Timer('msys_draw.gen_unstruct_mesh()', 'import msys_draw')
        t  = tt.timeit(number=1)
        print '[1;34mMechSys[0m: time spent on generation and drawing = [1;31m',t,'[0m [1;32mseconds[0m'

    # ----------------------------------------------------------------------------------- FEA 

    # add nodes boundary
    elif evt==EVT_FEA_ADDNBRY:
        obj = di.get_obj()
        if obj!=None:
            nbrys = di.get_nbrys (obj)
            di.set_nbry (obj, len(nbrys), '0.0','0.0','0.0', 'uy', '0.0')
            Blender.Window.QRedrawAll()

    # delete all nodes boundary
    elif evt==EVT_FEA_DELALLNBRY:
        obj = di.get_obj()
        if obj!=None:
            di.del_all_nbrys(obj)
            Blender.Window.QRedrawAll()

    # add edges boundary
    elif evt==EVT_FEA_ADDEBRY:
        obj = di.get_obj()
        if obj!=None:
            ebrys = di.get_ebrys (obj)
            di.set_ebry (obj, len(ebrys), '-10', 'uy', '0.0')
            Blender.Window.QRedrawAll()

    # delete all edges boundary
    elif evt==EVT_FEA_DELALLEBRY:
        obj = di.get_obj()
        if obj!=None:
            di.del_all_ebrys(obj)
            Blender.Window.QRedrawAll()

    # add faces boundary
    elif evt==EVT_FEA_ADDFBRY:
        obj = di.get_obj()
        if obj!=None:
            fbrys = di.get_fbrys (obj)
            di.set_fbry (obj, len(fbrys), '-100', 'uy', '0.0')
            Blender.Window.QRedrawAll()

    # delete all faces boundary
    elif evt==EVT_FEA_DELALLFBRY:
        obj = di.get_obj()
        if obj!=None:
            di.del_all_fbrys(obj)
            Blender.Window.QRedrawAll()

    # add elems attributes
    elif evt==EVT_FEA_ADDEATT:
        obj = di.get_obj()
        if obj!=None:
            eatts = di.get_eatts (obj)
            di.set_eatt (obj, len(eatts), '-1', 'Quad4PStrain', 'LinElastic', 'E=207 nu=0.3', 'Sx=0 Sy=0 Sz=0 Sxy=0')
            Blender.Window.QRedrawAll()

    # delete all elems attributes
    elif evt==EVT_FEA_DELALLEATT:
        obj = di.get_obj()
        if obj!=None:
            di.del_all_eatts(obj)
            Blender.Window.QRedrawAll()

    # run a FE simulation
    elif evt==EVT_FEA_RUN:
        edm, obj, msh = di.get_msh()
        fem.run_fea (obj)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    elif evt==EVT_FEA_SCRIPT:
        edm, obj, msh = di.get_msh()
        fem.gen_script (obj)
        Blender.Window.Redraw(Blender.Window.Types.TEXT)


# ================================================================================= Callbacks

# ---------------------------------- CAD

def setx_callback(evt,val):
    dict = di.load_dict()
    dict['newpoint_x'] = val

def sety_callback(evt,val):
    dict = di.load_dict()
    dict['newpoint_y'] = val

def setz_callback(evt,val):
    dict = di.load_dict()
    dict['newpoint_z'] = val

def fillet_radius_callback(evt,val):
    dict = di.load_dict()
    dict['fillet_radius'] = val

def fillet_steps_callback(evt,val):
    dict = di.load_dict()
    dict['fillet_steps'] = val


# ---------------------------------- Mesh

def show_props_callback(evt,val):
    dict = di.load_dict()
    dict['show_props'] = val
    Blender.Window.QRedrawAll()

def show_axes_callback(evt,val):
    dict = di.load_dict()
    dict['show_axes'] = val
    Blender.Window.QRedrawAll()

def show_v_ids_callback(evt,val):
    dict = di.load_dict()
    dict['show_v_ids'] = val
    Blender.Window.QRedrawAll()

def show_e_ids_callback(evt,val):
    dict = di.load_dict()
    dict['show_e_ids'] = val
    Blender.Window.QRedrawAll()

def show_f_ids_callback(evt,val):
    dict = di.load_dict()
    dict['show_f_ids'] = val
    Blender.Window.QRedrawAll()

def show_etags_callback(evt,val):
    dict = di.load_dict()
    dict['show_etags'] = val
    Blender.Window.QRedrawAll()

def show_ftags_callback(evt,val):
    dict = di.load_dict()
    dict['show_ftags'] = val
    Blender.Window.QRedrawAll()

def show_elems_callback(evt,val):
    dict = di.load_dict()
    dict['show_elems'] = val
    Blender.Window.QRedrawAll()

def ndivx_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_ndiv (o, 'x', val)
    Blender.Window.QRedrawAll()

def ndivy_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_ndiv (o, 'y', val)
    Blender.Window.QRedrawAll()

def ndivz_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_ndiv (o, 'z', val)
    Blender.Window.QRedrawAll()

def acoefx_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_acoef (o, 'x', val)
    Blender.Window.QRedrawAll()

def acoefy_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_acoef (o, 'y', val)
    Blender.Window.QRedrawAll()

def acoefz_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_acoef (o, 'z', val)
    Blender.Window.QRedrawAll()

def nonlinx_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_nonlin (o, 'x', val)
    Blender.Window.QRedrawAll()

def nonliny_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_nonlin (o, 'y', val)
    Blender.Window.QRedrawAll()

def nonlinz_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_nonlin (o, 'z', val)
    Blender.Window.QRedrawAll()


# ---------------------------------- etag, ftag, btag

def etag_callback(evt,val):
    dict = di.load_dict()
    dict['newetag'] = val

def ftag_callback(evt,val):
    dict = di.load_dict()
    dict['newftag'] = val

def btag_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_btag (o, val)
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
        di.set_reg_tag (obj, evt-1000, val)

def regs_maxarea_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_maxarea (obj, evt-1000, val)

def regs_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_x (obj, evt-1000, val)
        Blender.Window.QRedrawAll()

def regs_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_y (obj, evt-1000, val)
        Blender.Window.QRedrawAll()

def regs_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_z (obj, evt-1000, val)
        Blender.Window.QRedrawAll()

def regs_delonereg_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_reg (obj, evt-1000)
        Blender.Window.QRedrawAll()


# ---------------------------------- Holes

def holes_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_x (obj, evt-2000, val)
        Blender.Window.QRedrawAll()

def holes_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_y (obj, evt-2000, val)
        Blender.Window.QRedrawAll()

def holes_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_z (obj, evt-2000, val)
        Blender.Window.QRedrawAll()

def holes_delonehole_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_hol (obj, evt-2000)
        Blender.Window.QRedrawAll()


# ---------------------------------- nbrys

def nodebry_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_x (obj, evt-3000, val)

def nodebry_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_y (obj, evt-3000, val)

def nodebry_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_z (obj, evt-3000, val)

def nodebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_key (obj, evt-3000, val)

def nodebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_val (obj, evt-3000, val)

def nodebry_delonenbry_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_nbry (obj,evt-3000)
        Blender.Window.QRedrawAll()


# ---------------------------------- ebrys

def edgebry_settag_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_ebry_tag (obj, evt-4000, str(val))

def edgebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_ebry_key (obj, evt-4000, val)

def edgebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_ebry_val (obj, evt-4000, val)

def edgebry_deloneebry_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_ebry (obj,evt-4000)
        Blender.Window.QRedrawAll()


# ---------------------------------- fbrys

def facebry_settag_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_fbry_tag (obj, evt-5000, str(val))

def facebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_fbry_key (obj, evt-5000, val)

def facebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_fbry_val (obj, evt-5000, val)

def facebry_delonefbry_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_fbry (obj,evt-5000)
        Blender.Window.QRedrawAll()


# ---------------------------------- eatts

def elematt_settag_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_eatt_tag (obj, evt-6000, str(val))

def elematt_settype_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_eatt_type (obj, evt-6000, val)

def elematt_setmodel_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_eatt_model (obj, evt-6000, val)

def elematt_setprms_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_eatt_prms (obj, evt-6000, val)

def elematt_setinis_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_eatt_inis (obj, evt-6000, val)

def elematt_deloneeatt_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_eatt (obj,evt-6000)
        Blender.Window.QRedrawAll()


# ======================================================================================= GUI

# Draw GUI
def gui():
    # load dictionary
    d = di.load_dict()

    # Data from current selected object
    btag     = -1
    ndivx    = 2
    ndivy    = 2
    ndivz    = 1
    acoefx   = '0'
    acoefy   = '0'
    acoefz   = '0'
    nonlinx  = 0
    nonliny  = 0
    nonlinz  = 0
    minangle = '-1.0'
    maxarea  = '-1.0'
    regs     = []
    hols     = []
    nbrys    = []
    ebrys    = []
    fbrys    = []
    eatts    = []
    edm, obj, msh = di.get_msh (False)
    if obj!=None:
        btag     = di.get_btag     (obj)
        ndivx    = di.get_ndiv     (obj,'x')
        ndivy    = di.get_ndiv     (obj,'y')
        ndivz    = di.get_ndiv     (obj,'z')
        acoefx   = di.get_acoef    (obj,'x')
        acoefy   = di.get_acoef    (obj,'y')
        acoefz   = di.get_acoef    (obj,'z')
        nonlinx  = di.get_nonlin   (obj,'x')
        nonliny  = di.get_nonlin   (obj,'y')
        nonlinz  = di.get_nonlin   (obj,'z')
        minangle = di.get_minangle (obj)
        maxarea  = di.get_maxarea  (obj)
        regs     = di.get_regs     (obj)
        hols     = di.get_hols     (obj)
        nbrys    = di.get_nbrys    (obj)
        ebrys    = di.get_ebrys    (obj)
        fbrys    = di.get_fbrys    (obj)
        eatts    = di.get_eatts    (obj)
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
    h_cad       = 4*rh+ggy
    h_struct    = 6*rh+ggy+sgy
    h_unst_regs = (1+len(regs))*rh+sgy
    h_unst_hols = (1+len(hols))*rh+sgy
    h_unstruct  = ggy+3*sgy+2*rh + rh+h_unst_regs + rh+h_unst_hols
    h_mesh      = 2*ggy+4*rh+sgy + rh+h_struct + rh+h_unstruct
    h_fea_nbrys = (1+len(nbrys))*rh+sgy
    h_fea_ebrys = (1+len(ebrys))*rh+sgy
    h_fea_fbrys = (1+len(fbrys))*rh+sgy
    h_fea_eatts = (1+len(eatts))*rh+sgy
    h_fea       = ggy+3*rh + rh+h_fea_nbrys + rh+h_fea_ebrys + rh+h_fea_fbrys + rh+h_fea_eatts

    # Background color
    BGL.glClearColor (0.531, 0.543, 0.614, 0.0)
    BGL.glClear      (Blender.BGL.GL_COLOR_BUFFER_BIT)
    
    # CAD
    dx  = 2*gx+50;
    row = hei-gy-d['inirow']
    h   = h_cad
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('CAD')
    BGL.glColor3f     (0.85, 0.85, 0.85)
    BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Points:')
    Draw.String       ('X=',        EVT_NONE      , dx    , row, 80, rh, d['newpoint_x'],128,'Set the X value of the new point to be added',setx_callback)
    Draw.String       ('Y=',        EVT_NONE      , dx+ 80, row, 80, rh, d['newpoint_y'],128,'Set the Y value of the new point to be added',sety_callback)
    Draw.String       ('Z=',        EVT_NONE      , dx+160, row, 80, rh, d['newpoint_z'],128,'Set the Z value of the new point to be added',setz_callback)
    Draw.PushButton   ('Add x y z', EVT_CAD_ADDXYZ, dx+240, row, 80, rh, 'Add point from coordinates'); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Edges:');
    Draw.String       ('Radius=', EVT_NONE      , dx     , row , 120 , rh , d['fillet_radius'],  128,'Set radius for fillet operation',fillet_radius_callback)
    Draw.Number       ('Steps=',  EVT_NONE      , dx+120 , row , 120 , rh , d['fillet_steps' ],1,100,'Set steps for fillet operation' ,fillet_steps_callback)
    Draw.PushButton   ('Fillet',  EVT_CAD_FILLET, dx+240 , row ,  80 , rh , 'Create a fillet between two edges'); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Edges:');
    Draw.PushButton   ('Break edge',        EVT_CAD_EDGEBREAK , dx,      row ,  120, rh , 'Break an edge at a previously selected point')
    Draw.PushButton   ('Break at middle',   EVT_CAD_EDGEBREAKM, dx+120 , row ,  120, rh , 'Break an edge at its middle point')
    Draw.PushButton   ('Edge intersection', EVT_CAD_EDGEINTERS, dx+240 , row ,  120, rh , 'Find the intersection (smaller distance) between two edges'); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('File:')
    Draw.PushButton   ('Read Points', EVT_CAD_FPOINT,  dx,     row, 120, rh, 'Add points by reading a list from file')
    Draw.PushButton   ('Read Spline', EVT_CAD_FSPLINE, dx+120, row, 120, rh, 'Add a spline by reading its points from file')

    # Mesh
    dx   = 2*gx+55;
    row -= ggy
    h    = h_mesh
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('Mesh')
    BGL.glColor3f     (0.85, 0.85, 0.85)
    BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Show:')
    Draw.Toggle       ('Props', EVT_NONE, dx    , row, 60, rh, d['show_props'], 'Show mesh properties'      , show_props_callback)
    Draw.Toggle       ('V IDs', EVT_NONE, dx+ 60, row, 60, rh, d['show_v_ids'], 'Show vertex IDs'           , show_v_ids_callback)
    Draw.Toggle       ('E IDs', EVT_NONE, dx+120, row, 60, rh, d['show_e_ids'], 'Show edge IDs'             , show_e_ids_callback)
    Draw.Toggle       ('F IDs', EVT_NONE, dx+180, row, 60, rh, d['show_f_ids'], 'Show face IDs'             , show_f_ids_callback); row -= rh
    Draw.Toggle       ('Axes' , EVT_NONE, dx    , row, 60, rh, d['show_axes'],  'Show local system axes'    , show_axes_callback )
    Draw.Toggle       ('ETags', EVT_NONE, dx+ 60, row, 60, rh, d['show_etags'], 'Show edge tags'            , show_etags_callback)
    Draw.Toggle       ('FTags', EVT_NONE, dx+120, row, 60, rh, d['show_ftags'], 'Show face tags'            , show_ftags_callback)
    Draw.Toggle       ('Elems', EVT_NONE, dx+180, row, 60, rh, d['show_elems'], 'Show elements information' , show_elems_callback); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Set tags:')
    Draw.Number       ('',     EVT_NONE,          dx,     row, 80, rh, d['newetag'], -1000, 0,'New edge tag',etag_callback)
    Draw.PushButton   ('Edge', EVT_MESH_SETETAG,  dx+80,  row, 80, rh,                        'Set edges tag (0 => remove tag)')
    Draw.Number       ('',     EVT_NONE,          dx+160, row, 80, rh, d['newftag'], -1000, 0,'New face tag',ftag_callback)
    Draw.PushButton   ('Face', EVT_MESH_SETFTAG,  dx+240, row, 80, rh,                        'Set faces tag (0 => remove tag)'); row -= rh+sgy
    Draw.PushButton   ('Delete all properties', EVT_MESH_DELPROP, ggx, row, 200, rh , 'Delete all properties')

    # Mesh -- structured
    row -= gy
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
    BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
    BGL.glColor3f     (0.0, 0.0, 0.0)
    Draw.PushButton   ('Refresh', EVT_REFRESH, ggx+120, row+2, 60, rh-4, 'Refresh GUI'); row -= rh+gy
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Local Axes:')
    Draw.PushButton   ('Set x', EVT_MESH_SETX, dx    , row, 60, rh , 'Set local x axis (one edge must be previously selected)')
    Draw.PushButton   ('Set y', EVT_MESH_SETY, dx+ 60, row, 60, rh , 'Set local y axis (one edge must be previously selected)')
    Draw.PushButton   ('Set z', EVT_MESH_SETZ, dx+120, row, 60, rh , 'Set local z axis (one edge must be previously selected)'); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Divisions:')
    Draw.Number       ('nX=',    EVT_NONE, dx,     row, 80, rh, ndivx,  1, 1000,'Set number of divisions along x (local axis)',ndivx_callback)
    Draw.Number       ('nY=',    EVT_NONE, dx+ 80, row, 80, rh, ndivy,  1, 1000,'Set number of divisions along y (local axis)',ndivy_callback)
    Draw.Number       ('nZ=',    EVT_NONE, dx+160, row, 80, rh, ndivz,  1, 1000,'Set number of divisions along z (local axis)',ndivz_callback); row -= rh
    Draw.String       ('aX=',    EVT_NONE, dx    , row, 80, rh, acoefx, 32,'Set the a coeficient of the divisions function for the x direction (0 => equal size divisions)',acoefx_callback)
    Draw.String       ('aY=',    EVT_NONE, dx+ 80, row, 80, rh, acoefy, 32,'Set the a coeficient of the divisions function for the x direction (0 => equal size divisions)',acoefy_callback)
    Draw.String       ('aZ=',    EVT_NONE, dx+160, row, 80, rh, acoefz, 32,'Set the a coeficient of the divisions function for the x direction (0 => equal size divisions)',acoefz_callback); row -= rh
    Draw.Toggle       ('NonLinX',EVT_NONE, dx,     row, 80, rh, nonlinx, 'Set non-linear divisions along x', nonlinx_callback)
    Draw.Toggle       ('NonLinY',EVT_NONE, dx+ 80, row, 80, rh, nonliny, 'Set non-linear divisions along y', nonliny_callback)
    Draw.Toggle       ('NonLinZ',EVT_NONE, dx+160, row, 80, rh, nonlinz, 'Set non-linear divisions along z', nonlinz_callback); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Block:')
    Draw.Number       ('Tag=', EVT_NONE, dx, row, 100, rh, btag, -1000, 0,'Set the block tag (to be replicated to all elements)',btag_callback); row -= rh+sgy
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
    BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
    BGL.glColor3f     (0.0, 0.0, 0.0)
    Draw.PushButton   ('Refresh', EVT_REFRESH, ggx+120, row+2, 60, rh-4, 'Refresh GUI'); row -= rh+gy
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
    Draw.PushButton   ('Add',        EVT_MESH_ADDREG,     ggx+120, row+2, 60, rh-4, 'Add region')
    Draw.PushButton   ('Delete all', EVT_MESH_DELALLREGS, ggx+200, row+2, 80, rh-4, 'Delete all regions'); row -= rh
    BGL.glRasterPos2i (ggx, row+5); Draw.Text('     Tag        Max area            X                  Y                  Z'); row -= rh
    for i, reg in enumerate(regs):
        Draw.Number     ('',    1000+i, ggx    , row, 60, rh, int(reg[0]), -1000, 0,'Region tag',                  regs_tag_callback)
        Draw.String     ('',    1000+i, ggx+ 60, row, 80, rh,     reg[1],    128,   'Max area (-1 => use default)',regs_maxarea_callback)
        Draw.String     ('',    1000+i, ggx+140, row, 80, rh,     reg[2],    128,   'X of the rege',               regs_setx_callback)
        Draw.String     ('',    1000+i, ggx+220, row, 80, rh,     reg[3],    128,   'Y of the rege',               regs_sety_callback)
        Draw.String     ('',    1000+i, ggx+300, row, 80, rh,     reg[4],    128,   'Z of the rege',               regs_setz_callback)
        Draw.PushButton ('Del', 1000+i, ggx+380, row, 40, rh,                       'Delete this row',             regs_delonereg_callback); row -= rh

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
    Draw.PushButton   ('Add',        EVT_MESH_ADDHOLE,     ggx+120, row+2, 60, rh-4, 'Add hole')
    Draw.PushButton   ('Delete all', EVT_MESH_DELALLHOLES, ggx+200, row+2, 80, rh-4, 'Delete all holes'); row -= rh
    BGL.glRasterPos2i (ggx, row+5); Draw.Text('         X                  Y                  Z'); row -= rh
    for i, hol in enumerate(hols):
        Draw.String     ('',    2000+i, ggx    , row, 100, rh, hol[0],128,'X of the hole',holes_setx_callback)
        Draw.String     ('',    2000+i, ggx+100, row, 100, rh, hol[1],128,'Y of the hole',holes_sety_callback)
        Draw.String     ('',    2000+i, ggx+200, row, 100, rh, hol[2],128,'Z of the hole',holes_setz_callback)
        Draw.PushButton ('Del', 2000+i, ggx+300, row,  40, rh, 'Delete this row', holes_delonehole_callback); row -= rh

    # Mesh -- unstructured -- main
    row -= 2*sgy
    ggx -= gx
    Draw.PushButton ('Generate (triangles/tetrahedrons)', EVT_MESH_GENUNSTRU ,  ggx, row, 300, rh, 'Generated unstructured mesh')

    # FEA
    dx   = 2*gx+55;
    ggx  = 2*gx
    row -= ggy+gy
    h    = h_fea
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('FEA')
    BGL.glColor3f     (0.8, 0.8, 0.8)
    BGL.glRecti       (gx, row, wid-gx, row-h)

    # FEA -- nodes bry
    ggx += gx
    row -= gy
    h    = h_fea_nbrys
    BGL.glColor3f     (0.431, 0.443, 0.514)
    BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('Nodes boundaries')
    BGL.glColor3f     (0.72, 0.72, 0.8)
    BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
    BGL.glColor3f     (0.0, 0.0, 0.0)
    Draw.PushButton   ('Add',        EVT_FEA_ADDNBRY,    ggx+120, row+2, 60, rh-4, 'Add nodes boundary')
    Draw.PushButton   ('Delete all', EVT_FEA_DELALLNBRY, ggx+200, row+2, 80, rh-4, 'Delete all nodes boundary'); row -= rh
    BGL.glRasterPos2i (ggx, row+5); Draw.Text('         X                  Y                  Z           Key        Value'); row -= rh
    for i, nb in enumerate(nbrys):
        Draw.String     ('',    3000+i, ggx    , row, 80, rh, nb[0],128,'X of the node with boundary condition',                                      nodebry_setx_callback)
        Draw.String     ('',    3000+i, ggx+ 80, row, 80, rh, nb[1],128,'Y of the node with boundary condition',                                      nodebry_sety_callback)
        Draw.String     ('',    3000+i, ggx+160, row, 80, rh, nb[2],128,'Z of the node with boundary condition',                                      nodebry_setz_callback)
        Draw.String     ('',    3000+i, ggx+240, row, 40, rh, nb[3],  2,'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', nodebry_setkey_callback)
        Draw.String     ('',    3000+i, ggx+280, row, 80, rh, nb[4],128,'Value of essential/natural boundary condition',                              nodebry_setval_callback)
        Draw.PushButton ('Del', 3000+i, ggx+360, row, 40, rh, 'Delete this row',                                                                      nodebry_delonenbry_callback); row -= rh

    # FEA -- edges bry
    h = h_fea_ebrys
    BGL.glColor3f     (0.431, 0.443, 0.514)
    BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('Edges boundaries')
    BGL.glColor3f     (0.72, 0.72, 0.8)
    BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
    BGL.glColor3f     (0.0, 0.0, 0.0)
    Draw.PushButton   ('Add',        EVT_FEA_ADDEBRY,    ggx+120, row+2, 60, rh-4, 'Add edges boundary')
    Draw.PushButton   ('Delete all', EVT_FEA_DELALLEBRY, ggx+200, row+2, 80, rh-4, 'Delete all edges boundary'); row -= rh
    BGL.glRasterPos2i (ggx, row+5); Draw.Text('        Tag         Key        Value'); row -= rh
    for i, eb in enumerate(ebrys):
        Draw.Number     ('',    4000+i, ggx    , row, 80, rh, int(eb[0]),-1000,0,'Set the edge tag',                                                          edgebry_settag_callback)
        Draw.String     ('',    4000+i, ggx+ 80, row, 40, rh, eb[1],  2,         'Key such as ux, uy, fx, fz corresponding to the essential/natural variable',edgebry_setkey_callback)
        Draw.String     ('',    4000+i, ggx+120, row, 80, rh, eb[2],128,         'Value of essential/natural boundary condition',                             edgebry_setval_callback)
        Draw.PushButton ('Del', 4000+i, ggx+200, row, 40, rh, 'Delete this row',                                                                              edgebry_deloneebry_callback); row -= rh

    # FEA -- faces bry
    h = h_fea_fbrys
    BGL.glColor3f     (0.431, 0.443, 0.514)
    BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('Faces boundaries')
    BGL.glColor3f     (0.72, 0.72, 0.8)
    BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
    BGL.glColor3f     (0.0, 0.0, 0.0)
    Draw.PushButton   ('Add',        EVT_FEA_ADDFBRY,    ggx+120, row+2, 60, rh-4, 'Add faces boundary')
    Draw.PushButton   ('Delete all', EVT_FEA_DELALLFBRY, ggx+200, row+2, 80, rh-4, 'Delete all faces boundary'); row -= rh
    BGL.glRasterPos2i (ggx, row+5); Draw.Text('        Tag         Key        Value'); row -= rh
    for i, fb in enumerate(fbrys):
        Draw.Number     ('',    5000+i, ggx    , row, 80, rh, int(fb[0]),-1000,0,'Set the face tag',                                                           facebry_settag_callback)
        Draw.String     ('',    5000+i, ggx+ 80, row, 40, rh, fb[1],  2,         'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', facebry_setkey_callback)
        Draw.String     ('',    5000+i, ggx+120, row, 80, rh, fb[2],128,         'Value of essential/natural boundary condition',                              facebry_setval_callback)
        Draw.PushButton ('Del', 5000+i, ggx+200, row, 40, rh, 'Delete this row',                                                                               facebry_delonefbry_callback); row -= rh

    # FEA -- elems bry
    h = h_fea_eatts
    BGL.glColor3f     (0.431, 0.443, 0.514)
    BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('Elements attributes')
    BGL.glColor3f     (0.72, 0.72, 0.8)
    BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
    BGL.glColor3f     (0.0, 0.0, 0.0)
    Draw.PushButton   ('Add',        EVT_FEA_ADDEATT,    ggx+120, row+2, 60, rh-4, 'Add elems attributes')
    Draw.PushButton   ('Delete all', EVT_FEA_DELALLEATT, ggx+200, row+2, 80, rh-4, 'Delete all elems attributes'); row -= rh
    BGL.glRasterPos2i (ggx, row+5); Draw.Text('      Tag           Type                Model           Parameters        Initial Vals'); row -= rh
    for i, ea in enumerate(eatts):
        Draw.Number     ('',    6000+i, ggx    , row,  60, rh, int(ea[0]),-1000,0,'Set the element tag',                       elematt_settag_callback)
        Draw.String     ('',    6000+i, ggx+ 60, row, 100, rh, ea[1], 20,         'Element type: ex.: Quad4PStrain',           elematt_settype_callback)
        Draw.String     ('',    6000+i, ggx+160, row,  80, rh, ea[2], 20,         'Constitutive model: ex.: LinElastic',       elematt_setmodel_callback)
        Draw.String     ('',    6000+i, ggx+240, row, 100, rh, ea[3],128,         'Parameters: ex.: E=207_nu=0.3',             elematt_setprms_callback)
        Draw.String     ('',    6000+i, ggx+340, row,  80, rh, ea[4],128,         'Initial values: ex.: Sx=0 Sy=0 Sz=0 Sxy=0', elematt_setinis_callback)
        Draw.PushButton ('Del', 6000+i, ggx+420, row,  40, rh, 'Delete this row',                                              elematt_deloneeatt_callback); row -= rh

    # FEA -- main
    row -= 2*sgy
    ggx -= gx
    Draw.PushButton ('Run analysis',    EVT_FEA_RUN,    ggx,     row, 150, rh, 'Run a FE analysis directly (without script)')
    Draw.PushButton ('Generate script', EVT_FEA_SCRIPT, ggx+150, row, 150, rh, 'Generate script for FEA'); row -= rh


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
