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

# Constants
EVT_NONE    = 0 # used for buttons with callbacks
EVT_REFRESH = 1 # refresh all windows
# CAD
EVT_CAD_ADDXYZ     =  2 # type in values to define point(s)
EVT_CAD_FILLET     =  3 # 3D fillet
EVT_CAD_EDGEBREAK  =  4 # break edge
EVT_CAD_EDGEBREAKM =  5 # break edge at middle point
EVT_CAD_FPOINT     =  6 # read points from file
EVT_CAD_FSPLINE    =  7 # create a spline from points in a file
# Mesh
EVT_MESH_SETETAG =  8 # set edges tag
EVT_MESH_SETFTAG =  9 # set faces tag
EVT_MESH_DELPROP = 10 # delete all properties
# Mesh -- structured
EVT_MESH_SETX     = 11 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_SETY     = 12 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_SETZ     = 13 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_GENSTRU  = 14 # generate structured mesh using MechSys module
# FEA
EVT_FEA_ADDNBRY    = 15 # add nodes boundary
EVT_FEA_DELALLNBRY = 16 # delete all nodes boundary
EVT_FEA_ADDFBRY    = 17 # add faces boundary
EVT_FEA_DELALLFBRY = 18 # delete all faces boundary
EVT_FEA_ADDEATT    = 19 # add element attributes
EVT_FEA_DELALLEATT = 20 # delete all elements attributes
EVT_FEA_RUN        = 21 # run a FE simulation
EVT_FEA_SCRIPT     = 22 # generate script for FEA 


# Get current active mesh (will exit edit mode)
# Returns:
#          edm: edit mode == True or False
#          obj: current object
#          msh: current mesh
def get_mesh(with_error=True):
    edm = Blender.Window.EditMode()
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        return edm, obj, msh
    else:
        if with_error: Blender.Draw.PupMenu('ERROR|Please, select an object of type Mesh')
        else: return 0,0,0


# Handle input events
def event(evt, val):
    if evt==Draw.QKEY and not val: Draw.Exit()
    if evt==Draw.ESCKEY: Draw.Redraw(1)


# Handle button events
def button_event(evt):

    # load dictionary
    dict = di.load_dict()

    # refresh all windows
    if evt==EVT_REFRESH: Blender.Window.QRedrawAll()

    # =============================================================================== CAD

    # add vertices
    if evt==EVT_CAD_ADDXYZ:
        x = dict['newpoint_x']
        y = dict['newpoint_y']
        z = dict['newpoint_z']
        dr.add_point (Vector(x,y,z))

    # fillet two edges
    elif evt==EVT_CAD_FILLET: dr.fillet(dict['fillet_radius'],dict['fillet_steps'])

    # break an edge
    elif evt==EVT_CAD_EDGEBREAK: dr.break_edge()

    # break an edge at middle point
    elif evt==EVT_CAD_EDGEBREAKM: dr.break_edge(True)

    # read vertices from file
    elif evt==EVT_CAD_FPOINT: Blender.Window.FileSelector(dr.add_points_from_file, "Read X Y Z cols")

    # read spline from file
    elif evt==EVT_CAD_FSPLINE: Blender.Window.FileSelector(dr.add_spline_from_file, "Read X Y Z cols")

    # =============================================================================== Mesh

    # set edges tag
    elif evt==EVT_MESH_SETETAG:
        edm, obj, msh = get_mesh()
        tag = Blender.Draw.Create(-10)
        block = []
        block.append(("Tag = ", tag, -10000, 0))
        if Blender.Draw.PupBlock("Set edges tag", block)==1:
            for id in msh.edges.selected():
                di.set_tag (obj, 'edge', id, tag.val)
            Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set faces tag
    elif evt==EVT_MESH_SETFTAG:
        edm, obj, msh = get_mesh()
        tag = Blender.Draw.Create(-100)
        block = []
        block.append(("Tag = ", tag, -10000, 0))
        if Blender.Draw.PupBlock("Set faces tag", block)==1:
            for id in msh.faces.selected():
                di.set_tag (obj, 'face', id, tag.val)
            Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # delete all properties
    elif evt==EVT_MESH_DELPROP:
        scn = bpy.data.scenes.active
        obs = scn.objects.selected
        for o in obs:
            for p in o.getAllProperties():
                print o.name, p.name, p.data
                o.removeProperty(p)
            Blender.Window.QRedrawAll()

    # ------------------------------------------------------------------------------- Structured

    # set local x-axis
    elif evt==EVT_MESH_SETX:
        edm, obj, msh = get_mesh()
        if len(msh.edges.selected())==1:
            di.set_local_axis (obj, 'x', msh.edges.selected()[0])
            Blender.Window.QRedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select only one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set local y-axis
    elif evt==EVT_MESH_SETY:
        edm, obj, msh = get_mesh()
        if len(msh.edges.selected())==1:
            di.set_local_axis (obj, 'y', msh.edges.selected()[0])
            Blender.Window.QRedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select only one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set local z-axis
    elif evt==EVT_MESH_SETZ:
        edm, obj, msh = get_mesh()
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

    # =============================================================================== FEA

    # add nodes boundary
    elif evt==EVT_FEA_ADDNBRY:
        scn = bpy.data.scenes.active
        obj = scn.objects.active
        if obj!=None:
            nbrys = di.get_nbrys (obj)
            di.set_nbry (obj, len(nbrys), 0,0,0, 'uy', 0.0)
            Blender.Window.QRedrawAll()

    # delete all nodes boundary
    elif evt==EVT_FEA_DELALLNBRY:
        scn = bpy.data.scenes.active
        obj = scn.objects.active
        if obj!=None:
            di.del_all_nbrys(obj)
            Blender.Window.QRedrawAll()

    # add faces boundary
    elif evt==EVT_FEA_ADDFBRY:
        scn = bpy.data.scenes.active
        obj = scn.objects.active
        if obj!=None:
            fbrys = di.get_fbrys (obj)
            di.set_fbry (obj, len(fbrys), -10, 'uy', 0.0)
            Blender.Window.QRedrawAll()

    # delete all faces boundary
    elif evt==EVT_FEA_DELALLFBRY:
        scn = bpy.data.scenes.active
        obj = scn.objects.active
        if obj!=None:
            di.del_all_fbrys(obj)
            Blender.Window.QRedrawAll()

    # add elems attributes
    elif evt==EVT_FEA_ADDEATT:
        scn = bpy.data.scenes.active
        obj = scn.objects.active
        if obj!=None:
            eatts = di.get_eatts (obj)
            di.set_eatt (obj, len(eatts), -10, 'Quad4PStrain', 'LinElastic', 'E=207 nu=0.3', 'Sx=0 Sy=0 Sz=0 Sxy=0')
            Blender.Window.QRedrawAll()

    # delete all elems attributes
    elif evt==EVT_FEA_DELALLEATT:
        scn = bpy.data.scenes.active
        obj = scn.objects.active
        if obj!=None:
            di.del_all_eatts(obj)
            Blender.Window.QRedrawAll()

    # run a FE simulation
    elif evt==EVT_FEA_RUN:
        edm, obj, msh = get_mesh()
        nbrys = di.get_nbrys (obj)
        fbrys = di.get_fbrys (obj)
        eatts = di.get_eatts (obj)
        fem.run_fea (obj, nbrys, fbrys, eatts)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    elif evt==EVT_FEA_SCRIPT:
        edm, obj, msh = get_mesh()
        nbrys = di.get_nbrys (obj)
        fbrys = di.get_fbrys (obj)
        eatts = di.get_eatts (obj)
        fem.gen_script (obj, nbrys, fbrys, eatts)
        Blender.Window.Redraw(Blender.Window.Types.TEXT)


# CAD
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

# Mesh
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

def show_ele_tags_callback(evt,val):
    dict = di.load_dict()
    dict['show_ele_tags'] = val
    Blender.Window.QRedrawAll()

def ndivx_callback(evt,val):
    edm, obj, msh = get_mesh()
    di.set_ndiv (obj, 'x', val)
    if edm: Blender.Window.EditMode(1)
    Blender.Window.QRedrawAll()

def ndivy_callback(evt,val):
    edm, obj, msh = get_mesh()
    di.set_ndiv (obj, 'y', val)
    if edm: Blender.Window.EditMode(1)
    Blender.Window.QRedrawAll()

def ndivz_callback(evt,val):
    edm, obj, msh = get_mesh()
    di.set_ndiv (obj, 'z', val)
    if edm: Blender.Window.EditMode(1)
    Blender.Window.QRedrawAll()

def btag_callback(evt,val):
    edm, obj, msh = get_mesh()
    di.set_btag (obj, val)
    if edm: Blender.Window.EditMode(1)
    Blender.Window.QRedrawAll()


# FEA ------------------------------- Nodes boundaries
def nodebry_setx_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_nbry_x (obj, evt-1000, val)

def nodebry_sety_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_nbry_y (obj, evt-1000, val)

def nodebry_setz_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_nbry_z (obj, evt-1000, val)

def nodebry_setkey_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_nbry_key (obj, evt-1000, val)

def nodebry_setval_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_nbry_val (obj, evt-1000, val)

def nodebry_delonenbry_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.del_nbry (obj,evt-1000)
        Blender.Window.QRedrawAll()

# FEA ------------------------------- Faces boundaries
def facebry_settag_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_fbry_tag (obj, evt-10000, val)

def facebry_setkey_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_fbry_key (obj, evt-10000, val)

def facebry_setval_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_fbry_val (obj, evt-10000, val)

def facebry_delonefbry_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.del_fbry (obj,evt-10000)
        Blender.Window.QRedrawAll()

# FEA ------------------------------- Elements attributes
def elematt_settag_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_eatt_tag (obj, evt-10000, val)

def elematt_settype_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_eatt_type (obj, evt-10000, val)

def elematt_setmodel_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_eatt_model (obj, evt-10000, val)

def elematt_setprms_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_eatt_prms (obj, evt-10000, val)

def elematt_setinis_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.set_eatt_inis (obj, evt-10000, val)

def elematt_deloneeatt_callback(evt,val):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None:
        di.del_eatt (obj,evt-10000)
        Blender.Window.QRedrawAll()


# Draw GUI
def gui():
    # load dictionary
    d = di.load_dict()

    # Data from current selected object
    edm, obj, msh = get_mesh (False)
    if obj>0:
        btag   = di.get_btag  (obj)
        ndivx  = di.get_ndiv  (obj,'x')
        ndivy  = di.get_ndiv  (obj,'y')
        ndivz  = di.get_ndiv  (obj,'z')
        nbrys  = di.get_nbrys (obj)
        fbrys  = di.get_fbrys (obj)
        eatts  = di.get_eatts (obj)
        if edm: Blender.Window.EditMode(1)
    else:
        btag   = -1
        ndivx  = 2
        ndivy  = 2
        ndivz  = 1
        nbrys  = []
        fbrys  = []
        eatts  = []

    # Width and Height
    wid, hei = Blender.Window.GetAreaSize()

    # Col and Row to add entities
    gx  = 10
    gy  = 10
    ggx = 2*gx
    ggy = 2*gy
    col = 7
    row = 0
    cw  = 100  # column width
    rh  = 20   # row height
    ch  = rh+2 # caption height

    # Background color
    BGL.glClearColor (0.431, 0.443, 0.514, 0.0)
    BGL.glClear      (Blender.BGL.GL_COLOR_BUFFER_BIT)
    
    # CAD
    dx  = 2*gx+50;
    row = hei-gy
    h   = 4*rh+ggy
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('CAD')
    BGL.glColor3f     (0.8, 0.8, 0.8)
    BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Points:')
    Draw.Number       ('X=',        EVT_NONE      , dx    , row, 80, rh, d['newpoint_x'],-1.0e+7,1.0e+7,'Set the X value of the new point to be added',setx_callback)
    Draw.Number       ('Y=',        EVT_NONE      , dx+ 80, row, 80, rh, d['newpoint_y'],-1.0e+7,1.0e+7,'Set the Y value of the new point to be added',sety_callback)
    Draw.Number       ('Z=',        EVT_NONE      , dx+160, row, 80, rh, d['newpoint_z'],-1.0e+7,1.0e+7,'Set the Z value of the new point to be added',setz_callback)
    Draw.PushButton   ('Add x y z', EVT_CAD_ADDXYZ, dx+240, row, 80, rh, 'Add point from coordinates'); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Edges:');
    Draw.Number       ('Radius', EVT_NONE      , dx     , row , 120 , rh , d['fillet_radius'],0.0,1.0e+9,'Set radius for fillet operation',fillet_radius_callback)
    Draw.Number       ('Steps' , EVT_NONE      , dx+120 , row , 120 , rh , d['fillet_steps' ],1  ,100   ,'Set steps for fillet operation' ,fillet_steps_callback)
    Draw.PushButton   ('Fillet', EVT_CAD_FILLET, dx+240 , row ,  80 , rh , 'Create a fillet between two edges'); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Edges:');
    Draw.PushButton   ('Break edge',      EVT_CAD_EDGEBREAK , dx,      row ,  120, rh , 'Break an edge at a previously selected point')
    Draw.PushButton   ('Break at middle', EVT_CAD_EDGEBREAKM, dx+120 , row ,  120, rh , 'Break an edge at its middle point'); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('File:')
    Draw.PushButton   ('Read Points', EVT_CAD_FPOINT,  dx,     row, 120, rh, 'Add points by reading a list from file')
    Draw.PushButton   ('Read Spline', EVT_CAD_FSPLINE, dx+120, row, 120, rh, 'Add a spline by reading its points from file')

    # Mesh
    dx   = 2*gx+55;
    row -= ggy
    h    = 8*rh+ggy+5*gy
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('Mesh')
    BGL.glColor3f     (0.8, 0.8, 0.8)
    BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Show:')
    Draw.Toggle       ('Props', EVT_NONE, dx    , row, 60, rh, d['show_props'],    'Show mesh properties'  , show_props_callback)
    Draw.Toggle       ('Axes' , EVT_NONE, dx+ 60, row, 60, rh, d['show_axes'],     'Show local system axes', show_axes_callback )
    Draw.Toggle       ('V IDs', EVT_NONE, dx+120, row, 60, rh, d['show_v_ids'],    'Show vertex IDs'       , show_v_ids_callback)
    Draw.Toggle       ('E IDs', EVT_NONE, dx+180, row, 60, rh, d['show_e_ids'],    'Show edge IDs'         , show_e_ids_callback)
    Draw.Toggle       ('F IDs', EVT_NONE, dx+240, row, 60, rh, d['show_f_ids'],    'Show face IDs'         , show_f_ids_callback)
    Draw.Toggle       ('Ele T', EVT_NONE, dx+300, row, 60, rh, d['show_ele_tags'], 'Show elements tags'    , show_ele_tags_callback); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Set tags:')
    Draw.PushButton   ('Edge', EVT_MESH_SETETAG,  dx,    row, 80, rh, 'Set edges tag')
    Draw.PushButton   ('Face', EVT_MESH_SETFTAG,  dx+80, row, 80, rh, 'Set faces tag'); row -= rh+gy
    Draw.PushButton   ('Delete all properties', EVT_MESH_DELPROP, ggx, row, 200, rh , 'Delete all properties')
    
    # Mesh -- structured
    dx   = 2*gx+70;
    ggx += gx
    dx  += gx
    row -= gy
    h    = 3*rh+ggy+3*gy
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
    Draw.Number       ('nX=',  EVT_NONE, dx,     row, 80, rh, ndivx, 1, 10000,'Set number of divisions along x (local axis)',ndivx_callback)
    Draw.Number       ('nY=',  EVT_NONE, dx+ 80, row, 80, rh, ndivy, 1, 10000,'Set number of divisions along y (local axis)',ndivy_callback)
    Draw.Number       ('nZ=',  EVT_NONE, dx+160, row, 80, rh, ndivz, 1, 10000,'Set number of divisions along z (local axis)',ndivz_callback); row -= rh
    BGL.glRasterPos2i (ggx, row+4)
    Draw.Text         ('Block:')
    Draw.Number       ('Tag=',     EVT_NONE,          dx,  row, 100, rh, btag, -1000, 0,'Set the block tag (to be replicated to all elements)',btag_callback); row -= rh+gy
    Draw.PushButton   ('Generate', EVT_MESH_GENSTRU,  ggx, row, 200, rh, 'Generated structured mesh')

    # FEA
    dx   = 2*gx+55;
    ggx  = 2*gx
    row -= ggy+gy
    h    = (6+len(nbrys)+len(fbrys)+len(eatts))*rh+rh+8*gy
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
    h    = (1+len(nbrys))*rh+gy
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
    BGL.glRasterPos2i (ggx+28,    row+5); Draw.Text('X')
    BGL.glRasterPos2i (ggx+28+60, row+5); Draw.Text('Y')
    BGL.glRasterPos2i (ggx+28+120, row+5); Draw.Text('Z')
    BGL.glRasterPos2i (ggx+10+180, row+5); Draw.Text('Key')
    BGL.glRasterPos2i (ggx+23+220, row+5); Draw.Text('Value'); row -= rh
    for i, nb in enumerate(nbrys):
        Draw.Number     ('',    1000+i, ggx    , row, 60, rh, nb[0],-1.0e+7,1.0e+7,'X of the node with boundary condition',nodebry_setx_callback)
        Draw.Number     ('',    1000+i, ggx+ 60, row, 60, rh, nb[1],-1.0e+7,1.0e+7,'Y of the node with boundary condition',nodebry_sety_callback)
        Draw.Number     ('',    1000+i, ggx+120, row, 60, rh, nb[2],-1.0e+7,1.0e+7,'Z of the node with boundary condition',nodebry_setz_callback)
        Draw.String     ('',    1000+i, ggx+180, row, 40, rh, nb[3],2,'Key such as ux, uy, fx, fz corresponding to the essential/natural variable',nodebry_setkey_callback)
        Draw.Number     ('',    1000+i, ggx+220, row, 80, rh, nb[4],-1.0e+7,1.0e+7,'Value of essential/natural boundary condition',nodebry_setval_callback)
        Draw.PushButton ('Del', 1000+i, ggx+300, row, 40, rh, 'Delete this row', nodebry_delonenbry_callback); row -= rh

    # FEA -- faces bry
    h = (1+len(fbrys))*rh+gy
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
    BGL.glRasterPos2i (ggx+28,     row+5); Draw.Text('Tag')
    BGL.glRasterPos2i (ggx+10+60,  row+5); Draw.Text('Key')
    BGL.glRasterPos2i (ggx+23+100, row+5); Draw.Text('Value'); row -= rh
    for i, fb in enumerate(fbrys):
        Draw.Number     ('',    10000+i, ggx    , row, 60, rh, fb[0],-10000,0,'Set the face tag',facebry_settag_callback)
        Draw.String     ('',    10000+i, ggx+ 60, row, 40, rh, fb[1],2,'Key such as ux, uy, fx, fz corresponding to the essential/natural variable',facebry_setkey_callback)
        Draw.Number     ('',    10000+i, ggx+100, row, 80, rh, fb[2],-1.0e+7,1.0e+7,'Value of essential/natural boundary condition',facebry_setval_callback)
        Draw.PushButton ('Del', 10000+i, ggx+180, row, 40, rh, 'Delete this row', facebry_delonefbry_callback); row -= rh

    # FEA -- elems bry
    h = (1+len(eatts))*rh+gy
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
    BGL.glRasterPos2i (ggx+28,     row+5); Draw.Text('Tag')
    BGL.glRasterPos2i (ggx+10+60,  row+5); Draw.Text('Key')
    BGL.glRasterPos2i (ggx+23+100, row+5); Draw.Text('Value'); row -= rh
    for i, ea in enumerate(eatts):
        Draw.Number     ('',    10000+i, ggx    , row, 60, rh, ea[0],-10000,0,'Set the element tag',                elematt_settag_callback)
        Draw.String     ('',    10000+i, ggx+ 60, row, 60, rh, ea[1],20,'Element type: ex.: Quad4PStrain',          elematt_settype_callback)
        Draw.String     ('',    10000+i, ggx+120, row, 60, rh, ea[2],20,'Constitutive model: ex.: LinElastic',      elematt_setmodel_callback)
        Draw.String     ('',    10000+i, ggx+180, row, 60, rh, ea[3],64,'Parameters: ex.: E=207_nu=0.3',            elematt_setprms_callback)
        Draw.String     ('',    10000+i, ggx+240, row, 60, rh, ea[4],64,'Initial values: ex.: Sx=0 Sy=0 Sz=0 Sxy=0',elematt_setinis_callback)
        Draw.PushButton ('Del', 10000+i, ggx+300, row, 40, rh, 'Delete this row', elematt_deloneeatt_callback); row -= rh

    # FEA -- main
    ggx -= gx
    row -= ggy
    Draw.PushButton ('Run analysis',    EVT_FEA_RUN,    ggx,     row, 120, rh, 'Run a FE analysis directly (without script)')
    Draw.PushButton ('Generate script', EVT_FEA_SCRIPT, ggx+120, row, 120, rh, 'Generate script for FEA'); row -= rh


# Register GUI
Draw.Register (gui, event, button_event)

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
