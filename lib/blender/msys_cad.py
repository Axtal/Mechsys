#!BPY

"""
Name: 'MechSys'
Blender: 2.45
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

# Constants
EVT_NONE       =  0 # used for buttons with callbacks
EVT_VPOINT     =  1 # type in values to define point(s)
EVT_FPOINT     =  2 # read points from file
EVT_FSPLINE    =  3 # create a spline from points in a file
EVT_FILLET     =  4 # 3D fillet
EVT_BREAK      =  5 # break edge
EVT_READ_TRI   =  6 # read 2D mesh from Triangle
EVT_GEN_STRU   =  7 # generate structured mesh using MechSys module
EVT_SET_X      =  8 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_SET_Y      =  9 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_SET_Z      = 10 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_SET_NDIV   = 11 # set number of divisons along x, y and z
EVT_ETAG       = 12 # set edges tag
EVT_FTAG       = 13 # set faces tag
EVT_DEL_PROP   = 14 # delete all properties

# Get current active mesh (will exit edit mode)
# Returns:
#          edm: edit mode == True or False
#          obj: current object
#          msh: current mesh
def get_mesh():
    edm = Blender.Window.EditMode()
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        return edm, obj, msh
    else: Blender.Draw.PupMenu('ERROR|Please, select an object of type Mesh')

# Handle input events
def event(evt, val):
    if evt==Draw.QKEY and not val: Draw.Exit()
    if evt==Draw.ESCKEY: Draw.Redraw(1)

# Handle button events
def button_event(evt):

    # load dictionary
    dict = di.load_dict()

    # add vertices
    if evt==EVT_VPOINT:
        x    = Blender.Draw.Create(0.000)
        y    = Blender.Draw.Create(0.000)
        z    = Blender.Draw.Create(0.000)
        more = Blender.Draw.Create(1)
        block = []
        block.append("Coordinates")
        block.append(("X = ", x, -1.0e+7, 1.0e+7))
        block.append(("Y = ", y, -1.0e+7, 1.0e+7))
        block.append(("Z = ", z, -1.0e+7, 1.0e+7))
        block.append(("Add More", more, "Add more points?"))
        while more==1 and Blender.Draw.PupBlock("Type in the Coordinates", block)==1:
            dr.add_point(Vector(x.val,y.val,z.val))

    # read vertices from file
    elif evt==EVT_FPOINT: Blender.Window.FileSelector(dr.add_points_from_file, "Read X Y Z cols")

    # read spline from file
    elif evt==EVT_FSPLINE: Blender.Window.FileSelector(dr.add_spline_from_file, "Read X Y Z cols")

    # fillet two edges
    elif evt==EVT_FILLET: dr.fillet(dict['fillet_radius'],dict['fillet_steps'])

    # break an edge
    elif evt==EVT_BREAK: dr.break_edge()

    # read 2d mesh from 'Triangle'
    elif evt==EVT_READ_TRI: Blender.Window.FileSelector(dr.read_2d_mesh, "Read 2D mesh")

    # set local x-axis
    elif evt==EVT_SET_X:
        edm, obj, msh = get_mesh()
        if len(msh.edges.selected())==1:
            di.set_local_axis (obj, 'x', msh.edges.selected()[0])
            Blender.Window.QRedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select only one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set local y-axis
    elif evt==EVT_SET_Y:
        edm, obj, msh = get_mesh()
        if len(msh.edges.selected())==1:
            di.set_local_axis (obj, 'y', msh.edges.selected()[0])
            Blender.Window.QRedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select only one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set local z-axis
    elif evt==EVT_SET_Z:
        edm, obj, msh = get_mesh()
        if len(msh.edges.selected())==1:
            di.set_local_axis (obj, 'z', msh.edges.selected()[0])
            Blender.Window.QRedrawAll()
        else: Blender.Draw.PupMenu('ERROR|Please, select onlz one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # define number of divisions
    elif evt==EVT_SET_NDIV:
        scn   = bpy.data.scenes.active
        obj   = scn.objects.active
        ndivx = Blender.Draw.Create(di.get_ndiv(obj,'x'))
        ndivy = Blender.Draw.Create(di.get_ndiv(obj,'y'))
        ndivz = Blender.Draw.Create(di.get_ndiv(obj,'z'))
        block = []
        block.append(("nDivX = ", ndivx, 1, 100000))
        block.append(("nDivY = ", ndivy, 1, 100000))
        block.append(("nDivZ = ", ndivz, 1, 100000))
        if Blender.Draw.PupBlock("Set number of divisions", block)==1:
            scn = bpy.data.scenes.active
            obs = scn.objects.selected
            for o in obs:
                di.set_ndiv (o, 'x', ndivx.val)
                di.set_ndiv (o, 'y', ndivy.val)
                di.set_ndiv (o, 'z', ndivz.val)
            Blender.Window.QRedrawAll()

    # set edges tag
    elif evt==EVT_ETAG:
        edm, obj, msh = get_mesh()
        tag = Blender.Draw.Create(-10)
        block = []
        block.append(("Tag = ", tag, -10000, 0))
        if Blender.Draw.PupBlock("Set edges tag", block)==1:
            for id in msh.edges.selected():
                di.set_tags (obj, 'edge', id, tag.val)
            Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # set faces tag
    elif evt==EVT_FTAG:
        edm, obj, msh = get_mesh()
        tag = Blender.Draw.Create(-100)
        block = []
        block.append(("Tag = ", tag, -10000, 0))
        if Blender.Draw.PupBlock("Set faces tag", block)==1:
            for id in msh.faces.selected():
                di.set_tags (obj, 'face', id, tag.val)
            Blender.Window.QRedrawAll()
        if edm: Blender.Window.EditMode(1) # return to EditMode

    # generate structured mesh via MechSys
    elif evt==EVT_GEN_STRU:
        tt = timeit.Timer('msys_draw.gen_struct_mesh()', 'import msys_draw')
        t  = tt.timeit(number=1)
        print '[1;34mMechSys[0m: time spent on generation and drawing = [1;31m',t,'[0m [1;32mseconds[0m'

    # delete all properties
    elif evt==EVT_DEL_PROP:
        scn = bpy.data.scenes.active
        obs = scn.objects.selected
        for o in obs:
            for p in o.getAllProperties():
                o.removeProperty(p)
            Blender.Window.QRedrawAll()


def fillet_radius_callback(evt,val):
    dict = di.load_dict()
    dict['fillet_radius'] = val

def fillet_steps_callback(evt,val):
    dict = di.load_dict()
    dict['fillet_steps'] = val

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

def show_ndivs_callback(evt,val):
    dict = di.load_dict()
    dict['show_ndivs'] = val
    Blender.Window.QRedrawAll()


# Draw GUI
def gui():
    # load dictionary
    d = di.load_dict()

    # Width and Height
    wid, hei = Blender.Window.GetAreaSize()

    # Col and Row to add entities
    col = 7
    row = 0
    cw  = 100  # column width
    rh  = 20   # row height

    # Background color
    Blender.BGL.glClearColor (0.431, 0.443, 0.514, 0.0)
    Blender.BGL.glClear      (Blender.BGL.GL_COLOR_BUFFER_BIT)

    # Buttons
    Draw.PushButton("Point(s) x y z"    , EVT_VPOINT   ,   0 , row , 100 , rh , "Type in values to define point(s)")
    Draw.PushButton("Point(s) [file]"   , EVT_FPOINT   , 100 , row , 100 , rh , "Read point list from file")          ; row += rh
    Draw.PushButton("Spline [file]"     , EVT_FSPLINE  ,   0 , row , 100 , rh , "Add a spline with points from file") ; row += rh
    Draw.PushButton("Fillet edge"       , EVT_FILLET   ,   0 , row , 100 , rh , "Fillet two edges")
    Draw.Number    ("Radius"            , EVT_NONE     , 100 , row , 100 , rh , d['fillet_radius'],0.0,1.0e+9,"Set radius for fillet operation",fillet_radius_callback)
    Draw.Number    ("Steps"             , EVT_NONE     , 200 , row , 100 , rh , d['fillet_steps' ],1  ,100   ,"Set steps for fillet operation" ,fillet_steps_callback) ; row += rh
    Draw.PushButton("Break edge"        , EVT_BREAK    ,   0 , row , 100 , rh , "Break an edge")                      ; row += rh
    Draw.PushButton("Read from triangle", EVT_READ_TRI ,   0 , row , 120 , rh , "Read 2D mesh from Triangle")

    # Mesh properties
    BGL.glColor3f     (0.4,0.4,0.4)
    BGL.glRecti       (0, hei, wid, hei/3)
    BGL.glColor3f     (1,1,1)
    BGL.glRasterPos2i (3,hei-14)
    Draw.Text         ('Mesh properties:'); row = hei-rh
    Draw.Toggle       ('Show props'     , EVT_NONE      , 120 , row , 100 , rh, d['show_props'], 'Show mesh properties'    , show_props_callback)
    Draw.Toggle       ('Show axes'      , EVT_NONE      , 220 , row , 100 , rh, d['show_axes'],  'Show local system axes'  , show_axes_callback ) ; row -= rh
    Draw.Toggle       ('Show V IDs'     , EVT_NONE      ,   0 , row ,  80 , rh, d['show_v_ids'], 'Show vertex IDs'         , show_v_ids_callback) ;
    Draw.Toggle       ('Show E IDs'     , EVT_NONE      ,  80 , row ,  80 , rh, d['show_e_ids'], 'Show edge IDs'           , show_e_ids_callback) ;
    Draw.Toggle       ('Show F IDs'     , EVT_NONE      , 160 , row ,  80 , rh, d['show_f_ids'], 'Show face IDs'           , show_f_ids_callback) ;
    Draw.Toggle       ('nDivs'          , EVT_NONE      , 240 , row ,  80 , rh, d['show_ndivs'], 'Show number of divisions', show_ndivs_callback) ; row -= rh
    Draw.PushButton   ('Set local x'    , EVT_SET_X     ,   0 , row ,  80 , rh , 'Set local x axis (needs one edge selected)')
    Draw.PushButton   ('Set local y'    , EVT_SET_Y     ,  80 , row ,  80 , rh , 'Set local y axis (needs one edge selected)')
    Draw.PushButton   ('Set local z'    , EVT_SET_Z     , 160 , row ,  80 , rh , 'Set local z axis (needs one edge selected)')
    Draw.PushButton   ('Set nDivs'      , EVT_SET_NDIV  , 240 , row ,  80 , rh , 'Set set number of divisions along x, y, and z') ; row -= rh
    Draw.PushButton   ('Set E tag'      , EVT_ETAG      ,   0 , row ,  80 , rh , 'Set edges tag')
    Draw.PushButton   ('Set F tag'      , EVT_FTAG      ,  80 , row ,  80 , rh , 'Set faces tag')              ; row -= rh
    Draw.PushButton   ('Gen structured' , EVT_GEN_STRU  ,   0 , row , 100 , rh , 'Generated structured mesh')
    Draw.PushButton   ('Del properties' , EVT_DEL_PROP  , 100 , row , 100 , rh , 'Delete all properties')      ; row -= rh

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
