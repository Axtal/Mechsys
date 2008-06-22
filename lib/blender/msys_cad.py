#!BPY

"""
Name: 'MechSysCAD'
Blender: 2.45
Group: 'Themes'
Tooltip: 'CAD in Blender'
"""
__author__  = "Dorival Pedroso"
__version__ = "06.05.08"
__bpydoc__="""\
TODO
"""

# Modules
import Blender
from   Blender import Draw
from   Blender.Mathutils import Vector
import msys_draw as dr
import msys_dict as di

# Constants
EVT_VPOINT   = 1 # type in values to define point(s)
EVT_FPOINT   = 2 # read points from file
EVT_FSPLINE  = 3 # create a spline from points in a file
EVT_FILLET   = 4 # 3D fillet
EVT_BREAK    = 5 # break edge
EVT_READ_2DM = 6 # read 2D mesh from Triangle

# Handle input events
def event(evt, val):
    if evt==Draw.QKEY and not val: Draw.Exit()
    if evt==Draw.ESCKEY: Draw.Redraw(1)

# Handle button events
def button_event(evt):
    dict = di.load_dict()
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
    elif evt==EVT_FPOINT:   Blender.Window.FileSelector(dr.add_points_from_file, "Read X Y Z cols")
    elif evt==EVT_FSPLINE:  Blender.Window.FileSelector(dr.add_spline_from_file, "Read X Y Z cols")
    elif evt==EVT_FILLET:   dr.fillet(dict['fillet_radius'],dict['fillet_steps'])
    elif evt==EVT_BREAK:    dr.break_edge()
    elif evt==EVT_READ_2DM: Blender.Window.FileSelector(dr.read_2d_mesh, "Read 2D mesh")

def fillet_radius_callback(evt,val):
    dict = di.load_dict()
    dict['fillet_radius'] = val

def fillet_steps_callback(evt,val):
    dict = di.load_dict()
    dict['fillet_steps'] = val

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
    Draw.PushButton("Point(s) x y z" , EVT_VPOINT   ,   0 , row , 100 , rh , "Type in values to define point(s)")
    Draw.PushButton("Point(s) [file]", EVT_FPOINT   , 100 , row , 100 , rh , "Read point list from file")          ; row=row+rh
    Draw.PushButton("Spline [file]"  , EVT_FSPLINE  ,   0 , row , 100 , rh , "Add a spline with points from file") ; row=row+rh
    Draw.PushButton("Fillet edge"    , EVT_FILLET   ,   0 , row , 100 , rh , "Fillet two edges")                   ;
    Draw.Number    ("Radius"         , 1000         , 100 , row , 100 , rh , d['fillet_radius'],0.0,1.0e+9,"Set radius for fillet operation",fillet_radius_callback)
    Draw.Number    ("Steps"          , 1000         , 200 , row , 100 , rh , d['fillet_steps' ],1  ,100   ,"Set steps for fillet operation" ,fillet_steps_callback) ; row=row+rh
    Draw.PushButton("Break edge"     , EVT_BREAK    ,   0 , row , 100 , rh , "Break an edge")                      ; row=row+rh
    Draw.PushButton("Read 2D mesh"   , EVT_READ_2DM ,   0 , row , 100 , rh , "Read 2D mesh from Triangle")         ; row=row+rh

Draw.Register (gui, event, button_event)
