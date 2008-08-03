# SPACEHANDLER.VIEW3D.EVENT

import Blender
from   Blender import Draw
import bpy

evt = Blender.event
if evt==Draw.RIGHTMOUSE or evt==Draw.AKEY:
    Blender.Window.QRedrawAll()
    evt = None # tell that we want to return the event
else: evt = None # tell that we want to return the event

# set Blender.event to None so grabbed events won't be further processed by Blender
if evt: Blender.event = None
