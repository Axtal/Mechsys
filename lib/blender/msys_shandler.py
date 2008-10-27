# SPACEHANDLER.VIEW3D.EVENT

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
