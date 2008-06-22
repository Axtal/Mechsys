#!BPY
""" Registration info for Blender menus:
Name: 'TetGen export (.poly)...'
Blender: 241
Group: 'Export'
Tooltip: 'Export active object to TetGen .poly file (.poly)'
"""

__author__  = ("Dorival Pedroso")
__version__ = "2007/03/31"
__bpydoc__  = """\
This script exports to TetGen (.poly) format

You must select at least a mesh object. This mesh object may have
material applied. Each material color will be used to define a MARK.

You can select text3d objects in addition to ONE mesh object.
These text objects will be used to define the holes and 
attributes(regions).

The text for holes must start with "hol_". For example: hol_1
defines the hole number 1. The location of the text is the position
of the hole.

The text for attributes/regions must start with "att_". For example:
att_10 defines the region/attribute number 10. The ATTMARK will be
equal to -10, in this case. The location of the text is the position
of the attribute.

"""

# Copyright (C) 2007: Dorival Pedroso
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

import Blender
try:
    from os.path import exists, join
    pytinst = 1
except:
    print "No Python installed, for full features install Python (http://www.python.org/)."
    pytinst = 0

class TetGenExport:

    def __init__(self, filename, objects):
        self.filename    = filename
        self.file        = open(filename, "w")
        self.usecomments = False
        self.prec        = 3
        self.mesh        = []
        self.mats        = []
        self.nodes       = []
        self.faces       = []
        self.holes       = []
        self.atts        = []
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        for object in objects:
            data     = object.getData()
            objName  = object.getName()
            objType  = object.getType()
            dataType = type(data)
            self.file.write("# objName=%s, objType=%s, dataType=%s\n"% (objName, objType, dataType))
            if dataType==Blender.Types.MeshType or dataType==Blender.Types.NMeshType:
                mworld     = object.matrixWorld
                self.mesh  = data
                self.mesh.transform(mworld)
                self.mats  = self.mesh.materials
                self.nodes = self.mesh.verts
                self.faces = self.mesh.faces
            elif dataType==Blender.Types.Text3dType or dataType==Blender.Types.TextType:
                arr = data.getText().split("_")
                if len(arr)>0:
                    key = arr[0]
                    val = ""
                    des = ""
                    loc = object.matrixWorld[3][0:3]
                    if len(arr)>1:
                        val = arr[1]
                    if key=="hol":
                        hole = [ val, loc ]
                        self.holes.append(hole)
                    elif key=="att":
                        if len(arr)>2:
                            des = arr[2]
                        att = [ val, loc, des ]
                        self.atts.append(att)
            else:
                Blender.Draw.PupMenu('Stopped: Object has an invalid type.')
        if edm: Blender.Window.EditMode(1)

    def rgb2key(self,rgbcol,base=100):
        r = round(rgbcol[0]*base,0)
        g = round(rgbcol[1]*base,0)
        b = round(rgbcol[2]*base,0)
        return int(r*base*base + g*base + b)

    def export(self):

        # Header
        bfile = Blender.sys.expandpath(Blender.Get('filename'))
        self.file.write("# filename = %s\n" % Blender.sys.basename(bfile))
        self.file.write("# generator = Blender %s\n" % Blender.Get('version'))
        self.file.write("\n")

        if len(self.faces) == 0:
            return

        # Marks
        for mat in self.mats:
            self.file.write("#FMARK -%s %s\n" % (self.rgb2key(mat.rgbCol), mat.name))
        self.file.write("\n")
        for att in self.atts:
            self.file.write("#EMARK -%s %s\n" % (att[0], att[2]))
        self.file.write("\n")

        # Part 1 - nodes
        self.file.write("# Part 1 - nodes\n")
        self.file.write("%s 3 0 0  # nNodes 3D nAtts withBryMrks\n" % (len(self.nodes)))
        idx = 1
        for node in self.nodes:
            self.file.write("  %d  %s %s %s\n" % (idx,round(node[0],self.prec), round(node[1],self.prec), round(node[2],self.prec)))
            idx = idx + 1
        self.file.write("\n")

        # Part 2 - facets lists
        self.file.write("# Part 2 - facets lists\n")
        self.file.write("%s  1  # nFacets withBryMrks\n" % (len(self.faces)))
        idx = 1
        for face in self.faces:
            if face.mat<len(self.mats):
                key  = -self.rgb2key(self.mats[face.mat].rgbCol)
                desc = self.mats[face.mat].name
            else:
                key  = -1
                desc = "none"
            ncorners = len(face)
            cordStr  = ""
            for i in range(ncorners):
                indx    = self.nodes.index(face[i])+1
                cordStr = cordStr + "%s " % indx
            if self.usecomments:
                self.file.write("  #--------------- %s\n" % (desc))
            if self.usecomments or idx==1:
                self.file.write("  1 0 %s  # nPolygons nHoles bryMrk (%s)\n" % (key, desc))
            else:
                self.file.write("  1 0 %s\n" % (key))
            self.file.write("  %s  %s\n" % (ncorners, cordStr))
            idx = idx + 1
        self.file.write("\n")

        # Part 3 - volume holes
        self.file.write("# Part 3 - volume holes\n")
        self.file.write("%s  # nHoles\n" % (len(self.holes)))
        idx = 1
        for hole in self.holes:
            self.file.write("  %s  %s %s %s\n" % (idx,round(hole[1][0],self.prec),round(hole[1][1],self.prec),round(hole[1][2],self.prec)))
            idx = idx + 1
        self.file.write("\n")

        # Part 4 - region attributes
        self.file.write("# Part 4 - region attributes\n")
        self.file.write("%s  # nAtts\n" % (len(self.atts)))
        idx = 1
        for att in self.atts:
            self.file.write("  %s  %s %s %s  -%s  -1\n" % (idx,round(att[1][0],self.prec),round(att[1][1],self.prec),round(att[1][2],self.prec),att[0]))
            idx = idx + 1

##########################################################

def select_file(filename):
    if exists(filename):
        result = Blender.Draw.PupMenu("File Already Exists, Overwrite?%t|Yes%x1|No%x0")
        if(result != 1):
            return

    if not filename.endswith('.poly'):
        filename += '.poly'

    objects = Blender.Object.GetSelected()
    tet_export = TetGenExport(filename, objects)
    tet_export.export()

#########################################################

OBJS = Blender.Object.GetSelected()

if not OBJS:
    Blender.Draw.PupMenu('Stopped: Please, select objects.')
else:
    if pytinst==1:
        Blender.Window.FileSelector(select_file,"Export TetGen",Blender.sys.makename(ext='.poly'))
