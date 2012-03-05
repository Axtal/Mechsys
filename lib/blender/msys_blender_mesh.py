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

import math, subprocess
import time
import Blender
import bpy
import msys_blender_dict as di
import mechsys as ms


def print_timing(func):
    def wrapper(*arg):
        t1  = time.time ()
        res = func      (*arg)
        t2  = time.time ()
        print '[1;34mMechSys[0m: %s took [1;31m%f[0m [1;32mseconds[0m' % (func.func_name, (t2-t1))
        return res
    return wrapper


class MeshData:
    def __init__(self):
        # get active object
        self.edm, self.obj, ori_msh = di.get_msh()

        # transform vertices coordinates
        self.msh = ori_msh.copy()
        self.msh.transform (self.obj.matrix) # transform mesh to global coordinates

        # 3D mesh? Quadratic elements ?
        self.is3d = self.obj.properties['is3d'] if self.obj.properties.has_key('is3d') else False
        self.iso2 = self.obj.properties['iso2'] if self.obj.properties.has_key('iso2') else False
        self.ndim = 3 if self.is3d else 2

        # check if all z values are zero in 2D mesh
        if not self.is3d:
            for v in self.msh.verts:
                if abs(v.co[2])>1.0e-7: raise Exception('In 2D meshes, the z coordinate of all points must be zero (z=%g is invalid)'%v.co[2])

        # max area and min angle
        self.maxA   = self.obj.properties['maxarea'] if self.obj.properties.has_key('maxarea') else -1.0
        self.minAng = self.obj.properties['minang']  if self.obj.properties.has_key('minang')  else -1.0

        # number of regions and holes
        self.nregs = len(self.obj.properties['regs']) if self.obj.properties.has_key('regs') else 0
        self.nhols = len(self.obj.properties['hols']) if self.obj.properties.has_key('hols') else 0

        # number of vertices and segments
        self.nverts = len(self.msh.verts)
        self.nedges = len(self.msh.edges)
        self.nfaces = len(self.msh.faces)

        # get vtags
        self.vtags = {}
        if self.obj.properties.has_key('vtags'):
            for k, v in self.obj.properties['vtags'].iteritems():
                vid = int(k)
                self.vtags[vid] = v

        # get etags
        self.etags = {}
        if self.obj.properties.has_key('etags'):
            for k, v in self.obj.properties['etags'].iteritems():
                eid = int(k)
                self.etags[(self.msh.edges[eid].v1.index, self.msh.edges[eid].v2.index)] = v

        # get ftags
        self.ftags = {}
        if self.obj.properties.has_key('ftags'):
            for k, v in self.obj.properties['ftags'].iteritems():
                vids = [int(vid) for vid in k.split('_')]
                vids.sort()
                self.ftags[(vids[0],vids[1],vids[2])] = v

        # get blocks
        if self.obj.properties.has_key('blks'):

            neds2nverts = {4:4, 8:8, 12:8, 24:20}
            neds2edges  = {4:{(0,1):0, (1,2):1, (2,3):2, (0,3):3}, # quad4
                           8:{(0,4):0, (1,4):0, (1,5):1, (2,5):1, (2,6):2, (3,6):2, (3,7):3, (0,7):3}} # quad8
            neds2faces  = {12:{(0,3,4,7):0,(1,2,5,6):1,(0,1,4,5):2,(2,3,6,7):3,(0,1,2,3):4,(4,5,6,7):5}, # hex8
                           24:{(0,3,4,7, 11,15,16,19):0,
                               (1,2,5,6,  9,16,17,18):1,
                               (0,1,4,5,  8,12,16,17):2,
                               (2,3,6,7, 10,14,18,19):3,
                               (0,1,2,3,  8, 9,10,11):4,
                               (4,5,6,7, 12,13,14,15):5}} # hex20

            self.blks = []
            for k, v in self.obj.properties['blks'].iteritems(): # for each block

                # block data
                blk         = {}
                blk['ndim'] = self.ndim
                blk['tag']  = int(v[0])
                blk['nx']   = int(v[8])
                blk['ny']   = int(v[9])
                blk['nz']   = int(v[10])
                nedges      = int(v[17]) # number of edges

                # map: global vertices ids to local vertices ids
                gids2lids = di.block_local_ids (self.obj, v)
                if not gids2lids: raise Exception('Please define all blocks (local axes, divisions, etc.) first')

                # vertices
                blk['V'] = []
                for i in range(neds2nverts[nedges]):
                    if self.is3d: blk['V'].append ([-1, 0.0, 0.0, 0.0])
                    else:         blk['V'].append ([-1, 0.0, 0.0])
                for gvid, lvid in gids2lids.iteritems():
                    if self.is3d:
                        blk['V'][lvid][1] = self.msh.verts[gvid].co[0]
                        blk['V'][lvid][2] = self.msh.verts[gvid].co[1]
                        blk['V'][lvid][3] = self.msh.verts[gvid].co[2]
                    else:
                        blk['V'][lvid][1] = self.msh.verts[gvid].co[0]
                        blk['V'][lvid][2] = self.msh.verts[gvid].co[1]

                # find brytags
                if self.is3d: blk['brytags'] = [0 for i in range(6)]
                else:         blk['brytags'] = [0 for i in range(4)]
                if not self.is3d and len(self.etags)>0: # edge tags
                    for gedpair, etag in self.etags.iteritems():                    # gedpair=global edge pair
                        if gedpair[0] in gids2lids and gedpair[1] in gids2lids:     # if this edge is in the block
                            ledpair = (gids2lids[gedpair[0]],gids2lids[gedpair[1]])     # local edge pair
                            if ledpair[0]>ledpair[1]: ledpair = (ledpair[1],ledpair[0]) # sort local edge pair
                            loceid = neds2edges[nedges][ledpair]                        # local edge id
                            blk['brytags'][loceid] = etag

                if self.is3d and len(self.ftags)>0: # edge tags
                    for gvidtri, ftag in self.ftags.iteritems(): # gvidtup=global vertex triple
                        lvidtri = [gids2lids[gvidtri[0]],gids2lids[gvidtri[1]],gids2lids[gvidtri[2]]] # local vertex triplet
                        locfid  = -1
                        for lvids, lfid in neds2faces[nedges].iteritems():
                            if lvidtri[0] in lvids and lvidtri[1] in lvids and lvidtri[2] in lvids:
                                locfid = lfid
                                break
                        if locfid<0: raise Exception('msys_mesh.py: __internal error__: local face id (locfid) not found')
                        blk['brytags'][locfid] = ftag

                # add blk to blks
                self.blks.append (blk)

    def __del__(self):
        if self.edm: Blender.Window.EditMode(1)

    def get_face_tag(self, msh_face):
        tag  = 0
        vids = [v.index for v in msh_face.verts]
        vids.sort()
        vids = (vids[0],vids[1],vids[2])
        if vids in self.ftags: tag = self.ftags[vids]
        return tag


# =========================================================================== Linear mesh

@print_timing
def gen_frame_mesh(gen_script=False,txt=None,cpp=False,with_headers=True):
    m = MeshData()
    if gen_script:
        key = '.cpp' if cpp else '.py'
        if txt==None: txt = Blender.Text.New(m.obj.name+'_frame'+key)
        if cpp: # C++ script

            # header
            if with_headers:
                txt.write ('// MechSys\n')
                txt.write ('#include <mechsys/fem/fem.h>\n')
                txt.write ('\nint main(int argc, char **argv) try\n')
                txt.write ('{\n')

            # set
            txt.write ('    Mesh::Generic mesh(/*NDim*/%d);\n' % (m.ndim))
            txt.write ('    mesh.SetSize (%d/*vertices*/, %d/*cells*/);\n' % (m.nverts, m.nedges))

            # vertices
            for v in m.msh.verts:
                tag = 0
                if v.index in m.vtags: tag = m.vtags[v.index]
                if m.is3d: txt.write ('    mesh.SetVert (%4d,  %4d,  %6e, %6e, %6e);\n' % (v.index, tag, v.co[0], v.co[1], v.co[2]))
                else:      txt.write ('    mesh.SetVert (%4d,  %4d,  %6e, %6e);\n'      % (v.index, tag, v.co[0], v.co[1]))

            # cells
            for e in m.msh.edges:
                key = (e.v1.index, e.v2.index)
                tag = 0
                if key in m.etags: tag = m.etags[key]
                else: raise Exception('For frames, all edges==cells must have an edge tag')
                txt.write ('    mesh.SetCell (%4d, %4d, Array<int>(%4d, %4d));\n' % (e.index, tag, e.v1.index, e.v2.index))

            # bottom
            if with_headers:
                txt.write ('}\n')
                txt.write ('MECHSYS_CATCH\n')

        else: # Python script
            raise Exception('Frames: python script is not implemented yet')

    else: # run
        raise Exception('Frames: running of simulation through Blender is not implemented yet')


# =========================================================================== Structured mesh

@print_timing
def gen_struct_mesh (gen_script=False,txt=None,cpp=False,with_headers=True):
    m = MeshData()
    if not m.obj.properties.has_key('blks'): raise Exception("Please add blocks first")
    if gen_script:
        key = '.cpp' if cpp else '.py'
        if txt==None: txt = Blender.Text.New(m.obj.name+'_smesh'+key)
        if cpp: # C++ script

            # header
            if with_headers:
                txt.write ('// MechSys\n')
                txt.write ('#include <mechsys/fem/fem.h>\n')
                txt.write ('\nint main(int argc, char **argv) try\n')
                txt.write ('{\n')

            # blocks
            txt.write ('    Array<Mesh::Block> blks(%d);\n'%len(m.blks))
            for i, blk in enumerate(m.blks):
                txt.write ('    blks[%d].Set (/*NDim*/%d, /*Tag*/%d, /*NVert*/%d,\n'%(i, m.ndim, blk['tag'], len(blk['V'])))
                for v in blk['V']:
                    if m.is3d: txt.write ('                  %6.1f, %13.6e, %13.6e, %13.6e,\n'%(v[0],v[1],v[2],v[3]))
                    else:      txt.write ('                  %6.1f, %13.6e, %13.6e,\n'        %(v[0],v[1],v[2]))
                txt.write ('                  ')
                nb = len(blk['brytags'])
                for j, b in enumerate(blk['brytags']):
                    if j==nb-1: txt.write ('%6.1f);\n'%b)
                    else:       txt.write ('%6.1f, '%b)
                txt.write ('    blks[%d].SetNx (%d);\n'%(i,blk['nx']))
                txt.write ('    blks[%d].SetNy (%d);\n'%(i,blk['ny']))
                if m.is3d: txt.write ('    blks[%d].SetNz (%d);\n'%(i,blk['nz']))

            str_o2 = 'true' if m.iso2 else 'false'
            txt.write ('    Mesh::Structured mesh(/*NDim*/%d);\n' % (m.ndim))
            txt.write ('    mesh.Generate (blks, /*O2*/'+str_o2+');\n')
            if m.is3d: txt.write ('    mesh.WriteVTU (\"%s\");\n' % (m.obj.name+'_smesh'))
            else:      txt.write ('    mesh.WriteMPY (\"%s\", true);\n' % (m.obj.name+'_smesh'))

            # bottom
            if with_headers:
                txt.write ('}\n')
                txt.write ('MECHSYS_CATCH\n')

        else: # Python script
            str_o2 = 'True' if m.iso2 else 'False'
            if with_headers: txt.write ('from mechsys import *\n\n')
            txt.write ('blks = '+m.blks.__str__()+'\n')
            txt.write ('mesh = Structured(%d) # ndim\n' % (m.ndim))
            txt.write ('mesh.Generate (blks, '+str_o2+') # blks, o2?\n')
            if m.is3d: txt.write ('mesh.WriteVTU (\"%s\") # FileKey\n' % (m.obj.name+'_smesh'))
            else:      txt.write ('mesh.WriteMPY (\"%s\", True) # FileKey, WithTags\n' % (m.obj.name+'_smesh'))

    else: # run
        mesh = ms.Structured (m.ndim)
        mesh.Generate (m.blks, m.iso2)
        if m.is3d: mesh.WriteVTU (m.obj.name+'_smesh')
        else:      mesh.WriteMPY (m.obj.name+'_smesh', True) # FileKey, WithTags
        add_mesh (m.obj, mesh, 'struct')
        if with_headers: # for FEM
            mesh.WriteVTU (m.obj.name+'_smesh')
        return mesh


# ========================================================================= Unstructured mesh

@print_timing
def gen_unstruct_mesh (gen_script=False,txt=None,cpp=False,with_headers=True):
    m = MeshData()
    if gen_script:
        key = '.cpp' if cpp else '.py'
        if txt==None: txt = Blender.Text.New(m.obj.name+'_umesh'+key)
        if cpp: # C++ script

            # header
            if with_headers:
                txt.write ('// MechSys\n')
                txt.write ('#include <mechsys/fem/fem.h>\n')
                txt.write ('\nint main(int argc, char **argv) try\n')
                txt.write ('{\n')

            # points, regions, and holes
            if m.is3d: num_eorf, txt_eorf = m.nfaces, 'faces'
            else:      num_eorf, txt_eorf = m.nedges, 'segments'
            txt.write ('    Mesh::Unstructured mesh(/*NDim*/%d);\n' % (m.ndim))
            txt.write ('    mesh.Set (%d/*points*/, %d/*%s*/, %d/*regions*/, %d/*holes*/);\n' % (m.nverts, num_eorf, txt_eorf, m.nregs, m.nhols))
            lin = ''
            if m.nregs>0:
                idx = 0
                for k, v in m.obj.properties['regs'].iteritems():
                    if m.is3d: lin += '    mesh.SetReg (%4d, %4d,  %8e,  %6e, %6e, %6e);\n' % (idx, v[0], v[1], v[2], v[3], v[4])
                    else:      lin += '    mesh.SetReg (%4d, %4d,  %8e,  %6e, %6e);\n'      % (idx, v[0], v[1], v[2], v[3])
                    idx += 1
            if m.nhols>0:
                idx = 0
                for k, v in m.obj.properties['hols'].iteritems():
                    if m.is3d: lin += '    mesh.SetHol (%4d, %6e, %6e, %6e,\n' % (idx, v[0], v[1], v[2])
                    else:      lin += '    mesh.SetHol (%4d, %6e, %6e,\n'      % (idx, v[0], v[1])
                    idx += 1
            for v in m.msh.verts:
                tag = 0
                if v.index in m.vtags: tag = m.vtags[v.index]
                if m.is3d: lin += '    mesh.SetPnt (%4d,  %4d,  %6e, %6e, %6e);\n' % (v.index, tag, v.co[0], v.co[1], v.co[2])
                else:      lin += '    mesh.SetPnt (%4d,  %4d,  %6e, %6e);\n'      % (v.index, tag, v.co[0], v.co[1])
            txt.write (lin)
            if m.is3d:
                for idf, f in enumerate(m.msh.faces):
                    lin = ''
                    nvs = len(f.verts)
                    lin += '    mesh.SetFac (%4d, %4d, Array<int>('%(idf,m.get_face_tag(f))
                    for idv, v in enumerate(f.verts):
                        lin += '%d' % v.index
                        if idv==nvs-1: lin += '));\n'
                        else:          lin += ', '
                    txt.write (lin)
            else:
                for e in m.msh.edges:
                    key = (e.v1.index, e.v2.index)
                    tag = 0
                    if key in m.etags: tag = m.etags[key]
                    txt.write ('    mesh.SetSeg (%4d, %4d, %4d, %4d);\n' % (e.index, tag, e.v1.index, e.v2.index))
            str_o2 = 'true' if m.iso2 else 'false'
            txt.write ('    mesh.Generate (/*O2*/%s, /*GlobalMaxArea*/%s, /*Quiet*/true, /*MinAng*/%s);\n' % (str_o2,str(m.maxA),str(m.minAng)))
            if m.is3d: txt.write ('    mesh.WriteVTU (\"%s\");\n' % (m.obj.name+'_umesh'))
            else:      txt.write ('    mesh.WriteMPY (\"%s\", /*WithTags*/true);\n' % (m.obj.name+'_umesh'))

            # bottom
            if with_headers:
                txt.write ('}\n')
                txt.write ('MECHSYS_CATCH\n')

        else: # Python script
            if with_headers: txt.write ('from mechsys import *\n\n')
            txt.write ('mesh = Unstructured(%d)\n' % (m.ndim))
            lin = 'mesh.Set ({\'pts\':['
            for v in m.msh.verts:
                tag = 0
                if v.index in m.vtags: tag = m.vtags[v.index]
                if m.is3d: lin += '[%4d, %6e, %6e, %6e]' % (tag, v.co[0], v.co[1], v.co[2])
                else:      lin += '[%4d, %6e, %6e]'      % (tag, v.co[0], v.co[1])
                if v.index==m.nverts-1: lin += '],\n'
                else:                   lin += ',\n                  '
            lin += '           \'rgs\':['
            if m.nregs>0:
                idx = 0
                for k, v in m.obj.properties['regs'].iteritems():
                    if m.is3d: lin += '[%4d, %8e, %6e, %6e, %6e]' % (v[0], v[1], v[2], v[3], v[4])
                    else:      lin += '[%4d, %8e, %6e, %6e]'      % (v[0], v[1], v[2], v[3])
                    if idx==m.nregs-1: lin += '],\n'
                    else:              lin += ',\n                  '
                    idx += 1
            else: lin += '],\n'
            lin += '           \'hls\':['
            if m.nhols>0:
                idx = 0
                for k, v in m.obj.properties['hols'].iteritems():
                    if m.is3d: lin += '[%6e, %6e, %6e]' % (v[0], v[1], v[2])
                    else:      lin += '[%6e, %6e]'      % (v[0], v[1])
                    if idx==m.nhols-1: lin += '],\n'
                    else:              lin += ',\n                  '
                    idx += 1
            else: lin += '],\n'
            if m.is3d:
                lin += '           \'con\':['
                for idf, f in enumerate(m.msh.faces):
                    lin += '[%4d, [' % m.get_face_tag(f)
                    nvs = len(f.verts)
                    for idv, v in enumerate(f.verts):
                        lin += '%d' % v.index
                        if idv==nvs-1:
                            if idf==m.nfaces-1: lin += ']]]})\n'
                            else:               lin += ']],\n                  '
                        else:                   lin += ', '
            else:
                lin += '           \'con\':['
                idx = 0
                for e in m.msh.edges:
                    key = (e.v1.index, e.v2.index)
                    tag = 0
                    if key in m.etags: tag = m.etags[key]
                    lin += '[%4d, %4d, %4d]' % (tag, e.v1.index, e.v2.index)
                    if idx==m.nedges-1: lin += ']})\n'
                    else:               lin += ',\n                  '
                    idx += 1
            str_o2 = 'True' if m.iso2 else 'False'
            txt.write (lin)
            txt.write ('mesh.Generate (%s, %s, True, %s) # O2, GlobalMaxArea\n' % (str_o2,str(m.maxA),str(m.minAng)))
            if m.is3d: txt.write ('mesh.WriteVTU (\"%s\") # FileKey\n' % (m.obj.name+'_umesh'))
            else:      txt.write ('mesh.WriteMPY (\"%s\", True) # FileKey, WithTags\n' % (m.obj.name+'_umesh'))

    else: # run
        dat = {'pts':[], 'rgs':[], 'hls':[], 'con':[]}
        for v in m.msh.verts:
            tag = 0
            if v.index in m.vtags: tag = m.vtags[v.index]
            if m.is3d: dat['pts'].append ([tag, v.co[0], v.co[1], v.co[2]])
            else:      dat['pts'].append ([tag, v.co[0], v.co[1]])
        if m.nregs>0:
            for k, v in m.obj.properties['regs'].iteritems():
                if m.is3d: dat['rgs'].append ([int(v[0]), v[1], v[2], v[3], v[4]])
                else:      dat['rgs'].append ([int(v[0]), v[1], v[2], v[3]])
        if m.nhols>0:
            for k, v in m.obj.properties['hols'].iteritems():
                if m.is3d: dat['hls'].append ([v[0], v[1], v[2]])
                else:      dat['hls'].append ([v[0], v[1]])
        if m.is3d:
            dat['con'] = []
            for f in m.msh.faces:
                dat['con'].append ([m.get_face_tag(f), [v.index for v in f.verts]])
        else:
            for e in m.msh.edges:
                key = (e.v1.index, e.v2.index)
                tag = 0
                if key in m.etags: tag = m.etags[key]
                dat['con'].append ([tag, e.v1.index, e.v2.index])
        mesh = ms.Unstructured (m.ndim)
        mesh.Set      (dat)
        mesh.Generate (m.iso2, m.maxA, True, m.minAng)
        if m.is3d: mesh.WriteVTU (m.obj.name+'_umesh')
        else:      mesh.WriteMPY (m.obj.name+'_umesh', True) # FileKey, WithTags
        add_mesh (m.obj, mesh, 'unstruct')
        if with_headers: # for FEM
            mesh.WriteVTU (m.obj.name+'_umesh')
        return mesh

# ====================================================================================== Draw

@print_timing
def add_mesh(obj, mesh, msh_type):
    # delete old mesh
    if obj.properties.has_key('msh_name'):
        msh_name = obj.properties['msh_name']
        try:
            old_msh  = bpy.data.objects[msh_name]
            scn      = bpy.data.scenes.active
            scn.objects.unlink (old_msh)
        except: pass
    else: msh_name = obj.name+'_msh'

    # add new object/mesh to Blender
    scn                        = bpy.data.scenes.active
    new_msh                    = bpy.data.meshes.new      (msh_name)
    new_obj                    = scn.objects.new (new_msh, msh_name)
    new_obj.restrictSelect     = True
    obj.properties['msh_name'] = new_obj.name
    verts                      = []
    edges                      = []
    mesh.GetVertsEdges   (verts, edges)
    new_msh.verts.extend (verts)
    new_msh.edges.extend (edges)
    new_obj.select       (0)
    obj.makeParent       ([new_obj])

    # set mesh type
    obj.properties['msh_type'] = msh_type

    # redraw
    Blender.Window.QRedrawAll()


def paraview(unstructured=False):
    obj = di.get_obj()
    if obj.properties.has_key('msh_name'):
        if unstructured: fn = obj.name+'_umesh.vtu'
        else:            fn = obj.name+'_smesh.vtu'
        if not Blender.sys.exists(fn): raise Exception('File <'+fn+'> does not exist (please, generate mesh first)')
        try: pid = subprocess.Popen(['paraview', '--data='+fn]).pid
        except:
            Blender.Window.WaitCursor(0)
            raise Exception('Paraview is not available, please install it first')
    else: raise Exception('Please generate mesh first')
