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

import math
import time
import Blender
import bpy
import msys_dict as di
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

        # max area
        self.maxA = self.obj.properties['maxarea'] if self.obj.properties.has_key('maxarea') else -1.0

        # number of regions and holes
        self.nregs = len(self.obj.properties['regs']) if self.obj.properties.has_key('regs') else 0
        self.nhols = len(self.obj.properties['hols']) if self.obj.properties.has_key('hols') else 0

        # number of vertices and segments
        self.nverts = len(self.msh.verts)
        self.nedges = len(self.msh.edges)

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
                           8:[(0,1), (1,2), (2,3), (3,0)],    # quad8
                           12:[],   # hex8
                           24:[]}   # hex20

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
                if len(self.etags)>0: # edge tags
                    for gedpair, etag in self.etags.iteritems():                    # gedpair=global edge pair
                        ledpair = (gids2lids[gedpair[0]],gids2lids[gedpair[1]])     # local edge pair
                        if ledpair[0]>ledpair[1]: ledpair = (ledpair[1],ledpair[0]) # sort local edge pair
                        leid = neds2edges[nedges][ledpair]                          # local edge id
                        blk['brytags'][leid] = etag

                if len(self.ftags)>0: # edge tags
                    print self.ftags

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
def gen_frame_mesh(txt=None,cpp=False):
    # get selected object and mesh
    edm, obj, msh = di.get_msh()

    # 3D mesh?
    is3d = obj.properties['is3d'] if obj.properties.has_key('is3d') else False

    # transform vertices coordinates
    ori = [v for v in msh.verts] # create a copy in local coordinates
    msh.transform (obj.matrix)   # transform mesh to global coordinates

    # MechSys mesh
    mesh = ms.mesh_generic(is3d)

    # etags
    etags = {}
    if obj.properties.has_key('etags'):
        for k, v in obj.properties['etags'].iteritems():
            eid = int(k)
            etags[(msh.edges[eid].v1.index, msh.edges[eid].v2.index)] = v

    if txt!=None: # generate script

        if cpp:
            if is3d: txt.write ('Mesh::Generic mesh(true); // true=>3D\n')
            else:    txt.write ('Mesh::Generic mesh(false); // false=>2D\n')

            # vertices
            info = '// VertIdx, OnBry==True, X, Y, Z'
            txt.write('mesh.SetNVerts    (%d); // set number of vertices\n'%len(msh.verts))
            for i, v in enumerate(msh.verts):
                if is3d: txt.write ('mesh.SetVert      (%d, true, %s,%s,%s); %s\n' % (i, str(v.co[0]), str(v.co[1]), str(v.co[2]), info))
                else:    txt.write ('mesh.SetVert      (%d, true, %s,%s); %s\n'    % (i, str(v.co[0]), str(v.co[1]),               info))
                info = ''

            # elements
            inf1 = '// ElemIdx, ETag, OnBry==True, VTK_LINE==3'
            inf2 = '// ElemIdx, Vert1, Vert2'
            inf3 = '// ElemIdx, LocalEdgeID==0, Tag'
            txt.write('mesh.SetNElems    (%d); // set number of elements\n'%len(msh.edges))
            for i, e in enumerate(msh.edges):
                key = (e.v1.index, e.v2.index)
                if not key in etags: raise Exception('All edges must have a edge tag')
                etag = etags[key]
                txt.write ('mesh.SetElem     (%d,%d,true,VTK_LINE); %s\n' % (i, etag,       inf1)) # 3 == VTK_LINE
                txt.write ('mesh.SetElemCon  (%d,0,%d); %s\n'             % (i, e.v1.index, inf2)) # 0 == local node index
                txt.write ('mesh.SetElemCon  (%d,1,%d); %s\n'             % (i, e.v2.index, inf2)) # 1 == local node index
                txt.write ('mesh.SetElemEtag (%d,0,%d); %s\n'             % (i, etags[key], inf3)) # 0 == local edge index
                inf1 = inf2 = inf3 = ''

        else:
            if is3d: txt.write ('mesh = ms.mesh_generic(True) # True=>3D\n')
            else:    txt.write ('mesh = ms.mesh_generic(False) # False=>2D\n')

            # vertices
            info = '# VertIdx, OnBry==True, X, Y, Z'
            txt.write ('mesh.set_nverts    (%d) # set number of vertices\n' % len(msh.verts))
            for i, v in enumerate(msh.verts):
                if is3d: txt.write ('mesh.set_vert      (%d, True, %s,%s,%s) %s\n' % (i, str(v.co[0]), str(v.co[1]), str(v.co[2]), info))
                else:    txt.write ('mesh.set_vert      (%d, True, %s,%s) %s\n'    % (i, str(v.co[0]), str(v.co[1]),               info))
                info = ''

            # elements
            inf1 = '# ElemIdx, ETag, OnBry==True, VTK_LINE==3'
            inf2 = '# ElemIdx, Vert1, Vert2'
            inf3 = '# ElemIdx, LocalEdgeID==0, Tag'
            txt.write('mesh.set_nelems    (%d) # set number of elements\n'%len(msh.edges))
            for i, e in enumerate(msh.edges):
                key = (e.v1.index, e.v2.index)
                if not key in etags: raise Exception('All edges must have a edge tag')
                etag = etags[key]
                txt.write ('mesh.set_elem      (%d,%d,True,3) %s\n' % (i, etag,       inf1)) # 3 == VTK_LINE
                txt.write ('mesh.set_elem_con  (%d,0,%d) %s\n'      % (i, e.v1.index, inf2)) # 0 == local node index
                txt.write ('mesh.set_elem_con  (%d,1,%d) %s\n'      % (i, e.v2.index, inf2)) # 1 == local node index
                txt.write ('mesh.set_elem_etag (%d,0,%d) %s\n'      % (i, etags[key], inf3)) # 0 == local edge index
                inf1 = inf2 = inf3 = ''

    else: # run

        mesh = ms.mesh_generic(is3d)

        # vertices
        mesh.set_nverts (len(msh.verts))
        for i, v in enumerate(msh.verts):
            if is3d: mesh.set_vert (i, True, v.co[0], v.co[1], v.co[2])
            else:    mesh.set_vert (i, True, v.co[0], v.co[1])

        # elements
        mesh.set_nelems (len(msh.edges))
        for i, e in enumerate(msh.edges):
            key = (e.v1.index, e.v2.index)
            if not key in etags: raise Exception('All edges must have a edge tag')
            etag = etags[key]
            mesh.set_elem      (i,etag,True,3)    # 3 == VTK_LINE
            mesh.set_elem_con  (i, 0, e.v1.index) # 0 == local node index
            mesh.set_elem_con  (i, 1, e.v2.index) # 1 == local node index
            mesh.set_elem_etag (i, 0, etags[key]) # 0 == local edge index

    # Restore local coordinates
    msh.verts = ori
    if edm: Blender.Window.EditMode(1)

    # generate mesh
    if txt==None: return mesh


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
            txt.write ('    Mesh::Unstructured mesh(/*NDim*/%d);\n' % (m.ndim))
            txt.write ('    mesh.Set (%d/*points*/, %d/*segments*/, %d/*regions*/, %d/*holes*/,\n' % (m.nverts, m.nedges, m.nregs, m.nhols))
            lin = ''
            for v in m.msh.verts:
                tag = 0
                if v.index in m.vtags: tag = m.vtags[v.index]
                if m.is3d: lin += '    %4d.,  %4d.,  %6e, %6e, %6e,\n' % (v.index, tag, v.co[0], v.co[1], v.co[2])
                else:      lin += '    %4d.,  %4d.,  %6e, %6e,\n'      % (v.index, tag, v.co[0], v.co[1])
            if m.nregs>0:
                for k, v in m.obj.properties['regs'].iteritems():
                    if m.is3d: lin += '            %4d.,  %6e, %6e, %6e, %8e,\n' % (v[0], v[2], v[3], v[4], v[1])
                    else:      lin += '            %4d.,  %6e, %6e, %8e,\n'      % (v[0], v[2], v[3], v[1])
            if m.nhols>0:
                for k, v in m.obj.properties['hols'].iteritems():
                    if m.is3d: lin += '                 %6e, %6e, %6e,\n' % (v[0], v[1], v[2])
                    else:      lin += '                 %6e, %6e,\n'      % (v[0], v[1])
            txt.write (lin[:len(lin)-2])
            txt.write ('    );\n')
            for e in m.msh.edges:
                key = (e.v1.index, e.v2.index)
                tag = 0
                if key in m.etags: tag = m.etags[key]
                txt.write ('    mesh.SetSeg (%4d, %4d, %4d, %4d);\n' % (e.index, tag, e.v1.index, e.v2.index))
            str_o2 = 'true' if m.iso2 else 'false'
            txt.write ('    mesh.Generate (/*O2*/%s, /*GlobalMaxArea*/%s);\n' % (str_o2,str(m.maxA)))
            if m.is3d: txt.write ('    mesh.WriteVTU (\"%s\");\n' % (m.obj.name+'_umesh'))
            else:      txt.write ('    mesh.WriteMPY (\"%s\", /*WithTags*/true);\n' % (m.obj.name+'_umesh'))

            # bottom
            if with_headers:
                txt.write ('}\n')
                txt.write ('MECHSYS_CATCH\n')

        else: # Python script
            if with_headers: txt.write ('from mechsys import *\n\n')
            txt.write ('mesh = Unstructured(%d)\n' % (m.ndim))
            lin = 'mesh.Set ({\'P\':['
            for v in m.msh.verts:
                tag = 0
                if v.index in m.vtags: tag = m.vtags[v.index]
                if m.is3d: lin += '[%4d, %6e, %6e, %6e]' % (tag, v.co[0], v.co[1], v.co[2])
                else:      lin += '[%4d, %6e, %6e]'      % (tag, v.co[0], v.co[1])
                if v.index==m.nverts-1: lin += '],\n'
                else:                   lin += ',\n                '
            lin += '           \'R\':['
            if m.nregs>0:
                idx = 0
                for k, v in m.obj.properties['regs'].iteritems():
                    if m.is3d: lin += '[%4d, %6e, %6e, %6e, %6e]' % (v[0], v[2], v[3], v[4], v[1])
                    else:      lin += '[%4d, %6e, %6e, %6e]'      % (v[0], v[2], v[3], v[1])
                    if idx==m.nregs-1: lin += '],\n'
                    else:              lin += ',\n                '
                    idx += 1
            else: lin += '],\n'
            lin += '           \'H\':['
            if m.nhols>0:
                idx = 0
                for k, v in m.obj.properties['hols'].iteritems():
                    if m.is3d: lin += '[%6e, %6e, %6e]' % (v[0], v[1], v[2])
                    else:      lin += '[%6e, %6e]'      % (v[0], v[1])
                    if idx==m.nhols-1: lin += '],\n'
                    else:              lin += ',\n                '
                    idx += 1
            else: lin += '],\n'
            if m.is3d:
                lin += '           \'F\':['
                nfs = len(m.msh.faces)
                for idf, f in enumerate(m.msh.faces):
                    lin += '[%4d, [[' % m.get_face_tag(f)
                    nvs = len(f.verts)
                    for idv, v in enumerate(f.verts):
                        lin += '%d' % v.index
                        if idv==nvs-1:
                            if idf==nfs-1: lin += ']]]]})\n'
                            else:          lin += ']]],\n                '
                        else:          lin += ', '
            else:
                lin += '           \'S\':['
                idx = 0
                for e in m.msh.edges:
                    key = (e.v1.index, e.v2.index)
                    tag = 0
                    if key in m.etags: tag = m.etags[key]
                    lin += '[%4d, %4d, %4d]' % (tag, e.v1.index, e.v2.index)
                    if idx==m.nedges-1: lin += ']})\n'
                    else:               lin += ',\n                '
                    idx += 1
            str_o2 = 'True' if m.iso2 else 'False'
            txt.write (lin)
            txt.write ('mesh.Generate (%s, %s) # O2, GlobalMaxArea\n' % (str_o2,str(m.maxA)))
            if m.is3d: txt.write ('mesh.WriteVTU (\"%s\") # FileKey\n' % (m.obj.name+'_umesh'))
            else:      txt.write ('mesh.WriteMPY (\"%s\", True) # FileKey, WithTags\n' % (m.obj.name+'_umesh'))

    else: # run
        dat = {'P':[], 'R':[], 'H':[]}
        if m.is3d: dat['F'] = []
        else:      dat['S'] = []
        for v in m.msh.verts:
            tag = 0
            if v.index in m.vtags: tag = m.vtags[v.index]
            if m.is3d: dat['P'].append ([tag, v.co[0], v.co[1], v.co[2]])
            else:      dat['P'].append ([tag, v.co[0], v.co[1]])
        if m.nregs>0:
            for k, v in m.obj.properties['regs'].iteritems():
                if m.is3d: dat['R'].append ([int(v[0]), v[2], v[3], v[4], v[1]])
                else:      dat['R'].append ([int(v[0]), v[2], v[3], v[1]])
        if m.nhols>0:
            for k, v in m.obj.properties['hols'].iteritems():
                if m.is3d: dat['H'].append ([v[0], v[1], v[2]])
                else:      dat['H'].append ([v[0], v[1]])
        if m.is3d:
            dat['F'] = []
            for f in m.msh.faces:
                dat['F'].append ([m.get_face_tag(f), [[v.index for v in f.verts]]])
        else:
            for e in m.msh.edges:
                key = (e.v1.index, e.v2.index)
                tag = 0
                if key in m.etags: tag = m.etags[key]
                dat['S'].append ([tag, e.v1.index, e.v2.index])
        mesh = ms.Unstructured (m.ndim)
        mesh.Set      (dat)
        mesh.Generate (m.iso2, m.maxA)
        if m.is3d: mesh.WriteVTU (m.obj.name+'_umesh')
        else:      mesh.WriteMPY (m.obj.name+'_umesh', True) # FileKey, WithTags
        add_mesh (m.obj, mesh, 'unstruct')
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
