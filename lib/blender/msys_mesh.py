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
from   mechsys import *


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
        self.edm, self.obj, self.msh = di.get_msh()

        # transform vertices coordinates
        self.ori = [v for v in self.msh.verts] # create a copy in local coordinates
        self.msh.transform (self.obj.matrix)   # transform mesh to global coordinates

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

    def __del__(self):
        # restore local coordinates
        self.msh.verts = self.ori
        if self.edm: Blender.Window.EditMode(1)


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
def gen_struct_mesh(gen_script=False,txt=None,cpp=False):
    # get active object
    edm, obj, msh = di.get_msh()

    # check for blocks
    if not obj.properties.has_key('blks'):
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, define blocks first')
    nblks = len(obj.properties['blks'])
    if nblks<1:
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, define blocks first')

    # 3D mesh? Quadratic elements ?
    is3d = obj.properties['is3d'] if obj.properties.has_key('is3d') else False
    iso2 = obj.properties['iso2'] if obj.properties.has_key('iso2') else False

    # transform vertices coordinates
    ori = [v for v in msh.verts] # create a copy in local coordinates
    msh.transform (obj.matrix)   # transform mesh to global coordinates

    # create text for script
    if gen_script:

        if cpp:
            if txt==None:
                txt = Blender.Text.New(obj.name+'_smesh')
            txt.write ('	Array<Mesh::Block> bks(%d);\n\n' % nblks)

        else:
            if txt==None:
                txt = Blender.Text.New(obj.name+'_smesh')
                txt.write ('import Blender, bpy\n')
                txt.write ('import mechsys   as ms\n')
                txt.write ('import msys_mesh as me\n\n')
            txt.write ('bks = [ms.mesh_block() for i in range(%d)]\n\n' % nblks)

    # generate list with blocks and face colors
    ib  = 0
    bks = []
    for k, v in obj.properties['blks'].iteritems():

        # origin and local system
        origin, xp, yp, zp = int(v[4]), int(v[5]), int(v[6]), int(v[7])
        if origin<0:
            msh.verts = ori # restore local coordinates
            if edm: Blender.Window.EditMode(1)
            raise Exception('Please, assign local axes to this block (%d:%d)'%(int(k),int(v[0])))
        if is3d and zp<0:
            msh.verts = ori # restore local coordinates
            if edm: Blender.Window.EditMode(1)
            raise Exception('Please, define the Z-axis of this (3D) block (%d:%d)'%(int(k),int(v[0])))

        # divisions and weights
        nx,   ny,   nz   = int(v[ 8]), int(v[ 9]), int(v[10])
        ax,   ay,   az   =     v[11],      v[12],      v[13]
        linx, liny, linz = int(v[14]), int(v[15]), int(v[16])

        # vertices and edges
        verts = []
        edges = []
        eids  = []
        for i in range(18,18+int(v[17])): eids.append(int(v[i]))
        for e in eids:
            v1 = msh.edges[e].v1.index
            v2 = msh.edges[e].v2.index
            found_v1 = False
            found_v2 = False
            for ve in verts:
                if ve[0]==v1: found_v1 = True
                if ve[0]==v2: found_v2 = True
            if is3d:
                if not found_v1: verts.append ((v1, msh.verts[v1].co[0], msh.verts[v1].co[1], msh.verts[v1].co[2]))
                if not found_v2: verts.append ((v2, msh.verts[v2].co[0], msh.verts[v2].co[1], msh.verts[v2].co[2]))
            else:
                if not found_v1: verts.append ((v1, msh.verts[v1].co[0], msh.verts[v1].co[1]))
                if not found_v2: verts.append ((v2, msh.verts[v2].co[0], msh.verts[v2].co[1]))
            edges.append((v1,v2))

        # etags
        etags = []
        if obj.properties.has_key('etags'):
            for m, n in obj.properties['etags'].iteritems():
                e = int(m)
                if e in eids: etags.append((msh.edges[e].v1.index, msh.edges[e].v2.index, n))

        # ftags
        ftags = []
        veids = [x[0] for x in verts]
        if is3d and obj.properties.has_key('ftags'):
            for m, n in obj.properties['ftags'].iteritems():
                vids = [int(vid) for vid in m.split('_')]
                face_is_in_block = True
                for vid in vids:
                    if not vid in veids:
                        face_is_in_block = False
                        break
                if face_is_in_block: ftags.append ((vids[0],vids[1],vids[2],vids[3],n))

        # new block
        if gen_script:

            if cpp: # C++ script

                txt.write ('	// Block # %d ----------------------------------------------------------------\n'%ib)
                txt.write ('	Mesh::Verts_T ve%d(%d);\n' % (ib,len(verts)))
                txt.write ('	Mesh::Edges_T ed%d(%d);\n' % (ib,len(edges)))
                if len(etags)>0:          txt.write ('	Mesh::ETags_T et%d(%d);\n' % (ib,len(etags)))
                if is3d and len(ftags)>0: txt.write ('	Mesh::FTags_T ft%d(%d);\n' % (ib,len(ftags)))
                txt.write ('	ve%d = ' % (ib))
                for i, ve in enumerate(verts):
                    if is3d:
                        if i<(len(verts)-1): txt.write ('T(%d,%s,%s,%s), '  % (ve[0],str(ve[1]),str(ve[2]),str(ve[3])))
                        else:                txt.write ('T(%d,%s,%s,%s);\n' % (ve[0],str(ve[1]),str(ve[2]),str(ve[3])))
                    else:
                        if i<(len(verts)-1): txt.write ('T(%d,%s,%s,0.0), '  % (ve[0],str(ve[1]),str(ve[2])))
                        else:                txt.write ('T(%d,%s,%s,0.0);\n' % (ve[0],str(ve[1]),str(ve[2])))
                txt.write ('	ed%d = ' % (ib))
                for i, ed in enumerate(edges):
                    if i<(len(edges)-1): txt.write ('T(%d,%d), '  % (ed[0],ed[1]))
                    else:                txt.write ('T(%d,%d);\n' % (ed[0],ed[1]))
                if len(etags)>0:
                    txt.write ('	et%d = ' % (ib))
                    for i, et in enumerate(etags):
                        if i<(len(etags)-1): txt.write ('T(%d,%d,%d), '  % (et[0],et[1],et[2]))
                        else:                txt.write ('T(%d,%d,%d);\n' % (et[0],et[1],et[2]))
                if len(ftags)>0:
                    txt.write ('	ft%d = ' % (ib))
                    for i, ft in enumerate(ftags):
                        if i<(len(ftags)-1): txt.write ('T(%d,%d,%d,%d,%d), '  % (ft[0],ft[1],ft[2],ft[3],ft[4]))
                        else:                txt.write ('T(%d,%d,%d,%d,%d);\n' % (ft[0],ft[1],ft[2],ft[3],ft[4]))
                txt.write ('	bks[%d].Set (%d, ve%d, ed%d, ' % (ib,v[0],ib,ib))
                if len(etags)>0: txt.write ('&et%d, ' % ib)
                else:            txt.write ('NULL, ')
                if len(ftags)>0: txt.write ('&ft%d, ' % ib)
                else:            txt.write ('NULL, ')
                txt.write ('%d, %d, %d, %d);\n' % (origin,xp,yp,zp))
                txt.write ('	bks[%d].SetNx (%d,%s,%d);\n' % (ib,nx,str(ax),linx))
                txt.write ('	bks[%d].SetNy (%d,%s,%d);\n' % (ib,ny,str(ay),liny))
                if is3d: txt.write ('	bks[%d].SetNz (%d,%s,%d);\n' % (ib,nz,az,linz))
                txt.write ('\n')

            else: # Python script

                txt.write ('# Block # %d ----------------------------------------------------------------\n'%ib)
                txt.write ('bks[%d].set (%d, # Tag\n' % (ib,v[0]))
                txt.write ('            '+str(verts)+', # Vertices\n')
                txt.write ('            '+str(edges)+', # Edges\n')
                txt.write ('            '+str(etags)+', # Edge Tags\n')
                txt.write ('            '+str(ftags)+', # Face Tags\n')
                txt.write ('            %d, %d, %d, %d) # Origin, X-Plus, Y-Plus, Z-Plus\n' % (origin,xp,yp,zp))
                txt.write ('bks[%d].set_nx (%d,%s,%d)\n' % (ib,nx,str(ax),linx))
                txt.write ('bks[%d].set_ny (%d,%s,%d)\n' % (ib,ny,str(ay),liny))
                if is3d: txt.write ('bks[%d].set_nz (%d,%s,%d)\n' % (ib,nz,str(az),linz))
                txt.write ('\n')

        else: # Run

            bks.append     (ms.mesh_block())
            bks[-1].set    (int(v[0]), verts, edges, etags, ftags, origin, xp, yp, zp)
            bks[-1].set_nx (nx,ax,linx)
            bks[-1].set_ny (ny,ay,liny)
            bks[-1].set_nz (nz,az,linz)

        ib += 1

    # Restore local coordinates
    msh.verts = ori
    if edm: Blender.Window.EditMode(1)

    # generate mesh
    if gen_script:
        if cpp:
            txt.write('	// Generate\n')
            if is3d: txt.write('	Mesh::Structured mesh(true); // true=>3D\n')
            else:    txt.write('	Mesh::Structured mesh(false); // false=>2D\n')
            if iso2: txt.write('	mesh.SetO2     ();\n')
            txt.write('	mesh.SetBlocks (bks);\n')
            txt.write('	mesh.Generate  (true); // true=>WithInfo\n')
        else:
            txt.write('# Generate\n')
            if is3d: txt.write('mesh = ms.mesh_structured(True) # True=>3D\n')
            else:    txt.write('mesh = ms.mesh_structured(False) # False=>2D\n')
            if iso2: txt.write('mesh.set_o2     (True) # True=>Quadratic elements\n')
            txt.write('mesh.set_blocks (bks)\n')
            txt.write('mesh.generate   (True) # True=>WithInfo\n')
        return txt
    else:
        if len(bks)>0:
            mesh = ms.mesh_structured (is3d)
            if iso2: mesh.set_o2      (True)
            mesh.set_blocks (bks)
            mesh.generate   (True)
            print
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
            txt.write ('    mesh.WriteMPY (\"%s\", /*WithTags*/true);\n' % (m.obj.name+'_umesh'))

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
                pass
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
            txt.write ('mesh.WriteMPY (\"%s\", True) # FileKey, WithTags\n' % (m.obj.name+'_umesh'))

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
        for e in m.msh.edges:
            key = (e.v1.index, e.v2.index)
            tag = 0
            if key in m.etags: tag = m.etags[key]
            dat['S'].append ([tag, e.v1.index, e.v2.index])
        mesh = Unstructured (m.ndim)
        mesh.Set      (dat)
        mesh.Generate (m.iso2, m.maxA)
        add_mesh      (m.obj, mesh, 'unstruct')
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
