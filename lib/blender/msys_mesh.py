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
import mechsys as ms
import msys_dict as di


def print_timing(func):
    def wrapper(*arg):
        t1  = time.time ()
        res = func      (*arg)
        t2  = time.time ()
        print '[1;34mMechSys[0m: %s took [1;31m%f[0m [1;32mseconds[0m' % (func.func_name, (t2-t1))
        return res
    return wrapper

# =========================================================================== Linear mesh

@print_timing
def gen_frame_mesh(txt=None,cpp=False):
    if txt!=None: Blender.Window.WaitCursor(1)

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
                if is3d: txt.write ('mesh.SetVert      (%d, true, %g,%g,%g); %s\n' % (i, v.co[0], v.co[1], v.co[2], info))
                else:    txt.write ('mesh.SetVert      (%d, true, %g,%g); %s\n'    % (i, v.co[0], v.co[1],          info))
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
                if is3d: txt.write ('mesh.set_vert      (%d, True, %g,%g,%g) %s\n' % (i, v.co[0], v.co[1], v.co[2], info))
                else:    txt.write ('mesh.set_vert      (%d, True, %g,%g) %s\n'    % (i, v.co[0], v.co[1],          info))
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
        Blender.Window.WaitCursor(0)

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
def gen_struct_mesh(gen_script=False,txt=None,show_cursor=False,cpp=False):
    if show_cursor: Blender.Window.WaitCursor(1)

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
                        if i<(len(verts)-1): txt.write ('T(%d,%g,%g,%g), '  % (ve[0],ve[1],ve[2],ve[3]))
                        else:                txt.write ('T(%d,%g,%g,%g);\n' % (ve[0],ve[1],ve[2],ve[3]))
                    else:
                        if i<(len(verts)-1): txt.write ('T(%d,%g,%g,0.0), '  % (ve[0],ve[1],ve[2]))
                        else:                txt.write ('T(%d,%g,%g,0.0);\n' % (ve[0],ve[1],ve[2]))
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
                txt.write ('	bks[%d].SetNx (%d,%g,%d);\n' % (ib,nx,ax,linx))
                txt.write ('	bks[%d].SetNy (%d,%g,%d);\n' % (ib,ny,ay,liny))
                if is3d: txt.write ('	bks[%d].SetNz (%d,%g,%d);\n' % (ib,nz,az,linz))
                txt.write ('\n')

            else: # Python script

                txt.write ('# Block # %d ----------------------------------------------------------------\n'%ib)
                txt.write ('bks[%d].set (%d, # Tag\n' % (ib,v[0]))
                txt.write ('            '+str(verts)+', # Vertices\n')
                txt.write ('            '+str(edges)+', # Edges\n')
                txt.write ('            '+str(etags)+', # Edge Tags\n')
                txt.write ('            '+str(ftags)+', # Face Tags\n')
                txt.write ('            %d, %d, %d, %d) # Origin, X-Plus, Y-Plus, Z-Plus\n' % (origin,xp,yp,zp))
                txt.write ('bks[%d].set_nx (%d,%g,%d)\n' % (ib,nx,ax,linx))
                txt.write ('bks[%d].set_ny (%d,%g,%d)\n' % (ib,ny,ay,liny))
                if is3d: txt.write ('bks[%d].set_nz (%d,%g,%d)\n' % (ib,nz,az,linz))
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
        Blender.Window.WaitCursor(0)
        return txt
    else:
        if len(bks)>0:
            mesh = ms.mesh_structured (is3d)
            if iso2: mesh.set_o2      (True)
            mesh.set_blocks (bks)
            mesh.generate   (True)
            print
            if show_cursor: Blender.Window.WaitCursor(0)
            return mesh


# ========================================================================= Unstructured mesh

@print_timing
def gen_unstruct_mesh(gen_script=False,txt=None,show_cursor=True,cpp=False):
    if show_cursor: Blender.Window.WaitCursor(1)

    # get active object
    edm, obj, msh = di.get_msh()

    # 3D mesh? Quadratic elements ?
    is3d = obj.properties['is3d'] if obj.properties.has_key('is3d') else False
    iso2 = obj.properties['iso2'] if obj.properties.has_key('iso2') else False

    # transform vertices coordinates
    ori = [v for v in msh.verts] # create a copy in local coordinates
    msh.transform (obj.matrix)   # transform mesh to global coordinates

    # number of regions and holes
    nregs = len(obj.properties['regs']) if obj.properties.has_key('regs') else 0
    nhols = len(obj.properties['hols']) if obj.properties.has_key('hols') else 0

    # get etags
    etags = {}
    if obj.properties.has_key('etags'):
        for k, v in obj.properties['etags'].iteritems():
            eid = int(k)
            etags[(msh.edges[eid].v1.index, msh.edges[eid].v2.index)] = v

    # number of vertices and segments
    nverts = len(msh.verts)
    nedges = len(msh.edges)

    if gen_script:

        if cpp: # C++ script

            # unstructured mesh instance
            if txt==None:
                txt = Blender.Text.New(obj.name+'_umesh')
            if is3d: txt.write ('	Mesh::Unstructured mesh(true); // true=>3D\n')
            else:    txt.write ('	Mesh::Unstructured mesh(false); // false=>2D\n')
            if iso2: txt.write ('	mesh.SetO2(); // Quadratic elements\n')

            # set polygon
            txt.write ('	mesh.SetPolySize    (%d,%d,%d,%d);' % (nverts, nedges, nregs, nhols)+' // nVerts, nSegments, nRegs, nHols\n')

            # set vertices and edges
            info = ' // VertIdx, X, Y, Z'
            for v in msh.verts:
                txt.write ('	mesh.SetPolyPoint   (%d,%g,%g,%g);' % (v.index, v.co[0], v.co[1], v.co[2])+info+'\n')
                info = ''
            info = ' // SegmentIdx, VertIdx1, VertIdx2, Tag'
            for e in msh.edges:
                key = (e.v1.index, e.v2.index)
                if key in etags: txt.write ('	mesh.SetPolySegment (%d,%d,%d,%d);' % (e.index, e.v1.index, e.v2.index, etags[key])+info+'\n')
                else:            txt.write ('	mesh.SetPolySegment (%d,%d,%d);'    % (e.index, e.v1.index, e.v2.index)+info+'\n')
                info = ''

            # set regions and holes
            info = ' // RegIdx, Tag, MaxArea, X, Y, Z'
            if obj.properties.has_key('regs'):
                for k, v in obj.properties['regs'].iteritems():
                    txt.write ('	mesh.SetPolyRegion  (%d,%d,%g,%g,%g,%g);' % (int(k), int(v[0]), v[1], v[2], v[3], v[4])+info+'\n')
                    info = ''
            info = ' // HolIdx, X, Y, Z'
            if obj.properties.has_key('hols'):
                for k, v in obj.properties['hols'].iteritems():
                    txt.write ('	mesh.SetPolyHole    (%d,%g,%g,%g);' % (int(k), v[0], v[1], v[2])+info+'\n')
                    info = ''

        else: # Python script

            # unstructured mesh instance
            if txt==None:
                txt = Blender.Text.New(obj.name+'_umesh')
                txt.write ('import Blender, bpy\n')
                txt.write ('import mechsys   as ms\n')
                txt.write ('import msys_mesh as me\n\n')
            if is3d: txt.write ('mesh = ms.mesh_unstructured(True) # True=>3D\n')
            else:    txt.write ('mesh = ms.mesh_unstructured(False) # False=>2D\n')
            if iso2: txt.write ('mesh.set_o2           (True) # True=>Quadratic elements\n')

            # set polygon
            txt.write ('mesh.set_poly_size    (%d,%d,%d,%d)' % (nverts, nedges, nregs, nhols)+' # nVerts, nSegments, nRegs, nHols\n')

            # set vertices and edges
            info = ' # VertIdx, X, Y, Z'
            for v in msh.verts:
                txt.write ('mesh.set_poly_point   (%d,%g,%g,%g)' % (v.index, v.co[0], v.co[1], v.co[2])+info+'\n')
                info = ''
            info = ' # SegmentIdx, VertIdx1, VertIdx2, Tag'
            for e in msh.edges:
                key = (e.v1.index, e.v2.index)
                if key in etags: txt.write ('mesh.set_poly_segment (%d,%d,%d,%d)' % (e.index, e.v1.index, e.v2.index, etags[key])+info+'\n')
                else:            txt.write ('mesh.set_poly_segment (%d,%d,%d)'    % (e.index, e.v1.index, e.v2.index)+info+'\n')
                info = ''

            # set regions and holes
            info = ' # RegIdx, Tag, MaxArea, X, Y, Z'
            if obj.properties.has_key('regs'):
                for k, v in obj.properties['regs'].iteritems():
                    txt.write ('mesh.set_poly_region  (%d,%d,%g,%g,%g,%g)' % (int(k), int(v[0]), v[1], v[2], v[3], v[4])+info+'\n')
                    info = ''
            info = ' # HolIdx, X, Y, Z'
            if obj.properties.has_key('hols'):
                for k, v in obj.properties['hols'].iteritems():
                    txt.write ('mesh.set_poly_hole    (%d,%g,%g,%g)' % (int(k), v[0], v[1], v[2])+info+'\n')
                    info = ''

    else: # Run

        # unstructured mesh instance
        mesh = ms.mesh_unstructured (is3d)
        if iso2: mesh.set_o2        (True)

        # set polygon
        mesh.set_poly_size (nverts, nedges, nregs, nhols)

        # set vertices and edges
        for v in msh.verts:
            mesh.set_poly_point (v.index, v.co[0], v.co[1], v.co[2])
        for e in msh.edges:
            key = (e.v1.index, e.v2.index)
            if key in etags: mesh.set_poly_segment (e.index, e.v1.index, e.v2.index, etags[key])
            else:            mesh.set_poly_segment (e.index, e.v1.index, e.v2.index)

        # set regions and holes
        if obj.properties.has_key('regs'):
            for k, v in obj.properties['regs'].iteritems():
                mesh.set_poly_region (int(k), int(v[0]), v[1], v[2], v[3], v[4])
        if obj.properties.has_key('hols'):
            for k, v in obj.properties['hols'].iteritems():
                mesh.set_poly_hole (int(k), v[0], v[1], v[2])

    # Restore local coordinates
    msh.verts = ori
    if edm: Blender.Window.EditMode(1)

    # generate mesh
    maxa = obj.properties['maxarea'] if obj.properties.has_key('maxarea') else -1.0
    mina = obj.properties['minang']  if obj.properties.has_key('minang')  else -1.0
    if gen_script:
        if cpp:
            if maxa>0: txt.write ('	mesh.SetMaxAreaGlobal  (%g);\n' % (maxa))
            if mina>0: txt.write ('	mesh.SetMinAngleGlobal (%g);\n' % (mina))
            txt.write ('	mesh.Generate       (true); // true=>WithInfo\n')
        else:
            if maxa>0: txt.write ('mesh.set_max_area_global  (%g)\n' % (maxa))
            if mina>0: txt.write ('mesh.set_min_angle_global (%g)\n' % (mina))
            txt.write ('mesh.generate         (True) # True=>WithInfo\n')
        Blender.Window.WaitCursor(0)
        return txt
    else:
        if maxa>0: mesh.set_max_area_global  (maxa)
        if mina>0: mesh.set_min_angle_global (mina)
        mesh.generate (True)
        print
        if show_cursor: Blender.Window.WaitCursor(0)
        return mesh














@print_timing
def gen_unstruct_mesh_new(gen_script=False,txt=None,show_cursor=True,cpp=False):
    if show_cursor: Blender.Window.WaitCursor(1)

    # get active object
    edm, obj, msh = di.get_msh()

    # transform vertices coordinates
    ori = [v for v in msh.verts] # create a copy in local coordinates
    msh.transform (obj.matrix)   # transform mesh to global coordinates

    # 3D mesh? Quadratic elements ?
    is3d = obj.properties['is3d'] if obj.properties.has_key('is3d') else False
    iso2 = obj.properties['iso2'] if obj.properties.has_key('iso2') else False
    ndim = 3 if is3d else 2

    # max area
    maxA = obj.properties['maxarea'] if obj.properties.has_key('maxarea') else -1.0

    # number of regions and holes
    nregs = len(obj.properties['regs']) if obj.properties.has_key('regs') else 0
    nhols = len(obj.properties['hols']) if obj.properties.has_key('hols') else 0

    # number of vertices and segments
    nverts = len(msh.verts)
    nedges = len(msh.edges)

    # get vtags
    vtags = {}
    if obj.properties.has_key('vtags'):
        for k, v in obj.properties['vtags'].iteritems():
            vid = int(k)
            vtags[vid] = v

    # get etags
    etags = {}
    if obj.properties.has_key('etags'):
        for k, v in obj.properties['etags'].iteritems():
            eid = int(k)
            etags[(msh.edges[eid].v1.index, msh.edges[eid].v2.index)] = v

    if gen_script:

        if cpp: # C++ script

            if txt==None: txt = Blender.Text.New(obj.name+'_umesh')
            txt.write ('Mesh::Unstructured mesh(/*NDim*/%d);\n' % (ndim))
            txt.write ('mesh.Set (%d, %d, %d, %d,' % (nverts, nedges, nregs, nhols)+'                // nPoints, nSegments, nRegions, nHoles\n')
            lin = ''
            for v in msh.verts:
                tag = 0
                if v.index in vtags: tag = vtags[v.index]
                lin += ('%4d.,  %4d.,  %6e, %6e, %6e,    \n' % (v.index, tag, v.co[0], v.co[1], v.co[2]))
            if nregs>0:
                for k, v in obj.properties['regs'].iteritems():
                    lin += ('        %4d.,  %6e, %6e, %6e, %8e,\n' % (v[0], v[2], v[3], v[4], v[1]))
            if nhols>0:
                for k, v in obj.properties['hols'].iteritems():
                    lin += ('             %6e, %6e, %6e,    \n' % (v[0], v[1], v[2]))
            txt.write (lin[:len(lin)-2])
            txt.write (');\n')

            # segments
            for e in msh.edges:
                key = (e.v1.index, e.v2.index)
                tag = 0
                if key in etags: tag = etags[key]
                txt.write ('mesh.SetSeg (%d, %d, %4d, %4d);\n' % (e.index, tag, e.v1.index, e.v2.index))

            # generate
            str_o2 = 'true' if iso2 else 'false'
            txt.write ('mesh.Generate (/*O2*/%s, /*GlobalMaxArea*/%g);\n' % (str_o2,maxA))

            Blender.Window.WaitCursor(0)
            return txt

    else: raise Exception ('gen_unstruct_mesh_new: method not implemented yet')

    # restore local coordinates
    msh.verts = ori
    if edm: Blender.Window.EditMode(1)







# ====================================================================================== Draw

@print_timing
def set_elems(obj, elems):
    obj.properties['elems'] = elems
    for sid in obj.properties['stages']:
        stg = 'stg_'+sid
        if not obj.properties[stg].has_key('eatts'): obj.properties[stg]['eatts'] = {} 
        temp = {}
        id   = 0
        while obj.properties[stg]['eatts'].has_key(str(id)): id += 1
        for k, v in elems.iteritems():
            tag = int(v[0])
            vtk = int(v[1])
            if not temp.has_key(tag):
                temp[tag] = True
                old_id    = ''
                for m, n in obj.properties[stg]['eatts'].iteritems(): # find old tag
                    if int(n[0])==tag:
                        old_id = m
                        break
                if old_id=='': # add new
                    tid = di.props_push_new('texts', 'gam=20') # returns text_id  == tid
                    obj.properties[stg]['eatts'][str(id)]    = di.new_eatt_props()
                    obj.properties[stg]['eatts'][str(id)][0] = tag
                    obj.properties[stg]['eatts'][str(id)][1] = di.key('vtk2ety')[vtk]
                    obj.properties[stg]['eatts'][str(id)][2] = di.key('vtk2pty')[vtk]
                    obj.properties[stg]['eatts'][str(id)][7] = tid
                    id += 1
                else:
                    obj.properties[stg]['eatts'][old_id][1] = di.key('vtk2ety')[vtk]
                    obj.properties[stg]['eatts'][old_id][2] = di.key('vtk2pty')[vtk]

@print_timing
def set_ebrys(obj):
    if not di.key('mshsetfem'): return
    # set edges boundaries
    if obj.properties.has_key('etags'):
        for sid in obj.properties['stages']:
            stg = 'stg_'+sid
            if not obj.properties[stg].has_key('ebrys'): obj.properties[stg]['ebrys'] = {} 
            temp = {}
            id   = 0
            while obj.properties[stg]['ebrys'].has_key(str(id)): id += 1
            for k, v in obj.properties['etags'].iteritems():
                tag = int(v)
                if not temp.has_key(tag):
                    temp[tag] = True
                    old_id    = ''
                    for m, n in obj.properties[stg]['ebrys'].iteritems(): # find old tag
                        if int(n[0])==tag:
                            old_id = m
                            break
                    if old_id=='': # add new
                        obj.properties[stg]['ebrys'][str(id)]    = di.new_ebry_props()
                        obj.properties[stg]['ebrys'][str(id)][0] = tag
                        id += 1

@print_timing
def set_fbrys(obj):
    if not di.key('mshsetfem'): return
    # 3D mesh?
    is3d = obj.properties['is3d'] if obj.properties.has_key('is3d') else False
    # set faces boundaries
    if is3d and obj.properties.has_key('ftags'):
        for sid in obj.properties['stages']:
            stg = 'stg_'+sid
            if not obj.properties[stg].has_key('fbrys'): obj.properties[stg]['fbrys'] = {} 
            temp = {}
            id   = 0
            while obj.properties[stg]['fbrys'].has_key(str(id)): id += 1
            for k, v in obj.properties['ftags'].iteritems():
                tag = int(v)
                if not temp.has_key(tag):
                    temp[tag] = True
                    old_id    = ''
                    for m, n in obj.properties[stg]['fbrys'].iteritems(): # find old tag
                        if int(n[0])==tag:
                            old_id = m
                            break
                    if old_id=='': # add new
                        obj.properties[stg]['fbrys'][str(id)]    = di.new_fbry_props()
                        obj.properties[stg]['fbrys'][str(id)][0] = tag
                        id += 1

@print_timing
def add_mesh(obj, mesh, mesh_type):
    # delete old mesh
    if obj.properties.has_key('msh_name'):
        msh_name = obj.properties['msh_name']
        try:
            old_msh  = bpy.data.objects[msh_name]
            scn      = bpy.data.scenes.active
            scn.objects.unlink (old_msh)
        except: pass
    else: msh_name = obj.name+'_msh'

    # delete results
    if obj.properties.has_key('res'): obj.properties.pop('res')

    # add new object/mesh to Blender
    scn                        = bpy.data.scenes.active
    new_msh                    = bpy.data.meshes.new      (msh_name)
    new_obj                    = scn.objects.new (new_msh, msh_name)
    new_obj.restrictSelect     = True
    obj.properties['msh_name'] = new_obj.name
    verts                      = []
    edges                      = []
    mesh.get_verts       (verts)
    mesh.get_edges       (edges)
    new_msh.verts.extend (verts)
    new_msh.edges.extend (edges)
    new_obj.select       (0)
    obj.makeParent       ([new_obj])

    # new stage
    if not obj.properties.has_key('stages'): di.props_push_new_stage() # add to the current active object

    # set elements, edges brys, and faces brys
    elems = {}
    mesh.get_elems (elems)
    set_elems      (obj, elems)
    set_ebrys      (obj)
    set_fbrys      (obj)

    # set mesh type
    obj.properties['mesh_type'] = mesh_type

    # redraw
    Blender.Window.QRedrawAll()
