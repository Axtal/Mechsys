# Modules
import Blender
import bpy
import os
import math
import mechsys as ms
import msys_dict as di

def set_geometry(obj, msh, edge_brys):
    # Read a Blender object, with tags, and set a MechSys/FEM geometry object
    #
    # Example:
    #           edge_brys = [[-10, -20], ['uy', 'fy'], [0.0, -1]]

    # FE geometry structure
    g = ms.geom(2)

    # set nodes
    nn = len (msh.verts)
    g.set_nnodes (nn)
    for v in msh.verts:
        g.set_node (v.index, v.co[0], v.co[1])

    # set elements
    ne = len (msh.faces)
    g.set_nelems (ne)
    for f in msh.faces:
        if len(f.verts)==3: g.set_elem(f.index, 'Tri3PStrain',  1)
        else              : g.set_elem(f.index, 'Quad4PStrain', 1)
        for i in range(len(f.verts)):
            g.ele(f.index).connect (i, g.nod(f.verts[i].index))

    # set brys
    for f in msh.faces:
        eds            = f.edge_keys
        eds_global_ids = [msh.findEdges(vs[0], vs[1]) for vs in eds]
        eds_tags_list  = di.get_tags_list (obj, 'edge', eds_global_ids)
        ntags          = len(eds_tags_list)
        if ntags>0:
            for j in range(ntags):
                tag = eds_tags_list[j]
                if tag<0:
                    if tag in edge_brys[0]:
                        idx = edge_brys[0].index(tag)
                        g.ele(f.index).bry(edge_brys[1][idx], edge_brys[2][idx], j)

    # return FE geometry structure
    return g

def fill_mesh2d(obj, msh):

    # MechSys::Mesh::Generic
    mg = ms.mesh_generic()

    # transform mesh to global coordinates
    ori = msh.verts[:] # create a copy in local coordinates
    msh.transform (obj.matrix)

    # Vertices
    res = di.get_verts_on_edges_with_tags (obj, msh)
    mg.set_nverts (len(msh.verts))
    for i, v in enumerate(msh.verts):
        if i in res: onbry = True
        else:        onbry = False
        mg.set_vert (i, onbry, v.co[0], v.co[1]) # ,OnBry, Dupl, etc


    # Elements
    mg.set_nelems (len(msh.faces))
    for i, f in enumerate(msh.faces):
        eds            = f.edge_keys
        eds_global_ids = [msh.findEdges(vs[0], vs[1]) for vs in eds]
        eds_tags_list  = di.get_tags_list (obj, 'edge', eds_global_ids)
        if len(eds_tags_list)>0: onbry = True
        else:                    onbry = False
        mg.set_elem (i, -1, onbry, [v.index for v in f.verts], eds_tags_list) # Tag, OnBry, ETags, etc

    #mg.write_vtu('toto.vtu') # does NOT work because Mesh::General does not know the type of element (default is Lin2)

    # restore local coordinates
    msh.verts = ori

    return mg

def run_fea2d(obj, msh, nbrys, fbrys, eatts):
    # set cursor
    Blender.Window.WaitCursor(1)

    # mesh
    m2d = fill_mesh2d (obj, msh)

    # set geometry
    g2d = ms.geom(2)
    ms.set_geom (m2d, nbrys, fbrys, eatts, g2d)

    # solve
    sol = ms.solver('ForwardEuler')
    sol.set_geom(g2d).set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
    sol.solve()

    # output
    bfn = Blender.sys.expandpath (Blender.Get('filename'))
    key = Blender.sys.basename   (Blender.sys.splitext(bfn)[0])
    ms.write_vtu_equilib(g2d, key+'.vtu')
    print '[1;34mMechSys[0m: file <'+ key+'.vtu> generated'

    os.popen('paraview --data='+key+'.vtu')

    # restore cursor
    Blender.Window.WaitCursor(0)
