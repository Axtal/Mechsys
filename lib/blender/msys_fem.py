# Modules
import Blender
import bpy
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

def run_fea(obj, msh, edge_brys, prms, inis):
    # geometry
    g = set_geometry (obj, msh, edge_brys)

    # parameters and initial values
    for i in range(g.nelems()):
        g.ele(i).set_model('LinElastic', 'E=1000 nu=0.25', 'Sx=0.0')

    # solve
    sol = ms.solver('ForwardEuler')
    sol.set_geom(g).set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
    sol.solve()

    # output
    bfn = Blender.sys.expandpath (Blender.Get('filename'))
    key = Blender.sys.basename   (Blender.sys.splitext(bfn)[0])
    ms.write_vtu_equilib(g, key+'.vtu')
    print '[1;34mMechSys[0m: file <'+ key+'.vtu> generated'
