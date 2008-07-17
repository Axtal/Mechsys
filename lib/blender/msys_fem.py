# Modules
import Blender
import bpy
import os
import math
import mechsys as ms
import msys_dict as di

def fill_mesh(obj):
    if obj!=None and obj.type=='Mesh':
        # Mesh
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)

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
        VTK_TRIANGLE   =  5
        VTK_QUAD       =  9
        VTK_TETRA      = 10
        VTK_HEXAHEDRON = 12
        mg.set_nelems (len(msh.faces))
        for i, f in enumerate(msh.faces):
            eds            = f.edge_keys
            eds_global_ids = [msh.findEdges(vs[0], vs[1]) for vs in eds]
            eds_tags_list  = di.get_tags_list (obj, 'edge', eds_global_ids)
            if len(eds_tags_list)>0: onbry = True
            else:                    onbry = False
            mg.set_elem (i, -1, onbry, VTK_QUAD, [v.index for v in f.verts], eds_tags_list) # Tag, OnBry, ETags, etc

        # Write VTU
        bfn = Blender.sys.expandpath (Blender.Get('filename'))
        key = Blender.sys.basename   (Blender.sys.splitext(bfn)[0])
        mg.write_vtu (key+'_mesh.vtu')
        print '[1;34mMechSys[0m: file <'+ key+'_mesh.vtu> generated'

        # restore local coordinates
        msh.verts = ori

        # end
        if edm: Blender.Window.EditMode(1)
        return mg


def run_fea(obj, nbrys, fbrys, eatts):
    # set cursor
    Blender.Window.WaitCursor(1)

    # mesh
    m = fill_mesh (obj)

    # set geometry
    g = ms.geom(2)
    ms.set_geom (m, nbrys, fbrys, eatts, g)

    # solve
    sol = ms.solver('ForwardEuler')
    sol.set_geom(g).set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
    sol.solve()

    # output
    bfn = Blender.sys.expandpath (Blender.Get('filename'))
    key = Blender.sys.basename   (Blender.sys.splitext(bfn)[0])
    ms.write_vtu_equilib(g, key+'.vtu')
    print '[1;34mMechSys[0m: file <'+ key+'.vtu> generated'

    os.popen('paraview --data='+key+'.vtu')

    # restore cursor
    Blender.Window.WaitCursor(0)


def gen_script(obj, nbrys, fbrys, eatts):
    # text
    bfn = Blender.sys.expandpath (Blender.Get('filename'))
    key = Blender.sys.basename   (Blender.sys.splitext(bfn)[0])
    txt = Blender.Text.New(obj.name+'_script')
    txt.write ('import bpy\n')
    txt.write ('import mechsys\n')
    txt.write ('import msys_fem\n')
    txt.write ('obj   = bpy.data.objects["'+obj.name+'"]\n')
    txt.write ('mesh  = msys_fem.fill_mesh (obj)\n')
    txt.write ('nbrys = '+nbrys.__str__()+'\n')
    txt.write ('fbrys = '+fbrys.__str__()+'\n')
    txt.write ('eatts = '+eatts.__str__()+'\n')
    txt.write ('geom  = mechsys.geom(2)\n')
    txt.write ('mechsys.set_geom (mesh, nbrys, fbrys, eatts, geom)\n')
    txt.write ('sol   = mechsys.solver("ForwardEuler")\n')
    txt.write ('sol.set_geom(geom)\n')
    txt.write ('sol.set_lin_sol("LA").set_num_div(1).set_delta_time(0.0)\n')
    txt.write ('sol.solve()\n')
    txt.write ('mechsys.write_vtu_equilib(geom, "'+key+'.vtu")\n')
