# Modules

import os
import math
import pickle
import Blender
import bpy
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

        # Transform mesh to global coordinates
        ori = msh.verts[:] # create a copy in local coordinates
        msh.transform (obj.matrix)

        # Vertices
        res = di.get_verts_on_edges_with_tags (obj, msh)
        mg.set_nverts (len(msh.verts))
        for i, v in enumerate(msh.verts):
            if i in res: onbry = True
            else:        onbry = False
            mg.set_vert (i, onbry, v.co[0], v.co[1]) # ,OnBry, Dupl, etc

        # Read elements properties
        nelems = obj.properties['nelems']
        elems  = pickle.load (open(Blender.sys.makename(ext='_elems_'+obj.name+'.pickle'),'r'))

        # Elements
        mg.set_nelems (nelems)
        start = 0
        for i in range(nelems):
            mg.set_elem (i, elems['tags'][i], elems['onbs'][i], elems['vtks'][i], elems['cons'][i], elems['etags'][i])

        # Restore local coordinates
        msh.verts = ori

        # End
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
    fn = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
    ms.write_vtu_equilib (g, fn)
    print '[1;34mMechSys[0m: file <'+fn+'> generated'

    # call ParaView
    os.popen('paraview --data='+fn)

    # restore cursor
    Blender.Window.WaitCursor(0)


def gen_script(obj, nbrys, fbrys, eatts):
    # text
    fn  = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
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
    txt.write ('mechsys.write_vtu_equilib(geom, "'+fn+'")\n')
