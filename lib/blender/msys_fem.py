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
        if obj.properties['is3d']: mg.set_3d()

        # Transform mesh to global coordinates
        ori = msh.verts[:] # create a copy in local coordinates
        msh.transform (obj.matrix)

        # Vertices
        mg.set_nverts (len(msh.verts))
        for i, v in enumerate(msh.verts):
            onbry = True if (i in obj.properties['verts_bry']) else False
            if mg.is_3d(): mg.set_vert (i, onbry, v.co[0], v.co[1], v.co[2])
            else:          mg.set_vert (i, onbry, v.co[0], v.co[1])

        # Restore local coordinates
        msh.verts = ori

        # Elements
        nelems = obj.properties['nelems']
        mg.set_nelems (nelems)
        for i in range(nelems):
            # element
            mg.set_elem (i, obj.properties['elems']['tags'][i],
                            obj.properties['elems']['onbs'][i],
                            obj.properties['elems']['vtks'][i])

            # connectivities
            for j in range(len(obj.properties['elems']['cons'] [str(i)])):
                mg.set_elem_con (i, j, obj.properties['elems']['cons'] [str(i)][j])

            # edge tags
            for j in range(len(obj.properties['elems']['etags'][str(i)])):
                mg.set_elem_etag (i, j, obj.properties['elems']['etags'][str(i)][j])

            # face tags
            if obj.properties['is3d']:
                for j in range(len(obj.properties['elems']['ftags'][str(i)])):
                    mg.set_elem_ftag (i, j, obj.properties['elems']['ftags'][str(i)][j])

        # End
        if edm: Blender.Window.EditMode(1)
        return mg


def run_fea(obj):
    # set cursor
    Blender.Window.WaitCursor(1)

    # mesh
    m = fill_mesh (obj)

    # boundary conditions
    nbrys = di.get_nbrys_numeric (obj)
    ebrys = di.get_ebrys_numeric (obj)
    fbrys = di.get_fbrys_numeric (obj)
    eatts = di.get_eatts_numeric (obj)

    # set geometry
    g = ms.geom(3) if m.is_3d() else ms.geom(2)
    ms.set_nodes_elems (m, eatts, g)
    ms.set_brys        (m, nbrys, ebrys, fbrys, g)

    # solve
    sol = ms.solver('ForwardEuler')
    sol.set_geom(g).set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
    sol.solve()

    # output
    fn = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
    ms.out_vtu (g, fn)
    print '[1;34mMechSys[0m: file <'+fn+'> generated'

    # call ParaView
    os.popen('paraview --data='+fn)

    # restore cursor
    Blender.Window.WaitCursor(0)


def gen_script(obj):
    # boundary conditions
    nbrys = di.get_nbrys_numeric (obj)
    ebrys = di.get_ebrys_numeric (obj)
    fbrys = di.get_fbrys_numeric (obj)
    eatts = di.get_eatts_numeric (obj)

    # text
    ndim = 3 if obj.properties['is3d'] else 2
    fn   = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
    txt  = Blender.Text.New(obj.name+'_script')
    txt.write ('import bpy\n')
    txt.write ('import mechsys\n')
    txt.write ('import msys_fem\n')
    txt.write ('obj   = bpy.data.objects["'+obj.name+'"]\n')
    txt.write ('mesh  = msys_fem.fill_mesh (obj)\n')
    txt.write ('nbrys = '+nbrys.__str__()+'\n')
    txt.write ('ebrys = '+ebrys.__str__()+'\n')
    txt.write ('fbrys = '+fbrys.__str__()+'\n')
    txt.write ('eatts = '+eatts.__str__()+'\n')
    txt.write ('geom  = mechsys.geom('+str(ndim)+')\n')
    txt.write ('mechsys.set_nodes_elems (mesh, eatts, geom)\n')
    txt.write ('mechsys.set_brys        (mesh, nbrys, ebrys, fbrys, geom)\n')
    txt.write ('sol   = mechsys.solver("ForwardEuler")\n')
    txt.write ('sol.set_geom(geom)\n')
    txt.write ('sol.set_lin_sol("LA").set_num_div(1).set_delta_time(0.0)\n')
    txt.write ('sol.solve()\n')
    txt.write ('mechsys.out_vtu(geom, "'+fn+'")\n')
