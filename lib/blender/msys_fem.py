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

import subprocess
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

    else: raise Exception('Object must be of type Mesh')


def set_geo(obj,nbrys,ebrys,fbrys,eatts):
    if (len(eatts)<1): raise Exception('Element attributes MUST be defined in order to run a simulation')
    ndim = 3 if obj.properties['is3d'] else 2
    mesh = fill_mesh   (obj)
    geom = ms.geom     (ndim)
    ms.set_nodes_elems (mesh, eatts, geom)
    ms.set_brys        (mesh, nbrys, ebrys, fbrys, geom)
    return geom


def set_geo_linele(obj,nbrys,eatts):
    if (len(eatts)<1): raise Exception('Element attributes MUST be defined in order to run a simulation')
    if obj!=None and obj.type=='Mesh':
        # Check
        if len(obj.getAllProperties())<1: raise Exception('Please, set rods (edges) IDs first')

        # Mesh
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)

        # Transform mesh to global coordinates
        ori = msh.verts[:] # create a copy in local coordinates
        msh.transform (obj.matrix)

        # is 3d ?
        is_3d = False
        for i, v in enumerate(msh.verts):
            if i>0:
                if abs(v.co[2]-msh.verts[0].co[2])>1.0e-4:
                    is_3d = True
                    break
        ndim = 3 if is_3d else 2

        # FEM geometry
        g = ms.geom (ndim)

        # Vertices
        g.set_nnodes (len(msh.verts))
        for i, v in enumerate(msh.verts):
            g.set_node (i, v.co[0], v.co[1])

        # Restore local coordinates
        msh.verts = ori

        # Edge tags
        #etags = di.get_edgetags_ (obj)
        #def get_edgetags_(obj):
            # Out:   {id1:tag1, id2:tag2, ... num edges with tags}
            #res = {}
            #for p in obj.getAllProperties():
                #if p.name[:4]=='edge':
                    #eid      = int(p.name[5:])
                    #res[eid] = p.data
        #return res 
        if len(etags)!=len(msh.edges): raise Exception('All rods must have an edge tag')

        # Elements
        g.set_nelems (len(msh.edges))
        for i, e in enumerate(msh.edges):
            tag   = etags[i]
            found = False
            for ea in eatts:
                if ea[0] == tag:
                    g.set_elem         (i, ea[1], 1)
                    g.ele(i).connect   (0, g.nod(e.v1.index)).connect(1, g.nod(e.v2.index))
                    g.ele(i).set_model (ea[2],ea[3],ea[4])
                    found = True
                    break
            if not found: raise Exception('Tag = %d MUST be defined in elements attributes'%tag)

        # Boundary conditions
        for i in range(g.nnodes()):
            x = g.nod(i).x()
            y = g.nod(i).y()
            for nb in nbrys:
                d = math.sqrt((x-nb[0])**2+(y-nb[1])**2)
                if d<1.0e-5:
                    g.nod(i).bry(nb[3],nb[4])

        # End
        if edm: Blender.Window.EditMode(1)
        return g

    else: raise Exception('Object must be of type Mesh')


def store_scalars(geo,obj):
    obj.properties['scalars'] = {}
    keys    = ['ux','uy','uz','fx','fy','fz','u','q']
    vals    = [[]   for i in range(len(keys))]
    has_key = [True for i in range(len(keys))]
    for i, key in enumerate(keys):
        try:    val = geo.nod(0).val(key)
        except: has_key[i] = False
    for i in range(geo.nnodes()):
        for j, key in enumerate(keys):
            if has_key[j]: vals[j].append(geo.nod(i).val(key))
    for i, key in enumerate(keys):
        if has_key[i]: obj.properties['scalars'][key] = vals[i]


def get_brys_atts(obj):
    d = di.load_dict()

    # nbrys
    nbrys = []
    if obj.properties.has_key('nbrys'):
        for k, v in obj.properties['nbrys'].iteritems():
            nbrys.append([v[0], v[1], v[2], d['dofvars'][int(v[3])], v[4]])

    # nbrysID
    nbrysID = []
    if obj.properties.has_key('nbrysID'):
        for k, v in obj.properties['nbrysID'].iteritems():
            nbrysID.append([int(v[0]), d['dofvars'][int(v[1])], v[2]])

    # ebrys
    ebrys = []
    if obj.properties.has_key('ebrys'):
        for k, v in obj.properties['ebrys'].iteritems():
            ebrys.append([int(k), d['dofvars'][int(v[0])], v[1]])

    # fbrys
    fbrys = []
    if obj.properties.has_key('fbrys'):
        for k, v in obj.properties['fbrys'].iteritems():
            fbrys.append([int(k), d['dofvars'][int(v[0])], v[1]])

    # eatts
    eatts = []
    if obj.properties.has_key('eatts'):
        for k, v in obj.properties['eatts'].iteritems():
            r = v.split()
            eatts.append([int(k), d['etypes'][int(r[0])], d['models'][int(r[1])], r[2].replace('_',' '), r[3].replace('_',' ')])

    return nbrys, nbrysID, ebrys, fbrys, eatts


def run_analysis():
    obj = di.get_obj()

    # set cursor
    Blender.Window.WaitCursor(1)

    # check for properties
    linele = False
    try:    is_3d  = obj.properties['is3d']
    except: linele = True # assuming linear structures (trusses, beams)

    # boundary conditions & properties
    nbrys, nbrysID, ebrys, fbrys, eatts = get_brys_atts (obj)
    if (len(eatts)<1): raise Exception('Element attributes MUST be defined in order to run a simulation')

    # problem geometry
    if linele: geo = set_geo_linele (obj,nbrys,eatts)
    else:      geo = set_geo        (obj,nbrys,ebrys,fbrys,eatts)

    # nodes boundary conditions
    for nb in nbrysID:
        geo.nod(nb[0]).bry(nb[1],nb[2])

    # solution
    sol = ms.solver('ForwardEuler')
    sol.set_geom(geo)
    #sol.set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
    sol.solve()

    # save results in object
    store_scalars(geo,obj)

    # output
    fn = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
    ms.out_vtu (geo, fn)
    print '[1;34mMechSys[0m: file <'+fn+'> generated'

    # redraw
    Blender.Window.QRedrawAll()

    # restore cursor
    Blender.Window.WaitCursor(0)


def gen_script():
    obj = di.get_obj()

    # set cursor
    Blender.Window.WaitCursor(1)

    # check for properties
    linele = False
    try:    is_3d  = obj.properties['is3d']
    except: linele = True # assuming linear structures (trusses, beams)

    # boundary conditions
    nbrys, nbrysID, ebrys, fbrys, eatts = get_brys_atts (obj)

    # text
    fn  = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
    txt = Blender.Text.New(obj.name+'_script')
    txt.write ('import Blender\n')
    txt.write ('import bpy\n')
    txt.write ('import mechsys\n')
    txt.write ('import msys_fem as mf\n')
    txt.write ('\n# Show running cursor\n')
    txt.write ('Blender.Window.WaitCursor(1)\n')
    txt.write ('\n# Boundary conditions & properties\n')
    txt.write ('nbrys = '+nbrys.__str__()+'\n')
    txt.write ('ebrys = '+ebrys.__str__()+'\n')
    txt.write ('fbrys = '+fbrys.__str__()+'\n')
    txt.write ('eatts = '+eatts.__str__()+'\n')
    txt.write ('\n# Problem geometry\n')
    txt.write ('obj = bpy.data.objects["'+obj.name+'"]\n')
    if linele: txt.write ('geo = mf.set_geo_linele(obj,nbrys,eatts)\n')
    else:      txt.write ('geo = mf.set_geo(obj,nbrys,ebrys,fbrys,eatts)\n')
    txt.write ('\n# nodes boundary conditions\n')
    for nb in nbrysID:
        txt.write('geo.nod('+str(nb[0])+').bry("'+nb[1]+'",'+str(nb[2])+')\n')
    txt.write ('\n# Solution\n')
    txt.write ('sol = mechsys.solver("ForwardEuler")\n')
    txt.write ('sol.set_geom(geo)\n')
    #txt.write ('sol.set_lin_sol("LA").set_num_div(1).set_delta_time(0.0)\n')
    txt.write ('sol.solve()\n')
    txt.write ('\n# Store results in object\n')
    txt.write ('mf.store_scalars(geo,obj)\n')
    txt.write ('\n# Output\n')
    txt.write ('mechsys.out_vtu(geo, "'+fn+'")\n')
    txt.write ('\n# Hide running cursor\n')
    txt.write ('Blender.Window.WaitCursor(0)\n')

    # restore cursor
    Blender.Window.WaitCursor(0)
    Blender.Window.Redraw(Blender.Window.Types.TEXT)


def paraview(obj):
    fn = Blender.sys.makename (ext='_FEA_'+obj.name+'.vtu')
    if Blender.sys.exists(fn):
        Blender.Window.WaitCursor(1)
        try: pid = subprocess.Popen(['paraview', '--data='+fn]).pid
        except:
            Blender.Window.WaitCursor(0)
            raise Exception('Paraview is not available, please install it first')
        Blender.Window.WaitCursor(0)
    else: raise Exception('File <'+fn+'> does not exist (please, run analysis first)')
