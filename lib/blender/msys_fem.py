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
import Blender
from   Blender.Mathutils import Vector
import bpy
import mechsys   as ms
import msys_dict as di
import msys_mesh as me


def save_stage_mats_info():
    obj = di.get_obj()
    xyz = obj.getLocation('worldspace')
    t3d = Blender.Text3d.New(obj.name+'_saved')
    t3d.setText("SAVED STAGES INFORMATION FROM %s"%obj.name)
    scn = bpy.data.scenes.active
    new_obj = scn.objects.new (t3d)
    new_obj.makeDisplayList   () # rebuild the display list for this object
    new_obj.setLocation       (xyz)
    if obj.properties.has_key('stages'):
        new_obj.properties['stages'] = obj.properties['stages']
        new_obj.properties['texts']  = obj.properties['texts']
        new_obj.properties['mats']   = obj.properties['mats']
        for k, v in obj.properties.iteritems():
            if k[:4]=='stg_':
                new_obj.properties[k] = obj.properties[k]
    Blender.Window.QRedrawAll()


def read_stage_mats_info():
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    if not len(obs)==2: raise Exception('Please, select one Text3D object and one Mesh object first')
    if obs[0].type=='Mesh' and obs[1].type=='Text':
        src = obs[1]
        des = obs[0]
    elif obs[0].type=='Text' and obs[1].type=='Mesh':
        src = obs[0]
        des = obs[1]
    else: raise Exception('Please, select one Text3D object and one Mesh object first')
    msg = 'Confirm copying stages/materials information from ('+src.name+') to ('+des.name+') ?%t|Yes'
    res = Blender.Draw.PupMenu(msg)
    if res>0:
        if src.properties.has_key('stages'):
            des.properties['stages'] = src.properties['stages']
            des.properties['texts']  = src.properties['texts']
            des.properties['mats']   = src.properties['mats']
            for k, v in src.properties.iteritems():
                if k[:4]=='stg_':
                    des.properties[k] = src.properties[k]
    Blender.Window.QRedrawAll()


def get_mats(obj):
    d = di.load_dict()
    mats = {}
    if obj.properties.has_key('mats'):
        for k, v in obj.properties['mats'].iteritems():
            desc = obj.properties['texts'][str(int(v[11]))]
            if int(v[0])==0: # LinElastic
                mats[int(k)] = [d['mdl'][int(v[0])], 'E=%g nu=%g' % (v[1],v[2]), desc]
            elif int(v[0])==1: # LinDiffusion
                mats[int(k)] = [d['mdl'][int(v[0])], 'k=%g' % (v[3]), desc]
            elif int(v[0])==2: # CamClay
                mats[int(k)] = [d['mdl'][int(v[0])], 'lam=%g kap=%g phics=%g G=%g v=%g' % (v[4],v[5],v[6],v[7],v[8]), desc]
            elif int(v[0])==3: # BeamElastic
                mats[int(k)] = [d['mdl'][int(v[0])], 'E=%g A=%g Izz=%g' % (v[1],v[9],v[10]), desc]
    return mats


def get_eatts(obj,mats,stg):
    d     = di.load_dict()
    eatts = []
    if obj.properties[stg].has_key('eatts'):
        for k, v in obj.properties[stg]['eatts'].iteritems():
            tag   = int(v[0])
            inis  = 'ZERO'
            matID = int(v[2])
            prop  = obj.properties['texts'][str(int(v[3]))]
            act   = True if int(v[4]) else False
            if mats.has_key(matID):
                mdl     = mats[matID][0]
                prms    = mats[matID][1]
                matdesc = mats[matID][2]
            else:
                mdl     = 'LinElastic'
                prms    = 'E=200 nu=0.2'
                matdesc = '__no material__'
            eatts.append ([tag, d['ety'][int(v[1])], mdl, prms, inis, prop, act, matdesc])
    return eatts


def get_act_deact(obj,stg):
    d = di.load_dict()
    activate   = {}
    deactivate = {}
    if obj.properties[stg].has_key('eatts'):
        for k, v in obj.properties[stg]['eatts'].iteritems():
            tag             = int(v[0])
            activate  [tag] = int(v[5])
            deactivate[tag] = int(v[6])
    return activate, deactivate


def get_brys(obj,stg):
    d = di.load_dict()

    # nbrys
    nbrys = []
    if obj.properties[stg].has_key('nbrys'):
        for k, v in obj.properties[stg]['nbrys'].iteritems():
            nbrys.append([v[0], v[1], v[2], d['dfv'][int(v[3])], v[4]])

    # nbsID
    nbsID = []
    if obj.properties[stg].has_key('nbsID'):
        for k, v in obj.properties[stg]['nbsID'].iteritems():
            nbsID.append([int(v[0]), d['dfv'][int(v[1])], v[2]])

    # ebrys
    ebrys = []
    if obj.properties[stg].has_key('ebrys'):
        for k, v in obj.properties[stg]['ebrys'].iteritems():
            ebrys.append([int(v[0]), d['dfv'][int(v[1])], v[2]])

    # fbrys
    fbrys = []
    if obj.properties[stg].has_key('fbrys'):
        for k, v in obj.properties[stg]['fbrys'].iteritems():
            fbrys.append([int(v[0]), d['dfv'][int(v[1])], v[2]])

    return nbrys, nbsID, ebrys, fbrys


def get_mesh(obj, txt=None):
    # mesh
    edm = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0)
    msh = obj.getData(mesh=1)

    # check mesh
    if not obj.properties.has_key('nelems'): raise Exception('Please, generate mesh first')

    # transform mesh to global coordinates
    ori = msh.verts[:] # create a copy in local coordinates
    msh.transform (obj.matrix)

    # check for frame mesh
    frame = False
    if obj.properties.has_key('frame'): frame = obj.properties['frame']

    # check ndim
    is3d = obj.properties['3dmesh']
    ndim = 3 if is3d else 2

    if txt!=None: # generate script
        # mesh
        txt.write('\n# MechSys mesh structure\n')
        if is3d: txt.write('mesh = ms.mesh_generic(True)\n')
        else:    txt.write('mesh = ms.mesh_generic(False)\n')

        # vertices
        txt.write('\n# vertices\n')
        txt.write('mesh.set_nverts (%d)\n' % len(msh.verts))
        for i, v in enumerate(msh.verts):
            onbry = True if (i in obj.properties['verts_bry']) else False
            if is3d: txt.write('mesh.set_vert (%d,%d,%f,%f,%f)\n' % (i, onbry, v.co[0], v.co[1], v.co[2]))
            else:    txt.write('mesh.set_vert (%d,%d,%f,%f)\n'    % (i, onbry, v.co[0], v.co[1]))

        # elements
        txt.write('\n# elements\n')
        nelems = obj.properties['nelems']
        txt.write('mesh.set_nelems (%d)\n' % nelems)
        for i in range(nelems):
            txt.write('mesh.set_elem (%d,%d,%d,%d)\n' %(i, obj.properties['elems']['tags'][i], obj.properties['elems']['onbs'][i], obj.properties['elems']['vtks'][i]))

        # connectivities
        txt.write('\n# connectivities\n')
        for i in range(nelems):
            for j in range(len(obj.properties['elems']['cons'] [str(i)])):
                txt.write('mesh.set_elem_con (%d,%d,%d)\n' % (i, j, obj.properties['elems']['cons'] [str(i)][j]))

        # edge tags
        txt.write('\n# edge tags\n')
        for i in range(nelems):
            for j in range(len(obj.properties['elems']['etags'][str(i)])):
                tag = obj.properties['elems']['etags'][str(i)][j]
                if tag<0: txt.write('mesh.set_elem_etag (%d,%d,%d)\n' % (i, j, tag))

        # face tags
        if is3d:
            txt.write('\n# face tags\n')
            for i in range(nelems):
                for j in range(len(obj.properties['elems']['ftags'][str(i)])):
                    txt.write('mesh.set_elem_ftag (%d,%d,%d)\n' % (i, j, obj.properties['elems']['ftags'][str(i)][j]))

    else:
        # mesh structure
        mesh = ms.mesh_generic(is3d)

        # vertices
        mesh.set_nverts (len(msh.verts))
        for i, v in enumerate(msh.verts):
            onbry = True if (i in obj.properties['verts_bry']) else False
            if is3d: mesh.set_vert (i, onbry, v.co[0], v.co[1], v.co[2])
            else:    mesh.set_vert (i, onbry, v.co[0], v.co[1])

        # elements
        nelems = obj.properties['nelems']
        mesh.set_nelems (nelems)
        for i in range(nelems):
            # element
            mesh.set_elem (i, obj.properties['elems']['tags'][i],
                            obj.properties['elems']['onbs'][i],
                            obj.properties['elems']['vtks'][i])

            # connectivities
            for j in range(len(obj.properties['elems']['cons'] [str(i)])):
                mesh.set_elem_con (i, j, obj.properties['elems']['cons'] [str(i)][j])

            # edge tags
            for j in range(len(obj.properties['elems']['etags'][str(i)])):
                mesh.set_elem_etag (i, j, obj.properties['elems']['etags'][str(i)][j])

            # face tags
            if is3d:
                for j in range(len(obj.properties['elems']['ftags'][str(i)])):
                    mesh.set_elem_ftag (i, j, obj.properties['elems']['ftags'][str(i)][j])

    # end
    msh.verts = ori # restore local coordinates
    if edm: Blender.Window.EditMode(1)
    if txt==None: return mesh


def run_analysis(gen_script=False):
    # object
    obj = di.get_obj()
    Blender.Window.WaitCursor(1)

    # check number of stages
    if not obj.properties.has_key('stages'): raise Exception('Please, add stages first')
    nstages = len(obj.properties['stages'])

    # first stage
    for k, v in obj.properties['stages'].iteritems():
        if v[0]==1: stg = 'stg_'+k

    # materials and first eatts
    mats   = get_mats  (obj)
    eatts1 = get_eatts (obj,mats,stg)

    # check for frame mesh
    frame = False
    if obj.properties.has_key('frame'): frame = obj.properties['frame']

    # ndim
    is3d = obj.properties['3dmesh']
    ndim = 3 if is3d else 2

    if gen_script:
        # create new script
        txt = Blender.Text.New(obj.name+'_fem')

        # import libraries
        if not di.key('fullsc'):
            txt.write ('import Blender, bpy\n')
            txt.write ('import msys_fem as mf\n')
        txt.write     ('import mechsys  as ms\n')

        # change cursor
        if not di.key('fullsc'):
            txt.write ('\n# Show running cursor\n')
            txt.write ('Blender.Window.WaitCursor(1)\n')

        # get/set mesh
        if not di.key('fullsc'):
            txt.write ('\n# Get mesh\n')
            txt.write ('obj  = bpy.data.objects["'+obj.name+'"]\n')
            txt.write ('mesh = mf.get_mesh(obj)\n')
        else:
            txt.write ('\n# Set mesh\n')
            get_mesh  (obj, txt)

        # allocate geometry
        txt.write ('\n# Geometry\n')
        txt.write ('geo = ms.geom (%d)\n' % ndim)

        # element attributes
        neatt = len(eatts1)
        txt.write ('\n# Element attributes\n')
        for i, ea in enumerate(eatts1):
            if i==0:       estr  = 'eatts = ['
            else:          estr  = '         '
            if ea[6]:      estr += '[%d, "%s", "%s", "%s", "%s", "%s", True ]'#,  # %s\n'
            else:          estr += '[%d, "%s", "%s", "%s", "%s", "%s", False]'#,  # %s\n'
            if i==neatt-1: estr += ']  # %s\n'
            else:          estr += ',  # %s\n'
            txt.write (estr % (ea[0],ea[1],ea[2],ea[3],ea[4],ea[5],ea[7]))

        # set nodes and elements
        txt.write           ('\n# Set nodes and elements\n')
        if frame: txt.write ('ms.set_nodes_elems (mesh, eatts, geo, 1.0e-5, True)\n')
        else:     txt.write ('ms.set_nodes_elems (mesh, eatts, geo)\n')

        # allocate solver
        txt.write ('\n# Solver\n')
        txt.write ('sol = ms.solver ("ForwardEuler")\n')
        txt.write ('sol.set_geom    (geo)\n')

        # allocate output and open collection
        txt.write ('\n# Open collection for output\n')
        txt.write ('out = ms.output     ()\n')
        txt.write ('out.open_collection ("'+obj.name+'")\n')

        # solve each stage
        for num in range(1,nstages+1):

            # find stage info
            for k, v in obj.properties['stages'].iteritems():
                if int(v[0])==num:
                    stg   = 'stg_'+k
                    desc  = obj.properties['texts'][str(int(v[1]))]
                    abf   = True if int(v[2]) else False # apply body forces ?
                    cdi   = True if int(v[3]) else False # clear displacements ?
                    ndiv  = int(v[4])
                    dtime = v[5]
                    break

            # activate and deactivate elements
            txt.write ('\n# Stage # %d --------------------------------------------------------------\n'%num)
            act, deact = get_act_deact (obj,stg)
            for k, v in act.iteritems():
                if v: txt.write ('geo.activate (%d)\n'%(k))
            for k, v in deact.iteritems():
                if v: txt.write ('geo.deactivate (%d)\n'%(k))

            # boundary conditions
            nbrys, nbsID, ebrys, fbrys = get_brys  (obj,stg)
            txt.write ('nbrys = '+nbrys.__str__()+'\n')
            txt.write ('ebrys = '+ebrys.__str__()+'\n')
            txt.write ('fbrys = '+fbrys.__str__()+'\n')
            txt.write ('ms.set_brys             (mesh, nbrys, ebrys, fbrys, geo)\n')
            for nb in nbsID:
                txt.write ('geo.nod('+str(nb[0])+').bry("'+nb[1]+'",'+str(nb[2])+')\n')

            # apply body forces
            if abf: txt.write ('geo.apply_body_forces   ()\n')

            # solve
            txt.write ('sol.solve_with_info     (%d, %g, %d, "%s\\n")\n'%(ndiv,dtime,num,desc))

            # clear displacements
            if cdi: txt.write ('geo.clear_displacements ()\n')

            # output
            txt.write ('out.vtu                 (geo, sol.time())\n')

            # save results
            if not di.key('fullsc'):
                txt.write ('mf.save_results         (geo,obj, %d)\n'%num)

        # close collection
        txt.write ('\n# Close collection\n')
        txt.write ('out.close_collection()\n')

        # change cursor
        if not di.key('fullsc'):
            txt.write ('\n# Hide running cursor\n')
            txt.write ('Blender.Window.WaitCursor(0)\n')

    else:
        # get/set mesh
        mesh = get_mesh(obj)

        # allocate geometry
        geo = ms.geom (ndim)

        # element attributes
        eatts = []
        for ea in eatts1: eatts.append(ea[:7])

        # set nodes and elements
        if frame: ms.set_nodes_elems (mesh, eatts, geo, 1.0e-5, True)
        else:     ms.set_nodes_elems (mesh, eatts, geo)

        # allocate solver
        sol = ms.solver ("ForwardEuler")
        sol.set_geom    (geo)

        # allocate output and open collection
        out = ms.output     ()
        out.open_collection (obj.name)

        # solve each stage
        for num in range(1,nstages+1):

            # find stage info
            for k, v in obj.properties['stages'].iteritems():
                if int(v[0])==num:
                    stg   = 'stg_'+k
                    desc  = obj.properties['texts'][str(int(v[1]))]
                    abf   = True if int(v[2]) else False # apply body forces ?
                    cdi   = True if int(v[3]) else False # clear displacements ?
                    ndiv  = int(v[4])
                    dtime = v[5]
                    break

            # activate and deactivate elements
            act, deact = get_act_deact (obj,stg)
            for k, v in act.iteritems():
                if v: geo.activate(k)
            for k, v in deact.iteritems():
                if v: geo.deactivate(k)

            # boundary conditions
            nbrys, nbsID, ebrys, fbrys = get_brys  (obj,stg)
            ms.set_brys (mesh, nbrys, ebrys, fbrys, geo)
            for nb in nbsID: geo.nod(nb[0]).bry(nb[1],nb[2])

            # apply body forces
            if abf: geo.apply_body_forces()

            # solve
            sol.solve_with_info (ndiv,dtime,num,desc+'\n')

            # clear displacements
            if cdi: geo.clear_displacements()

            # output
            out.vtu (geo, sol.time())

            # save results
            save_results (geo,obj, num)

        # close collection
        out.close_collection()

    # end
    #Blender.Run(txt.name)
    Blender.Window.WaitCursor(0)


def save_results(geo,obj, stage_num):
    # dictionary
    if not obj.properties.has_key('res'): obj.properties['res'] = {}
    s                        = str(stage_num)
    obj.properties['res'][s] = {}

    # geometry limits
    sp, ep = [], []
    geo.bounds_3d (sp, ep)
    v1    = Vector(sp)
    v2    = Vector(ep)
    delta = v2-v1
    obj.properties['res'][s]['length'] = delta.length

    # save results at nodes
    for k, key in di.key('dfv').iteritems():
        vals = []
        for i in range(geo.nnodes()):
            try:    vals.append(geo.nod(i).val(key))
            except: vals.append(0.0)
        obj.properties['res'][s][key] = vals

    # save extra output
    obj.properties['res'][s]['extra'] = {}
    for i in range(geo.nelems()):
        if geo.ele(i).has_extra() and geo.ele(i).is_active():
            ide = str(i)
            obj.properties['res'][s]['extra'][ide] = {}
            geo.ele(i).calc_dep_vars()
            co, no, va = dict(), [], dict()
            geo.ele(i).out_extra (co, no, va)
            if co.has_key('X'):
                if len(co['X'])>0:
                    obj.properties['res'][s]['extra'][ide]['coords'] = co
                    obj.properties['res'][s]['extra'][ide]['normal'] = no
                    obj.properties['res'][s]['extra'][ide]['values'] = {}
                    for k, v in va.iteritems():
                        obj.properties['res'][s]['extra'][ide]['values'][k] = v
                        # max value
                        key  = 'max_'+k
                        maxv = max([abs(val) for val in v])
                        if obj.properties['res'][s].has_key(key):
                            if maxv>obj.properties['res'][s][key][0]: obj.properties['res'][s][key] = [maxv, i, geo.ele(i).nod(0).x(), geo.ele(i).nod(0).y()]
                        else: obj.properties['res'][s][key] = [maxv, i, geo.ele(i).nod(0).x(), geo.ele(i).nod(0).y()]

    Blender.Window.QRedrawAll()


def paraview():
    obj = di.get_obj()
    fn  = obj.name+'.pvd'
    if Blender.sys.exists(fn):
        Blender.Window.WaitCursor(1)
        try: pid = subprocess.Popen(['paraview', '--data='+fn]).pid
        except:
            Blender.Window.WaitCursor(0)
            raise Exception('Paraview is not available, please install it first')
        Blender.Window.WaitCursor(0)
    else: raise Exception('File <'+fn+'> does not exist (please, run analysis first)')
