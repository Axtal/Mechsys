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


def get_brys_atts(obj):
    d = di.load_dict()

    # nbrys
    nbrys = []
    if obj.properties.has_key('nbrys'):
        for k, v in obj.properties['nbrys'].iteritems():
            nbrys.append([v[0], v[1], v[2], d['dfv'][int(v[3])], v[4]])

    # nbsID
    nbsID = []
    if obj.properties.has_key('nbsID'):
        for k, v in obj.properties['nbsID'].iteritems():
            nbsID.append([int(v[0]), d['dfv'][int(v[1])], v[2]])

    # ebrys
    ebrys = []
    if obj.properties.has_key('ebrys'):
        for k, v in obj.properties['ebrys'].iteritems():
            ebrys.append([int(v[0]), d['dfv'][int(v[1])], v[2]])

    # fbrys
    fbrys = []
    if obj.properties.has_key('fbrys'):
        for k, v in obj.properties['fbrys'].iteritems():
            fbrys.append([int(v[0]), d['dfv'][int(v[1])], v[2]])

    # mats
    mats = {}
    if obj.properties.has_key('mats'):
        for k, v in obj.properties['mats'].iteritems():
            if int(v[0])==0: # LinElastic
                mats[int(k)] = [d['mdl'][int(v[0])], 'E=%f nu=%f' % (v[1],v[2])]
            elif int(v[0])==1: # LinDiffusion
                mats[int(k)] = [d['mdl'][int(v[0])], 'k=%f' % (v[3])]
            elif int(v[0])==2: # CamClay
                mats[int(k)] = [d['mdl'][int(v[0])], 'lam=%f kap=%f phics=%f G=%f v=%f' % (v[4],v[5],v[6],v[7],v[8])]
            elif int(v[0])==3: # BeamElastic
                mats[int(k)] = [d['mdl'][int(v[0])], 'E=%f A=%f Izz=%f' % (v[1],v[9],v[10])]

    # eatts
    eatts = []
    if obj.properties.has_key('eatts'):
        for k, v in obj.properties['eatts'].iteritems():
            inis  = 'ZERO'
            matID = int(v[2])
            if mats.has_key(matID):
                mdl  = mats[matID][0]
                prms = mats[matID][1]
            else:
                mdl  = 'LinElastic'
                prms = 'E=200 nu=0.2'
            eatts.append([int(v[0]), d['ety'][int(v[1])], mdl, prms, inis, '', True])

    return nbrys, nbsID, ebrys, fbrys, eatts


def set_geo(obj,nbrys,ebrys,fbrys,eatts, gen_script=False,txt=None):
    # mesh
    edm = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0)
    msh = obj.getData(mesh=1)

    # transform mesh to global coordinates
    ori = msh.verts[:] # create a copy in local coordinates
    msh.transform (obj.matrix)

    # check for frame mesh
    frame = False
    if obj.properties.has_key('frame'): frame = obj.properties['frame']

    # set geometry
    is3d = obj.properties['3dmesh']
    ndim = 3 if is3d else 2
    if gen_script: txt.write('geo = ms.geom (%d)\n' % ndim)
    else: geo = ms.geom (ndim)

    if gen_script:
        txt.write('\n# mesh structure\n')
        if is3d: txt.write('mg = ms.mesh_generic(True)\n')
        else:    txt.write('mg = ms.mesh_generic(False)\n')

        txt.write('\n# vertices\n')
        txt.write('mg.set_nverts (%d)\n' % len(msh.verts))
        for i, v in enumerate(msh.verts):
            onbry = True if (i in obj.properties['verts_bry']) else False
            if is3d: txt.write('mg.set_vert (%d,%d,%f,%f,%f)\n' % (i, onbry, v.co[0], v.co[1], v.co[2]))
            else:    txt.write('mg.set_vert (%d,%d,%f,%f)\n'    % (i, onbry, v.co[0], v.co[1]))

        txt.write('\n# elements\n')
        nelems = obj.properties['nelems']
        txt.write('mg.set_nelems (%d)\n' % nelems)
        for i in range(nelems):
            txt.write('mg.set_elem (%d,%d,%d,%d)\n' %(i, obj.properties['elems']['tags'][i], obj.properties['elems']['onbs'][i], obj.properties['elems']['vtks'][i]))

        txt.write('\n# connectivities\n')
        for i in range(nelems):
            for j in range(len(obj.properties['elems']['cons'] [str(i)])):
                txt.write('mg.set_elem_con (%d,%d,%d)\n' % (i, j, obj.properties['elems']['cons'] [str(i)][j]))

        txt.write('\n# edge tags\n')
        for i in range(nelems):
            for j in range(len(obj.properties['elems']['etags'][str(i)])):
                tag = obj.properties['elems']['etags'][str(i)][j]
                if tag<0: txt.write('mg.set_elem_etag (%d,%d,%d)\n' % (i, j, tag))

        if is3d:
            txt.write('\n# face tags\n')
            for i in range(nelems):
                for j in range(len(obj.properties['elems']['ftags'][str(i)])):
                    txt.write('mg.set_elem_ftag (%d,%d,%d)\n' % (i, j, obj.properties['elems']['ftags'][str(i)][j]))

        txt.write('\n# set geometry\n')
        if frame: txt.write('ms.set_nodes_elems (mg, eatts, geo, 1.0e-5, True)\n')
        else:     txt.write('ms.set_nodes_elems (mg, eatts, geo, 1.0e-5, False)\n')
        txt.write('ms.set_brys        (mg, nbrys, ebrys, fbrys, geo)\n')

    else:
        # mesh structure
        mg = ms.mesh_generic(is3d)

        # vertices
        mg.set_nverts (len(msh.verts))
        for i, v in enumerate(msh.verts):
            onbry = True if (i in obj.properties['verts_bry']) else False
            if is3d: mg.set_vert (i, onbry, v.co[0], v.co[1], v.co[2])
            else:    mg.set_vert (i, onbry, v.co[0], v.co[1])

        # elements
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
            if is3d:
                for j in range(len(obj.properties['elems']['ftags'][str(i)])):
                    mg.set_elem_ftag (i, j, obj.properties['elems']['ftags'][str(i)][j])

        # set geometry
        ms.set_nodes_elems (mg, eatts, geo, 1.0e-5, frame)
        ms.set_brys        (mg, nbrys, ebrys, fbrys, geo)

    # end
    msh.verts = ori # restore local coordinates
    if edm: Blender.Window.EditMode(1)
    if not gen_script: return geo


def save_results(geo,obj):
    # check what variables are available
    obj.properties['res']        = {}
    obj.properties['res']['l2g'] = {}             # dof vars map: local to global
    dfvmnu                       = 'DOF Vars %t|' # dof vars menu
    i                            = 0
    for k, v in di.key('dfv').iteritems():
        try:
            val = geo.nod(0).val(v)
            obj.properties['res']['l2g'][str(i)] = k
            dfvmnu += v+' %x'+str(i+1)+'|'
            i += 1
        except: pass
    obj.properties['res']['dfvmnu'] = dfvmnu

    # save values in object
    for k, v in obj.properties['res']['l2g'].iteritems():
        key  = di.key('dfv')[v]
        vals = []
        for i in range(geo.nnodes()):
            try:    vals.append(geo.nod(i).val(key))
            except: vals.append(0.0)
        obj.properties['res'][key] = vals

    # limits
    sp, ep = [], []
    geo.bounds_3d (sp, ep)
    v1    = Vector(sp)
    v2    = Vector(ep)
    delta = v2-v1
    obj.properties['res']['length'] = delta.length

    # save extra output
    obj.properties['res']['extra'] = {}
    for i in range(geo.nelems()):
        if geo.ele(i).has_extra():
            ide = str(i)
            obj.properties['res']['extra'][ide] = {}
            geo.ele(i).calc_dep_vars()
            co, no, va = dict(), [], dict()
            geo.ele(i).out_extra (co, no, va)
            if co.has_key('X'):
                if len(co['X'])>0:
                    obj.properties['res']['extra'][ide]['coords'] = co
                    obj.properties['res']['extra'][ide]['normal'] = no
                    obj.properties['res']['extra'][ide]['values'] = {}
                    for k, v in va.iteritems():
                        obj.properties['res']['extra'][ide]['values'][k] = v
                        # max value
                        key  = 'max_'+k
                        maxv = max([abs(val) for val in v])
                        if obj.properties['res'].has_key(key):
                            if maxv>obj.properties['res'][key]: obj.properties['res'][key] = maxv
                        else: obj.properties['res'][key] = maxv

def run_analysis(gen_script=False):
    # get object and set cursor
    obj = di.get_obj()
    Blender.Window.WaitCursor(1)

    # boundary conditions & properties
    nbrys, nbsID, ebrys, fbrys, eatts = get_brys_atts (obj)

    if gen_script:
        txt = Blender.Text.New(obj.name+'_fem')
        txt.write ('import Blender, bpy\n')
        txt.write ('import mechsys  as ms\n')
        txt.write ('import msys_fem as mf\n')
        txt.write ('\n# Show running cursor\n')
        txt.write ('Blender.Window.WaitCursor(1)\n')
        txt.write ('\n# Boundary conditions & properties\n')
        txt.write ('nbrys = '+nbrys.__str__()+'\n')
        txt.write ('ebrys = '+ebrys.__str__()+'\n')
        txt.write ('fbrys = '+fbrys.__str__()+'\n')
        txt.write ('eatts = '+eatts.__str__()+'\n')
        txt.write ('\n# FEM data\n')
        if di.key('fullsc'):
            geo = set_geo (obj,nbrys,ebrys,fbrys,eatts, True,txt)
        else:
            txt.write ('obj = bpy.data.objects["'+obj.name+'"]\n')
            txt.write ('geo = mf.set_geo(obj,nbrys,ebrys,fbrys,eatts)\n')
        if len(nbsID)>0: txt.write ('\n# nodes boundary conditions\n')
        for nb in nbsID:
            txt.write('geo.nod('+str(nb[0])+').bry("'+nb[1]+'",'+str(nb[2])+')\n')
        txt.write ('\n# Solution\n')
        txt.write ('sol = ms.solver     ("ForwardEuler")\n')
        txt.write ('sol.set_geom        (geo)\n')
        txt.write ('sol.solve_with_info ()\n')
        if not di.key('fullsc'):
            txt.write ('\n# Save results in object\n')
            txt.write ('mf.save_results(geo,obj)\n')
        txt.write ('\n# Output\n')
        txt.write ('out = ms.output()\n')
        txt.write ('out.vtu(geo, \''+obj.name+'_FEM.vtu\')\n')
        txt.write ('\n# Hide running cursor\n')
        txt.write ('Blender.Window.WaitCursor(0)\n')
    else:
        if (len(eatts)<1):
            if edm: Blender.Window.EditMode(1)
            raise Exception('Element attributes must be defined before running the simulation')

        # FEM data
        geo = set_geo (obj,nbrys,ebrys,fbrys,eatts)

        # nodes boundary conditions
        for nb in nbsID:
            geo.nod (nb[0]).bry(nb[1],nb[2])

        # solution
        sol = ms.solver     ('ForwardEuler')
        sol.set_geom        (geo)
        sol.solve_with_info ()

        # save results in object
        save_results (geo,obj)

        # output
        fn = obj.name+'_FEM.vtu'
        out = ms.output()
        out.vtu (geo, fn)
        print '[1;34mMechSys[0m: file <'+fn+'> generated'

    # redraw and restore cursor
    Blender.Window.WaitCursor(0)
    Blender.Window.QRedrawAll()


def paraview():
    obj = di.get_obj()
    fn  = obj.name+'_FEM.vtu'
    if Blender.sys.exists(fn):
        Blender.Window.WaitCursor(1)
        try: pid = subprocess.Popen(['paraview', '--data='+fn]).pid
        except:
            Blender.Window.WaitCursor(0)
            raise Exception('Paraview is not available, please install it first')
        Blender.Window.WaitCursor(0)
    else: raise Exception('File <'+fn+'> does not exist (please, run analysis first)')
