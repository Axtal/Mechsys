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
            elif int(v[0])==4: # BiotElastic
                mats[int(k)] = [d['mdl'][int(v[0])], 'E=%g nu=%g k=%g gw=%g' % (v[1],v[2],v[3],v[12]), desc]
            elif int(v[0])==5: # Reinforcement
                mats[int(k)] = [d['mdl'][int(v[0])], 'E=%g Ar=%g At=%g ks=%g c=%g phi=%g' % (v[1],v[13],v[14],v[15],v[16],v[17]), desc]
            elif int(v[0])==6: # SpringElastic
                mats[int(k)] = [d['mdl'][int(v[0])], 'ks=%g' % (v[15]), desc]
            elif int(v[0])==7: # RodElastic
                mats[int(k)] = [d['mdl'][int(v[0])], 'E=%g A=%g' % (v[1],v[9]), desc]
    return mats


def get_eatts(obj,mats,stg):
    d     = di.load_dict()
    eatts = []
    if obj.properties[stg].has_key('eatts'):
        for k, v in obj.properties[stg]['eatts'].iteritems():
            tag   = int(v[0])
            gty   = d['gty'][int(v[1])]
            pty   = d['pty'][int(v[2])]
            inis  = 'ZERO'
            matID = int(v[3])
            prop  = obj.properties['texts'][str(int(v[7]))]
            act   = True if int(v[4]) else False
            if mats.has_key(matID):
                mdl     = mats[matID][0]
                prms    = mats[matID][1]
                matdesc = mats[matID][2]
            else:
                mdl     = ''
                prms    = ''
                matdesc = '__no material__'
            eatts.append ([tag, gty, pty, mdl, prms, inis, prop, act, matdesc])
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


def get_reinforcements(obj,is3d):
    reinfs = {}
    if obj.properties.has_key('reinfs'):
        for k, v in obj.properties['reinfs'].iteritems():
            tag = int(v[0])
            if is3d: reinfs[(v[1],v[2],v[3], v[4],v[5],v[6])] = tag
            else:    reinfs[(v[1],v[2],      v[4],v[5]     )] = tag
    return reinfs


def get_linear_elems(obj):
    lines= {}
    if obj.properties.has_key('lines'):
        for k, v in obj.properties['lines'].iteritems():
            lines[(v[1],v[2])] = v[0]
    return lines


def run_analysis(gen_script=False):
    Blender.Window.WaitCursor(1)

    # get active object
    obj = di.get_obj()

    # check mesh type
    if obj.properties.has_key('mesh_type'): mesh_type = obj.properties['mesh_type']
    else: raise Exception('Please, generate mesh first or set "Frame mesh" toggle')

    # check number of stages
    if not obj.properties.has_key('stages'): raise Exception('Please, add stages first')
    nstages = len(obj.properties['stages'])

    # first stage
    sid, stg = di.find_stage (obj, 1)

    # materials and first eatts
    mats  = get_mats  (obj)
    eatts = get_eatts (obj, mats, stg)
    if len(eatts)<1: raise Exception('Please, define element attributes first')

    # ndim
    is3d = obj.properties['is3d'] if obj.properties.has_key('is3d') else False
    ndim = 3 if is3d else 2
    print '[1;34mMechSys[0m: Running [1;31m%dD[0m simulation'%ndim

    # reinforcements
    reinfs    = get_reinforcements (obj, is3d)
    has_reinf = len(reinfs)>0

    # linear elements
    lines     = get_linear_elems (obj)
    has_lines = len(lines)>0

    if gen_script:

        # create new script
        txt = Blender.Text.New(obj.name+'_fem')

        if di.key('fem_cpp'): # C++ script
            # includes
            txt.write ('// Std Lib\n')
            txt.write ('#include <iostream>\n\n')
            txt.write ('// MechSys\n')
            txt.write ('#include "mechsys.h"\n')
            txt.write ('#include "util/exception.h"\n\n')
            txt.write ('#define T boost::make_tuple\n\n')

            # main
            txt.write ('int main(int argc, char **argv) try\n')
            txt.write ('{\n\n')

            # mesh
            txt.write ('	////////////////////////////////////////////////////////////////////////////////////// Mesh /////\n\n')

            if   mesh_type=='struct':   me.gen_struct_mesh   (True, txt, False, True)
            elif mesh_type=='unstruct': me.gen_unstruct_mesh (True, txt, False, True)
            elif mesh_type=='frame':    me.gen_frame_mesh    (      txt, False, True)
            txt.write ('\n')

            txt.write ('	////////////////////////////////////////////////////////////////////////////////////// FEM //////\n\n')

            # data and solver
            txt.write ('	// Data and Solver\n')
            txt.write ('	FEM::Data   dat (%d);\n'        % ndim)
            txt.write ('	FEM::Solver sol (dat, "%s");\n' % obj.name)

            # element attributes
            neatt = len(eatts)
            txt.write ('\n	// Element attributes\n')
            txt.write ('	FEM::EAtts_T eatts(%d);\n'%len(eatts))
            for i, ea in enumerate(eatts):
                if i==0:       estr  = '	eatts = '
                else:          estr  = '	        '
                if ea[7]:      estr += 'T(%d, "%s", "%s", "%s", "%s", "%s", "%s", true )'#,  # %s\n'
                else:          estr += 'T(%d, "%s", "%s", "%s", "%s", "%s", "%s", false)'#,  # %s\n'
                if i==neatt-1: estr += ';  // %s\n'
                else:          estr += ',  // %s\n'
                txt.write (estr % (ea[0],ea[1],ea[2],ea[3],ea[4],ea[5],ea[6],ea[8]))

            # set geometry: nodes and elements
            txt.write ('\n	// Set nodes and elements (geometry)\n')
            if mesh_type=='frame': txt.write ('	dat.SetOnlyFrame  (); // frame (beam/truss) mesh only\n')
            txt.write ('	dat.SetNodesElems (&mesh, &eatts);\n')

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
                        dtime =     v[5]
                        act   = int(v[6])
                        break

                # run only if stage is active
                if act:

                    # activate and deactivate elements
                    txt.write ('\n	// Stage # %d --------------------------------------------------------------\n'%num)
                    elem_act, elem_deact = get_act_deact (obj,stg)
                    for k, v in elem_act.iteritems():
                        if v: txt.write ('	dat.Activate (%d);\n'%(k))
                    for k, v in elem_deact.iteritems():
                        if v: txt.write ('	dat.Deactivate (%d);\n'%(k))

                    # boundary conditions
                    nbrys, nbsID, ebrys, fbrys = get_brys (obj,stg)
                    if len(nbrys)>0: txt.write ('	FEM::NBrys_T nbrys(%d);\n'%len(nbrys))
                    if len(ebrys)>0: txt.write ('	FEM::EBrys_T ebrys(%d);\n'%len(ebrys))
                    if len(fbrys)>0: txt.write ('	FEM::FBrys_T fbrys(%d);\n'%len(fbrys))
                    if len(nbrys)>0:
                        txt.write ('	nbrys = ')
                        for i, nb in enumerate(nbrys):
                            if i<(len(nbrys)-1): txt.write('T(%g,%g,%g,"%s",%g), ' %(nb[0],nb[1],nb[2],nb[3],nb[4]))
                            else:                txt.write('T(%g,%g,%g,"%s",%g);\n'%(nb[0],nb[1],nb[2],nb[3],nb[4]))
                    if len(ebrys)>0:
                        txt.write ('	ebrys = ')
                        for i, eb in enumerate(ebrys):
                            if i<(len(ebrys)-1): txt.write('T(%d,"%s",%g), ' %(eb[0],eb[1],eb[2]))
                            else:                txt.write('T(%d,"%s",%g);\n'%(eb[0],eb[1],eb[2]))
                    if len(fbrys)>0:
                        txt.write ('	fbrys = ')
                        for i, fb in enumerate(fbrys):
                            if i<(len(fbrys)-1): txt.write('T(%d,"%s",%g), ' %(fb[0],fb[1],fb[2]))
                            else:                txt.write('T(%d,"%s",%g);\n'%(fb[0],fb[1],fb[2]))
                    txt.write ('	dat.SetBrys       (&mesh, ')
                    if len(nbrys)>0: txt.write ('&nbrys, ')
                    else:            txt.write ('NULL, ')
                    if len(ebrys)>0: txt.write ('&ebrys, ')
                    else:            txt.write ('NULL, ')
                    if len(fbrys)>0: txt.write ('&fbrys);\n')
                    else:            txt.write ('NULL);\n')
                    for nb in nbsID:
                        txt.write ('	dat.Nod('+str(nb[0])+')->Bry   ("'+nb[1]+'",'+str(nb[2])+');\n')

                    # apply body forces
                    if abf: txt.write ('	dat.AddVolForces   ();\n')

                    # solve
                    txt.write ('	sol.SolveWithInfo (%d, %g, %d, "%s\\n");\n'%(ndiv,dtime,num,desc))

                    # clear displacements
                    if cdi: txt.write ('	dat.ClearDisp ();\n')

            # main
            txt.write ('\n}\n')
            txt.write ('catch (Exception * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }\n')
            txt.write ('catch (char const * m) { std::cout<<"Fatal: "<<m<<std::endl;  exit(1); }\n')
            txt.write ('catch (...) { std::cout << "Some exception (...) ocurred\\n"; } \n')

        else: # Python script

            # import libraries
            if not di.key('fullsc'):
                txt.write ('import Blender, bpy\n')
                txt.write ('import msys_mesh as me\n')
                txt.write ('import msys_fem  as mf\n')
                txt.write ('import mechsys   as ms\n')
            else: txt.write ('import mechsys as ms\n')

            # change cursor
            if not di.key('fullsc'):
                txt.write ('\n# Show running cursor\n')
                txt.write ('Blender.Window.WaitCursor(1)\n')

            # mesh
            txt.write ('\n###################################################################################### Mesh #####\n')

            # generate mesh
            txt.write ('\n# Mesh generation\n')
            if not di.key('fullsc'):
                txt.write ('obj  = bpy.data.objects["'+obj.name+'"]\n')
                if   mesh_type=='struct':   txt.write ('mesh = me.gen_struct_mesh()\n')
                elif mesh_type=='unstruct': txt.write ('mesh = me.gen_unstruct_mesh()\n')
                elif mesh_type=='frame':    txt.write ('mesh = me.gen_frame_mesh()\n')
            else:
                if   mesh_type=='struct':   me.gen_struct_mesh   (True, txt)
                elif mesh_type=='unstruct': me.gen_unstruct_mesh (True, txt)
                elif mesh_type=='frame':    me.gen_frame_mesh    (      txt)

            # FEM
            txt.write ('\n###################################################################################### FEM ######\n')

            # data and solver
            txt.write ('\n# Data and Solver\n')
            txt.write ('dat = ms.data   (%d)\n'        % ndim)
            txt.write ('sol = ms.solver (dat, "%s")\n' % obj.name)

            # element attributes
            neatt = len(eatts)
            txt.write ('\n# Element attributes\n')
            for i, ea in enumerate(eatts):
                if i==0:       estr  = 'eatts = ['
                else:          estr  = '         '
                if ea[7]:      estr += '[%d, "%s", "%s", "%s", "%s", "%s", "%s", True ]'#,  # %s\n'
                else:          estr += '[%d, "%s", "%s", "%s", "%s", "%s", "%s", False]'#,  # %s\n'
                if i==neatt-1: estr += ']  # %s\n'
                else:          estr += ',  # %s\n'
                txt.write (estr % (ea[0],ea[1],ea[2],ea[3],ea[4],ea[5],ea[6],ea[8]))

            # set geometry: nodes and elements
            txt.write ('\n# Set nodes and elements (geometry)\n')
            if mesh_type=='frame': txt.write ('dat.set_only_frame  () # frame (beam/truss) mesh only\n')
            txt.write ('dat.set_nodes_elems (mesh, eatts)\n')

            # add reinforcements
            if has_reinf:
                txt.write ('\n# Add reinforcements\n')
                txt.write ('reinfs = '+reinfs.__str__()+'\n')
                txt.write ('dat.add_reinfs (reinfs, eatts)\n')

            # add linear elements
            if has_lines:
                txt.write ('\n# Add linear elements\n')
                txt.write ('lines = '+lines.__str__()+'\n')
                txt.write ('dat.add_lin_elems (lines, eatts)\n')

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
                        dtime =     v[5]
                        act   = int(v[6])
                        break

                # run only if stage is active
                if act:

                    # activate and deactivate elements
                    txt.write ('\n# Stage # %d --------------------------------------------------------------\n'%num)
                    elem_act, elem_deact = get_act_deact (obj,stg)
                    for k, v in elem_act.iteritems():
                        if v: txt.write ('dat.activate (%d)\n'%(k))
                    for k, v in elem_deact.iteritems():
                        if v: txt.write ('dat.deactivate (%d)\n'%(k))

                    # boundary conditions
                    nbrys, nbsID, ebrys, fbrys = get_brys (obj,stg)
                    if len(nbrys)>0: txt.write ('nbrys = '+nbrys.__str__()+'\n')
                    if len(ebrys)>0: txt.write ('ebrys = '+ebrys.__str__()+'\n')
                    if len(fbrys)>0: txt.write ('fbrys = '+fbrys.__str__()+'\n')
                    txt.write ('dat.set_brys         (mesh, ')
                    if len(nbrys)>0: txt.write ('nbrys, ')
                    else:            txt.write ('[], ')
                    if len(ebrys)>0: txt.write ('ebrys, ')
                    else:            txt.write ('[], ')
                    if len(fbrys)>0: txt.write ('fbrys)\n')
                    else:            txt.write ('[])\n')
                    for nb in nbsID:
                        txt.write ('dat.nod('+str(nb[0])+').bry       ("'+nb[1]+'",'+str(nb[2])+')\n')

                    # apply body forces
                    if abf: txt.write ('dat.add_vol_forces   ()\n')

                    # solve
                    txt.write ('sol.solve_with_info  (%d, %g, %d, "%s\\n")\n'%(ndiv,dtime,num,desc))

                    # clear displacements
                    if cdi: txt.write ('dat.clear_disp       ()\n')

                    # save results
                    if not di.key('fullsc'): txt.write ('mf.save_results      (sol, dat, obj, %d)\n'%num)

            # change cursor
            if not di.key('fullsc'):
                txt.write ('\n# Hide running cursor\n')
                txt.write ('Blender.Window.WaitCursor(0)\n')

    else:
        # get/set mesh
        if   mesh_type=='struct':   mesh = me.gen_struct_mesh   ()
        elif mesh_type=='unstruct': mesh = me.gen_unstruct_mesh ()
        elif mesh_type=='frame':    mesh = me.gen_frame_mesh    ()

        # data and solver
        dat = ms.data   (ndim)
        sol = ms.solver (dat, obj.name)

        # element attributes
        elem_atts = []
        for ea in eatts: elem_atts.append(ea[:8])

        # set geometry: nodes and elements
        if mesh_type=='frame': dat.set_only_frame()
        dat.set_nodes_elems (mesh, elem_atts)

        # add reinforcements
        if has_reinf: dat.add_reinfs (reinfs, elem_atts)

        # add linear elements
        if has_lines: dat.add_lin_elems (lines, elem_atts)

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
                    dtime =     v[5]
                    act   = int(v[6])
                    break

            # run only if stage is active
            if act:

                # activate and deactivate elements
                elem_act, elem_deact = get_act_deact (obj,stg)
                for k, v in elem_act.iteritems():
                    if v: dat.activate(k)
                for k, v in elem_deact.iteritems():
                    if v: dat.deactivate(k)

                # boundary conditions
                nbrys, nbsID, ebrys, fbrys = get_brys (obj,stg)
                dat.set_brys (mesh, nbrys, ebrys, fbrys)
                for nb in nbsID: dat.nod(nb[0]).bry(nb[1],nb[2])

                # apply body forces
                if abf: dat.apply_body_forces()

                # solve
                sol.solve_with_info (ndiv,dtime,num,desc+'\n')

                # clear displacements
                if cdi: dat.clear_displacements()

                # save results
                save_results (sol, dat, obj, num)

    # end
    #Blender.Run(txt.name)
    Blender.Window.WaitCursor(0)


def save_results(sol, dat, obj, stage_num):
    # dictionary
    s = str(stage_num)
    if not obj.properties.has_key('res'): obj.properties['res'] = {}
    obj.properties['res'][s] = {}

    # menu with labels
    obj.properties['res'][s]['idx2lbl'] = {} # map label index to label key
    lbs = []
    sol.out().get_labels (lbs)
    menu = 'Labels %t|'
    for i, l in enumerate(lbs):
        obj.properties['res'][s]['idx2lbl'][str(i)] = l
        menu += l + ' %x' + str(i+1) + '|'
    obj.properties['res'][s]['lblmnu'] = menu

    # save results at nodes
    for l in lbs:
        vals = []
        for i in range(dat.nnodes()): vals.append (sol.out().val(i, l))
        obj.properties['res'][s][l] = vals

    # save extra output
    obj.properties['res'][s]['extra'] = {}
    for i in range(dat.nelems()):
        if dat.ele(i).has_extra() and dat.ele(i).is_active():
            ide = str(i)
            obj.properties['res'][s]['extra'][ide] = {}
            dat.ele(i).calc_dep_vars()
            co, no, va = dict(), [], dict()
            dat.ele(i).out_extra (co, no, va)
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
                            if maxv>obj.properties['res'][s][key][0]: obj.properties['res'][s][key] = [maxv, i, dat.ele(i).nod(0).x(), dat.ele(i).nod(0).y()]
                        else: obj.properties['res'][s][key] = [maxv, i, dat.ele(i).nod(0).x(), dat.ele(i).nod(0).y()]

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
