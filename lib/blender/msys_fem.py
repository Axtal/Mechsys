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

from   multiprocessing import Process
import subprocess, math
import Blender, bpy
from   Blender.Mathutils import Vector
import mechsys   as ms
import msys_dict as di
import msys_mesh as me


def get_act_deact(obj, stg):
    d = di.load_dict()
    activate   = {}
    deactivate = {}
    if obj.properties[stg].has_key('eatts'):
        for k, v in obj.properties[stg]['eatts'].iteritems():
            tag             = int(v[0])
            activate  [tag] = int(v[7])
            deactivate[tag] = int(v[8])
    return activate, deactivate


def get_prps_mdls(obj, pty, stg, script=False, cpp=False):
    d    = di.load_dict()
    prps = {}
    mdls = {}

    for k, v in obj.properties[stg]['eatts'].iteritems():
        # prps
        tag   = int(v[0])                      # ex: -1
        geom  = d['pty2geom'][pty][int(v[2])]  # ex: Hex8
        gty   = d['pty2gty'] [pty][int(v[3])]  # ex: psa
        matID = str(int(v[4]))                 # mat id
        act   = True if int(v[6]) else False   # active ?
        prps[tag] = {'prob':d['pty'][pty], 'geom':geom, gty:1.0}
        if gty=='pse': prps[tag]['h'] = v[5]   # thickness (for pse only)

        # mdls
        m    = obj.properties['mats'][matID]
        name = d['mdl'][int(m[0])]
        desc = obj.properties['texts'][str(int(m[1]))]
        if name=='Rod':
            mdls[tag] = {'name':name, gty:1.0, 'E':m[2], 'A':m[3]}
        elif name=='Beam':
            mdls[tag] = {'name':name, gty:1.0, 'E':m[2], 'A':m[3], 'Izz':m[4]}
        elif name=='LinElastic':
            mdls[tag] = {'name':name, gty:1.0, 'E':m[2], 'nu':m[5]}
        elif name=='ElastoPlastic':
            mdls[tag] = {'name':name, gty:1.0, 'E':m[2], 'nu':m[5], 'fc':d['fc'][m[6]], 'sY':m[7], 'cu':m[8]}
        elif name=='CamClay':
            mdls[tag] = {'name':name, gty:1.0, 'lam':m[9], 'kap':m[10], 'nu':m[5], 'phi':m[11]}
        elif name=='LinFlow':
            mdls[tag] = {'name':name, gty:1.0, 'k':m[12]}

    if script and cpp:
        # prps
        for tag, prp in prps.iteritems():
            keys, vals, count, num = '"', '', 1, len(prp)
            for k, v in prp.iteritems():
                keys += (k+'"' if count==num else k+' ')
                if   k=='prob': vals += ('PROB("'+v+'")' if count==num else 'PROB("'+v+'"), ')
                elif k=='geom': vals += ('GEOM("'+v+'")' if count==num else 'GEOM("'+v+'"), ')
                else:           vals += (str(v)          if count==num else str(v)+', ')
                count += 1
            prp['CPPKEYS'] = keys
            prp['CPPVALS'] = vals

        # mdls
        for tag, mdl in mdls.iteritems():
            keys, vals, count, num = '"', '', 1, len(mdl)
            for k, v in mdl.iteritems():
                keys += (k+'"' if count==num else k+' ')
                if k=='name': vals += ('MODEL("'+v+'")' if count==num else 'MODEL("'+v+'"), ')
                else:         vals += (str(v)           if count==num else str(v) + ', ')
                count += 1
            mdl['CPPKEYS'] = keys
            mdl['CPPVALS'] = vals

    return prps, mdls


def get_brys(obj, stg, script=False, cpp=False):
    d   = di.load_dict()
    pty = obj.properties['pty']
    nbrys, ebrys, fbrys = {}, {}, {}

    if obj.properties[stg].has_key('nbrys'):
        for k, v in obj.properties[stg]['nbrys'].iteritems():
            nbrys[int(v[0])] = {d['pty2Ndfv'][pty][int(v[1])] : v[2]}

    if obj.properties[stg].has_key('ebrys'):
        for k, v in obj.properties[stg]['ebrys'].iteritems():
            ebrys[int(v[0])] = {d['pty2Fdfv'][pty][int(v[1])] : v[2]}

    if obj.properties[stg].has_key('fbrys'):
        for k, v in obj.properties[stg]['fbrys'].iteritems():
            fbrys[int(v[0])] = {d['pty2Fdfv'][pty][int(v[1])] : v[2]}

    if script and cpp:
        for tag, nbc in nbrys.iteritems():
            keys, vals, count, num = '"', '', 1, len(nbc)
            for k, v in nbc.iteritems():
                keys  += (k+'"'  if count==num else k+' ')
                vals  += (str(v) if count==num else str(v)+', ')
                count += 1
            nbc['CPPKEYS'] = keys
            nbc['CPPVALS'] = vals

        for tag, ebc in ebrys.iteritems():
            keys, vals, count, num = '"', '', 1, len(ebc)
            for k, v in ebc.iteritems():
                keys  += (k+'"'  if count==num else k+' ')
                vals  += (str(v) if count==num else str(v)+', ')
                count += 1
            ebc['CPPKEYS'] = keys
            ebc['CPPVALS'] = vals

        for tag, fbc in fbrys.iteritems():
            keys, vals, count, num = '"', '', 1, len(fbc)
            for k, v in fbc.iteritems():
                keys  += (k+'"'  if count==num else k+' ')
                vals  += (str(v) if count==num else str(v)+', ')
                count += 1
            fbc['CPPKEYS'] = keys
            fbc['CPPVALS'] = vals

    return nbrys, ebrys, fbrys


def get_stage_info(obj, num):
    stage = {}
    for k, v in obj.properties['stages'].iteritems():
        if int(v[0])==num:
            stage['stg']   = 'stg_'+k
            stage['desc']  = obj.properties['texts'][str(int(v[1]))]
            stage['abf']   = True if int(v[2]) else False # apply body forces ?
            stage['cdi']   = True if int(v[3]) else False # clear displacements ?
            stage['ndiv']  = int(v[4])
            stage['act']   = int(v[5])
            stage['type']  = int(v[6])
            if stage['type']==2: # dynamics
                stage['tf']     = float(v[7])
                stage['dt']     = float(v[8])
                stage['dtOut']  = float(v[9])
                stage['outVTU'] = True if int(v[31]) else False
    return stage


class FEMData:
    def __init__(self, script=False):
        # get active object
        self.obj = di.get_obj()
        if not self.obj.properties.has_key('mats'):     raise Exception('Please define materials first')
        if not self.obj.properties.has_key('pty'):      raise Exception('Please define the problem type first')
        if not self.obj.properties.has_key('msh_type'): raise Exception('Please generate the mesh first or set the "Frame mesh" toggle')
        if not self.obj.properties.has_key('stages'):   raise Exception('Please add stages first')
        if not self.obj.properties.has_key('is3d'):     self.obj.properties['is3d'] = False

        # input data
        self.filekey  = self.obj.name+'_fem'                    # file key
        self.pty      = self.obj.properties['pty']              # problem type
        self.str_pty  = di.key('pty')[self.pty]                 # string for problem type, "Equilib", "Flow", ...
        self.msh_type = self.obj.properties['msh_type']         # mesh type
        self.nstages  = len(self.obj.properties['stages'])      # number of stages
        self.ndim     = 3 if self.obj.properties['is3d'] else 2 # space dimension

        # get first stage
        stage_id, first_stg = di.find_stage (self.obj, 1)
        if not self.obj.properties[first_stg].has_key('eatts'): raise Exception('Please define element attributes/properties first')

        # materials and prps
        self.cpp = True if di.key('fem_cpp') else False
        self.prps, self.mdls = get_prps_mdls (self.obj, self.pty, first_stg, script, self.cpp)

        # nodes and elements for output
        self.nout, self.eout = [], []
        if self.obj.properties.has_key('nout'): self.nout = [int(n) for n in self.obj.properties['nout'].split(' ')]
        if self.obj.properties.has_key('eout'): self.eout = [int(e) for e in self.obj.properties['eout'].split(' ')]

def gen_script():
    dat = FEMData(True)
    if dat.cpp: # C++ script ===================================================================================================
        txt = Blender.Text.New(dat.filekey+'.cpp')

        # headers
        txt.write ('// MechSys\n')
        txt.write ('#include <mechsys/fem/fem.h>\n')
        txt.write ('\nusing FEM::PROB;\n')
        txt.write ('using FEM::GEOM;\n')
        txt.write ('\nint main(int argc, char **argv) try\n')
        txt.write ('{\n')

        # mesh
        txt.write ('    // mesh\n')
        if   dat.msh_type=='struct':   me.gen_struct_mesh   (True, txt, True, False) # gen_script, txt, cpp, with_headers
        elif dat.msh_type=='unstruct': me.gen_unstruct_mesh (True, txt, True, False) # gen_script, txt, cpp, with_headers
        elif dat.msh_type=='frame':    me.gen_frame_mesh    (      txt, True, False) # gen_script, txt, cpp, with_headers

        # elements properties
        txt.write ('\n    // elements properties\n')
        txt.write ('    Dict prps;\n')
        for tag, prp in dat.prps.iteritems(): txt.write ('    prps.Set ('+str(tag)+', '+prp['CPPKEYS']+', '+prp['CPPVALS']+');\n')

        # models
        txt.write ('\n    // models\n')
        txt.write ('    Dict mdls;\n')
        for tag, mdl in dat.mdls.iteritems(): txt.write ('    mdls.Set ('+str(tag)+', '+mdl['CPPKEYS']+', '+mdl['CPPVALS']+');\n')

        # initial values
        txt.write ('\n    // initial values\n')
        txt.write ('    Dict inis;\n')
        for tag, mdl in dat.mdls.iteritems():
            if   dat.str_pty=='Equilib': txt.write ('    inis.Set ('+str(tag)+', "sx sy sz sxy", 0.0, 0.0, 0.0, 0.0);\n')
            elif dat.str_pty=='Flow':    txt.write ('    inis.Set ('+str(tag)+', "vx vy", 0.0, 0.0);\n')

        # domain
        txt.write ('\n    // domain\n')
        txt.write ('    FEM::Domain dom(mesh, prps, mdls, inis);\n')

        # solver
        txt.write ('\n    // solver\n')
        txt.write ('    FEM::Solver sol(dom);\n')

        # solve each stage
        for num in range(1,dat.nstages+1):

            # find stage info
            stage = get_stage_info (dat.obj, num)

            # run only if stage is active
            if stage['act']:
                # stage #
                txt.write ('\n    // stage # %d ==============================================================\n'%num)

                # boundary conditions
                txt.write ('\n    // boundary conditions\n')
                nbrys, ebrys, fbrys = get_brys (dat.obj, stage['stg'], True, True) # obj, stg, script, cpp
                if num==1: txt.write ('    Dict bcs;\n')
                for tag, nbc in nbrys.iteritems(): txt.write ('    bcs.Set ('+str(tag)+', '+nbc['CPPKEYS']+', '+nbc['CPPVALS']+');\n')
                for tag, ebc in ebrys.iteritems(): txt.write ('    bcs.Set ('+str(tag)+', '+ebc['CPPKEYS']+', '+ebc['CPPVALS']+');\n')
                for tag, fbc in fbrys.iteritems(): txt.write ('    bcs.Set ('+str(tag)+', '+fbc['CPPKEYS']+', '+fbc['CPPVALS']+');\n')
                txt.write ('    dom.SetBCs (bcs);\n')

                # apply body forces
                if stage['abf']:
                    txt.write ('\n    // apply body forces\n')
                    txt.write ('    dom.Gravity ();\n')

                # solve
                txt.write ('\n    // solve\n')
                if stage['type']==0: # equilib/steady
                    txt.write ('    sol.Solve (%d);\n'%(stage['ndiv']))
                elif stage['type']==2: # dynamics
                    extra = ', "'+dat.filekey+'"' if stage['outVTU'] else ''
                    txt.write ('    sol.DynSolve (%g, %g, %g%s); // tf,dt,dtOut,filekey\n' % (stage['tf'],stage['dt'],stage['dtOut'],extra))

                # clear displacements
                if stage['cdi']:
                    txt.write ('\n    // clear displacements\n')
                    txt.write ('    SDPair uvs;\n')
                    if dat.ndim==2: txt.write ('    uvs.Set ("ux uy", 0.0,0.0);\n')
                    else:           txt.write ('    uvs.Set ("ux uy uz", 0.0,0.0,0.0);\n')
                    txt.write ('    dom.SetUVals (uvs);\n')

                # output
                txt.write ('\n    // output\n')
                txt.write ('    dom.WriteVTU ("'+dat.filekey+'_'+stage['stg']+'");\n')

        # footer
        txt.write ('\n    return 0;\n')
        txt.write ('}\n')
        txt.write ('MECHSYS_CATCH\n')

    else: # Python script ===================================================================================================
        txt = Blender.Text.New(dat.filekey+'.py')

        # headers
        txt.write ('# headers\n')
        txt.write ('from mechsys import *\n')

        # mesh
        txt.write ('\n# mesh\n')
        if   dat.msh_type=='struct':   me.gen_struct_mesh   (True, txt, False, False) # gen_script, txt, cpp, with_headers
        elif dat.msh_type=='unstruct': me.gen_unstruct_mesh (True, txt, False, False) # gen_script, txt, cpp, with_headers
        elif dat.msh_type=='frame':    me.gen_frame_mesh    (      txt, False, False) # gen_script, txt, cpp, with_headers

        # elements properties
        txt.write ('\n# elements properties\n')
        txt.write ('prps = Dict()\n')
        for tag, prp in dat.prps.iteritems():
            str_prp = prp.__str__().replace('\''+prp['prob']+'\'','PROB("'+prp['prob']+'")').replace('\''+prp['geom']+'\'','GEOM("'+prp['geom']+'")')
            txt.write ('prps.Set ('+str(tag)+', '+str_prp+')\n')

        # models
        txt.write ('\n# models\n')
        txt.write ('mdls = Dict()\n')
        for tag, mdl in dat.mdls.iteritems():
            str_mdl = mdl.__str__().replace('\''+mdl['name']+'\'','MODEL("'+mdl['name']+'")')
            txt.write ('mdls.Set ('+str(tag)+', '+str_mdl+')\n')

        # initial values
        txt.write ('\n# initial values\n')
        txt.write ('inis = Dict()\n')
        for tag, mdl in dat.mdls.iteritems():
            if   dat.str_pty=='Equilib': txt.write ('inis.Set ('+str(tag)+', {"sx":0.0, "sy":0.0, "sz":0.0, "sxy":0.0})\n')
            elif dat.str_pty=='Flow':    txt.write ('inis.Set ('+str(tag)+', {"vx":0.0, "vy":0.0})\n')

        # domain
        txt.write ('\n# domain\n')
        txt.write ('dom = FEM_Domain (mesh, prps, mdls, inis)\n')

        # solver
        txt.write ('\n# solver\n')
        txt.write ('sol = FEM_Solver (dom)\n')

        # solve each stage
        for num in range(1,dat.nstages+1):

            # find stage info
            stage = get_stage_info (dat.obj, num)

            # run only if stage is active
            if stage['act']:

                # stage
                txt.write ('\n# stage # %d ==============================================================\n'%num)

                # boundary conditions
                txt.write ('\n# boundary conditions\n')
                nbrys, ebrys, fbrys = get_brys (dat.obj, stage['stg'], True) # obj, stg, script
                if num==1: txt.write ('bcs = Dict()\n')
                for tag, nbc in nbrys.iteritems(): txt.write ('bcs.Set ('+str(tag)+', '+nbc.__str__()+')\n')
                for tag, ebc in ebrys.iteritems(): txt.write ('bcs.Set ('+str(tag)+', '+ebc.__str__()+')\n')
                for tag, fbc in fbrys.iteritems(): txt.write ('bcs.Set ('+str(tag)+', '+fbc.__str__()+')\n')
                txt.write ('dom.SetBCs (bcs)\n')

                # activate and deactivate elements
                elem_act, elem_deact = get_act_deact (dat.obj, stage['stg'])
                for k, v in elem_act.iteritems():
                    if v: txt.write ('dom.Activate (%d)\n'%(k))
                for k, v in elem_deact.iteritems():
                    if v: txt.write ('dom.Deactivate (%d)\n'%(k))

                # apply body forces
                if stage['abf']:
                    txt.write ('\n# apply body forces\n')
                    txt.write ('dom.Gravity ()\n')

                # solve
                txt.write ('\n# solve\n')
                if stage['type']==0: # equilib/steady
                    txt.write ('sol.Solve (%d)\n'%(stage['ndiv']))
                elif stage['type']==2: # dynamics
                    extra = ', "'+dat.filekey+'"' if stage['outVTU'] else ''
                    txt.write ('sol.DynSolve (%g, %g, %g%s) # tf,dt,dtOut,filekey\n' % (stage['tf'],stage['dt'],stage['dtOut'],extra))

                # clear displacements
                if stage['cdi']:
                    txt.write ('\n# clear displacements\n')
                    txt.write ('uvs = SDPair()\n')
                    if dat.ndim==2: txt.write ('uvs.Set ({"ux":0.0, "uy":0.0})\n')
                    else:           txt.write ('uvs.Set ({"ux":0.0, "uy":0.0, "uz":0.0})\n')
                    txt.write ('dom.SetUVals (uvs)\n')

                # output
                txt.write ('\n# output\n')
                txt.write ('dom.WriteVTU ("'+dat.filekey+'_'+stage['stg']+'")\n')


def run_simulation(running, fatal):
    try:
        dat = FEMData()

        # mesh
        if   dat.msh_type=='struct':   mesh = me.gen_struct_mesh   (False)
        elif dat.msh_type=='unstruct': mesh = me.gen_unstruct_mesh (False)
        elif dat.msh_type=='frame':    mesh = me.gen_frame_mesh    ()

        # elements properties
        fem_prps = ms.Dict()
        for tag, prp in dat.prps.iteritems():
            prp['prob'] = ms.PROB(prp['prob'])
            prp['geom'] = ms.GEOM(prp['geom'])
            fem_prps.Set (tag, prp)

        # models
        fem_mdls = ms.Dict()
        for tag, mdl in dat.mdls.iteritems():
            mdl['name'] = ms.MODEL(mdl['name'])
            fem_mdls.Set (tag, mdl)

        # initial values
        fem_inis = ms.Dict()
        for tag, mdl in dat.mdls.iteritems():
            if   dat.str_pty=='Equilib': fem_inis.Set (tag, {"sx":1.0, "sy":1.0, "sz":1.0, "sxy":1.0, "v0":1.5})
            elif dat.str_pty=='Flow':    fem_inis.Set (tag, {"vx":0.0, "vy":0.0})

        # domain
        dom = ms.FEM_Domain (mesh, fem_prps, fem_mdls, fem_inis)

        # nodes and elements for output
        if len(dat.nout)>0: dom.SetOutNods (dat.filekey, dat.nout)
        if len(dat.eout)>0: dom.SetOutEles (dat.filekey, dat.eout)

        # solver
        sol = ms.FEM_Solver (dom)
        sol.TolR = 1.0e-2

        # solve each stage
        bcs = ms.Dict()
        for num in range(1,dat.nstages+1):

            # find stage info
            stage = get_stage_info (dat.obj, num)

            # run only if stage is active
            if stage['act']:

                # boundary conditions
                nbrys, ebrys, fbrys = get_brys (dat.obj, stage['stg'])
                for tag, nbc in nbrys.iteritems(): bcs.Set (tag, nbc)
                for tag, ebc in ebrys.iteritems(): bcs.Set (tag, ebc)
                for tag, fbc in fbrys.iteritems(): bcs.Set (tag, fbc)
                dom.SetBCs (bcs)

                # activate and deactivate elements
                elem_act, elem_deact = get_act_deact (dat.obj, stage['stg'])
                for k, v in elem_act.iteritems():
                    if v: dom.Activate (k)
                for k, v in elem_deact.iteritems():
                    if v: dom.Deactivate (k)

                # apply body forces
                if stage['abf']: dom.Gravity ()

                # solve
                if stage['type']==0: # equilib/steady
                    sol.Solve (stage['ndiv'])
                elif stage['type']==2: # dynamics
                    if stage['outVTU']: sol.DynSolve (stage['tf'],stage['dt'],stage['dtOut'],dat.filekey)
                    else:               sol.DynSolve (stage['tf'],stage['dt'],stage['dtOut'])

                # clear displacements
                if stage['cdi']:
                    uvs = ms.SDPair()
                    if dat.ndim==2: uvs.Set ({"ux":0.0, "uy":0.0})
                    else:           uvs.Set ({"ux":0.0, "uy":0.0, "uz":0.0})
                    dom.SetUVals (uvs)

                # output
                dom.WriteVTU (dat.filekey+'_'+stage['stg'])

        # notify parent that we have finished
        running.value = 0

    except Exception, inst:
        print '[1;34mMechSys[0m: Error: '+'[1;31m'+inst.args[0]+'[0m'
        print '>>>>>>>>>>> exception caught by msys_fem.run_simulation <<<<<<<<<<<'
        running.value = 0
        fatal  .value = 1


def run():
    d = di.load_dict()
    if d['fem_running'].value: raise Exception('Another FEM simulation is already running')
    d['fem_running'].value = 1
    d['fem_fatal']  .value = 0
    d['fem_process'] = Process(target=run_simulation, args=(d['fem_running'],d['fem_fatal']))
    print "\n#############################  FEM: running simulation  #################################\n"
    Blender.Window.QRedrawAll()
    d['fem_process'].start()


def stop():
    d = di.load_dict()
    if d['fem_running'].value:
        d['fem_process'].terminate()
        d['fem_running'].value = 0
        print "\n#############################  FEM: simulation stoped  ##################################\n"
        Blender.Window.QRedrawAll()
    else: raise Exception('There is no FEM simulation running')


def paraview(dynamics=False):
    obj = di.get_obj()
    if obj.properties.has_key('stages'):
        if dynamics:
            fn = obj.name+'_fem_00000001.vtu'
            if not Blender.sys.exists(fn): raise Exception('File <'+fn+'> does not exist (please, run analysis first)')
            strkey = '_fem_..vtu'
        else:
            nstages = len(obj.properties['stages'])
            for i in range(nstages):
                fn = obj.name+'_fem_stg_%d.vtu'%i
                if not Blender.sys.exists(fn): raise Exception('File <'+fn+'> does not exist (please, run analysis first)')
            strkey = '_fem_stg_..vtu' if nstages>1 else '_fem_stg_0.vtu'
        try: pid = subprocess.Popen(['paraview', '--data='+obj.name+strkey]).pid
        except:
            Blender.Window.WaitCursor(0)
            raise Exception('Paraview is not available, please install it first')
    else: raise Exception('Please add stages first and run the simulation')

def plot_res(ele):
    obj = di.get_obj()
    if ele:
        if obj.properties.has_key('eout'):
            eout = [int(e) for e in obj.properties['eout'].split()][di.key('fem_eout_plt')]
            fn   = obj.name+"_fem_ele_%d_-1.res"%eout
            if not Blender.sys.exists(fn): raise Exception('File <'+fn+'> does not exist (please, set elements for output and run analysis first)')
            lin  = 'from msys_plotter import *\n'
            lin += 'plt = Plotter()\n'
            lin += 'plt.plot ("%s")\n'%fn
            lin += 'show()\n'
            f = open(obj.name+'_plot.py', 'w')
            f.write (lin)
            f.close ()
            try: pid = subprocess.Popen(['python', obj.name+'_plot.py']).pid
            except:
                Blender.Window.WaitCursor(0)
                raise Exception('FEM:plot_res: Python (Matplotlib) command failed')
        else: raise Exception('FEM: Calculation did not have any element for output')
    else:
        if obj.properties.has_key('nout'):
            nout = [int(n) for n in obj.properties['nout'].split()][di.key('fem_nout_plt')]
            fn   = obj.name+"_fem_nod_%d_-1.res"%nout
            if not Blender.sys.exists(fn): raise Exception('File <'+fn+'> does not exist (please, set nodes for output and run analysis first)')
            lin  = 'from msys_plotter import *\n'
            lin += 'plt = Plotter()\n'
            lin += 'plt.plot_node ("%s")\n'%fn
            lin += 'show()\n'
            f = open(obj.name+'_plot.py', 'w')
            f.write (lin)
            f.close ()
            try: pid = subprocess.Popen(['python', obj.name+'_plot.py']).pid
            except:
                Blender.Window.WaitCursor(0)
                raise Exception('FEM:plot_res: Python (Matplotlib) command failed')
        else: raise Exception('FEM: Calculation did not have any node for output')
