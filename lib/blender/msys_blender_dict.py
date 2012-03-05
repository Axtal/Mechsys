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

from   multiprocessing import Value
import Blender
import bpy
import string


# ================================================================================ Dictionary

def load_dict_for_scroll():
    dict = Blender.Registry.GetKey('MechSysDict')
    if not dict: dict, dict['gui_inirow'] = {}, 0
    return dict

def load_dict():
    dict = Blender.Registry.GetKey('MechSysDict')
    if not dict:
        dict                  = {}
        # GUI
        dict['gui_show_cad']  = False
        dict['gui_show_mesh'] = False
        dict['gui_show_mat']  = False
        dict['gui_show_fem']  = False
        dict['gui_show_res']  = False
        dict['gui_show_dem']  = False
        dict['gui_inirow']    = 0
        # SETTINGS
        dict['show_props']    = False
        dict['show_e_ids']    = False # edges
        dict['show_v_ids']    = False # vertices
        dict['show_blks']     = True
        dict['show_Bids']     = False
        dict['show_axes']     = True
        dict['show_regs']     = True
        dict['show_lins']     = True  # show linear elements
        dict['show_vtags']    = True
        dict['show_etags']    = True
        dict['show_ftags']    = True
        # CAD
        dict['cad_x']         = '0.0'
        dict['cad_y']         = '0.0'
        dict['cad_z']         = '0.0'
        dict['cad_rad']       = '0.0'
        dict['cad_stp']       = 10
        # MESH
        dict['newvtag']       = -100
        dict['newetag']       = -10
        dict['newftag']       = -10
        dict['hide_mesh']     = False
        dict['fmsh_cpp']      = True  # frame mesh: generate C++ script
        dict['smsh_cpp']      = False # structured mesh: generate C++ script
        dict['umsh_cpp']      = False # unstructured mesh: generate C++ script
        # FEM
        dict['fem_prob']      = 0
        dict['fem_stage']     = 0     # stage ID
        dict['fem_cpp']       = False # generate C++ script
        dict['fem_process']   = None
        dict['fem_running']   = Value('i',0)
        dict['fem_fatal']     = Value('i',0)
        dict['fem_nout_plt']  = 0
        dict['fem_eout_plt']  = 0
        # DEM
        dict['dem_Lx']         = 4.0
        dict['dem_Ly']         = 4.0
        dict['dem_Lz']         = 4.0
        dict['dem_Nx']         = 2
        dict['dem_Ny']         = 2
        dict['dem_Nz']         = 2
        dict['dem_R']          = 0.1
        dict['dem_rho']        = 1.0
        dict['dem_seed']       = 123
        dict['dem_prob']       = 1.0
        dict['dem_pkg']        = 1    # 0:Spheres, 1:HCP, 2:Voronoi, 3:Mesh
        dict['dem_res']        = 8
        dict['dem_draw_verts'] = True
        dict['dem_draw_edges'] = True
        dict['dem_Kn']         = 10000.0
        dict['dem_Kt']         = 5000.0
        dict['dem_Gn']         = 16.0
        dict['dem_Gt']         = 0.0
        dict['dem_mu']         = 0.4
        dict['dem_beta']       = 0.12
        dict['dem_eta']        = 1.0
        dict['dem_iso_pf']     = 5.0
        dict['dem_iso_timef']  = 50.0
        dict['dem_iso_dt']     = 0.001
        dict['dem_iso_dtout']  = 1.0
        dict['dem_iso_render'] = True
        dict['dem_ttt_comp']   = True
        dict['dem_ttt_ext']    = False
        dict['dem_ttt_pcte']   = False
        dict['dem_ttt_pf']     = 5.0
        dict['dem_ttt_qf']     = 0.0
        dict['dem_ttt_thf']    = 30.0
        dict['dem_ttt_timef']  = 200.0
        dict['dem_ttt_dt']     = 0.001
        dict['dem_ttt_dtout']  = 1.0
        dict['dem_ttt_pex']    = False
        dict['dem_ttt_pey']    = False
        dict['dem_ttt_pez']    = True
        dict['dem_ttt_exf']    = 0.0
        dict['dem_ttt_eyf']    = 0.0
        dict['dem_ttt_ezf']    = -0.2
        dict['dem_ttt_render'] = True
        dict['dem_cpp_script'] = False
        dict['dem_process']    = None
        dict['dem_running']    = Value('i',0)
        dict['dem_fatal']      = Value('i',0)

        # DEM packings
        dict['dem_pkgs']     = {0:'Spheres', 1:'Spheres HCP', 2:'Voronoi', 3:'From Mesh'}
        dict['dem_pkgs_mnu'] = 'Packing type %t|From Mesh %x4|Voronoi %x3|Spheres HCP %x2|Spheres %x1'

        # Models
        dict['mdl']      = {0:'Rod', 1:'Beam', 2:'LinElastic', 3:'ElastoPlastic', 4:'CamClay', 5:'LinFlow'}
        dict['mdlmnu']   = 'Models %t|LinFlow %x6|CamClay %x5|ElastoPlastic %x4|LinElastic %x3|Beam %x2|Rod %x1'
        dict['mdl2nrow'] = {0:1, 1:1, 2:1, 3:2, 4:2, 5:1}
        dict['fc']       = {0:'VM', 1:'DP', 2:'MC' }
        dict['fcmnu']    = 'Failure Criteria %t|Mohr-Coulomb %x3|Drucker-Prager %x2|von Mises %x1'

        # Problem types
        dict['pty'] = { 0:'Equilib', 1:'Flow' }
        dict['ptymnu'] = 'Problem Types %t|Flow %x2|Equilibrium %x1'

        # DOF vars
        dict['pty2Ndfv'] = { # problem type to Node DOF vars
                0:{0:'ux', 1:'uy', 2:'uz',  3:'fx', 4:'fy', 5:'fz',  6:'wz', 7:'mz'}, # Equilib 
                1:{0:'H',  1:'Q'} }                                                   # Flow

        dict['pty2Fdfv'] = { # problem type to Edge/Face DOF vars
                0:{0:'ux', 1:'uy', 2:'uz',  3:'qx', 4:'qy', 5:'qz',  6:'qn', 7:'qnl', 8:'qnr'}, # Equilib 
                1:{0:'conv', 1:'flux'} }                                                        # Flow

        dict['pty2Nmnu'] = {
                0:'N vars %t|mz %x8|wz %x7|fz %x6|fy %x5|fx %x4|uz %x3|uy %x2|ux %x1', # Equilib
                1:'N vars %t|Q %x2|H %x1' }                                            # Flow

        dict['pty2Fmnu'] = {
                0:'E/F vars %t|qnr %x9|qnl %x8|qn %x7|qz %x6|qy %x5|qx %x4|uz %x3|uy %x2|ux %x1', # Equilib
                1:'E/F vars %t|flux %x2|conv %x1' }                                               # Flow

        # Problem/geometry
        dict['pty2geom'] = { # problem type to geom type
                0:{0:'Tri3', 1:'Tri6', 2:'Quad4', 3:'Quad8', 4:'Tet4', 5:'Tet10', 6:'Hex8', 7:'Hex20', 8:'Rod', 9:'Beam' }, # Equilib
                1:{0:'Tri3', 1:'Tri6', 2:'Quad4', 3:'Quad8', 4:'Tet4', 5:'Tet10', 6:'Hex8', 7:'Hex20' } }                   # Flow

        dict['pty2gpmnu'] = {
                0:'Types %t|Beam %x10|Rod %x9|Hex20 %x8|Hex8 %x7|Tet10 %x6|Tet4 %x5|Quad8 %x4|Quad4 %x3|Tri6 %x2|Tri3 %x1', # Equilib
                1:'Types %t|Hex20 %x8|Hex8 %x7|Tet10 %x6|Tet4 %x5|Quad8 %x4|Quad4 %x3|Tri6 %x2|Tri3 %x1' }                  # Flow

        # Geometry types
        dict['pty2gty'] = {
                0:{0:'pse', 1:'psa', 2:'axs', 3:'d3d', 4:'fra'},
                1:{0:'d2d', 1:'d3d' } }

        dict['pty2gmnu'] = {
                0:'Geometry type %t|fra %x5|d3d %x4|axs %x3|psa %x2|pse %x1',
                1:'Geometry type %t|d3d %x2|d2d %x1' }

        # FEM constants
        dict['Scheme_t']     = { 0:'FE', 1:'ME', 2:'NR' }             # Steady time integration scheme: Forward-Euler, Modified-Euler, Newton-Rhapson
        dict['TScheme_t']    = { 0:'SS11' }                           # Transient time integration scheme: (Single step/O1/1st order)
        dict['DScheme_t']    = { 0:'SS22', 1:'GN22', 2:'GNHMCoup' }   # Dynamic time integration scheme: (Single step/O2/2nd order), (Generalized Newmark/O2/2nd order)
        dict['Damping_t']    = { 0:'None', 1:'Rayleigh', 2:'HMCoup' } # Damping type: none, Rayleigh type (C=alp*M+bet*K), HydroMechCoupling
        dict['schememnu']    = 'Solver scheme %t|Newton-Rhapson (NR) %x3|Modified-Euler (ME) %x2|Forward-Euler (FE) %x1'
        dict['dschememnu']   = 'Dynamic Solver scheme %t|SS22 %x2|Newmark (GN22) %x1'
        dict['stagetypemnu'] = 'Solver type %t|Dynamics %x3|Transient %x2|Quasi-static/Steady %x1'
        dict['damptymnu']    = 'Damping type %t|Rayleigh %x2|None %x1'

        Blender.Registry.SetKey('MechSysDict', dict)
        print '[1;34mMechSys[0m: dictionary created'
    return dict

def set_key(key,value):
    Blender.Registry.GetKey('MechSysDict')[key] = value
    Blender.Window.QRedrawAll()

def toggle_key(key):
    Blender.Registry.GetKey('MechSysDict')[key] = not Blender.Registry.GetKey('MechSysDict')[key]
    Blender.Window.QRedrawAll()

def key(key):
    try: return Blender.Registry.GetKey('MechSysDict')[key]
    except: return {}


# ============================================================================== Objects and Meshes

def get_msh_obj(obj,with_error=True):
    msh_obj = None
    if obj.properties.has_key('msh_type'): msh_type = obj.properties['msh_type']
    else:
        if with_error: raise Exception('Please, generate mesh first or set "Frame mesh" toggle')
        else: return msh_obj
    if msh_type=='frame': msh_obj = obj
    else:
        if obj.properties.has_key('msh_name'): msh_obj = bpy.data.objects[obj.properties['msh_name']]
        else:
            if with_error: raise Exception('Please, generate Mesh first')
    return msh_obj

def get_obj(with_error=True):
    # Get current active object
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if with_error:
        if obj==None: raise Exception('Please, select an object first')
    return obj

def get_msh(with_error=True):
    # Get current active object and mesh (will exit edit mode)
    edm = Blender.Window.EditMode()
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if with_error:
        if obj==None:        raise Exception('Please, select an object first')
        if obj.type!='Mesh': raise Exception('Please, select an object of MESH type')
    if edm: Blender.Window.EditMode(0)
    msh = obj.getData(mesh=1)
    return edm, obj, msh

def get_selected_edges(msh):
    # Output: the 'really' selected edges (those that have the 'sel' property == 1)
    # Note: msh.edges.selected() == (selected + others) returns all possible selected
    #       edges, including those that can be defined by conecting the 'real' selected ones,
    #       since it returns 'edges for which BOTH vertices are selected'
    sel = []
    for eid in msh.edges.selected():
        if msh.edges[eid].sel==True: sel.append(eid)
    return sel


# ============================================================================== New Properties

def new_blk_props():
    edm, obj, msh = get_msh()
    eids          = get_selected_edges(msh)
    neds          = len(eids)
    is3d          = obj.properties['is3d'] if obj.properties.has_key('is3d') else False
    # check
    if is3d:
        if not (neds==12 or neds==24):
            if edm: Blender.Window.EditMode(1)
            raise Exception('To set a 3D block, 12 or 24 edges must be selected (%d is invalid)'%neds)
    else:
        if not (neds==4 or neds==8):
            if edm: Blender.Window.EditMode(1)
            raise Exception('To set a 2D block, 4 or 8 edges must be selected (%d is invalid)'%neds)
    if edm: Blender.Window.EditMode(1)
    # new block properties
    blk = [  -1, #   0:  block tag
             -1, #   1:  edge ID: X-axis
             -1, #   2:  edge ID: Y-axis
             -1, #   3:  edge ID: Z-axis
             -1, #   4:  vertex ID: Origin (local axes)
             -1, #   5:  vertex ID: X-plus (local axes)
             -1, #   6:  vertex ID: Y-plus (local axes)
             -1, #   7:  vertex ID: Z-plus (local axes)
              2, #   8:  nx: number of divisions along X-axis
              2, #   9:  ny: number of divisions along Y-axis
              2, #  10:  nz: number of divisions along Z-axis
            0.0, #  11:  ax coefficient
            0.0, #  12:  ay coefficient
            0.0, #  13:  az coefficient
              0, #  14:  linx: nonlinear X divisions ?
              0, #  15:  liny: nonlinear Y divisions ?
              0, #  16:  linz: nonlinear Z divisions ?
           neds, #  17:  number of edges
             -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 ] # edges IDs
    blk[18:18+neds] = eids
    return blk

def new_reg_props():
    x, y, z = Blender.Window.GetCursorPos()
    return [-1, -1.0, x,y,z] # 0:tag, 1:maxarea, 2:x, 3:y, 4:z

def new_hol_props():
    x, y, z = Blender.Window.GetCursorPos()
    return [x,y,z] # 0:x, 1:y, 2:z

def new_lin_props(): return [-200, 0, 1] # tag, N0, N1

def new_mat_props():
    return [      0,    #   0:  ModelType
                 -1,    #   1:  idx -- index to material name/description in 'texts'
             1000.0,    #   2:  E   -- Young (LinElastic)
                1.0,    #   3:  A   -- Beam: Area
                1.0,    #   4:  Izz -- Beam: Inertia
               0.25,    #   5:  nu  -- Poisson (LinElastic)
                  0,    #   6:  fc  -- Failure criterion (ElastoPlastic)
                0.0,    #   7:  sY  -- Yield stress (ElastoPlastic)
               -1.0,    #   8:  cu  -- Undrained cohesion (ElastoPlastic)
                0.01,   #   9:  lam -- lambda (CamClay)
                0.001,  #  10:  kap -- kappa (CamClay)
               25.0,    #  11:  phi -- shear angle at CS (CamClay)
                1.0 ]   #  12:  k   -- diffusion coefficient (LinFlow)

def new_stage_props(): return [1,        #  0: number
                               -1,       #  1: idx_desc(in texts)
                               0,        #  2: apply_body_forces?
                               0,        #  3: clear_disps?
                               1,        #  4: ndiv
                               1,        #  5: active?
                               0,        #  6: type: 0=equilib/steady, 1=transient, 2=dynamics
                               1.0,      #  7: tf        Final time
                               0.01,     #  8: dt        time increment
                               0.1,      #  9: dtOut     time increment for output
                               1.0e-7,   # 10: TolR      Tolerance for the norm of residual
                               False,    # 11: CalcWork  Calc work done == twice the stored (elastic) strain energy ? 
                               1,        # 12: Scheme    Scheme: FE_t (Forward-Euler), ME_t (Modified-Euler), NE_t (Newton-Rhapson)
                               1,        # 13: nSS       FE and NR: number of substeps 
                               1.0e-5,   # 14: STOL      ME: 
                               1.0,      # 15: dTini     ME: 
                               0.1,      # 16: mMin      ME: 
                               10.0,     # 17: mMax      ME: 
                               2000,     # 18: MaxSS     ME: 
                               False,    # 19: CteTg     Constant tangent matrices (linear problems) => K will be calculated once 
                               False,    # 20: ModNR     Modified Newton-Rhapson ? 
                               20,       # 21: MaxIt     Max iterations (for Newton-Rhapson) 
                               0,        # 22: TScheme   Transient scheme 
                               2./3.,    # 23: Theta     Transient scheme constant 
                               0,        # 24: DScheme   Dynamic scheme 
                               0,        # 25: DampTy    Damping type 
                               0.005,    # 26: DampAm    Rayleigh damping Am coefficient (C = Am*M + Ak*K) 
                               0.5,      # 27: DampAk    Rayleigh damping Ak coefficient (C = Am*M + Ak*K) 
                               0.5,      # 28: DynTh1    Dynamic coefficient Theta 1 
                               0.5,      # 29: DynTh2    Dynamic coefficient Theta 2 
                               True,     # 30: DynCont   Continue dynamic simulation (after DySolve was called once) 
                               True]     # 31: Output VTU file each dtOut ?

def new_nbry_props():  return [-100, 0, 0.0]           # tag, key, val
def new_ebry_props():  return [-10,  0, 0.0]           # tag, key, val
def new_fbry_props():  return [-10,  0, 0.0]           # tag, key, val
def new_eatt_props():  return [ -1,    # 0: tag
                                 0,    # 1: empty
                                 0,    # 2: geom/prob: Beam, Quad4, Hex8, ...
                                 0,    # 3: gty: psa, pse, ...
                                -1,    # 4: material ID
                                 1.0,  # 5: h (thickness)
                                 1,    # 6: active
                                 0,    # 7: activate?
                                 0 ]   # 8: deactivate?


# ============================================================================== Object Properties

def props_set_with_tag(key,subkey,tag,props):
    obj = get_obj()
    if not obj.properties.has_key(key): obj.properties[key] = {}
    if tag==0: obj.properties[key].pop    (subkey)
    else:      obj.properties[key].update({subkey:props})
    if len(obj.properties[key])==0: obj.properties.pop(key)

def props_del_all_tags(key):
    obj = get_obj()
    if obj.properties.has_key(key):
        msg = 'Confirm delete all tags ('+key+')?%t|Yes'
        res = Blender.Draw.PupMenu(msg)
        if res>0:
            obj.properties.pop(key)
            Blender.Window.QRedrawAll()

def props_push_new(key,props,check=False,ic=0,fc=0):
    # ic(inicomp),fc(fincomp): initial and final indexes for comparision between props elements,
    # in order to check whether the property was added or not
    obj = get_obj()
    if not obj.properties.has_key(key): obj.properties[key] = {}
    if check:
        for k, v in obj.properties[key].iteritems():
            eds = []
            for i in range(ic,fc): eds.append(v[i])
            if props[ic:fc]==eds: raise Exception('Property was already added')
    id = 0
    while obj.properties[key].has_key(str(id)): id += 1
    obj.properties[key][str(id)] = props
    Blender.Window.QRedrawAll()
    return id

def props_set_item(key,id,item,val):
    obj = get_obj()
    obj.properties[key][str(id)][item] = val
    Blender.Window.QRedrawAll()

def props_set_text(key,id,txt):
    obj = get_obj()
    obj.properties[key][str(id)] = txt
    Blender.Window.QRedrawAll()

def props_del(key,id):
    msg = 'Confirm delete this item?%t|Yes'
    res = Blender.Draw.PupMenu(msg)
    if res>0:
        obj = get_obj()
        obj.properties[key].pop(str(id))
        if len(obj.properties[key])==0: obj.properties.pop(key)
        Blender.Window.QRedrawAll()

def props_set_val(key, val):
    obj = get_obj()
    obj.properties[key] = val
    Blender.Window.QRedrawAll()

def props_del_all(key):
    obj = get_obj()
    if not obj.properties.has_key(key): return
    msg = 'Confirm delete ALL?%t|Yes'
    res  = Blender.Draw.PupMenu(msg)
    if res>0:
        obj.properties.pop(key)
        Blender.Window.QRedrawAll()

# ------------------------------------------------------------------------------ Materials

def props_push_new_mat():
    obj = get_obj()
    mid = props_push_new('mats', new_mat_props()) # returns material_id == mid
    tid = props_push_new('texts', 'material')     # returns text_id     == tid
    props_set_item('mats',mid,1,tid)              # description
    Blender.Window.QRedrawAll()

def props_del_all_mats():
    obj = get_obj()
    msg = 'Confirm delete ALL materials?%t|Yes'
    res  = Blender.Draw.PupMenu(msg)
    if res>0:
        # delete description
        for mid, v in obj.properties['mats'].iteritems():
            tid = str(int(v[1])) # text_id
            obj.properties['texts'].pop(tid)
        if len(obj.properties['texts'])==0: obj.properties.pop('texts')
        # delete all materials
        obj.properties.pop('mats')
        Blender.Window.QRedrawAll()

# ------------------------------------------------------------------------------ Stages

def props_push_new_stage():
    obj = get_obj()
    sid = props_push_new('stages', new_stage_props()) # returns stage_id == sid
    tid = props_push_new('texts', 'simulation stage') # returns text_id  == tid
    stg = 'stg_'+str(sid)
    nst = len(obj.properties['stages']) # number of stages
    props_set_item('stages',sid,0,nst)  # num
    props_set_item('stages',sid,1,tid)  # description
    obj.properties[stg] = {}
    set_key('fem_stage', sid)
    if nst>1:
        # copy eatts from first stage
        for k, v in obj.properties['stages'].iteritems():
            if int(v[0])==1: # first stage
                obj.properties[stg] = obj.properties['stg_'+k]
                break
    else: obj.properties['pty'] = 0
    Blender.Window.QRedrawAll()

def props_push_new_fem(stage_ids,key,props):
    obj = get_obj()
    for sid in stage_ids:
        stg = 'stg_'+str(sid)
        if not obj.properties[stg].has_key(key): obj.properties[stg][key] = {}
        id = 0
        while obj.properties[stg][key].has_key(str(id)): id += 1
        obj.properties[stg][key][str(id)] = props
    Blender.Window.QRedrawAll()

def props_set_fem(stage_ids,key,id,item,val):
    obj = get_obj()
    for sid in stage_ids:
        stg = 'stg_'+str(sid)
        obj.properties[stg][key][str(id)][item] = val
    Blender.Window.QRedrawAll()

def props_set_fem_all_stg(key,id,item,val):
    obj  = get_obj()
    sids = [int(k) for k, v in obj.properties['stages'].iteritems()]
    props_set_fem (sids, key, id, item, val)

def props_del_fem(stage_ids,key,id):
    msg = 'Confirm delete this item?%t|Yes'
    res = Blender.Draw.PupMenu(msg)
    if res>0:
        for sid in stage_ids:
            stg = 'stg_'+str(sid)
            obj = get_obj()
            obj.properties[stg][key].pop(str(id))
            if len(obj.properties[stg][key])==0: obj.properties[stg].pop(key)
        Blender.Window.QRedrawAll()

def props_del_fem_all_stg(key,id):
    obj  = get_obj()
    sids = [int(k) for k, v in obj.properties['stages'].iteritems()]
    props_del_fem(sids, key, id)

def props_del_all_fem(stage_ids,key,with_confirmation=True):
    obj = get_obj()
    if with_confirmation:
        msg = 'Confirm delete ALL?%t|Yes'
        res  = Blender.Draw.PupMenu(msg)
    else: res = 1
    if res>0:
        for sid in stage_ids:
            stg = 'stg_'+str(sid)
            if obj.properties[stg].has_key(key):
                obj.properties[stg].pop(key)
        Blender.Window.QRedrawAll()

def props_del_stage():
    msg = 'Confirm delete this stage?%t|Yes'
    res = Blender.Draw.PupMenu(msg)
    if res>0:
        obj = get_obj()
        sid = str(Blender.Registry.GetKey('MechSysDict')['fem_stage'])
        tid = str(int(obj.properties['stages'][sid][1]))
        stg = 'stg_'+sid
        # delete description
        obj.properties['texts'].pop(tid)
        # delete stage
        num = int(obj.properties['stages'][sid][0])
        obj.properties['stages'].pop(sid)
        obj.properties.pop(stg)
        # reset other stages numbers and set current fem_stage
        for k, v in obj.properties['stages'].iteritems():
            other_num = int(obj.properties['stages'][k][0])
            if other_num>num:
                obj.properties['stages'][k][0] = other_num-1
            if int(obj.properties['stages'][k][0])==1: set_key('fem_stage', int(k))
        # delete 'stages' and 'texts'
        if len(obj.properties['stages'])==0: obj.properties.pop('stages')
        if len(obj.properties['texts'] )==0: obj.properties.pop('texts')
        # set current fem stage
        Blender.Window.QRedrawAll()

def props_del_all_stages(with_confirmation=True):
    obj = get_obj()
    if not obj.properties.has_key('stages'): return
    if with_confirmation:
        msg = 'Confirm delete ALL stages?%t|Yes'
        res = Blender.Draw.PupMenu(msg)
    else: res = 1
    if res>0:
        for sid, v in obj.properties['stages'].iteritems():
            # delete description
            obj.properties['texts'].pop(str(int(v[1]))) # v[1] == text ID
            # delete stage
            stg = 'stg_'+sid
            obj.properties.pop(stg)
        if len(obj.properties['texts'])==0: obj.properties.pop('texts')
        obj.properties.pop('stages')
        set_key('fem_stage', 0)
        Blender.Window.QRedrawAll()

def find_stage(obj,num):
    for k, v in obj.properties['stages'].iteritems():
        if int(v[0])==num:
            sid = k
            stg = 'stg_'+k
            return sid, stg

# ------------------------------------------------------------------------------ Blocks

def blk_set_local_system(obj,msh,id):
    # set local system: origin and vertices on positive directions
    key = str(id)
    iex = int(obj.properties['blks'][key][1])
    iey = int(obj.properties['blks'][key][2])
    if iex>=0 and iey>=0:
        ex     = msh.edges[iex]
        ey     = msh.edges[iey]
        origin = -1
        x_plus = -1
        y_plus = -1
        if ex.v1.index==ey.v1.index:
            origin = ex.v1.index
            x_plus = ex.v2.index
            y_plus = ey.v2.index
        elif ex.v1.index==ey.v2.index:
            origin = ex.v1.index
            x_plus = ex.v2.index
            y_plus = ey.v1.index
        elif ex.v2.index==ey.v1.index:
            origin = ex.v2.index
            x_plus = ex.v1.index
            y_plus = ey.v2.index
        elif ex.v2.index==ey.v2.index:
            origin = ex.v2.index
            x_plus = ex.v1.index
            y_plus = ey.v1.index
        obj.properties['blks'][key][4] = origin
        obj.properties['blks'][key][5] = x_plus
        obj.properties['blks'][key][6] = y_plus
        iez = int(obj.properties['blks'][key][3])
        if iez>=0:
            z_plus = -1
            ez = msh.edges[iez]
            if ez.v1.index==origin:
                z_plus = ez.v2.index
            elif ez.v2.index==origin:
                z_plus = ez.v1.index
            obj.properties['blks'][key][7] = z_plus

def blk_set_axis(id,item):
    edm, obj, msh = get_msh()
    sel           = get_selected_edges(msh)
    if len(sel)==1:
        key  = str(id)
        neds = int(obj.properties['blks'][key][17])
        eid  = sel[0]
        eds  = []
        for i in range(18,18+neds): eds.append(int(obj.properties['blks'][key][i]))
        if eid in eds:
            obj.properties['blks'][key][item] = eid
            blk_set_local_system(obj,msh,id)
        else:
            if edm: Blender.Window.EditMode(1)
            raise Exception('This edge # '+str(eid)+' is not part of this block')
        Blender.Window.QRedrawAll()
    else:
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, select one edge (only) in order to define the local axis')
    if edm: Blender.Window.EditMode(1)


# ====================================================================================== Util

def sarray_set_val(strarr,item,value):
    # Set value of item item in a StringArray of length len
    # Ex.: Input: strarr = 'A BB CCC DD_DD' (len==4)
    #             item   = 2
    #             value  = 'XX'
    # Output: 'A BB XX DD_DD'
    d       = strarr.split()
    d[item] = value
    for i, v in enumerate(d):
        if i==0: res  = v
        else:    res += ' '+v
    return res

def rgb2html(rgb):
    # convert (R, G, B) tuple to a html value RRGGBB
    return '%02x%02x%02x' % (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255))

def html2rgb(html):
    # convert RRGGBB or a (R, G, B) tuple
    if len(html)!=6: raise Exception('html2rgb: %s is not in RRGGBB format' % html)
    return tuple(float(int(html[i:i+2], 16)/255.0) for i in range(0, 6, 2))

def rgb2hex(rgb):
    # convert (R, G, B) tuple to a hex value 0x000000
    return int('0x%02x%02x%02x' % (int(rgb[0]*255), int(rgb[1]*255), int(rgb[2]*255)), 16)

def hex2rgb(hex):
    # convert a hex value 0x000000 to a (R, G, B) tuple
    return html2rgb('%06x' % hex)


#  2+       r_idx    = 0        # start right vertex index
#   |\      edge_ids = [0,1,2]  # indexes of all edges to search for r_idx
#  0| \1    v1_ids   = [1,0,0]  # v1 vertex indexes
#   |  \    v2_ids   = [2,2,1]  # v2 vertex indexes
#  1+---+0  eds      = [1,0,2]  # result: edge indexes
#     2     ids      = [2,1,0]  # result: vertex indexes
def sort_edges_and_verts(msh, edge_ids, r_idx):
    # loop over all connected edges
    eds = []
    ids = []
    while len(edge_ids)>0:
        v1_ids = [msh.edges[ie].v1.index for ie in edge_ids] # ie => index of an edge
        v2_ids = [msh.edges[ie].v2.index for ie in edge_ids]
        if r_idx in v1_ids:
            idx   = v1_ids.index(r_idx)
            r_idx = v2_ids[idx]
            eds.append   (edge_ids[idx])
            ids.append   (r_idx)
            edge_ids.pop (idx)
        elif r_idx in v2_ids:
            idx   = v2_ids.index(r_idx)
            r_idx = v1_ids[idx]
            eds.append   (edge_ids[idx])
            ids.append   (r_idx)
            edge_ids.pop (idx)
        else: break
    return eds, ids

def bounding_box(msh):
    minx, miny, minz = msh.verts[0].co[0], msh.verts[0].co[1], msh.verts[0].co[2]
    maxx, maxy, maxz = minx, miny, minz
    for v in msh.verts:
        if v.co[0]<minx: minx = v.co[0]
        if v.co[1]<miny: miny = v.co[1]
        if v.co[2]<minz: minz = v.co[2]
        if v.co[0]>maxx: maxx = v.co[0]
        if v.co[1]>maxy: maxy = v.co[1]
        if v.co[2]>maxz: maxz = v.co[2]
    p0 = Blender.Mathutils.Vector([minx,miny,minz])
    p1 = Blender.Mathutils.Vector([maxx,maxy,maxz])
    return p0, p1, (p1-p0).length

def find_node(obj):
    x, y, z = Blender.Window.GetCursorPos()
    pt  = Blender.Mathutils.Vector([x,y,z])
    msh = obj.getData (mesh=1)
    ori = [v for v in msh.verts] # create a copy before transforming to global coordinates
    msh.transform (obj.matrix) # transform to global coordinates
    nid = None # id of the vertex close to x,y,z
    p0, p1, maxlen = bounding_box(msh)
    for v in msh.verts:
        if (pt-v.co).length<0.05*maxlen:
            nid = v.index
            break
    msh.verts = ori # restore mesh to local coordinates
    if nid==None: raise Exception('Could not find any Node close to Cursor Position: x=%g, y=%g, z=%g'%(x,y,z))
    return nid


######################################################################## Shortest path

def find_shortest_path(graph, start, end, path=[]):
    path = path + [start]
    if start==end: return path
    if not graph.has_key(start): return None
    shortest = None
    for node in graph[start]:
        if node not in path:
            newpath = find_shortest_path (graph, node, end, path)
            if newpath:
                if not shortest or len(newpath)<len(shortest):
                    shortest = newpath
    return shortest


##################################################################### Block local IDs

def block_local_ids(obj, blk_dat):
    neds                     = int(blk_dat[17]) # number of edges
    xeid, yeid, zeid         = int(blk_dat[1]), int(blk_dat[2]), int(blk_dat[3])
    origin, xvid, yvid, zvid = int(blk_dat[4]), int(blk_dat[5]), int(blk_dat[6]), int(blk_dat[7])
    if xeid<0 or yeid<0 or origin<0 or xvid<0 or yvid<0: return None

    # build connectivity map: vertex id => edge ids
    blk_eids = []
    for i in range(neds): blk_eids.append (int(blk_dat[18+i]))
    msh = obj.getData (mesh=1)
    edges_to_remove = [xeid,yeid,zeid]
    vid2eids = {}
    for eid, ed in enumerate(msh.edges):
        if eid in blk_eids and not eid in edges_to_remove:
            v0, v1 = ed.key[0], ed.key[1]
            if v0 in vid2eids: vid2eids[v0].append (eid)
            else:              vid2eids[v0] = [eid]
            if v1 in vid2eids: vid2eids[v1].append (eid)
            else:              vid2eids[v1] = [eid]

    # find neighbours
    neigh = {} # map: vertex id => ids of vertex neighbours
    for vid, eids in vid2eids.iteritems():
        neigh[vid] = []
        for eid in eids:
            v0, v1 = msh.edges[eid].key[0], msh.edges[eid].key[1]
            if not v0==vid: neigh[vid].append (v0)
            if not v1==vid: neigh[vid].append (v1)

    # find local IDs
    loc_vids = {}
    if zeid<0 and zvid<0: # 2D

        if neds==4: # quad4
            neigh.pop(xvid)
            neigh.pop(yvid)
            v2 = neigh.keys()[0]
            loc_vids[origin] = 0
            loc_vids[xvid]   = 1
            loc_vids[v2]     = 2
            loc_vids[yvid]   = 3
            return loc_vids

        elif neds==8: # quad8
            loc_vids[origin] = 0
            loc_vids[xvid]   = 4
            loc_vids[yvid]   = 7
            res = find_shortest_path (neigh, xvid, yvid)
            loc_vids[res[1]] = 1
            loc_vids[res[2]] = 5
            loc_vids[res[3]] = 2
            loc_vids[res[4]] = 6
            loc_vids[res[5]] = 3
            return loc_vids

    else: # 3D

        if neds==12: # hex8
            v2 = find_shortest_path (neigh, xvid, yvid)[1]
            v5 = find_shortest_path (neigh, xvid, zvid)[1]
            v7 = find_shortest_path (neigh, yvid, zvid)[1]
            v6 = find_shortest_path (neigh, v5, v7)[1]
            loc_vids[origin]= 0
            loc_vids[xvid]  = 1
            loc_vids[v2]    = 2
            loc_vids[yvid]  = 3
            loc_vids[zvid]  = 4
            loc_vids[v5]    = 5
            loc_vids[v6]    = 6
            loc_vids[v7]    = 7
            return loc_vids

        elif neds==24: # hex20
            res = find_shortest_path (neigh, xvid, yvid)
            v2  = res[3]
            loc_vids[origin] = 0
            loc_vids[xvid]   = 8
            loc_vids[yvid]   = 11
            loc_vids[zvid]   = 16
            loc_vids[res[1]] = 1
            loc_vids[res[2]] = 9
            loc_vids[res[3]] = 2
            loc_vids[res[4]] = 10
            loc_vids[res[5]] = 3
            res = find_shortest_path (neigh, xvid, zvid)
            v5  = res[3]
            v4  = res[5]
            loc_vids[res[2]] = 17
            loc_vids[res[3]] = 5
            loc_vids[res[4]] = 12
            loc_vids[res[5]] = 4
            res = find_shortest_path (neigh, yvid, zvid)
            v7  = res[3]
            loc_vids[res[2]] = 19
            loc_vids[res[3]] = 7
            loc_vids[res[4]] = 15
            neigh.pop(v4) # remove v4
            res = find_shortest_path (neigh, v5, v7)
            v6  = res[2]
            loc_vids[res[1]] = 13
            loc_vids[res[2]] = 6
            loc_vids[res[3]] = 14
            res = find_shortest_path (neigh, v6, v2)
            loc_vids[res[1]] = 18
            return loc_vids
