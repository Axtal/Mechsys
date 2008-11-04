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
        dict['gui_show_set']  = True
        dict['gui_show_cad']  = True
        dict['gui_show_mesh'] = True
        dict['gui_show_fem']  = True
        dict['gui_show_res']  = True
        dict['gui_inirow']    = 0
        # SETTINGS
        dict['show_props']    = False
        dict['show_e_ids']    = True
        dict['show_v_ids']    = True
        dict['show_blks']     = True
        dict['show_axes']     = True
        dict['show_etags']    = True
        dict['show_ftags']    = True
        dict['show_elems']    = True
        dict['show_regs']     = True
        dict['show_opac']     = 0.3
        # CAD
        dict['cad_x']         = '0.0'
        dict['cad_y']         = '0.0'
        dict['cad_z']         = '0.0'
        dict['cad_rad']       = '0.0'
        dict['cad_stp']       = 10
        # MESH
        dict['newetag']       = [-10, 0]         # tag, type
        dict['newftag']       = [-100, 0x000080] # tag, colour
        # FEM
        dict['fullsc']        = False # generate full script (for FEA)
        # RESULTS
        dict['show_res']        = False
        dict['res_dfv']         = 0
        dict['res_show_scalar'] = False
        dict['res_warp_scale']  = '10'
        dict['res_show_warp']   = 1

        # DOF Vars
        dict['dfv']    = { 0:'ux', 1:'uy', 2:'uz', 3:'fx', 4:'fy', 5:'fz', 6:'u', 7:'q' }
        dict['dfvmnu'] = 'DOF Vars %t|q %x8|u %x7|fz %x6|fy %x5|fx %x4|uz %x3|uy %x2|ux %x1'

        # Element types
        dict['ety']  = {  0:'Hex8Equilib',    1:'Hex8Diffusion',
                          2:'Quad4PStrain',   3:'Quad4PStress',   4:'Quad4Diffusion',
                          5:'Tri3PStrain',    6:'Tri3PStress',    7:'Tri3Diffusion',
                          8:'Rod',
                          9:'Quad8PStrain',  10:'Quad8PStress',  11:'Quad8Diffusion',
                         12:'Tri6PStrain',   13:'Tri6PStress',   14:'Tri6Diffusion',
                         15:'Hex20Equilib',  16:'Hex20Diffusion',
                         17:'Tet4Equilib',   18:'Tet4Diffusion',
                         19:'Tet10Equilib',  20:'Tet10Diffusion' }

        dict['etymnu'] = 'Element Types %t|Rod %x9|Tri3Diffusion %x8|Tri3PStress %x7|Tri3PStrain %x6|Quad4Diffusion %x5|Quad4PStress %x4|Quad4PStrain %x3|Hex8Diffusion %x2|Hex8Equilib %x1'

        # Models
        dict['mdl']    = { 0:'LinElastic', 1:'LinDiffusion' }
        dict['mdlmnu'] = 'Constitutive Models %t|LinDiffusion %x2|LinElastic %x1'

        # VTK Cell Type (tentative mapping)
        dict['vtk2ety'] = {  5: 5,   # VTK_TRIANGLE             => Tri3PStrain
                             9: 2,   # VTK_QUAD                 => Quad4PStrain
                            10:17,   # VTK_TETRA                => Tet4Equilib
                            12: 0,   # VTK_HEXAHEDRON           => Hex8Equilib
                            22:12,   # VTK_QUADRATIC_TRIANGLE   => Tri6PStrain
                            23: 9,   # VTK_QUADRATIC_QUAD       => Quad8PStrain
                            24:19,   # VTK_QUADRATIC_TETRA      => Tet10Equilib
                            25:15 }  # VTK_QUADRATIC_HEXAHEDRON => Hex20Equilib

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
    return Blender.Registry.GetKey('MechSysDict')[key]


# ============================================================================== Objects and Meshes

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
    # check
    if obj.properties['3dmesh']:
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

def new_mat_props():
    return [   200.0,     #   0:  E -- Young
                 0.2,     #   1:  nu -- Poisson
                 0.0891,  #   2:  lam -- lambda
                 0.0196,  #   3:  kap -- kappa
                31.6,     #   4:  phics -- shear angle at CS
             18130.0 ]    #   5:  G -- shear modulus (kPa)

def new_ini_props():
    #                                       0  1  2    3   4   5
    return [ 0.0,0.0,0.0, 0.0,0.0,0.0,  #  Sx,Sy,Sz, Sxy,Syz,Szx
                              1.6910 ]  #  6:  v -- initial specific volume

def new_nbry_props(): return [0.0,0.0,0.0, 0, 0.0]             # x,y,z, ux, val
def new_nbID_props(): return [0, 0, 0.0]                       # ID, ux, val
def new_ebry_props(): return [-10, 0, 0.0]                     # tag, ux, val
def new_fbry_props(): return [-100, 0, 0.0, key('newftag')[1]] # tag, ux, val, colour
def new_eatt_props(): return [-1, 2, 0, -1, -1]                # tag ElemType Model MatID IniID


# ============================================================================== Object Properties

def props_set_with_tag(key,subkey,tag,props):
    obj = get_obj()
    if not obj.properties.has_key(key): obj.properties[key] = {}
    if tag==0: obj.properties[key].pop    (subkey)
    else:      obj.properties[key].update({subkey:props})
    if len(obj.properties[key])==0: obj.properties.pop(key)

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
    while obj.properties[key].has_key(str(id)):
        id += 1
    obj.properties[key][str(id)] = props
    Blender.Window.QRedrawAll()

def props_set_item(key,id,item,val):
    obj = get_obj()
    obj.properties[key][str(id)][item] = val
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
    msg = 'Confirm delete ALL?%t|Yes'
    res  = Blender.Draw.PupMenu(msg)
    if res>0:
        obj.properties.pop(key)
        Blender.Window.QRedrawAll()

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
        raise Exception('Please, select (only) one edge to define the local axis')
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






def get_cg(msh, vids, vtkcelltype):
    # In:   vids = verts_ids
    x, y, z = 0, 0, 0
    if   vtkcelltype== 5: # VTK_TRIANGLE
        x = (msh.verts[vids[0]].co[0]+msh.verts[vids[1]].co[0]+msh.verts[vids[2]].co[0])/3.0
        y = (msh.verts[vids[0]].co[1]+msh.verts[vids[1]].co[1]+msh.verts[vids[2]].co[1])/3.0
        z = (msh.verts[vids[0]].co[2]+msh.verts[vids[1]].co[2]+msh.verts[vids[2]].co[2])/3.0
    elif vtkcelltype== 9: # VTK_QUAD
        x = (msh.verts[vids[0]].co[0]+msh.verts[vids[2]].co[0])/2.0
        y = (msh.verts[vids[0]].co[1]+msh.verts[vids[2]].co[1])/2.0
        z = (msh.verts[vids[0]].co[2]+msh.verts[vids[2]].co[2])/2.0
    elif vtkcelltype==10: # VTK_TETRA
        pass
    elif vtkcelltype==12: # VTK_HEXAHEDRON
        x = (msh.verts[vids[0]].co[0]+msh.verts[vids[6]].co[0])/2.0
        y = (msh.verts[vids[0]].co[1]+msh.verts[vids[6]].co[1])/2.0
        z = (msh.verts[vids[0]].co[2]+msh.verts[vids[6]].co[2])/2.0
    elif vtkcelltype==22: # VTK_QUADRATIC_TRIANGLE
        pass
    elif vtkcelltype==23: # VTK_QUADRATIC_QUAD
        pass
    elif vtkcelltype==24: # VTK_QUADRATIC_TETRA
        pass
    elif vtkcelltype==25: # VTK_QUADRATIC_HEXAHEDRON
        pass
    return x, y, z


def get_poly_area(msh, vids):
    a = 0.0
    n = len(vids)
    for i in range(n):
        j  = (i+1) % n
        a += msh.verts[vids[i]].co[0] * msh.verts[vids[j]].co[1] - msh.verts[vids[j]].co[0] * msh.verts[vids[i]].co[1]
    return a/2.0


def get_poly_cg(msh, vids):
    x, y, z = 0, 0, 0
    n = len(vids)
    for i in range(n):
        j  = (i+1) % n
        p  =  msh.verts[vids[i]].co[0] * msh.verts[vids[j]].co[1] - msh.verts[vids[j]].co[0] * msh.verts[vids[i]].co[1]
        x += (msh.verts[vids[i]].co[0] + msh.verts[vids[j]].co[0]) * p 
        y += (msh.verts[vids[i]].co[1] + msh.verts[vids[j]].co[1]) * p 
        z +=  msh.verts[vids[i]].co[2]
    area = get_poly_area (msh, vids)
    x /= 6.0*area
    y /= 6.0*area
    z /= n
    return x, y, z


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


def get_file_key():
    bfn = Blender.sys.expandpath (Blender.Get('filename'))
    key = Blender.sys.basename   (Blender.sys.splitext(bfn)[0])
    return key
