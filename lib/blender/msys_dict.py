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
    if not dict:
        dict               = {}
        dict['gui_inirow'] = 0
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
        dict['newetag']       = [-10, 0]
        dict['newftag']       = [-100, 0x000080]
        # FEM
        dict['fem_fullsc']    = False
        dict['fem_struct']    = False
        dict['fem_unstru']    = True
        # RESULTS
        dict['res_show']        = False
        dict['res_scalar']      = 0
        dict['res_show_scalar'] = False
        dict['res_warp_scale']  = '10'
        dict['res_show_warp']   = 1

        dict['dofvars'] = { 0:'ux', 1:'uy', 2:'uz', 3:'fx', 4:'fy', 5:'fz', 6:'u', 7:'q' }
        dict['dofvars_menu'] = 'DOF Vars %t|q %x8|u %x7|fz %x6|fy %x5|fx %x4|uz %x3|uy %x2|ux %x1'

        dict['etypes'] = { 0:'Hex8Equilib', 1:'Hex8Diffusion', 2:'Quad4PStrain', 3:'Quad4PStress', 4:'Quad4Diffusion', 5:'Tri3PStrain', 6:'Tri3PStress', 7:'Tri3Diffusion', 8:'Rod' }
        dict['etypes_menu'] = 'Element Types %t|Rod %x9|Tri3Diffusion %x8|Tri3PStress %x7|Tri3PStrain %x6|Quad4Diffusion %x5|Quad4PStress %x4|Quad4PStrain %x3|Hex8Diffusion %x2|Hex8Equilib %x1'

        dict['models'] = { 0:'LinElastic', 1:'LinDiffusion' }
        dict['models_menu'] = 'Constitutive Models %t|LinDiffusion %x2|LinElastic %x1'

        Blender.Registry.SetKey('MechSysDict', dict)
        print '[1;34mMechSys[0m: dictionary created'
    return dict

def set_key(key,value):
    Blender.Registry.GetKey('MechSysDict')[key] = value

def set_key_and_redraw(key,value):
    Blender.Registry.GetKey('MechSysDict')[key] = value
    Blender.Window.QRedrawAll()

def toggle_key_and_redraw(key):
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
