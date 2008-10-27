# Modules
import Blender
import bpy
import string


# ================================================================================ Dictionary

def load_dict():
    dict = Blender.Registry.GetKey('MechSysDict')
    if not dict:
        dict                  = {}
        dict['gui_show_set']  = 1
        dict['gui_show_cad']  = 1
        dict['gui_show_mesh'] = 1
        dict['gui_show_fem']  = 1
        dict['gui_show_res']  = 1
        dict['inirow']        = 0
        dict['newpoint_x']    = '0.0'
        dict['newpoint_y']    = '0.0'
        dict['newpoint_z']    = '0.0'
        dict['newetag']       = [-10, 0]
        dict['newftag']       = [-100, 0x000080]
        dict['fillet_radius'] = '0.0'
        dict['fillet_steps']  = 10
        dict['show_props']    = 0
        dict['show_v_ids']    = 1
        dict['show_e_ids']    = 1
        dict['show_f_ids']    = 0
        dict['show_etags']    = 1
        dict['show_ftags']    = 1
        dict['show_tags_txt'] = 1
        dict['ftags_opac']    = 0.3
        dict['show_elems']    = 1
        dict['show_axes']     = 1
        dict['show_results']  = 0
        dict['show_scalar']   = 0
        dict['show_warp']     = 1
        dict['scalar_key']    = 'uy'
        dict['warp_scale']    = '10'

        dict['dofvars'] = { 0:'ux', 1:'uy', 2:'uz', 3:'fx', 4:'fy', 5:'fz', 6:'u', 7:'q' }
        dict['dofvars_menu'] = 'DOF Vars %t|q %x8|u %x7|fz %x6|fy %x5|fx %x4|uz %x3|uy %x2|ux %x1'

        dict['etypes'] = { 0:'Hex8Equilib', 1:'Hex8Diffusion', 2:'Quad4PStrain', 3:'Quad4PStress', 4:'Quad4Diffusion', 5:'Tri3PStrain', 6:'Tri3PStress', 7:'Tri3Diffusion', 8:'Rod' }
        dict['etypes_menu'] = 'Element Types %t|Rod %x9|Tri3Diffusion %x8|Tri3PStress %x7|Tri3PStrain %x6|Quad4Diffusion %x5|Quad4PStress %x4|Quad4PStrain %x3|Hex8Diffusion %x2|Hex8Equilib %x1'

        dict['models'] = { 0:'LinElastic', 1:'LinDiffusion' }
        dict['models_menu'] = 'Constitutive Models %t|LinDiffusion %x2|LinElastic %x1'

        Blender.Registry.SetKey('MechSysDict', dict)
        print '[1;34mMechSys[0m: dictionary created'
    return dict


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


def get_key():
    bfn = Blender.sys.expandpath (Blender.Get('filename'))
    key = Blender.sys.basename   (Blender.sys.splitext(bfn)[0])
    return key


# Get current active object
def get_obj():
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    return obj


# Get current active mesh (will exit edit mode)
# Returns:
#          edm: edit mode == True or False
#          obj: current object
#          msh: current mesh
def get_msh(with_error=True):
    edm = Blender.Window.EditMode()
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        return edm, obj, msh
    else: raise Exception('The selected object must be of MESH type')


# ================================================================================ Properties

def set_int_property(obj, key, value):
    try:
        prop      = obj.getProperty (key)
        prop.data = value
    except:
        obj.addProperty (key, value, 'INT')

def set_str_property(obj, key, value):
    try:
        prop      = obj.getProperty (key)
        prop.data = value
    except:
        obj.addProperty (key, value, 'STRING')

def set_float_property(obj, key, value):
    try:
        prop      = obj.getProperty (key)
        prop.data = value
    except:
        obj.addProperty (key, value, 'FLOAT')


# ====================================================== Local system and number of divisions

def set_local_axis(obj, key, id):
    # In:
    #      key = 'x', 'y', or 'z'
    #      id  = ID of the edge corresponding to the local x, y and z axes
    set_int_property (obj, key+'_axis', id)

    # set local system: origin and vertices on positive directions
    ix = get_local_axis(obj, 'x')
    iy = get_local_axis(obj, 'y')
    if ix>-1 and iy>-1:
        msh = obj.getData(mesh=1)
        ex  = msh.edges[ix]
        ey  = msh.edges[iy]
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
        else:
            raise Exception('local x-y axes must share the same origin vertex (obj=%s)' % obj.name)
        set_int_property (obj, 'origin', origin)
        set_int_property (obj, 'x_plus', x_plus)
        set_int_property (obj, 'y_plus', y_plus)
        iz = get_local_axis(obj, 'z')
        if iz>-1:
            ez = msh.edges[iz]
            if ez.v1.index==origin:
                z_plus = ez.v2.index
            elif ez.v2.index==origin:
                z_plus = ez.v1.index
            else:
                raise Exception('local x-y-z axes must share the same origin vertex (obj=%s)' % obj.name)
            set_int_property (obj, 'z_plus', z_plus)

def get_local_axis(obj, key):
    # In:
    #      key = 'x', 'y', or 'z'
    # Out:
    #      ID of the edge corresponding to the local x, y and z axes
    try:    id = obj.getProperty(key+'_axis').data
    except: id = -1
    return id

def get_local_system(obj):
    try:    origin = obj.getProperty('origin').data
    except: origin = -1
    try:    x_plus = obj.getProperty('x_plus').data
    except: x_plus = -1
    try:    y_plus = obj.getProperty('y_plus').data
    except: y_plus = -1
    try:    z_plus = obj.getProperty('z_plus').data
    except: z_plus = -1
    return origin, x_plus, y_plus, z_plus


def set_ndiv(obj, key, ndiv):
    # In:
    #      key  = 'x', 'y', or 'z'
    #      ndiv = number of divisions along 'x', 'y', or 'z'
    set_int_property (obj, 'ndiv_'+key, ndiv)

def get_ndiv(obj, key):
    # In:
    #      key = 'x', 'y', or 'z'
    # Out:
    #      ndiv = number of divisions along 'x', 'y', or 'z'
    try:    ndiv = obj.getProperty('ndiv_'+key).data
    except: ndiv = -1
    return ndiv


def set_acoef(obj, key, acoef):
    # In:
    #      key  = 'x', 'y', or 'z'
    #      acoef = a coefficient for the divisions function along 'x', 'y', or 'z'
    set_str_property (obj, 'acoef_'+key, acoef)

def get_acoef(obj, key):
    # In:
    #      key = 'x', 'y', or 'z'
    # Out:
    #      acoef = a coefficient for the divisions function along 'x', 'y', or 'z'
    try:    acoef = obj.getProperty('acoef_'+key).data
    except: acoef = '0'
    return acoef


def set_nonlin(obj, key, nonlin):
    # In:
    #      key    = 'x', 'y', or 'z'
    #      nonlin = 0 or 1
    set_int_property (obj, 'nonlin_'+key, nonlin)

def get_nonlin(obj, key):
    # In:
    #      key = 'x', 'y', or 'z'
    # Out:
    #      nonlin = 0 or 1
    try:    nonlin = obj.getProperty('nonlin_'+key).data
    except: nonlin = 0
    return nonlin


# ================================================================================ Block tags

def set_btag(obj, tag):
    set_int_property (obj, 'btag', tag)

def get_btag(obj):
    try: return obj.getProperty('btag').data
    except:
        set_btag(obj, -1)
        return -1


# ====================================================================== Minangle and maxarea

def set_minangle(obj, value):
    set_str_property (obj, 'minangle', value)

def get_minangle(obj):
    try:    value = obj.getProperty('minangle').data
    except: value = '-1.0'
    return value

def set_maxarea(obj, value):
    set_str_property (obj, 'maxarea', value)

def get_maxarea(obj):
    try:    value = obj.getProperty('maxarea').data
    except: value = '-1.0'
    return value


# =================================================================================== Regions

def set_reg(obj, id, tag, maxarea, x, y, z):
    try:
        prop      = obj.getProperty ('reg_'+str(id))
        prop.data = str(tag)+' '+maxarea+' '+x+' '+y+' '+z
    except:
        obj.addProperty ('reg_'+str(id), str(tag)+' '+maxarea+' '+x+' '+y+' '+z, 'STRING')

def set_reg_tag(obj, id, tag): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = str(tag)+' '+d[1]+' '+d[2]+' '+d[3]+' '+d[4]

def set_reg_maxarea(obj, id, maxarea): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+maxarea+' '+d[2]+' '+d[3]+' '+d[4]

def set_reg_x(obj, id, x): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+x+' '+d[3]+' '+d[4]

def set_reg_y(obj, id, y): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+d[2]+' '+y+' '+d[4]

def set_reg_z(obj, id, z): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+d[2]+' '+d[3]+' '+z

def get_regs(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:3]=='reg':
            d = p.data.split()
            res.append ([d[0], d[1], d[2], d[3], d[4]])
    return res 

def del_all_regs(obj):
    for p in obj.getAllProperties():
        if p.name[:3]=='reg': obj.removeProperty(p)

def del_reg(obj, id):
    regs = get_regs (obj)
    del_all_regs (obj)
    k = 0
    for i in range(len(regs)):
        if i!=id:
            set_reg (obj, k, regs[i][0], regs[i][1], regs[i][2], regs[i][3], regs[i][4])
            k += 1


# ===================================================================================== Holes

def set_hol(obj, id, x, y, z):
    try:
        prop      = obj.getProperty ('hol_'+str(id))
        prop.data = x+' '+y+' '+z
    except:
        obj.addProperty ('hol_'+str(id), x+' '+y+' '+z, 'STRING')

def set_hol_x(obj, id, x): # property must exist
    p = obj.getProperty ('hol_'+str(id))
    d = p.data.split()
    p.data = x+' '+d[1]+' '+d[2]

def set_hol_y(obj, id, y): # property must exist
    p = obj.getProperty ('hol_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+y+' '+d[2]

def set_hol_z(obj, id, z): # property must exist
    p = obj.getProperty ('hol_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+z

def get_hols(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:3]=='hol':
            d = p.data.split()
            res.append ([d[0], d[1], d[2]])
    return res 

def del_all_hols(obj):
    for p in obj.getAllProperties():
        if p.name[:3]=='hol': obj.removeProperty(p)

def del_hol(obj, id):
    hols = get_hols (obj)
    del_all_hols (obj)
    k = 0
    for i in range(len(hols)):
        if i!=id:
            set_hol (obj, k, hols[i][0], hols[i][1], hols[i][2])
            k += 1


# ========================================================================== Nodes boundaries

def set_nbry(obj, id, x, y, z, key, value):
    try:
        prop      = obj.getProperty ('nbry_'+str(id))
        prop.data = x+' '+y+' '+z+' '+key+' '+value
    except:
        obj.addProperty ('nbry_'+str(id), x+' '+y+' '+z+' '+key+' '+value, 'STRING')

def set_nbry_x(obj, id, x): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = x+' '+d[1]+' '+d[2]+' '+d[3]+' '+d[4]

def set_nbry_y(obj, id, y): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+y+' '+d[2]+' '+d[3]+' '+d[4]

def set_nbry_z(obj, id, z): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+z+' '+d[3]+' '+d[4]

def set_nbry_key(obj, id, key): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+d[2]+' '+key+' '+d[4]

def set_nbry_val(obj, id, val): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+d[2]+' '+d[3]+' '+val

def get_nbrys(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='nbry':
            d = p.data.split()
            res.append ([d[0], d[1], d[2], d[3], d[4]])
    return res 

def get_nbrys_numeric(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='nbry':
            d = p.data.split()
            res.append ([float(d[0]), float(d[1]), float(d[2]), d[3], float(d[4])])
    return res 

def del_all_nbrys(obj):
    for p in obj.getAllProperties():
        if p.name[:4]=='nbry': obj.removeProperty(p)

def del_nbry(obj, id):
    nbrys = get_nbrys (obj)
    del_all_nbrys (obj)
    k = 0
    for i in range(len(nbrys)):
        if i!=id:
            set_nbry (obj, k, nbrys[i][0], nbrys[i][1], nbrys[i][2], nbrys[i][3], nbrys[i][4])
            k += 1


# ========================================================== Nodes boundaries (given node ID)

def set_nbryid(obj, id, nid, key, value):
    try:
        prop      = obj.getProperty ('nbid_'+str(id))
        prop.data = nid+' '+key+' '+value
    except:
        obj.addProperty ('nbid_'+str(id), nid+' '+key+' '+value, 'STRING')

def set_nbryid_nid(obj, id, nid): # property must exist
    p = obj.getProperty ('nbid_'+str(id))
    d = p.data.split()
    p.data = nid+' '+d[1]+' '+d[2]

def set_nbryid_key(obj, id, key): # property must exist
    p = obj.getProperty ('nbid_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+key+' '+d[2]

def set_nbryid_val(obj, id, val): # property must exist
    p = obj.getProperty ('nbid_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+val

def get_nbryids(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='nbid':
            d = p.data.split()
            res.append ([d[0], d[1], d[2]])
    return res 

def get_nbryids_numeric(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='nbid':
            d = p.data.split()
            res.append ([int(d[0]), d[1], float(d[2])])
    return res 

def del_all_nbryids(obj):
    for p in obj.getAllProperties():
        if p.name[:4]=='nbid': obj.removeProperty(p)

def del_nbryid(obj, id):
    nbryids = get_nbryids (obj)
    del_all_nbryids (obj)
    k = 0
    for i in range(len(nbryids)):
        if i!=id:
            set_nbryid (obj, k, nbryids[i][0], nbryids[i][1], nbryids[i][2])
            k += 1


# ========================================================================== Edges boundaries

def get_ebrys(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='ebry':
            d = p.data.split()
            res.append ([d[0], d[1], d[2]])
    return res 

def get_ebrys_numeric(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='ebry':
            d = p.data.split()
            res.append ([int(d[0]), d[1], float(d[2])])
    return res 
