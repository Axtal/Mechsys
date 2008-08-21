# Modules
import Blender
import bpy
import string


# ================================================================================ Dictionary

def load_dict():
    dict = Blender.Registry.GetKey('MechSysDict')
    if not dict:
        dict                  = {}
        dict['inirow']        = 0
        dict['newpoint_x']    = '0.0'
        dict['newpoint_y']    = '0.0'
        dict['newpoint_z']    = '0.0'
        dict['newetag']       = -10
        dict['newftag']       = -100
        dict['fillet_radius'] = '0.0'
        dict['fillet_steps']  = 10
        dict['show_props']    = 0
        dict['show_v_ids']    = 1
        dict['show_e_ids']    = 1
        dict['show_f_ids']    = 0
        dict['show_etags']    = 1
        dict['show_ftags']    = 1
        dict['show_elems']    = 1
        dict['show_axes']     = 1
        Blender.Registry.SetKey('MechSysDict', dict)
        print '[1;34mMechSys[0m: dictionary created'
    return dict


# ====================================================================================== Util

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
    else:
        if with_error: Blender.Draw.PupMenu('ERROR|Please, select an object of type Mesh')
        else: return None,None,None


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
            Blender.Draw.PupMenu('ERROR|local x-y axes must share the same origin vertex (obj=%s)' % obj.name)
            return
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
                Blender.Draw.PupMenu('ERROR|local x-y-z axes must share the same origin vertex (obj=%s)' % obj.name)
                return
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
    try:    return obj.getProperty('btag').data
    except: return 0


# ====================================================================================== Tags

def set_etag(obj, id, tag):
    # In:
    #      id  = global ID of edge/face/elem with tag
    #      tag = the tag of edge/face/elem with tag
    set_int_property (obj, 'edge_'+str(id), tag)

def set_ftag(obj, edges_ids, tag):
    # In:
    #      edges_ids = global IDs of edges on selected face
    #      tag       = the tag of face
    eids = '_'.join([str(id) for id in edges_ids])
    if not obj.properties.has_key('ftags'): obj.properties['ftags'] = {}
    if tag==0: obj.properties['ftags'].pop(eids)
    else:      obj.properties['ftags'].update({eids:tag})

def get_tags(obj, key):
    # In:
    #      key = 'edge', 'face', or 'elem'
    # Out:
    #      [ [e/f/e_global_id, tag], ... num edges/face/elem with tags]
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]==key:
            res.append ([int(p.name[5:]), p.data])
    return res 

def get_etags(obj, msh):
    # Out:   {(v1,v2):tag1,  (v1,v2):tag2,  ... num edges with tags}
    res = {}
    for p in obj.getAllProperties():
        if p.name[:4]=='edge':
            eid = int(p.name[5:])
            res[(msh.edges[eid].v1.index, msh.edges[eid].v2.index)] = p.data
    return res

def get_ftags(obj, msh):
    # Out:   {(e1,e2,e3,e4,..):tag1,  (e1,e2,e3,e4..):tag2,  ... num faces with tags}
    res = {}
    if obj.properties.has_key('ftags'):
        for eids in obj.properties['ftags']:
            ids = [int(id) for id in eids.split('_')]
            res[tuple(ids)] = obj.properties['ftags'][eids]
    return res

def get_tags_list(obj, key, global_ids):
    # In:
    #      key        = 'edge' or 'face'
    #      global_ids = edge/face global ids sorted in _Blender_ way
    # Out:
    #      tags_list  = list with tags sorted in _MechSys_ way
    #
    # MechSys:                Blender (returned by 'face.edge_keys')
    #            3                            2
    #      +-----------+                +-----------+   
    #      |           |                |           |
    #      |           |                |           |   
    #    0 |           | 1    <<==    3 |           | 1
    #      |           |                |           |
    #      |           |                |           |   
    #      +-----------+                +-----------+   
    #            2                            0
    #
    # Examples:
    #   2D (key=='edge'):
    #                         3                           12 (-12)
    #                   +-----------+                +-----------+   
    #                   |           |                |           |
    #                   |           |                |           |   
    #                 0 |           | 1    <<==   13 |           | 11 (-11)
    #                   |           |           (-13)|           |
    #                   |           |                |           |   
    #                   +-----------+                +-----------+   
    #                         2                           10 (-10)
    #                In:  global_ids = [10, 11, 12, 13]
    #                Out: tags_list  = [-13, -11, -10, -12] (local/MechSys sorted)
    #
    ntags     = len(global_ids) # 4 for 2D or 12(edges) or 6(faces) for 3D
    tags      = get_tags (obj, key)
    tags_list = [0 for i in range(ntags)]
    if ntags==4: map = [2, 1, 3, 0] # map Blender face edge local ID to MechSys edge local ID
    for t in tags:
        if t[0] in global_ids and t[1]<0:
            local_id            = map[global_ids.index(t[0])]
            tags_list[local_id] = t[1]
    if sum(tags_list)==0: tags_list=[]
    return tags_list


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

def set_reg(obj, id, maxarea, x, y, z):
    try:
        prop      = obj.getProperty ('reg_'+str(id))
        prop.data = maxarea+' '+x+' '+y+' '+z
    except:
        obj.addProperty ('reg_'+str(id), maxarea+' '+x+' '+y+' '+z, 'STRING')

def set_reg_maxarea(obj, id, maxarea): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = maxarea+' '+d[1]+' '+d[2]+' '+d[3]

def set_reg_x(obj, id, x): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+x+' '+d[2]+' '+d[3]

def set_reg_y(obj, id, y): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+y+' '+d[3]

def set_reg_z(obj, id, z): # property must exist
    p = obj.getProperty ('reg_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+d[2]+' '+z

def get_regs(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:3]=='reg':
            d = p.data.split()
            res.append ([d[0], d[1], d[2], d[3]])
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
            set_reg (obj, k, regs[i][0], regs[i][1], regs[i][2], regs[i][3])
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
    p.data = d[0]+' '+d[1]+' '+d[2]+' '+d[3]+' '+value

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


# ========================================================================== Edges boundaries

def set_ebry(obj, id, tag, key, value):
    try:
        prop      = obj.getProperty ('ebry_'+str(id))
        prop.data = tag+' '+key+' '+value
    except:
        obj.addProperty ('ebry_'+str(id), tag+' '+key+' '+value, 'STRING')

def set_ebry_tag(obj, id, tag): # property must exist
    p = obj.getProperty ('ebry_'+str(id))
    d = p.data.split()
    p.data = tag+' '+d[1]+' '+d[2]

def set_ebry_key(obj, id, key): # property must exist
    p = obj.getProperty ('ebry_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+key+' '+d[2]

def set_ebry_val(obj, id, val): # property must exist
    p = obj.getProperty ('ebry_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+val

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

def del_all_ebrys(obj):
    for p in obj.getAllProperties():
        if p.name[:4]=='ebry': obj.removeProperty(p)

def del_ebry(obj, id):
    ebrys = get_ebrys (obj)
    del_all_ebrys (obj)
    k = 0
    for i in range(len(ebrys)):
        if i!=id:
            set_ebry (obj, k, ebrys[i][0], ebrys[i][1], ebrys[i][2])
            k += 1


# ========================================================================== Faces boundaries

def set_fbry(obj, id, tag, key, value):
    try:
        prop      = obj.getProperty ('fbry_'+str(id))
        prop.data = tag+' '+key+' '+value
    except:
        obj.addProperty ('fbry_'+str(id), tag+' '+key+' '+value, 'STRING')

def set_fbry_tag(obj, id, tag): # property must exist
    p = obj.getProperty ('fbry_'+str(id))
    d = p.data.split()
    p.data = tag+' '+d[1]+' '+d[2]

def set_fbry_key(obj, id, key): # property must exist
    p = obj.getProperty ('fbry_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+key+' '+d[2]

def set_fbry_val(obj, id, val): # property must exist
    p = obj.getProperty ('fbry_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+val

def get_fbrys(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='fbry':
            d = p.data.split()
            res.append ([d[0], d[1], d[2]])
    return res 

def get_fbrys_numeric(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='fbry':
            d = p.data.split()
            res.append ([int(d[0]), d[1], float(d[2])])
    return res 

def del_all_fbrys(obj):
    for p in obj.getAllProperties():
        if p.name[:4]=='fbry': obj.removeProperty(p)

def del_fbry(obj, id):
    fbrys = get_fbrys (obj)
    del_all_fbrys (obj)
    k = 0
    for i in range(len(fbrys)):
        if i!=id:
            set_fbry (obj, k, fbrys[i][0], fbrys[i][1], fbrys[i][2])
            k += 1


# ======================================================================= Elements attributes

def set_eatt(obj, id, tag, type, model, prms, inis):
    prms = '_'.join(prms.split())
    inis = '_'.join(inis.split())
    try:
        prop      = obj.getProperty ('eatt_'+str(id))
        prop.data = tag+' '+type+' '+model+' '+prms+' '+inis
    except:
        obj.addProperty ('eatt_'+str(id), tag+' '+type+' '+model+' '+prms+' '+inis, 'STRING')

def set_eatt_tag(obj, id, tag): # property must exist
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = tag+' '+d[1]+' '+d[2]+' '+d[3]+' '+d[4]

def set_eatt_type(obj, id, type): # property must exist
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+type+' '+d[2]+' '+d[3]+' '+d[4]

def set_eatt_model(obj, id, model): # property must exist
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+model+' '+d[3]+' '+d[4]

def set_eatt_prms(obj, id, prms): # property must exist
    prms = '_'.join(prms.split())
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+d[2]+' '+prms+' '+d[4]

def set_eatt_inis(obj, id, inis): # property must exist
    inis = '_'.join(inis.split())
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = d[0]+' '+d[1]+' '+d[2]+' '+d[3]+' '+inis

def get_eatts(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='eatt':
            d = p.data.split()
            res.append ([d[0], d[1], d[2], d[3].replace('_',' '), d[4].replace('_',' ')])
    return res 

def get_eatts_numeric(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='eatt':
            d = p.data.split()
            res.append ([int(d[0]), d[1], d[2], d[3].replace('_',' '), d[4].replace('_',' ')])
    return res 

def del_all_eatts(obj):
    for p in obj.getAllProperties():
        if p.name[:4]=='eatt': obj.removeProperty(p)

def del_eatt(obj, id):
    eatts = get_eatts (obj)
    del_all_eatts (obj)
    k = 0
    for i in range(len(eatts)):
        if i!=id:
            set_eatt (obj, k, eatts[i][0], eatts[i][1], eatts[i][2], eatts[i][3], eatts[i][4])
            k += 1
