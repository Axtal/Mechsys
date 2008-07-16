# Modules
import Blender
import string

# Dictionary
def load_dict():
    dict = Blender.Registry.GetKey('MechSysDict')
    if not dict:
        dict                  = {}
        dict['newpoint_x']    = 0.0
        dict['newpoint_y']    = 0.0
        dict['newpoint_z']    = 0.0
        dict['fillet_radius'] = 0.0
        dict['fillet_steps']  = 10
        dict['show_props']    = 1
        dict['show_v_ids']    = 1
        dict['show_e_ids']    = 1
        dict['show_f_ids']    = 1
        dict['show_ele_tags'] = 0
        dict['show_axes']     = 1
        Blender.Registry.SetKey('MechSysDict', dict)
        print '[1;34mMechSys[0m: dictionary created'
    return dict


def set_int_property(obj, key, value):
    try:
        prop      = obj.getProperty (key)
        prop.data = value
    except:
        obj.addProperty (key, value, 'INT')


def set_btag(obj, tag):
    set_int_property (obj, 'btag', tag)


def get_btag(obj):
    tag = 0
    for p in obj.getAllProperties():
        if p.name[:4]=='btag':
            tag = p.data
    return tag


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


def set_tag(obj, key, id, tag):
    # In:
    #      key = 'edge', 'face', or 'elem'
    #      id  = global ID of edge/face/elem with tag
    #      tag = the tag of edge/face/elem with tag
    set_int_property (obj, key+'_'+str(id), tag)


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


def get_verts_on_edges_with_tags(obj, msh):
    etags = get_tags (obj, 'edge')
    verts = []
    for e in etags:
        idx1 = msh.edges[e[0]].v1.index
        idx2 = msh.edges[e[0]].v2.index
        if idx1 not in verts: verts.append(idx1)
        if idx2 not in verts: verts.append(idx2)
    return verts


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


# ==================================================== Nodes boundaries

def set_nbry(obj, id, x, y, z, key, value):
    try:
        prop      = obj.getProperty ('nbry_'+str(id))
        prop.data = '%8.4f %8.4f %8.4f %s %12.6f'%(x,y,z,key,value)
    except:
        obj.addProperty ('nbry_'+str(id), '%8.4f %8.4f %8.4f %s %12.6f'%(x,y,z,key,value), 'STRING')

def set_nbry_x(obj, id, val): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = '%8.4f %8.4f %8.4f %s %12.6f'%(val,float(d[1]),float(d[2]), d[3] ,float(d[4]))

def set_nbry_y(obj, id, val): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = '%8.4f %8.4f %8.4f %s %12.6f'%(float(d[0]),val,float(d[2]), d[3] ,float(d[4]))

def set_nbry_z(obj, id, val): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = '%8.4f %8.4f %8.4f %s %12.6f'%(float(d[0]),float(d[1]),val, d[3] ,float(d[4]))

def set_nbry_key(obj, id, val): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = '%8.4f %8.4f %8.4f %s %12.6f'%(float(d[0]),float(d[1]),float(d[2]), val ,float(d[4]))

def set_nbry_val(obj, id, val): # property must exist
    p = obj.getProperty ('nbry_'+str(id))
    d = p.data.split()
    p.data = '%8.4f %8.4f %8.4f %s %12.6f'%(float(d[0]),float(d[1]),float(d[2]), d[3] ,val)

def get_nbrys(obj):
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]=='nbry':
            d = p.data.split()
            res.append ([float(d[0]),float(d[1]),float(d[2]), d[3], float(d[4])])
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


# ==================================================== Faces boundaries

def set_fbry(obj, id, tag, key, value):
    try:
        prop      = obj.getProperty ('fbry_'+str(id))
        prop.data = '%d %s %12.6f'%(tag,key,value)
    except:
        obj.addProperty ('fbry_'+str(id), '%d %s %12.6f'%(tag,key,value), 'STRING')

def set_fbry_tag(obj, id, val): # property must exist
    p = obj.getProperty ('fbry_'+str(id))
    d = p.data.split()
    p.data = '%d %s %12.6f'%(val, d[1] ,float(d[2]))

def set_fbry_key(obj, id, val): # property must exist
    p = obj.getProperty ('fbry_'+str(id))
    d = p.data.split()
    p.data = '%d %s %12.6f'%(int(d[0]), val ,float(d[2]))

def set_fbry_val(obj, id, val): # property must exist
    p = obj.getProperty ('fbry_'+str(id))
    d = p.data.split()
    p.data = '%d %s %12.6f'%(int(d[0]), d[1] ,val)

def get_fbrys(obj):
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


# ==================================================== Elements attributes

def set_eatt(obj, id, tag, type, model, prms, inis):
    prms = '_'.join(prms.split())
    inis = '_'.join(inis.split())
    try:
        prop      = obj.getProperty ('eatt_'+str(id))
        prop.data = '%d %s %s %s %s'%(tag,type,model,prms,inis)
    except:
        obj.addProperty ('eatt_'+str(id), '%d %s %s %s %s'%(tag,type,model,prms,inis), 'STRING')

def set_eatt_tag(obj, id, val): # property must exist
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = '%d %s %s %s %s'%(val, d[1], d[2], d[3], d[4])

def set_eatt_type(obj, id, val): # property must exist
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = '%d %s %s %s %s'%(int(d[0]), val, d[2], d[3], d[4])

def set_eatt_model(obj, id, val): # property must exist
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = '%d %s %s %s %s'%(int(d[0]), d[1], val, d[3], d[4])

def set_eatt_prms(obj, id, val): # property must exist
    val = '_'.join(val.split())
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = '%d %s %s %s %s'%(int(d[0]), d[1], d[2], val, d[4])

def set_eatt_inis(obj, id, val): # property must exist
    val = '_'.join(val.split())
    p = obj.getProperty ('eatt_'+str(id))
    d = p.data.split()
    p.data = '%d %s %s %s %s'%(int(d[0]), d[1], d[2], d[3], val)

def get_eatts(obj):
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
