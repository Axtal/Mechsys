# Modules
import Blender
import string

# Dictionary
def load_dict():
    dict = Blender.Registry.GetKey('MechSysDict')
    if not dict:
        dict                  = {}
        dict['fillet_radius'] = 0.0
        dict['fillet_steps']  = 10
        dict['show_props']    = 1
        dict['show_v_ids']    = 1
        dict['show_e_ids']    = 1
        dict['show_f_ids']    = 1
        dict['show_ndivs']    = 1
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
    #      key = 'edge' or 'face'
    #      id  = global ID of edge or face with tag
    #      tag = the tag of edge or face with tag
    set_int_property (obj, key+'_'+str(id), tag)


def get_tags(obj, key):
    # In:
    #      key = 'edge' or 'face'
    # Out:
    #      [ [edge_global_id, tag], ... num edges with tags]
    # or
    #      [ [face_global_id, tag], ... num faces with tags]
    res = []
    for p in obj.getAllProperties():
        if p.name[:4]==key:
            res.append ([int(p.name[5:]), p.data])
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
