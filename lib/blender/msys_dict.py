# Modules
import Blender

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
        Blender.Registry.SetKey('MechSysDict', dict)
        print '[1;34mMechSys[0m: dictionary created'
    return dict

 # IDs of the edges corresponding to the local x, y and z axes
def get_local_axes_props(obj):
    props = {}
    keys  = ['edge_x','edge_y','edge_z']
    for k in keys:
        try: props[k] = obj.getProperty(k)
        except:
            obj.addProperty (k, -1, 'INT')
            props[k] = obj.getProperty(k)
    return props

# Number of divisions along x, y, and z local axes
def get_ndivs_props(obj):
    props = {}
    keys  = ['ndiv_x','ndiv_y','ndiv_z']
    for k in keys:
        try: props[k] = obj.getProperty(k)
        except:
            obj.addProperty (k, 2, 'INT')
            props[k] = obj.getProperty(k)
    return props

# IDs and marks of vertices on boundary
def get_v_bry_props(obj):
    props = {}
    keys  = ['v_bry_ids','v_bry_mks']
    for k in keys:
        try: props[k] = obj.getProperty(k)
        except:
            obj.addProperty (k, '', 'STRING')
            props[k] = obj.getProperty(k)
    return props

# IDs and marks of edges on boundary
def get_e_bry_props(obj):
    props = {}
    keys  = ['e_bry_ids','e_bry_mks']
    for k in keys:
        try: props[k] = obj.getProperty(k)
        except:
            obj.addProperty (k, '', 'STRING')
            props[k] = obj.getProperty(k)
    return props

# IDs and marks of faces on boundary
def get_f_bry_props(obj):
    props = {}
    keys  = ['f_bry_ids','f_bry_mks']
    for k in keys:
        try: props[k] = obj.getProperty(k)
        except:
            obj.addProperty (k, '', 'STRING')
            props[k] = obj.getProperty(k)
    return props
