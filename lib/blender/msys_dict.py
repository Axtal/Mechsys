# Modules
import Blender

# Dictionary
def load_dict():
    dict = Blender.Registry.GetKey('MechSysDict')
    if not dict:
        dict                   = {}
        dict['fillet_radius']  = 0.0
        dict['fillet_steps']   = 10
        dict['display_props']  = 1
        dict['vertex_ids']     = 1
        dict['edge_ids']       = 1
        dict['face_ids']       = 1
        dict['disp_ndivs']     = 1
        Blender.Registry.SetKey('MechSysDict', dict)
        print '[1;34mMechSys[0m: dictionary created'
    return dict

def get_all_props(obj):
    props = {}

    # ID of the edge corresponding to the local x axis
    try: props['edge_x'] = obj.getProperty('edge_x')
    except:
        obj.addProperty('edge_x',-1,'INT')
        props['edge_x'] = obj.getProperty('edge_x')

    # ID of the edge corresponding to the local y ayis
    try: props['edge_y'] = obj.getProperty('edge_y')
    except:
        obj.addProperty('edge_y',-1,'INT')
        props['edge_y'] = obj.getProperty('edge_y')

    # ID of the edge corresponding to the local z azis
    try: props['edge_z'] = obj.getProperty('edge_z')
    except:
        obj.addProperty('edge_z',-1,'INT')
        props['edge_z'] = obj.getProperty('edge_z')

    # number of divisions along x
    try: props['ndiv_x'] = obj.getProperty('ndiv_x')
    except:
        obj.addProperty('ndiv_x',2,'INT')
        props['ndiv_x'] = obj.getProperty('ndiv_x')

    # number of divisions along x
    try: props['ndiv_y'] = obj.getProperty('ndiv_y')
    except:
        obj.addProperty('ndiv_y',2,'INT')
        props['ndiv_y'] = obj.getProperty('ndiv_y')

    # number of divisions along x
    try: props['ndiv_z'] = obj.getProperty('ndiv_z')
    except:
        obj.addProperty('ndiv_z',2,'INT')
        props['ndiv_z'] = obj.getProperty('ndiv_z')

    return props
