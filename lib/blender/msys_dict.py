# Modules
import Blender

# Dictionary
def load_dict():
    dict = Blender.Registry.GetKey('MechSysDict')
    if not dict:
        dict                  = {}
        dict['fillet_radius'] = 0.0
        dict['fillet_steps']  = 10
        Blender.Registry.SetKey('MechSysDict', dict)
        print '[1;34mMechSysCAD[0m: dictionary created'
    return dict
