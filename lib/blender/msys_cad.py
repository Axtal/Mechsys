#!BPY

"""
Name: 'MechSys'
Blender: 2.46
Group: 'Themes'
Tooltip: 'CAD, Mesh generator, and FEM in Blender'
"""
__author__  = "Dorival Pedroso"
__version__ = "06.05.08"
__bpydoc__="""\
TODO
"""

# In Mac OS X, the scripts folder is:
#   /Applications/blender.app/Contents/MacOS/.blender/scripts/

# Modules
import timeit
import Blender
from   Blender import Draw, BGL
from   Blender.Mathutils import Vector
import bpy
import msys_draw as dr
import msys_dict as di
import msys_fem  as fem


# ================================================================================= Constants

EVT_DEL              = 10000 # used for passing IDs through events
EVT_NONE             =  0 # used for buttons with callbacks
EVT_REFRESH          =  1 # refresh all windows
# SETTINGS
EVT_SHOWHIDE_SET     =  2 # show/hide SET box
EVT_SET_DELPROPS     =  3 # delete all properties
# CAD
EVT_SHOWHIDE_CAD     =  4 # show/hide CAD box
EVT_CAD_ADDXYZ       =  5 # type in values to define point(s)
EVT_CAD_FILLET       =  6 # 3D fillet
EVT_CAD_EDGEBREAK    =  7 # break edge
EVT_CAD_EDGEBREAKM   =  8 # break edge at middle point
EVT_CAD_EDGEINTERS   =  9 # edge closest distance
EVT_CAD_FPOINT       = 10 # read points from file
EVT_CAD_FSPLINE      = 11 # create a spline from points in a file
# Mesh
EVT_SHOWHIDE_MESH    = 12 # show/hide MESH box
EVT_MESH_SETETAG     = 13 # set edges tag
EVT_MESH_SETFTAG     = 14 # set faces tag
# Mesh -- structured
EVT_MESH_SETX        = 15 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_SETY        = 16 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_SETZ        = 17 # set local x-y coordinates of a block (need two edges pre-selected)
EVT_MESH_GENSTRU     = 18 # generate structured mesh using MechSys module
# Mesh -- unstructured
EVT_MESH_GENUNSTRU   = 19 # generate structured mesh using MechSys module
EVT_MESH_ADDREG      = 20 # add region
EVT_MESH_DELALLREGS  = 21 # delete all regions
EVT_MESH_ADDHOLE     = 22 # add hole
EVT_MESH_DELALLHOLES = 23 # delete all holes
# FEM
EVT_SHOWHIDE_FEM     = 24 # show/hide FEM box
EVT_FEM_ADDNBRY      = 25 # add nodes boundary (given coordinates)
EVT_FEM_DELALLNBRY   = 26 # delete all nodes boundary
EVT_FEM_ADDNBRYID    = 27 # add nodes boundary (given nodes IDs)
EVT_FEM_DELALLNBRYID = 28 # delete all nodes boundary
EVT_FEM_ADDEBRY      = 29 # add edges boundary
EVT_FEM_DELALLEBRY   = 30 # delete all edges boundary
EVT_FEM_RUN          = 31 # run a FE simulation
EVT_FEM_SCRIPT       = 32 # generate script for FEM 
EVT_FEM_PARAVIEW     = 33 # view in ParaView
# Results
EVT_SHOWHIDE_RES     = 34 # show/hide results box


# ==================================================================================== Events

# Handle input events
def event(evt, val):
    if   evt==Draw.QKEY and not val: Draw.Exit()
    elif evt==Draw.ESCKEY: Draw.Redraw(1)
    elif evt==Draw.WHEELDOWNMOUSE or evt==Draw.DOWNARROWKEY:
        d = di.load_dict()
        d['inirow'] -= 40
        Blender.Window.QRedrawAll()
    elif evt==Draw.WHEELUPMOUSE or evt==Draw.UPARROWKEY:
        d = di.load_dict()
        if d['inirow']<0:
            d['inirow'] += 40
            Blender.Window.QRedrawAll()
    elif evt==Draw.PAGEDOWNKEY:
        d = di.load_dict()
        d['inirow'] -= 100
        Blender.Window.QRedrawAll()
    elif evt==Draw.PAGEUPKEY:
        d = di.load_dict()
        if d['inirow']<0:
            d['inirow'] += 100
            Blender.Window.QRedrawAll()


# Handle button events
def button_event(evt):
    # load dictionary
    dict = di.load_dict()

    # refresh all windows
    if evt==EVT_REFRESH: Blender.Window.QRedrawAll()

    if True:#try:
        # ----------------------------------------------------------------------------------- SETTINGS

        # show/hide CAD box
        if evt==EVT_SHOWHIDE_SET:
            dict['gui_show_set'] = 0 if dict['gui_show_set'] else 1
            Blender.Window.QRedrawAll()

        # delete all properties
        elif evt==EVT_SET_DELPROPS:
            scn = bpy.data.scenes.active
            obs = scn.objects.selected
            if len(obs)>0:
                message = 'Confirm delete ALL properties?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    for o in obs:
                        for k, v in o.properties.iteritems():
                            print '[1;34mMechSys[0m: >> Object = ', o.name, ' [1;35mDeleted property:[0m ', k, str(v)
                            o.properties.pop(k)
                        for p in o.getAllProperties():
                            print '[1;34mMechSys[0m: Object = ', o.name, ' Deleted property: ', p.name, p.data
                            o.removeProperty(p)
                Blender.Window.QRedrawAll()

        # ----------------------------------------------------------------------------------- CAD

        # show/hide CAD box
        elif evt==EVT_SHOWHIDE_CAD:
            dict['gui_show_cad'] = 0 if dict['gui_show_cad'] else 1
            Blender.Window.QRedrawAll()

        # add vertices
        elif evt==EVT_CAD_ADDXYZ:
            x = dict['newpoint_x']
            y = dict['newpoint_y']
            z = dict['newpoint_z']
            dr.add_point (float(x), float(y), float(z))

        # fillet two edges
        elif evt==EVT_CAD_FILLET: dr.fillet(float(dict['fillet_radius']),dict['fillet_steps'])

        # break an edge
        elif evt==EVT_CAD_EDGEBREAK: dr.break_edge()

        # break an edge at middle point
        elif evt==EVT_CAD_EDGEBREAKM: dr.break_edge(True)

        # edge closest distance
        elif evt==EVT_CAD_EDGEINTERS: dr.edge_intersect()

        # read vertices from file
        elif evt==EVT_CAD_FPOINT: Blender.Window.FileSelector(dr.add_points_from_file, "Read X Y Z cols")

        # read spline from file
        elif evt==EVT_CAD_FSPLINE: Blender.Window.FileSelector(dr.add_spline_from_file, "Read X Y Z cols")

        # ---------------------------------------------------------------------------------- Mesh

        # show/hide MESH box
        elif evt==EVT_SHOWHIDE_MESH:
            dict['gui_show_mesh'] = 0 if dict['gui_show_mesh'] else 1
            Blender.Window.QRedrawAll()

        # set edges tag
        elif evt==EVT_MESH_SETETAG:
            edm, obj, msh = di.get_msh()
            for id in msh.edges.selected():
                di.set_etag (obj, id, dict['newetag'])
            Blender.Window.QRedrawAll()
            if edm: Blender.Window.EditMode(1) # return to EditMode

        # set faces tag
        elif evt==EVT_MESH_SETFTAG:
            edm, obj, msh = di.get_msh()
            nedges = len(msh.edges.selected())
            if nedges==3 or nedges==6 or nedges==4 or nedges==8:
                eids = '_'.join([str(id) for id in msh.edges.selected()])
                if not obj.properties.has_key('ftags'): obj.properties['ftags'] = {}
                if dict['newftag'][0]==0: obj.properties['ftags'].pop(eids)
                else:                     obj.properties['ftags'].update({eids:dict['newftag']})
                if len(obj.properties['ftags'])==0: obj.properties.pop('ftags')
            else: raise Exception('To set a face tag (FTAG): 3, 4, 6, or 8 edges must be selected')
            Blender.Window.QRedrawAll()
            if edm: Blender.Window.EditMode(1) # return to EditMode

        # ---------------------------------------------------------------------------- Structured

        # set local x-axis
        elif evt==EVT_MESH_SETX:
            edm, obj, msh = di.get_msh()
            if len(msh.edges.selected())==1:
                di.set_local_axis (obj, 'x', msh.edges.selected()[0])
                Blender.Window.QRedrawAll()
            else: raise Exception('Please, select only one edge (obj=%s)' % obj.name)
            if edm: Blender.Window.EditMode(1) # return to EditMode

        # set local y-axis
        elif evt==EVT_MESH_SETY:
            edm, obj, msh = di.get_msh()
            if len(msh.edges.selected())==1:
                di.set_local_axis (obj, 'y', msh.edges.selected()[0])
                Blender.Window.QRedrawAll()
            else: raise Exception('Please, select only one edge (obj=%s)' % obj.name)
            if edm: Blender.Window.EditMode(1) # return to EditMode

        # set local z-axis
        elif evt==EVT_MESH_SETZ:
            edm, obj, msh = di.get_msh()
            if len(msh.edges.selected())==1:
                di.set_local_axis (obj, 'z', msh.edges.selected()[0])
                Blender.Window.QRedrawAll()
            else: raise Exception('Please, select only one edge (obj=%s)' % obj.name)
            if edm: Blender.Window.EditMode(1) # return to EditMode

        # generate structured mesh via MechSys
        elif evt==EVT_MESH_GENSTRU:
            tt = timeit.Timer('msys_draw.gen_struct_mesh()', 'import msys_draw')
            t  = tt.timeit(number=1)
            print '[1;34mMechSys[0m: time spent on generation and drawing = [1;31m',t,'[0m [1;32mseconds[0m'

        # -------------------------------------------------------------------------- Unstructured

        # add region
        elif evt==EVT_MESH_ADDREG:
            obj = di.get_obj()
            if obj!=None:
                x, y, z = Blender.Window.GetCursorPos()
                regs = di.get_regs (obj)
                di.set_reg (obj, len(regs), -1, '-1', str(x),str(y),str(z))
                Blender.Window.QRedrawAll()

        # delete all regions
        elif evt==EVT_MESH_DELALLREGS:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_regs(obj)
                    Blender.Window.QRedrawAll()

        # add hole
        elif evt==EVT_MESH_ADDHOLE:
            obj = di.get_obj()
            if obj!=None:
                x, y, z = Blender.Window.GetCursorPos()
                hols = di.get_hols (obj)
                di.set_hol (obj, len(hols), str(x),str(y),str(z))
                Blender.Window.QRedrawAll()

        # delete all holes
        elif evt==EVT_MESH_DELALLHOLES:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_hols(obj)
                    Blender.Window.QRedrawAll()

        # generate unstructured mesh via MechSys
        elif evt==EVT_MESH_GENUNSTRU:
            tt = timeit.Timer('msys_draw.gen_unstruct_mesh()', 'import msys_draw')
            t  = tt.timeit(number=1)
            print '[1;34mMechSys[0m: time spent on generation and drawing = [1;31m',t,'[0m [1;32mseconds[0m'

        # ----------------------------------------------------------------------------------- FEM 

        # show/hide FEM box
        elif evt==EVT_SHOWHIDE_FEM:
            dict['gui_show_fem'] = 0 if dict['gui_show_fem'] else 1
            Blender.Window.QRedrawAll()

        # add nodes boundary
        elif evt==EVT_FEM_ADDNBRY:
            obj = di.get_obj()
            if obj!=None:
                nbrys = di.get_nbrys (obj)
                di.set_nbry (obj, len(nbrys), '0.0','0.0','0.0', 'uy', '0.0')
                Blender.Window.QRedrawAll()

        # delete all nodes boundary
        elif evt==EVT_FEM_DELALLNBRY:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_nbrys(obj)
                    Blender.Window.QRedrawAll()

        # add nodes boundary (given tags)
        elif evt==EVT_FEM_ADDNBRYID:
            obj = di.get_obj()
            if obj!=None:
                nbryids = di.get_nbryids (obj)
                di.set_nbryid (obj, len(nbryids), '0', 'ux', '0.0')
                Blender.Window.QRedrawAll()

        # delete all nodes boundary (given tags)
        elif evt==EVT_FEM_DELALLNBRYID:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_nbryids(obj)
                    Blender.Window.QRedrawAll()

        # add edges boundary
        elif evt==EVT_FEM_ADDEBRY:
            obj = di.get_obj()
            if obj!=None:
                ebrys = di.get_ebrys (obj)
                di.set_ebry (obj, len(ebrys), '-10', 'uy', '0.0')
                Blender.Window.QRedrawAll()

        # delete all edges boundary
        elif evt==EVT_FEM_DELALLEBRY:
            obj = di.get_obj()
            if obj!=None:
                message = 'Confirm delete ALL?%t|Yes'
                result  = Blender.Draw.PupMenu(message)
                if result>0:
                    di.del_all_ebrys(obj)
                    Blender.Window.QRedrawAll()

        # run a FE simulation
        elif evt==EVT_FEM_RUN:
            obj = di.get_obj ()
            fem.run_analysis (obj)

        # generate FE script
        elif evt==EVT_FEM_SCRIPT:
            obj = di.get_obj ()
            fem.gen_script   (obj)
            Blender.Window.Redraw(Blender.Window.Types.TEXT)

        # view results in ParaView
        elif evt==EVT_FEM_PARAVIEW:
            obj = di.get_obj ()
            fem.paraview (obj)

        # ----------------------------------------------------------------------------------- RES 

        # show/hide RES box
        elif evt==EVT_SHOWHIDE_RES:
            dict['gui_show_res'] = 0 if dict['gui_show_res'] else 1
            Blender.Window.QRedrawAll()


    #except Exception, inst:
        #msg = inst.args[0]
        #print '[1;34mMechSys[0m: Error: '+'[1;31m'+msg+'[0m'
        #Blender.Draw.PupMenu('ERROR|'+msg)


# ================================================================================= Callbacks

# ---------------------------------- CAD

def setx_callback(evt,val):
    dict = di.load_dict()
    dict['newpoint_x'] = val

def sety_callback(evt,val):
    dict = di.load_dict()
    dict['newpoint_y'] = val

def setz_callback(evt,val):
    dict = di.load_dict()
    dict['newpoint_z'] = val

def fillet_radius_callback(evt,val):
    dict = di.load_dict()
    dict['fillet_radius'] = val

def fillet_steps_callback(evt,val):
    dict = di.load_dict()
    dict['fillet_steps'] = val


# ---------------------------------- Mesh

def show_props_callback(evt,val):
    dict = di.load_dict()
    dict['show_props'] = val
    Blender.Window.QRedrawAll()

def show_axes_callback(evt,val):
    dict = di.load_dict()
    dict['show_axes'] = val
    Blender.Window.QRedrawAll()

def show_v_ids_callback(evt,val):
    dict = di.load_dict()
    dict['show_v_ids'] = val
    Blender.Window.QRedrawAll()

def show_e_ids_callback(evt,val):
    dict = di.load_dict()
    dict['show_e_ids'] = val
    Blender.Window.QRedrawAll()

def show_f_ids_callback(evt,val):
    dict = di.load_dict()
    dict['show_f_ids'] = val
    Blender.Window.QRedrawAll()

def ftags_opac_callback(evt,val):
    dict = di.load_dict()
    dict['ftags_opac'] = val
    Blender.Window.QRedrawAll()

def show_etags_callback(evt,val):
    dict = di.load_dict()
    dict['show_etags'] = val
    Blender.Window.QRedrawAll()

def show_ftags_callback(evt,val):
    dict = di.load_dict()
    dict['show_ftags'] = val
    Blender.Window.QRedrawAll()

def show_tags_txt_callback(evt,val):
    dict = di.load_dict()
    dict['show_tags_txt'] = val
    Blender.Window.QRedrawAll()

def show_elems_callback(evt,val):
    dict = di.load_dict()
    dict['show_elems'] = val
    Blender.Window.QRedrawAll()

def ndivx_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_ndiv (o, 'x', val)
    Blender.Window.QRedrawAll()

def ndivy_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_ndiv (o, 'y', val)
    Blender.Window.QRedrawAll()

def ndivz_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_ndiv (o, 'z', val)
    Blender.Window.QRedrawAll()

def acoefx_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_acoef (o, 'x', val)
    Blender.Window.QRedrawAll()

def acoefy_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_acoef (o, 'y', val)
    Blender.Window.QRedrawAll()

def acoefz_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_acoef (o, 'z', val)
    Blender.Window.QRedrawAll()

def nonlinx_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_nonlin (o, 'x', val)
    Blender.Window.QRedrawAll()

def nonliny_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_nonlin (o, 'y', val)
    Blender.Window.QRedrawAll()

def nonlinz_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_nonlin (o, 'z', val)
    Blender.Window.QRedrawAll()


# ---------------------------------- etag, ftag, btag

def etag_callback(evt,val):
    dict = di.load_dict()
    dict['newetag'] = val

def ftag_callback(evt,val):
    dict = di.load_dict()
    dict['newftag'] = [val, dict['newftag'][1]]
    Blender.Window.QRedrawAll()

def fclr_callback(evt,val):
    dict = di.load_dict()
    dict['newftag'] = [dict['newftag'][0], di.rgb2hex(val)]
    Blender.Window.QRedrawAll()

def btag_callback(evt,val):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for o in obs: di.set_btag (o, val)
    Blender.Window.QRedrawAll()


# ---------------------------------- Unstructured

def minangle_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_minangle (obj, val)
        Blender.Window.QRedrawAll()

def maxarea_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_maxarea (obj, val)
        Blender.Window.QRedrawAll()


# ---------------------------------- Regions

def regs_tag_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_tag (obj, evt-EVT_DEL, val)

def regs_maxarea_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_maxarea (obj, evt-EVT_DEL, val)

def regs_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_x (obj, evt-EVT_DEL, val)
        Blender.Window.QRedrawAll()

def regs_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_y (obj, evt-EVT_DEL, val)
        Blender.Window.QRedrawAll()

def regs_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_reg_z (obj, evt-EVT_DEL, val)
        Blender.Window.QRedrawAll()

def regs_delonereg_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_reg (obj, evt-EVT_DEL)
        Blender.Window.QRedrawAll()


# ---------------------------------- Holes

def holes_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_x (obj, evt-EVT_DEL, val)
        Blender.Window.QRedrawAll()

def holes_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_y (obj, evt-EVT_DEL, val)
        Blender.Window.QRedrawAll()

def holes_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_hol_z (obj, evt-EVT_DEL, val)
        Blender.Window.QRedrawAll()

def holes_delonehole_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_hol (obj, evt-EVT_DEL)
        Blender.Window.QRedrawAll()


# ---------------------------------- nbrys

def nodebry_setx_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_x (obj, evt-EVT_DEL, val)

def nodebry_sety_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_y (obj, evt-EVT_DEL, val)

def nodebry_setz_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_z (obj, evt-EVT_DEL, val)

def nodebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_key (obj, evt-EVT_DEL, val)

def nodebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbry_val (obj, evt-EVT_DEL, val)

def nodebry_delonenbry_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_nbry (obj,evt-EVT_DEL)
        Blender.Window.QRedrawAll()


# ---------------------------------- nbrys (given nodes IDs)

def nodebryid_settag_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbryid_nid (obj, evt-EVT_DEL, str(val))

def nodebryid_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbryid_key (obj, evt-EVT_DEL, val)

def nodebryid_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_nbryid_val (obj, evt-EVT_DEL, val)

def nodebryid_deloneebry_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_nbryid (obj,evt-EVT_DEL)
        Blender.Window.QRedrawAll()


# ---------------------------------- ebrys

def edgebry_settag_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_ebry_tag (obj, evt-EVT_DEL, str(val))

def edgebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_ebry_key (obj, evt-EVT_DEL, val)

def edgebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.set_ebry_val (obj, evt-EVT_DEL, val)

def edgebry_deloneebry_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        di.del_ebry (obj,evt-EVT_DEL)
        Blender.Window.QRedrawAll()


# ---------------------------------- fbrys

def facebry_setkey_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_DEL)
        obj.properties['fbrys'][tag][0] = val-1 # DOFVar==uz,fx,...

def facebry_setval_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_DEL)
        obj.properties['fbrys'][tag][1] = float(val) # Val


# ---------------------------------- eatts

def elematt_settype_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_DEL)
        res = di.sarray_set_val(obj.properties['eatts'][tag], 0, str(val-1))
        obj.properties['eatts'].pop(tag)
        obj.properties['eatts'][tag] = res

def elematt_setmodel_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_DEL)
        res = di.sarray_set_val(obj.properties['eatts'][tag], 1, str(val-1))
        obj.properties['eatts'].pop(tag)
        obj.properties['eatts'][tag] = res

def elematt_setprms_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_DEL)
        val = '_'.join(val.split())
        res = di.sarray_set_val(obj.properties['eatts'][tag], 2, val)
        obj.properties['eatts'].pop(tag)
        obj.properties['eatts'][tag] = res

def elematt_setinis_callback(evt,val):
    obj = di.get_obj()
    if obj!=None:
        tag = str(evt-EVT_DEL)
        val = '_'.join(val.split())
        res = di.sarray_set_val(obj.properties['eatts'][tag], 3, val)
        obj.properties['eatts'].pop(tag)
        obj.properties['eatts'][tag] = res


# ---------------------------------- Results

def show_results_callback(evt,val):
    dict = di.load_dict()
    dict['show_results'] = val
    Blender.Window.QRedrawAll()

def set_scalar_callback(evt,val):
    dict = di.load_dict()
    dict['scalar_key'] = val
    Blender.Window.QRedrawAll()

def show_scalar_callback(evt,val):
    dict = di.load_dict()
    dict['show_scalar'] = val
    Blender.Window.QRedrawAll()

def set_warp_callback(evt,val):
    dict = di.load_dict()
    dict['warp_scale'] = val
    Blender.Window.QRedrawAll()

def show_warp_callback(evt,val):
    dict = di.load_dict()
    dict['show_warp'] = val
    Blender.Window.QRedrawAll()


# ======================================================================================= GUI

def draw_label(x,y,w,h,txt):
    BGL.glColor3f     (0.663, 0.663, 0.663)
    BGL.glRecti       (x, y, x+w, y+h)
    BGL.glColor3f     (0.0, 0.0, 0.0)
    BGL.glRasterPos2i (x+5, y+5)
    Draw.Text         (txt)

# Draw GUI
def gui():
    # load dictionary
    d = di.load_dict()

    # Get current selected mesh
    try: edm, obj, msh = di.get_msh (False)
    except Exception, inst:
        edm = 0
        obj = None
        msh = None

    # Data from current selected object
    btag     = -1
    ndivx    = 2
    ndivy    = 2
    ndivz    = 1
    acoefx   = '0'
    acoefy   = '0'
    acoefz   = '0'
    nonlinx  = 0
    nonliny  = 0
    nonlinz  = 0
    minangle = '-1.0'
    maxarea  = '-1.0'
    regs     = []
    hols     = []
    nbrys    = []
    nbryids  = []
    ebrys    = []
    fbrys    = {}
    eatts    = {}
    # Set default values
    if obj!=None:
        btag     = di.get_btag     (obj)
        ndivx    = di.get_ndiv     (obj,'x')
        ndivy    = di.get_ndiv     (obj,'y')
        ndivz    = di.get_ndiv     (obj,'z')
        acoefx   = di.get_acoef    (obj,'x')
        acoefy   = di.get_acoef    (obj,'y')
        acoefz   = di.get_acoef    (obj,'z')
        nonlinx  = di.get_nonlin   (obj,'x')
        nonliny  = di.get_nonlin   (obj,'y')
        nonlinz  = di.get_nonlin   (obj,'z')
        minangle = di.get_minangle (obj)
        maxarea  = di.get_maxarea  (obj)
        regs     = di.get_regs     (obj)
        hols     = di.get_hols     (obj)
        nbrys    = di.get_nbrys    (obj)
        nbryids  = di.get_nbryids  (obj)
        ebrys    = di.get_ebrys    (obj)
        fbrys    = obj.properties['fbrys'] if obj.properties.has_key('fbrys') else {}
        eatts    = obj.properties['eatts'] if obj.properties.has_key('eatts') else {}

    # restore EditMode
    if edm: Blender.Window.EditMode(1)

    # Width and Height
    wid, hei = Blender.Window.GetAreaSize()

    # Col and Row to add entities
    gx  = 15   # x gap
    gy  = 20   # y gap
    sgy = 10   # small gy
    ggx = 2*gx 
    ggy = 2*gy
    cw  = 100  # column width
    rh  = 20   # row height
    ch  = rh+2 # caption height

    # Heights of windows
    h_set         = 3*rh+sgy+ggy
    h_cad         = 4*rh+ggy
    h_struct      = 6*rh+ggy+sgy
    h_unst_regs   = (1+len(regs))*rh+sgy
    h_unst_hols   = (1+len(hols))*rh+sgy
    h_unstruct    = ggy+3*sgy+2*rh + rh+h_unst_regs + rh+h_unst_hols
    h_mesh        = 5*rh + rh+h_struct + rh+h_unstruct
    h_fea_nbrys   = (1+len(nbrys))*rh+sgy
    h_fea_nbryids = (1+len(nbryids))*rh+sgy
    h_fea_ebrys   = (1+len(ebrys))*rh+sgy
    h_fea_fbrys   = (1+len(fbrys))*rh+sgy
    h_fea_eatts   = (1+len(eatts))*rh+sgy
    h_fea         = sgy+ggy+3*rh + rh+h_fea_nbrys + rh+h_fea_nbryids + rh+h_fea_ebrys + rh+h_fea_fbrys + rh+h_fea_eatts
    h_res         = rh+ggy

    # Background color
    BGL.glClearColor (0.531, 0.543, 0.614, 0.0)
    BGL.glClear      (Blender.BGL.GL_COLOR_BUFFER_BIT)
    

    # SETTINGS ========================================================================================================================================

    dx  = 2*gx+50;
    row = hei-gy-d['inirow']
    h   = h_set
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('SETTINGS')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,      wid-ggx-80-70, row+3, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_SET, wid-ggx-80,    row+3, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_set']:
        BGL.glColor3f     (0.85, 0.85, 0.85)
        BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Show:')
        Draw.Toggle       ('Props',     EVT_NONE, dx    , row, 60, rh, d['show_props'],   'Show mesh properties'      , show_props_callback)
        Draw.Toggle       ('V IDs',     EVT_NONE, dx+ 60, row, 60, rh, d['show_v_ids'],   'Show vertex IDs'           , show_v_ids_callback)
        Draw.Toggle       ('E IDs',     EVT_NONE, dx+120, row, 60, rh, d['show_e_ids'],   'Show edge IDs'             , show_e_ids_callback)
        Draw.Toggle       ('F IDs',     EVT_NONE, dx+180, row, 60, rh, d['show_f_ids'],   'Show face IDs'             , show_f_ids_callback)
        Draw.Slider       ('',          EVT_NONE, dx+240, row, 80, rh, d['ftags_opac'],0.0,1.0,0,'Set opacitity to paint faces with tags', ftags_opac_callback); row -= rh
        Draw.Toggle       ('Axes' ,     EVT_NONE, dx    , row, 60, rh, d['show_axes'],    'Show local system axes'    , show_axes_callback )
        Draw.Toggle       ('ETags',     EVT_NONE, dx+ 60, row, 60, rh, d['show_etags'],   'Show edge tags'            , show_etags_callback)
        Draw.Toggle       ('FTags',     EVT_NONE, dx+120, row, 60, rh, d['show_ftags'],   'Show face tags'            , show_ftags_callback)
        Draw.Toggle       ('Elems',     EVT_NONE, dx+180, row, 60, rh, d['show_elems'],   'Show elements information' , show_elems_callback)
        Draw.Toggle       ('Tags text', EVT_NONE, dx+240, row, 80, rh, d['show_tags_txt'],'Show label (text) of tags' , show_tags_txt_callback); row -= rh+sgy
        Draw.PushButton   ('Delete all properties', EVT_SET_DELPROPS, ggx, row, 200, rh , 'Delete all properties')
        row -= ggy

    else: row -= gy


    # CAD ==============================================================================================================================================

    h = h_cad
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('CAD')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,      wid-ggx-80-70, row+2, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_CAD, wid-ggx-80,    row+2, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_cad']:
        BGL.glColor3f     (0.85, 0.85, 0.85)
        BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Points:')
        Draw.String       ('X=',        EVT_NONE      , dx    , row, 80, rh, d['newpoint_x'],128,'Set the X value of the new point to be added',setx_callback)
        Draw.String       ('Y=',        EVT_NONE      , dx+ 80, row, 80, rh, d['newpoint_y'],128,'Set the Y value of the new point to be added',sety_callback)
        Draw.String       ('Z=',        EVT_NONE      , dx+160, row, 80, rh, d['newpoint_z'],128,'Set the Z value of the new point to be added',setz_callback)
        Draw.PushButton   ('Add x y z', EVT_CAD_ADDXYZ, dx+240, row, 80, rh, 'Add point from coordinates'); row -= rh
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Edges:');
        Draw.String       ('Radius=', EVT_NONE      , dx     , row , 120 , rh , d['fillet_radius'],  128,'Set radius for fillet operation',fillet_radius_callback)
        Draw.Number       ('Steps=',  EVT_NONE      , dx+120 , row , 120 , rh , d['fillet_steps' ],1,100,'Set steps for fillet operation' ,fillet_steps_callback)
        Draw.PushButton   ('Fillet',  EVT_CAD_FILLET, dx+240 , row ,  80 , rh , 'Create a fillet between two edges'); row -= rh
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Edges:');
        Draw.PushButton   ('Break edge',        EVT_CAD_EDGEBREAK , dx,      row ,  120, rh , 'Break an edge at a previously selected point')
        Draw.PushButton   ('Break at middle',   EVT_CAD_EDGEBREAKM, dx+120 , row ,  120, rh , 'Break an edge at its middle point')
        Draw.PushButton   ('Edge intersection', EVT_CAD_EDGEINTERS, dx+240 , row ,  120, rh , 'Find the intersection (smaller distance) between two edges'); row -= rh
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('File:')
        Draw.PushButton   ('Read Points', EVT_CAD_FPOINT,  dx,     row, 120, rh, 'Add points by reading a list from file')
        Draw.PushButton   ('Read Spline', EVT_CAD_FSPLINE, dx+120, row, 120, rh, 'Add a spline by reading its points from file')
        dx   = 2*gx+55;
        row -= ggy

    else: row -= gy


    # Mesh ==============================================================================================================================================

    h = h_mesh
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('MESH')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,       wid-ggx-80-70, row+2, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_MESH, wid-ggx-80,    row+2, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_mesh']:
        nftag = d['newftag'][0]
        nfclr = d['newftag'][1]
        BGL.glColor3f     (0.85, 0.85, 0.85)
        BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Set tags:')
        Draw.Number       ('',     EVT_NONE,          dx,     row, 80, rh, d['newetag'], -100, 0, 'New edge tag',etag_callback)
        Draw.PushButton   ('Edge', EVT_MESH_SETETAG,  dx+80,  row, 80, rh,                        'Set edges tag (0 => remove tag)')
        Draw.Number       ('',     EVT_NONE,          dx+160, row, 80, rh, nftag,       -1000, 0, 'New face tag', ftag_callback)
        Draw.ColorPicker  (        EVT_NONE,          dx+240, row, 80, rh, di.hex2rgb(nfclr),     'Select color to paint tagged face',   fclr_callback)
        Draw.PushButton   ('Face', EVT_MESH_SETFTAG,  dx+320, row, 80, rh,                        'Set faces tag (0 => remove tag)'); row -= rh

        # Mesh -- structured
        dx   = 2*gx+70
        ggx += gx
        dx  += gx
        h    = h_struct
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Structured')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h); row -= rh+gy
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Local Axes:')
        Draw.PushButton   ('Set x', EVT_MESH_SETX, dx    , row, 60, rh , 'Set local x axis (one edge must be previously selected)')
        Draw.PushButton   ('Set y', EVT_MESH_SETY, dx+ 60, row, 60, rh , 'Set local y axis (one edge must be previously selected)')
        Draw.PushButton   ('Set z', EVT_MESH_SETZ, dx+120, row, 60, rh , 'Set local z axis (one edge must be previously selected)'); row -= rh
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Divisions:')
        Draw.Number       ('nX=',    EVT_NONE, dx,     row, 80, rh, ndivx,  1, 1000,'Set number of divisions along x (local axis)',ndivx_callback)
        Draw.Number       ('nY=',    EVT_NONE, dx+ 80, row, 80, rh, ndivy,  1, 1000,'Set number of divisions along y (local axis)',ndivy_callback)
        Draw.Number       ('nZ=',    EVT_NONE, dx+160, row, 80, rh, ndivz,  1, 1000,'Set number of divisions along z (local axis)',ndivz_callback); row -= rh
        Draw.String       ('aX=',    EVT_NONE, dx    , row, 80, rh, acoefx, 32,'Set the a coeficient of the divisions function for the x direction (0 => equal size divisions)',acoefx_callback)
        Draw.String       ('aY=',    EVT_NONE, dx+ 80, row, 80, rh, acoefy, 32,'Set the a coeficient of the divisions function for the x direction (0 => equal size divisions)',acoefy_callback)
        Draw.String       ('aZ=',    EVT_NONE, dx+160, row, 80, rh, acoefz, 32,'Set the a coeficient of the divisions function for the x direction (0 => equal size divisions)',acoefz_callback); row -= rh
        Draw.Toggle       ('NonLinX',EVT_NONE, dx,     row, 80, rh, nonlinx, 'Set non-linear divisions along x', nonlinx_callback)
        Draw.Toggle       ('NonLinY',EVT_NONE, dx+ 80, row, 80, rh, nonliny, 'Set non-linear divisions along y', nonliny_callback)
        Draw.Toggle       ('NonLinZ',EVT_NONE, dx+160, row, 80, rh, nonlinz, 'Set non-linear divisions along z', nonlinz_callback); row -= rh
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Block:')
        Draw.Number       ('Tag=', EVT_NONE, dx, row, 100, rh, btag, -100, 0,'Set the block tag (to be replicated to all elements)',btag_callback); row -= rh+sgy
        Draw.PushButton   ('Generate (quadrilaterals/hexahedrons)', EVT_MESH_GENSTRU,  ggx, row, 300, rh, 'Generated structured mesh')

        # Mesh -- unstructured
        row -= ggy
        dx   = 2*gx+70;
        dx  += gx
        h    = h_unstruct
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Unstructured')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h); row -= rh+gy
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Quality:')
        Draw.String       ('q=',  EVT_NONE, dx,     row, 100, rh, minangle, 128, 'Set the minimum angle between edges/faces (-1 => use default)',                       minangle_callback)
        Draw.String       ('a=',  EVT_NONE, dx+100, row, 100, rh, maxarea,  128, 'Set the maximum area/volume (uniform) for triangles/tetrahedrons (-1 => use default)',maxarea_callback)

        # Mesh -- unstructured -- regs
        row -= sgy
        ggx += gx
        h    = h_unst_regs
        BGL.glColor3f     (0.53, 0.54, 0.6)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Regions')
        BGL.glColor3f     (0.82, 0.82, 0.9)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_MESH_ADDREG,     wid-ggx-70-80, row+2, 60, rh-4, 'Add region')
        Draw.PushButton   ('Delete all', EVT_MESH_DELALLREGS, wid-ggx-80,    row+2, 80, rh-4, 'Delete all regions'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('     Tag        Max area            X                  Y                  Z'); row -= rh
        for i, reg in enumerate(regs):
            Draw.Number     ('',    EVT_DEL+i, ggx    , row, 60, rh, int(reg[0]), -100, 0,'Region tag',                  regs_tag_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+ 60, row, 80, rh,     reg[1],   128,   'Max area (-1 => use default)',regs_maxarea_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+140, row, 80, rh,     reg[2],   128,   'X of the rege',               regs_setx_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+220, row, 80, rh,     reg[3],   128,   'Y of the rege',               regs_sety_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+300, row, 80, rh,     reg[4],   128,   'Z of the rege',               regs_setz_callback)
            Draw.PushButton ('Del', EVT_DEL+i, ggx+380, row, 40, rh,                      'Delete this row',             regs_delonereg_callback); row -= rh

        # Mesh -- unstructured -- holes
        h = h_unst_hols
        BGL.glColor3f     (0.53, 0.54, 0.6)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Holes')
        BGL.glColor3f     (0.82, 0.82, 0.9)
        BGL.glRecti       (3*gx, row, wid-3*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_MESH_ADDHOLE,     wid-ggx-70-80, row+2, 60, rh-4, 'Add hole')
        Draw.PushButton   ('Delete all', EVT_MESH_DELALLHOLES, wid-ggx-80,    row+2, 80, rh-4, 'Delete all holes'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('         X                  Y                  Z'); row -= rh
        for i, hol in enumerate(hols):
            Draw.String     ('',    EVT_DEL+i, ggx    , row, 100, rh, hol[0],128,'X of the hole',holes_setx_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+100, row, 100, rh, hol[1],128,'Y of the hole',holes_sety_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+200, row, 100, rh, hol[2],128,'Z of the hole',holes_setz_callback)
            Draw.PushButton ('Del', EVT_DEL+i, ggx+300, row,  40, rh, 'Delete this row', holes_delonehole_callback); row -= rh

        # Mesh -- unstructured -- main
        row -= 2*sgy
        ggx -= gx
        Draw.PushButton ('Generate (triangles/tetrahedrons)', EVT_MESH_GENUNSTRU ,  ggx, row, 300, rh, 'Generated unstructured mesh')

        # gaps/increments
        dx   = 2*gx+55;
        ggx  = 2*gx
        row -= ggy+gy

    else: row -= gy


    # FEM ==============================================================================================================================================

    h = h_fea
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('FEM')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,      wid-ggx-80-70, row+2, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_FEM, wid-ggx-80,    row+2, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_fem']:
        BGL.glColor3f     (0.8, 0.8, 0.8)
        BGL.glRecti       (gx, row, wid-gx, row-h)

        # FEM -- nodes bry (given coordinates)
        ggx += gx
        row -= gy
        h    = h_fea_nbrys
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Nodes boundaries (given X-Y-Z coordinates)')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_FEM_ADDNBRY,    wid-ggx-70-80, row+2, 60, rh-4, 'Add nodes boundary (given X-Y-Z coordinates)')
        Draw.PushButton   ('Delete all', EVT_FEM_DELALLNBRY, wid-ggx-80,    row+2, 80, rh-4, 'Delete all nodes boundary (X-Y-Z)'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('         X                  Y                  Z           Key        Value'); row -= rh
        for i, nb in enumerate(nbrys):
            Draw.String     ('',    EVT_DEL+i, ggx    , row, 80, rh, nb[0],128,'X of the node with boundary condition',                                      nodebry_setx_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+ 80, row, 80, rh, nb[1],128,'Y of the node with boundary condition',                                      nodebry_sety_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+160, row, 80, rh, nb[2],128,'Z of the node with boundary condition',                                      nodebry_setz_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+240, row, 40, rh, nb[3],  2,'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', nodebry_setkey_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+280, row, 80, rh, nb[4],128,'Value of essential/natural boundary condition',                              nodebry_setval_callback)
            Draw.PushButton ('Del', EVT_DEL+i, ggx+360, row, 40, rh, 'Delete this row',                                                                      nodebry_delonenbry_callback); row -= rh

        # FEM -- nodes bry (given nodes ID)
        h = h_fea_nbryids
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Nodes boundaries (given nodes ID)')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_FEM_ADDNBRYID,    wid-ggx-70-80, row+2, 60, rh-4, 'Add nodes boundary (given node ID)')
        Draw.PushButton   ('Delete all', EVT_FEM_DELALLNBRYID, wid-ggx-80,    row+2, 80, rh-4, 'Delete all nodes boundary (node ID)'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('         ID         Key        Value'); row -= rh
        for i, nbid in enumerate(nbryids):
            Draw.Number     ('',    EVT_DEL+i, ggx    , row, 80, rh, int(nbid[0]),0,10000,'Set the node ID',                                                           nodebryid_settag_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+ 80, row, 40, rh, nbid[1],  2,         'Key such as ux, uy, fx, fz corresponding to the essential/natural variable',nodebryid_setkey_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+120, row, 80, rh, nbid[2],128,         'Value of essential/natural boundary condition',                             nodebryid_setval_callback)
            Draw.PushButton ('Del', EVT_DEL+i, ggx+200, row, 40, rh, 'Delete this row',                                                                               nodebryid_deloneebry_callback); row -= rh

        # FEM -- edges bry
        h = h_fea_ebrys
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Edges boundaries')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0)
        Draw.PushButton   ('Add',        EVT_FEM_ADDEBRY,    wid-ggx-70-80, row+2, 60, rh-4, 'Add edges boundary')
        Draw.PushButton   ('Delete all', EVT_FEM_DELALLEBRY, wid-ggx-80,    row+2, 80, rh-4, 'Delete all edges boundary'); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('        Tag         Key        Value'); row -= rh
        for i, eb in enumerate(ebrys):
            Draw.Number     ('',    EVT_DEL+i, ggx    , row, 80, rh, int(eb[0]),-1000,0,'Set the edge tag',                                                          edgebry_settag_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+ 80, row, 40, rh, eb[1],  2,         'Key such as ux, uy, fx, fz corresponding to the essential/natural variable',edgebry_setkey_callback)
            Draw.String     ('',    EVT_DEL+i, ggx+120, row, 80, rh, eb[2],128,         'Value of essential/natural boundary condition',                             edgebry_setval_callback)
            Draw.PushButton ('Del', EVT_DEL+i, ggx+200, row, 40, rh, 'Delete this row',                                                                              edgebry_deloneebry_callback); row -= rh

        # FEM -- faces bry
        h = h_fea_fbrys
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Faces boundaries')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('        Tag         Key        Value'); row -= rh
        for k, v in fbrys.iteritems():
            tag = int(k)
            clr = di.hex2rgb(v[2])
            draw_label        (ggx,row,40,rh, k)
            BGL.glColor3f     (clr[0], clr[1], clr[2])
            BGL.glRecti       (ggx+40, row, ggx+80, row+rh)
            Draw.Menu   (d['dofvars_menu'], EVT_DEL+tag, ggx+ 80, row, 40, rh, int(v[0])+1,    'Key such as ux, uy, fx, fz corresponding to the essential/natural variable', facebry_setkey_callback)
            Draw.String ('',                EVT_DEL+tag, ggx+120, row, 80, rh, str(v[1]), 128, 'Value of essential/natural boundary condition',                              facebry_setval_callback); row -= rh

        # FEM -- elems atts
        h = h_fea_eatts
        BGL.glColor3f     (0.431, 0.443, 0.514)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-rh); row -= rh
        BGL.glColor3f     (1.0, 1.0, 1.0)
        BGL.glRasterPos2i (ggx, row+5)
        Draw.Text         ('Elements attributes')
        BGL.glColor3f     (0.72, 0.72, 0.8)
        BGL.glRecti       (2*gx, row, wid-2*gx, row-h);
        BGL.glColor3f     (0.0, 0.0, 0.0); row -= rh
        BGL.glRasterPos2i (ggx, row+5); Draw.Text('      Tag           Type                Model           Parameters        Initial Vals'); row -= rh
        for k, v in eatts.iteritems():
            tag = int(k)
            r   = v.split()
            draw_label      (ggx,row,40,rh, k)
            Draw.Menu       (d['etypes_menu'], EVT_DEL+tag, ggx+ 40, row, 120, rh, int(r[0])+1,               'Element type: ex.: Quad4PStrain',           elematt_settype_callback)
            Draw.Menu       (d['models_menu'], EVT_DEL+tag, ggx+160, row, 100, rh, int(r[1])+1,               'Constitutive model: ex.: LinElastic',       elematt_setmodel_callback)
            Draw.String     ('',               EVT_DEL+tag, ggx+260, row, 100, rh, r[2].replace('_',' '),128, 'Parameters: ex.: E=200 nu=0.25',            elematt_setprms_callback)
            Draw.String     ('',               EVT_DEL+tag, ggx+360, row,  80, rh, r[3].replace('_',' '),128, 'Initial values: ex.: Sx=0 Sy=0 Sz=0 Sxy=0', elematt_setinis_callback); row -= rh

        # FEM -- main
        row -= 2*sgy
        ggx -= gx
        Draw.PushButton ('Run analysis',     EVT_FEM_RUN,      ggx,     row, 150, rh, 'Run a FE analysis directly (without script)')
        Draw.PushButton ('Generate script',  EVT_FEM_SCRIPT,   ggx+150, row, 150, rh, 'Generate script for FEM')
        Draw.PushButton ('View in ParaView', EVT_FEM_PARAVIEW, ggx+300, row, 150, rh, 'View results in ParaView'); row -= rh

        # increments
        row = row + rh - ggy

    else: row -= gy


    # Results ===========================================================================================================================================

    h = h_res
    BGL.glColor3f     (0.4, 0.4, 0.4)
    BGL.glRecti       (gx, row, wid-gx, row-ch); row -= ch
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (ggx, row+5)
    Draw.Text         ('RESULTS')
    Draw.PushButton   ('Refresh',   EVT_REFRESH,      wid-ggx-80-70, row+2, 60, rh-4, 'Refresh GUI')
    Draw.PushButton   ('Show/Hide', EVT_SHOWHIDE_RES, wid-ggx-80,    row+2, 80, rh-4, 'Show/Hide this box')
    if d['gui_show_res']:
        BGL.glColor3f     (0.85, 0.85, 0.85)
        BGL.glRecti       (gx, row, wid-gx, row-h); row -= gy+rh
        BGL.glColor3f     (0.0, 0.0, 0.0)
        BGL.glRasterPos2i (ggx, row+4)
        Draw.Text         ('Show:')
        Draw.Toggle       ('Results',EVT_NONE, dx    , row, 60, rh, d['show_results'],     'Show results'               , show_results_callback)
        Draw.String       ('Val='   ,EVT_NONE, dx+60 , row, 60, rh, d['scalar_key']  , 32, 'Set value to be shown'      , set_scalar_callback)
        Draw.Toggle       ('Scalar' ,EVT_NONE, dx+120, row, 60, rh, d['show_scalar'] ,     'Show scalar values'         , show_scalar_callback)
        Draw.String       ('M='     ,EVT_NONE, dx+180, row, 60, rh, d['warp_scale']  , 32, 'Set warp (deformed) scale'  , set_warp_callback)
        Draw.Toggle       ('Warp'   ,EVT_NONE, dx+240, row, 60, rh, d['show_warp']   ,     'Show warped (deformed) mesh', show_warp_callback)


# Register GUI
Draw.Register (gui, event, button_event)


# ================================================================================ ScriptLink

# Load script for View3D drawing
dict       = di.load_dict()
script_lnk = 'msys_3dlink.py'
script_dir = Blender.Get('scriptsdir')+Blender.sys.sep
if script_lnk in [t.name for t in Blender.Text.Get()]:
    Blender.Text.unlink(Blender.Text.Get(script_lnk))
print '[1;34mMechSys[0m: loading <'+script_dir+script_lnk+'>'
Blender.Text.Load(script_dir+script_lnk)


# Connect View3D drawing script
scn = bpy.data.scenes.active
scn.clearScriptLinks()
if scn.getScriptLinks('Redraw')==None: 
    scn.addScriptLink(script_lnk,'Redraw')
elif script_lnk not in scn.getScriptLinks('Redraw'): 
    scn.addScriptLink(script_lnk,'Redraw')


# ============================================================================== SpaceHandler

# Load script for View3D space handler
script_shandler = 'msys_shandler.py'
if script_shandler in [t.name for t in Blender.Text.Get()]:
    Blender.Text.unlink(Blender.Text.Get(script_shandler))
print '[1;34mMechSys[0m: loading <'+script_dir+script_shandler+'>'
Blender.Text.Load(script_dir+script_shandler)

# TODO: Connect Space Handler to View3D
