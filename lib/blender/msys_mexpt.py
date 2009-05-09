#!BPY

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

"""
Name: 'MeshExport'
Blender: 2.46
Group: 'Themes'
Tooltip: 'Mesh export'
"""
__author__  = "Dorival Pedroso"
__version__ = "1.0.0"
__bpydoc__  = """TODO"""

import os
from   Blender import *
from   bpy     import *

############################################### get object and mesh

# get active mesh object
def get_active_mesh():
    scn = data.scenes.active
    obj = scn.objects.active
    if obj.type=='Mesh':
        edm = Window.EditMode()
        if edm: Window.EditMode(0)
        msh = obj.getData(mesh=1)
        return obj, msh, edm
    else: return None, None, 0

# get selected edges
def get_selected_edges(msh):
    sel = []
    for eid in msh.edges.selected():
        if msh.edges[eid].sel==True: sel.append(eid)
    return sel

############################################################ Export

# export mesh to a file
def export_mesh(fn,outedg):
    obj, msh, edm = get_active_mesh()
    if msh!=None:

        # check
        for f in msh.faces:
            nv = len(f.verts)
            va = f.verts[0]
            for i in range(1,nv-1):
                vb = f.verts[i]
                vc = f.verts[i+1]
                a  = vb.co-va.co
                b  = vc.co-va.co
                c  = Mathutils.CrossVecs(a, b)
                # check if face is parallel to x-y a plane
                if abs(c[0])>0.00001 or abs(c[1])>0.00001:
                    if edm: Window.EditMode(1)
                    raise Exception('All faces must be parallel to x-y plane')
                # check connectivity
                if c[1]<0.0:
                    if edm: Window.EditMode(1)
                    raise Exception('Mesh connectivity is invalid (nodes numbering order is incorrect)')

        # open file
        fil = open(fn,'w')

        # vertices
        lin = 'V = ['
        nv1 = len(msh.verts)-1
        for i, v in enumerate(msh.verts):
            # tag
            tag = 0
            if obj.properties.has_key('vtags'):
                vid = str(v.index)
                if obj.properties['vtags'].has_key(vid):
                    tag = obj.properties['vtags'][vid]
            # line
            s1  = ''    if i==0   else '     '
            s2  = ']\n' if i==nv1 else ',\n'
            lin = '%s%s[%3d,%14.6e,%14.6e, %4d]%s' % (lin,s1,v.index,v.co[0],v.co[1],tag,s2)
        fil.write(lin)

        # edges
        if outedg:
            lin = '\nE = ['
            ne1 = len(msh.edges)-1
            for i, e in enumerate(msh.edges):
                # tag
                tag = 0
                if obj.properties.has_key('etags'):
                    eid = str(e.index)
                    if obj.properties['etags'].has_key(eid):
                        tag = obj.properties['etags'][eid]
                # line
                s1  = ''    if i==0   else '     '
                s2  = ']\n' if i==ne1 else ',\n'
                lin = '%s%s[%3d,%3d,%3d, %4d]%s' % (lin,s1,e.index,e.v1.index,e.v2.index,tag,s2)
            fil.write(lin)

        # solids
        lin = '\nS = ['
        nf1 = len(msh.faces)-1
        for i, f in enumerate(msh.faces):
            # tag
            tag = 0
            if obj.properties.has_key('stags'):
                sid = str(f.index)
                if obj.properties['stags'].has_key(sid):
                    tag = obj.properties['stags'][sid]
            # nodes
            ln1  = ' ('
            fnv1 = len(f.verts)-1
            for j, v in enumerate(f.verts):
                ss  = ')' if j==fnv1 else ','
                ln1 = '%s%3d%s' % (ln1,v.index,ss)
            # edges
            edgs = msh.findEdges(f.edge_keys)
            ln2  = ' ('
            fne1 = len(edgs)-1
            for j, e in enumerate(edgs):
                ss  = ')' if j==fne1 else ','
                if outedg: ln2 = '%s%3d%s' % (ln2,e,ss)
                else:
                    etag = 0
                    if obj.properties.has_key('etags'):
                        eid = str(e)
                        if obj.properties['etags'].has_key(eid):
                            etag = obj.properties['etags'][eid]
                    ln2 = '%s%4d%s' % (ln2,etag,ss)
            # line
            s1  = ''    if i==0   else '     '
            s2  = ']\n' if i==nf1 else ',\n'
            lin = '%s%s[%3d,%s,%s, %4d]%s' % (lin,s1,f.index,ln1,ln2,tag,s2)
        fil.write(lin)

        # close file
        fil.close()
        if edm: Window.EditMode(1)

######################################################## Dictionary

# dictionary to hold GUI data
def load_dict():
    d = Registry.GetKey('gui_dict')
    if not d:
        d            = {}
        d['show']    = True
        d['showvid'] = True
        d['showeid'] = True
        d['showsid'] = True
        d['showvta'] = True
        d['showeta'] = True
        d['showsta'] = True
        d['newvta']  = -100
        d['neweta']  = -10
        d['newsta']  = -1
        d['outedg']  = False
        Registry.SetKey('gui_dict', d)
    return d

def set_key(key,val):
    Registry.GetKey('gui_dict')[key] = val
    Window.QRedrawAll()

def get_key(key):
    return Registry.GetKey('gui_dict')[key]

######################################################## Properties

def set_tag(obj,key,id,tag):
    print id, tag
    sid = str(id)
    if not obj.properties.has_key(key): obj.properties[key] = {}
    if tag==0: obj.properties[key].pop    (sid)
    else:      obj.properties[key].update({sid:tag})
    if len(obj.properties[key])==0: obj.properties.pop(key)

######################################################### try_catch

def TC(func):
    def wrapper(*arg):
        try: func(*arg)
        except Exception, inst:
            msg = inst.args[0]
            print 'ERROR: '+msg
            Draw.PupMenu('ERROR|'+msg)
    return wrapper

############################################################### GUI

# constants
EVT_NONE  = 0 # none
EVT_SHOW  = 1 # show
EVT_SHWET = 2 # show edge tags
EVT_VTAG  = 3 # tag vertex
EVT_ETAG  = 4 # tag edge
EVT_STAG  = 5 # tag solid
EVT_EXPT  = 6 # export

# events
def event(evt, val):
    if   evt==Draw.QKEY and not val: Draw.Exit()
    elif evt==Draw.ESCKEY: Draw.Redraw(1)

# button events
@TC
def bevent(evt):
    if evt==EVT_VTAG: # ===========================================
        obj, msh, edm = get_active_mesh()
        if msh!=None:
            for vid in msh.verts.selected():
                set_tag(obj,'vtags',vid,get_key('newvta'))
            if edm: Window.EditMode(1)
            Window.QRedrawAll()
    elif evt==EVT_ETAG: # =========================================
        obj, msh, edm = get_active_mesh()
        if msh!=None:
            for eid in get_selected_edges(msh):
                set_tag(obj,'etags',eid,get_key('neweta'))
            if edm: Window.EditMode(1)
            Window.QRedrawAll()
    elif evt==EVT_STAG: # =========================================
        obj, msh, edm = get_active_mesh()
        if msh!=None:
            for fid in msh.faces.selected():
                set_tag(obj,'stags',fid,get_key('newsta'))
            if edm: Window.EditMode(1)
            Window.QRedrawAll()
    elif evt==EVT_EXPT: # =========================================
        Window.FileSelector(cb_expmsh, 'Export mesh', 'outmesh.py')

# callbacks
def cb_show   (evt,val): set_key('show',   val)
def cb_showvid(evt,val): set_key('showvid',val)
def cb_showeid(evt,val): set_key('showeid',val)
def cb_showsid(evt,val): set_key('showsid',val)
def cb_showvta(evt,val): set_key('showvta',val)
def cb_showeta(evt,val): set_key('showeta',val)
def cb_showsta(evt,val): set_key('showsta',val)
def cb_vta    (evt,val): set_key('newvta', val)
def cb_eta    (evt,val): set_key('neweta', val)
def cb_sta    (evt,val): set_key('newsta', val)
def cb_outedg (evt,val): set_key('outedg', val)
def cb_expmsh (fn):
    if os.path.exists(fn):
        msg = 'Overwrite file <'+fn+'> ?%t|Yes'
        res = Draw.PupMenu(msg)
        if res>0: export_mesh(fn,get_key('outedg'))
    else: export_mesh(fn,get_key('outedg'))

# window
def gui():
    # load dictionary
    d = load_dict()

    # constants
    rg  = 20 # row gap
    cg  = 14 # column gap
    rh  = 20 # row height
    srh = 10 # small row height

    # window size
    W,H = Window.GetAreaSize()
    r   = H-rh-rg
    c   = cg
    w   = W-2*c

    # background
    BGL.glClearColor (0.725, 0.725, 0.725, 0.0)
    BGL.glClear      (BGL.GL_COLOR_BUFFER_BIT)

    # draw SETTINGS box ===========================================
    h = 3*rh
    BGL.glColor3f     (0.573, 0.512, 0.585)
    BGL.glRecti       (c, r, c+w, r+rh)
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         ('SETTINGS')
    BGL.glColor3f     (0.656, 0.696, 0.617)
    BGL.glRecti       (c,r-h,c+w,r)

    # draw SETTINGS widgets
    c += cg
    r -= (rh+srh)
    Draw.Toggle ('Show',   EVT_NONE, c,r-rh, 60, 2*rh, d['show'],   'Show mesh properties', cb_show)
    Draw.Toggle ('V IDs',  EVT_NONE, c+ 60, r, 60, rh, d['showvid'], 'Show vertices IDs',   cb_showvid)
    Draw.Toggle ('E Ids',  EVT_NONE, c+120, r, 60, rh, d['showeid'], 'Show edges IDs',      cb_showeid)
    Draw.Toggle ('S IDs',  EVT_NONE, c+180, r, 60, rh, d['showsid'], 'Show solids IDs',     cb_showsid)
    r -= rh
    Draw.Toggle ('V Tags', EVT_NONE, c+ 60, r, 60, rh, d['showvta'], 'Show vertices tags',  cb_showvta)
    Draw.Toggle ('E Tags', EVT_NONE, c+120, r, 60, rh, d['showeta'], 'Show edges tags',     cb_showeta)
    Draw.Toggle ('S Tags', EVT_NONE, c+180, r, 60, rh, d['showsta'], 'Show solids tags',    cb_showsta)

    # draw MESH box ===============================================
    h  = 4*rh
    c -= cg
    r -= (srh+2*rh)
    BGL.glColor3f     (0.573, 0.512, 0.585)
    BGL.glRecti       (c, r, c+w, r+rh)
    BGL.glColor3f     (1.0, 1.0, 1.0)
    BGL.glRasterPos2i (c+5, r+5)
    Draw.Text         ('MESH')
    BGL.glColor3f     (0.656, 0.696, 0.617)
    BGL.glRecti       (c,r-h,c+w,r)

    # draw MESH widgets
    c += cg
    r -= (rh+srh)
    Draw.Number     ('',          EVT_NONE, c   , r, 60, rh, d['newvta'], -1000, 0, 'New vertex tag', cb_vta)
    Draw.PushButton ('Set V tag', EVT_VTAG, c+60, r, 60, rh,                        'Set vertices tag (0 => remove tag)')
    r -= rh
    Draw.Number     ('',          EVT_NONE, c   , r, 60, rh, d['neweta'], -1000, 0, 'New edge tag', cb_eta)
    Draw.PushButton ('Set E tag', EVT_ETAG, c+60, r, 60, rh,                        'Set edges tag (0 => remove tag)')
    r -= rh
    Draw.Number     ('',          EVT_NONE, c   , r, 60, rh, d['newsta'], -1000, 0, 'New solid tag', cb_sta)
    Draw.PushButton ('Set S tag', EVT_STAG, c+60, r, 60, rh,                        'Set solids tag (0 => remove tag)')

    r += rh
    Draw.Toggle     ('Out edges', EVT_NONE, c+120+cg,r, 80, rh, d['outedg'], 'Output edges when exporting ?', cb_outedg)
    r -= rh
    Draw.PushButton ('Export',    EVT_EXPT, c+120+cg,r, 80, rh, 'Export mesh')

# register window
Draw.Register (gui, event, bevent)

# load dictionary
d = load_dict()

# load script-link in case it is not available
lnk = 'msys_mex3dlink.py'
if not lnk in [t.name for t in Text.Get()]:
    sdir = Get('scriptsdir')+sys.sep
    print sdir+lnk
    Text.Load(sdir+lnk)

# connect View3D drawing script
scn = data.scenes.active
scn.clearScriptLinks()
if scn.getScriptLinks('Redraw')==None:
    scn.addScriptLink(lnk,'Redraw')
elif lnk not in scn.getScriptLinks('Redraw'):
    scn.addScriptLink(lnk,'Redraw')
