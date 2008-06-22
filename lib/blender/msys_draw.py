# Modules
import Blender
import bpy
import math
import msys_dict
import mechsys as ms

def add_point(xyz):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj==None:
        msh = bpy.data.meshes.new('points')
        obj = scn.objects.new(msh,'points')
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        msh.verts.extend(xyz)
        if edm: Blender.Window.EditMode(1)
        obj.select(1)
        Blender.Window.RedrawAll()
    else:
        Blender.Draw.PupMenu('ERROR%t|Select a Mesh object before calling this function')

def add_points_from_file(filename):
    Blender.Window.WaitCursor(1)
    edm  = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0) # exit edm
    key  = Blender.sys.basename(Blender.sys.splitext(filename)[0])
    file = open(filename, 'r')
    msh  = bpy.data.meshes.new('points')
    scn  = bpy.data.scenes.active
    obj  = scn.objects.new(msh,key)
    for line in file.readlines():
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif words[0]=='x' or words[0]=='X': pass
        else:
            x, y, z = float(words[0]), float(words[1]), float(words[2])
            msh.verts.extend(x,y,z)
    if edm: Blender.Window.EditMode(1) # enter edm
    Blender.Window.RedrawAll()
    Blender.Window.WaitCursor(0)

def add_spline_from_file(filename):
    Blender.Window.WaitCursor(1)
    edm  = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0) # exit edm
    key  = Blender.sys.basename(Blender.sys.splitext(filename)[0])
    file = open(filename, 'r')
    scn  = bpy.data.scenes.active
    spl  = bpy.data.curves.new(key,'Curve')
    obj  = scn.objects.new(spl,key)
    first = 1
    for line in file.readlines():
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif words[0]=='x' or words[0]=='X': pass
        else:
            point = Blender.BezTriple.New((float(words[0]), float(words[1]), float(words[2])))
            if first==1:
                spl.appendNurb(point)
                path  = spl[0]
                first = 0
            else:
                path.append(point)
    for point in path:
        point.handleTypes = [Blender.BezTriple.HandleTypes.AUTO, Blender.BezTriple.HandleTypes.AUTO]
    spl.update()
    if edm: Blender.Window.EditMode(1) # enter edm
    Blender.Window.RedrawAll()
    Blender.Window.WaitCursor(0)

def distance(p1,p2):
    d = p1 - p2
    return math.sqrt(math.pow(d[0],2)+math.pow(d[1],2)+math.pow(d[2],2))

def closest(pt,pointlist):
    idx = 0
    min = distance(pt,pointlist[idx])
    for i in range(1,len(pointlist)):
        d = distance(pt,pointlist[i])
        if d<min:
            min = d
            idx = i
    return min, idx

def arc_point(msh,cen,sp,ep,steps):
    dr1 = sp-cen
    dr2 = ep-cen
    axi = Blender.Mathutils.CrossVecs(dr2,dr1)
    ang = Blender.Mathutils.AngleBetweenVecs(dr2,dr1)
    rot = Blender.Mathutils.RotationMatrix(ang/steps, 3, 'r', axi)
    for i in range(steps-1):
        pt  = cen+rot*dr1
        msh.verts.extend(pt)
        msh.edges.extend(msh.verts[len(msh.verts)-2],msh.verts[len(msh.verts)-1])
        dr1 = pt-cen

def fillet(radius,steps):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        if len(msh.edges.selected())==2:
            # points (vertices)
            v1 = msh.edges[msh.edges.selected()[0]].v1
            v2 = msh.edges[msh.edges.selected()[0]].v2
            v3 = msh.edges[msh.edges.selected()[1]].v1
            v4 = msh.edges[msh.edges.selected()[1]].v2
            # intersection
            i1, i2 = Blender.Mathutils.LineIntersect(v1.co,v2.co,v3.co,v4.co)
            # fillet
            if i1!=i2: Blender.Draw.PupMenu('ERROR%t|These two edges do not intersect')
            else:
                #
                #                   ep  _,--* p2
                #                  _,-*'
                #     a=|ep-p1|,--' /
                #        _,--'     /radius  cen
                #    p1 *_--------|----------*-  d=|c-p1|
                #         `--,_    \
                #              `--,_\
                #                   `-*,_
                #                   sp   `--* p3
                di, i = closest(i1,[v1.co,v2.co])
                dj, j = closest(i1,[v3.co,v4.co])
                if i==0:
                    v1.co = i1
                    p1,p2 = v1,v2
                else:
                    v2.co = i1
                    p1,p2 = v2,v1
                if j==0:
                    msh.edges[msh.edges.selected()[1]].v1 = p1
                    if v3!=p1: msh.verts.delete(v3)
                    p3 = msh.edges[msh.edges.selected()[1]].v2
                else:
                    msh.edges[msh.edges.selected()[1]].v2 = p1
                    if v4!=p1: msh.verts.delete(v4)
                    p3 = msh.edges[msh.edges.selected()[1]].v1
                if radius>0.0:
                    # vectors along the edges
                    e1  = p2.co-p1.co; e1.normalize()
                    e2  = p3.co-p1.co; e2.normalize()
                    alp = Blender.Mathutils.AngleBetweenVecs(e1,e2)/2.0
                    a   = radius/math.tan(math.radians(alp))
                    d   = math.sqrt(math.pow(a,2)+math.pow(radius,2))
                    ep  = p1.co+a*e1
                    sp  = p1.co+a*e2
                    mid = Blender.Mathutils.MidpointVecs(sp,ep); mid-=p1.co; mid.normalize()
                    cen = p1.co+d*mid
                    # add points to the arch
                    msh.verts.extend(sp)
                    msh.edges.extend(p3,msh.verts[len(msh.verts)-1])
                    arc_point(msh,cen,sp,ep,steps)
                    msh.verts.extend(ep)
                    msh.edges.extend(msh.verts[len(msh.verts)-2],msh.verts[len(msh.verts)-1])
                    msh.edges.extend(msh.verts[len(msh.verts)-1],p2)
                    msh.verts.delete(p1)
                Blender.Window.RedrawAll()
        else: Blender.Draw.PupMenu('ERROR%t|Please, select (only) two edges')
        if edm: Blender.Window.EditMode(1)
    else: Blender.Draw.PupMenu('ERROR%t|Please, select a Mesh object before calling this function')

def break_edge():
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        if len(msh.edges.selected())==1:
            # points (vertices)
            v1 = msh.edges[msh.edges.selected()[0]].v1
            v2 = msh.edges[msh.edges.selected()[0]].v2
            # break-point
            bp = Blender.Window.GetCursorPos()
            if bp!=v1.co and bp!=v2.co:
                msh.verts.extend(bp[0],bp[1],bp[2])
                new = msh.verts[len(msh.verts)-1]
                msh.verts.extend(v2.co)
                msh.edges[msh.edges.selected()[0]].v2 = new
                msh.edges.extend(new,msh.verts[len(msh.verts)-1])
        else: Blender.Draw.PupMenu('ERROR%t|Please, select only one edge')
        if edm: Blender.Window.EditMode(1)
    else: Blender.Draw.PupMenu('ERROR%t|Please, select a Mesh object before calling this function')

def read_2d_mesh(filename):
    Blender.Window.WaitCursor(1)
    edm   = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0) # exit edm
    key   = Blender.sys.basename(Blender.sys.splitext(filename)[0])
    fnode = open(Blender.sys.splitext(filename)[0]+".node",'r')
    felem = open(Blender.sys.splitext(filename)[0]+".ele" ,'r')
    fedge = open(Blender.sys.splitext(filename)[0]+".edge",'r')
    fneig = open(Blender.sys.splitext(filename)[0]+".neigh",'r')
    fnod  = open(Blender.sys.splitext(filename)[0]+".out.node",'w')
    fele  = open(Blender.sys.splitext(filename)[0]+".out.ele" ,'w')
    fedg  = open(Blender.sys.splitext(filename)[0]+".out.face",'w')
    msh   = bpy.data.meshes.new('mesh2d')
    scn   = bpy.data.scenes.active
    obj   = scn.objects.new(msh,key)
    # read nodes
    idx = 0
    for line in fnode.readlines():
        idx  += 1
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif idx==1:
            fnod.write('%s\n' % (words[0]))
            print '[1;34mMechSysCAD[0m: reading ' + words[0] + ' nodes'
            nnod = int(words[0])
        else:
            # read
            x, y, z = float(words[1]), float(words[2]), 0.0
            msh.verts.extend(x,y,z)
            # write
            fnod.write('  %s   %s  %s   %s\n' % (words[0],words[1],words[2],words[3]))
    # read elements
    idx = 0
    elems = []
    for line in felem.readlines():
        idx  += 1
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif idx==1:
            fele.write('%s\n' % (words[0]))
            print '[1;34mMechSysCAD[0m: reading ' + words[0] + ' elements'
            nele = int(words[0])
        else:
            # read
            verts    = []
            allverts = []
            for i in words[1:4]:    verts.append(int(i)-1) # 1,2,3
            for i in words[1:7]: allverts.append(int(i)-1) # 1,2,3,4,5,6
            elems.append(verts)
            msh.faces.extend(verts)
            # write
            fele.write('  %s   6   %s %s %s %s %s %s\n' % (words[0],words[1],words[2],words[3],words[4],words[5],words[6]))
    print elems
    # read neighbours
    idx   = 0
    neigh = []
    for line in fneig.readlines():
        idx  += 1
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif idx==1: print '[1;34mMechSysCAD[0m: reading the neighbours of ' + words[0] + ' triangles'
        else:
            verts = []
            for i in words[1:4]: verts.append(int(i)-1) # 1,2,3
            neigh.append(verts)
    # read edges
    idx = 0
    for line in fedge.readlines():
        idx  += 1
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif idx==1:
            fedg.write('%s\n' % (words[0]))
            print '[1;34mMechSysCAD[0m: reading ' + words[0] + ' edges'
            nedg = int(words[0])
        else:
            # read
            verts = []
            for i in words[1:3]: verts.append(int(i)-1) # 1,2
            # write
            fedg.write('  %s   %s %s %s   %s %s\n' % (words[0],words[1],words[2],words[3],words[4],words[5]))
            # joints
            if words[4]=='-9999':
                # new nodes
                nnod += 1; fnod.write('  %d   %f  %f   %s\n' % (nnod,msh.verts[int(words[1])-1].co[0],msh.verts[int(words[1])-1].co[1],words[4]))
                nnod += 1; fnod.write('  %d   %f  %f   %s\n' % (nnod,msh.verts[int(words[2])-1].co[0],msh.verts[int(words[2])-1].co[1],words[4]))
                nnod += 1; fnod.write('  %d   %f  %f   %s\n' % (nnod,msh.verts[int(words[3])-1].co[0],msh.verts[int(words[3])-1].co[1],words[4]))
                # new elements
                nele += 1; fele.write('  %d   6   %s %d %d %d %s %s\n' % (nele,words[1],nnod-2,nnod,nnod-1,words[2],words[3]))
                # new edges
                nedg += 1; fedg.write('  %d   3   %s %s %s   %s %d\n' % (nedg,nnod-2,nnod-1,nnod,words[4],nele))
                # new connectivity
                own1 = int(words[5])-1
                for i in neigh[own1]:
                    if i>=0:
                        l = int(words[1])-1
                        r = int(words[2])-1
                        for j in elems[i]:
                            if j==l: print 'left'
                            if j==r: print 'right'
                        #if elems[i].count(l)+elems[i].count(r)==2:
                            #fele.write('  %s   6   %s %s %s %s %s %s\n' % (i+1,words[1],words[2],words[3],words[4],words[5],words[6]))
    if edm: Blender.Window.EditMode(1) # enter edm
    fnod.close()
    fele.close()
    fedg.close()
    Blender.Window.RedrawAll()
    Blender.Window.WaitCursor(0)

def gen_struct_mesh():
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    edm = Blender.Window.EditMode()
    bks = []
    for obj in obs:
        if obj!=None and obj.type=='Mesh':
            if edm: Blender.Window.EditMode(0)
            msh = obj.getData(mesh=1)
            if len(msh.verts)==8: # 2D
                c = [[],[]] # coordinates
                for v in msh.verts:
                    c[0].append(v.co[0])
                    c[1].append(v.co[1])
                print c
                #b = ms.mesh_block()
                #b.set ()
                #bks.append(ms.mesh_block())
            else: Blender.Draw.PupMenu('ERROR%t|Each block must have 8 nodes exactly')
    if edm: Blender.Window.EditMode(1)
