from mechsys import *

# 2D: structured
def tst1():
    blks = [{'ndim':2, 'tag':-1, 'nx':2, 'ny':3, 'brytags':[-10,-20,-30,-40],
             'V':[[-1,  0.0, 0.0],
                  [-2,  1.0, 0.0],
                  [-3,  1.0, 1.0],
                  [-4,  0.0, 1.0]]},
            {'ndim':2, 'tag':-2, 'nx':4, 'ny':3, 'brytags':[-11,-22,-33,-44],
             'V':[[-5,  1.0, 0.0],
                  [-6,  2.0, 0.0],
                  [-7,  2.0, 1.0],
                  [-8,  1.0, 1.0]]}]
    mesh = Structured(2)       # 2D
    mesh.Generate (blks, True) # O2
    mesh.WriteVTU ("mesh01_quad_py")
    print   " File <mesh01_quad_py.vtu> generated"

# 2D: structured
def tst2():
    mesh = Structured(2) # 2D
    mesh.GenQRing (True, 4, 1, 100., 200., 6)    # O2 Nx Ny r R Nb
    mesh.WriteMPY ("mesh01_quad_ring_py", False) # OnlyMesh
    mesh.WriteVTU ("mesh01_quad_ring_py", 0)     # VolSurfOrBoth
    print   " File <mesh01_quad_ring_py.vtu> generated"

# 3D: structured
def tst3():
    mesh = Structured(3)  # 3D
    mesh.GenBox   (True)  # O2
    mesh.WriteVTU ("mesh01_hex_box_py", 0) # VolSurfOrBoth
    print   " File <mesh01_hex_box_py.vtu> generated"

# 2D: unstructured
def tst4():
    mesh = Unstructured(2)                    # 2D
    mesh.Set ({'pts':[[ -1, 0.0, 0.0],        # tag, x, y <<<<<<<<<<<<<<<<< points
                      [ -2, 1.5, 0.0],        # tag, x, y
                      [ -3, 1.5, 1.5],        # tag, x, y
                      [ -4, 0.0, 1.5],        # tag, x, y
                      [  0, 0.5, 0.5],        # tag, x, y
                      [  0, 1.0, 0.5],        # tag, x, y
                      [  0, 1.0, 1.0],        # tag, x, y
                      [  0, 0.5, 1.0]],       # tag, x, y
               'rgs':[[ -1, -1.0, 0.2, 0.8]], # tag, max{area}, x, y <<<<<< regions
               'hls':[[     0.7, 0.7]],       #      x, y  <<<<<<<<<<<<<<<< holes
               'con':[[-10,  0, 1],           # tag, L, R  <<<<<<<<<<<<<<<< segments
                      [-20,  1, 2],           # tag, L, R
                      [-30,  2, 3],           # tag, L, R
                      [-40,  3, 0],           # tag, L, R
                      [  0,  4, 5],           # tag, L, R
                      [  0,  5, 6],           # tag, L, R
                      [  0,  6, 7],           # tag, L, R
                      [  0,  7, 4]]})         # tag, L, R
    mesh.Generate ()
    mesh.WriteVTU ("mesh01_tri_py")
    print   " File <mesh01_tri_py.vtu> generated"

# 3D: unstructured
def tst5():
    mesh = Unstructured(3)
    mesh.Set ({'pts':[[-1, 0.0, 0.0, 0.0],
                      [-2, 1.0, 0.0, 0.0],
                      [-3, 0.0, 1.0, 0.0],
                      [-4, 0.0, 0.0, 1.0]],
               'rgs':[[-1, -1.0, 0.1, 0.1, 0.1]],
               'hls':[],
               'con':[[-1, [0, 2, 3]],
                      [-2, [0, 3, 1]],
                      [-3, [0, 1, 2]],
                      [-4, [1, 2, 3]]]})
    mesh.Generate (True) # O2
    mesh.WriteVTU ("mesh01_1tet_py")
    print   " File <mesh01_1tet_py.vtu> generated"

# 3D: unstructured
def tst6():
    mesh = Unstructured(3)
    mesh.GenBox   (True, 0.1) # O2 MaxVolume
    mesh.WriteVTU ("mesh01_tet_box_py")
    print "   File <mesh01_tet_box_py.vtu> generated"

# 3D: unstructured
def tst7():
    mesh = Unstructured(3)                           # 3D
    mesh.Set ({'pts':[[-1,  0.0, 0.0, 0.0],          # id, vtag, x, y, z, <<<<<< points
                      [-2,  1.5, 0.0, 0.0],          # id, vtag, x, y, z,
                      [-3,  1.5, 1.5, 0.0],          # id, vtag, x, y, z,
                      [-4,  0.0, 1.5, 0.0],          # id, vtag, x, y, z,
                      [ 0,  0.0, 0.0, 1.5],          # id, vtag, x, y, z, <<<<<< points
                      [ 0,  1.5, 0.0, 1.5],          # id, vtag, x, y, z,
                      [ 0,  1.5, 1.5, 1.5],          # id, vtag, x, y, z,
                      [ 0,  0.0, 1.5, 1.5],          # id, vtag, x, y, z,
                      [ 0,  0.5, 0.5, 0.5],          # id, vtag, x, y, z,
                      [ 0,  1.0, 0.5, 0.5],          # id, vtag, x, y, z,
                      [ 0,  1.0, 1.0, 0.5],          # id, vtag, x, y, z,
                      [ 0,  0.5, 1.0, 0.5],          # id, vtag, x, y, z,
                      [ 0,  0.5, 0.5, 1.0],          # id, vtag, x, y, z,
                      [ 0,  1.0, 0.5, 1.0],          # id, vtag, x, y, z,
                      [ 0,  1.0, 1.0, 1.0],          # id, vtag, x, y, z,
                      [ 0,  0.5, 1.0, 1.0]],         # id, vtag, x, y, z,
               'rgs':[[-1,  -1.0, 0.2, 0.2, 0.2]],   #      tag, max{vol}, x, y, z <<<<<<< regions
               'hls':[[     0.7, 0.7, 0.7]],         #           x, y, z, <<<<<<< holes
               'con':[[-1, [ 0, 3, 7, 4]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [-2, [ 1, 2, 6, 5]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [-3, [ 0, 1, 5, 4]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [-4, [ 2, 3, 7, 6]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [-5, [ 0, 1, 2, 3]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [-6, [ 4, 5, 6, 7]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [ 0, [ 8,11,15,12]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [ 0, [ 9,10,14,13]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [ 0, [ 8, 9,13,12]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [ 0, [10,11,15,14]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [ 0, [ 8, 9,10,11]],           # id, ftag, npolygons,  npoints, point0,point1,point2,point3
                      [ 0, [12,13,14,15]]]})         # id, ftag, npolygons,  npoints, point0,point1,point2,point3
    mesh.Generate ()
    mesh.WritePLY ("mesh01_tet_hole_py")
    mesh.WriteVTU ("mesh01_tet_hole_py")
    print   " File <mesh01_tet_hole_py.vtu> generated"

# main
if __name__=='__main__':
    tst1()
    tst2()
    tst3()
    tst4()
    tst5()
    tst6()
    tst7()
