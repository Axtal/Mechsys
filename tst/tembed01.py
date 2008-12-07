#                       | P                                                   | P = -1.0
#        3            2 V              5                                      V
#     8> @-----------------------------@       -                              .
#        |            .'|'.            |       |                            .' '.
#        |          .'  |  '.          |       |                          .'     '.
#        |        .'    |    '.        |       |                        .'         '.
#        |      .'      |      '.      |       |  H                   .'             '.
#        |    .'        |        '.    |       |                    .'                 '.
#        |  .'         1|          '.  |       |                  .'                     '.
#       0|.'............|............'.|4      |                .'.........................'.
#     8>@------------------------------@       -               /_\                         /_\
#       /_\                           /_\                      ///                          o
#       ///                            o
#
#        <------------ W ------------->

import mechsys as ms

ndivy = 1
ndivx = 2*ndivy
H     = 1.0
W     = 2.0*H
is_o2 = False

#///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

# Block # 1
blks = [ms.mesh_block()]
blks[0].set_coords (-1,                                          # tag to be replicated to all generated elements inside this block
                   {0:(0.0,0.0), 1:(W,0.0), 2:(W,H), 3:(0.0,H)}, # vertices' coordinates
                   [(0,1), (1,2), (2,3), (0,3)],                 # connectivity
                   {(0,1):-30},                                  # edge tags
                   {},                                           # face tags
                   [1 for i in range(ndivx)],                    # weights x
                   [1 for i in range(ndivy)],                    # weights y
                   [],                                           # weights z
                   0, 1, 3, 0)                                   # origin, x-plus, y-plus, z-plus

# Generate
mesh = ms.mesh_structured (False)
if (is_o2): mesh.set_o2()
mesh.set_blocks (blks)
mesh.generate   (True)

#////////////////////////////////////////////////////////////////////////////////////////// FEM /////

# weight
P = 1.0

# Geometry
g = ms.geom(2)

# Nodes brys
nbrys = [[  0.0, 0.0, 0.0, "ux", 0.0],
         [  0.0, 0.0, 0.0, "uy", 0.0],
         [    W, 0.0, 0.0, "uy", 0.0],
         [W/2.0,   H, 0.0, "fy",  -P]]

# Elements attributes
if (is_o2): eatts = [[-1, "Quad8PStrain", "LinElastic", "E=1.0 nu=0.0", "ZERO", "", True]]
else:       eatts = [[-1, "Quad4PStrain", "LinElastic", "E=1.0 nu=0.0", "ZERO", "", True]]

# Set geometry: nodes, elements and attributes
ms.set_nodes_elems (mesh, eatts, g)

# Add reinforcements
if False:
    ms.add_reinf (0.0, 0.0, 0.0, 1.0, 1.0, 0.0, "E=1.0e+8 Ar=0.1 ks=1.0e+12", True, -10, g)
    ms.add_reinf (1.0, 1.0, 0.0, 2.0, 0.0, 0.0, "E=1.0e+8 Ar=0.1 ks=1.0e+12", True, -20, g)
    ms.add_reinf (0.0, 0.0, 0.0, 2.0, 0.0, 0.0, "E=1.0e+8 Ar=0.1 ks=1.0e+12", True, -30, g)
else:
    ratts = [[-10, "Reinforcement", "", "E=1.0e+8 Ar=0.1 ks=1.0e+12", "ZERO", "", True],
             [-20, "Reinforcement", "", "E=1.0e+8 Ar=0.1 ks=1.0e+12", "ZERO", "", True],
             [-30, "Reinforcement", "", "E=1.0e+8 Ar=0.1 ks=1.0e+12", "ZERO", "", True]]
    ms.add_reinfs (False,                                      # Is3D
                   {0:(0.0, 0.0), 1:(1.0, 1.0), 2:(2.0, 0.0)}, # vertices
                   {(0,1):-10, (1,2):-20, (0,2):-30},          # edges/connectivity
                   ratts,                                      # reinforcement attributes
                   g)                                          # geometry

# Set boundary conditions
ms.set_brys (mesh, nbrys, [], [], g)

# Solve
sol = ms.solver     ("ForwardEuler")
sol.set_geom        (g)
sol.solve_with_info ()

# output
out = ms.output ()
out.vtu         (g, "tembed01_py.vtu")

#//////////////////////////////////////////////////////////////////////////////////////// Check /////

# Test
err = []

for i in range(g.nelems()):
    tag = g.ele(i).tag()
    nn  = g.ele(i).nnodes()
    if   (tag==-10 and nn<=3): err.append(abs(g.ele(i).val(0, "Sa") - (-0.707106781*P*10)))
    elif (tag==-20 and nn<=3): err.append(abs(g.ele(i).val(0, "Sa") - (-0.707106781*P*10)))
    elif (tag==-30 and nn<=3): err.append(abs(g.ele(i).val(0, "Sa") - ( 0.500000000*P*10)))

# Error summary
min_err = min(err)
max_err = max(err)
print
print 'Errors: -----------------------------------------------------------------------'
print "      Min            Max"
print "      %8e   %8e" % (min_err,max_err)
print
