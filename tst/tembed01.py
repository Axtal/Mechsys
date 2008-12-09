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

P     = 1.0
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

# Data and Solver
dat = ms.data   (2)
sol = ms.solver (dat, "tembed01_py")

# Elements attributes
if (is_o2): elem_type = "Quad8PStrain"
else:       elem_type = "Quad4PStrain"
eatts = [[-1, elem_type, "LinElastic", "E=1.0 nu=0.0", "ZERO", "", True],
         [-10, "Reinforcement", "", "E=1.0e+8 Ar=0.1 ks=1.0e+12", "ZERO", "", True],
         [-20, "Reinforcement", "", "E=1.0e+8 Ar=0.1 ks=1.0e+12", "ZERO", "", True],
         [-30, "Reinforcement", "", "E=1.0e+8 Ar=0.1 ks=1.0e+12", "ZERO", "", True]]

# Set geometry: nodes and elements
dat.set_nodes_elems (mesh, eatts)

# Add reinforcements
if False:
    ms.add_reinf (0.0, 0.0, 0.0, 1.0, 1.0, 0.0, "E=1.0e+8 Ar=0.1 ks=1.0e+12", True, -10, dat)
    ms.add_reinf (1.0, 1.0, 0.0, 2.0, 0.0, 0.0, "E=1.0e+8 Ar=0.1 ks=1.0e+12", True, -20, dat)
    ms.add_reinf (0.0, 0.0, 0.0, 2.0, 0.0, 0.0, "E=1.0e+8 Ar=0.1 ks=1.0e+12", True, -30, dat)
else:
    dat.add_reinfs ({(0.0,0.0, 1.0,1.0):-10, (1.0,1.0, 2.0,0.0):-20, (0.0,0.0, 2.0,0.0):-30}, eatts)

# Stage # 1
nbrys = [[  0.0, 0.0, 0.0, "ux", 0.0],
         [  0.0, 0.0, 0.0, "uy", 0.0],
         [    W, 0.0, 0.0, "uy", 0.0],
         [W/2.0,   H, 0.0, "fy",  -P]]
dat.set_brys        (mesh, nbrys, [], [])
sol.solve_with_info ()

#//////////////////////////////////////////////////////////////////////////////////////// Check /////

# Test
err = []

for i in range(dat.nelems()):
    tag = dat.ele(i).tag()
    nn  = dat.ele(i).nnodes()
    if   (tag==-10 and nn<=3): err.append(abs(dat.ele(i).val(0, "Sa") - (-0.707106781*P*10)))
    elif (tag==-20 and nn<=3): err.append(abs(dat.ele(i).val(0, "Sa") - (-0.707106781*P*10)))
    elif (tag==-30 and nn<=3): err.append(abs(dat.ele(i).val(0, "Sa") - ( 0.500000000*P*10)))

# Error summary
min_err = min(err)
max_err = max(err)
print
print 'Errors: -----------------------------------------------------------------------'
print "      Min            Max"
print "      %8e   %8e" % (min_err,max_err)
print
