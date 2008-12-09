import mechsys as ms

#           12
#  o|\ +-----------+ /|o
#  o|/ |  embank 2 | \|o
#      |-----------|
#      |  embank 1 |  12
#      |-----------|
#      |           |
#      +-----------+
#     /_\         /_\
#     o o         o o
#

# Constants
E     = 5000.0  # Young
nu    = 0.3     # Poisson
gw    = 10.0    # gamma w
k     = 1.0e-2  # permeability
ndivx = 10      # number of divisions along x and y (for each block)
ndivy =  4      # number of divisions along x and y (for each block)

#///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

bks = []
bks.append(ms.mesh_block())
bks[-1].set_coords(-1,
                   {0: (0.0, 0.0, 0.0), 1: (12.0, 0.0, 0.0), 2: (12.0, 4.0, 0.0), 3: (0.0, 4.0, 0.0)},
                   [(0, 1), (1, 2), (2, 3), (0, 3)],
                   {(1, 2): -10, (0, 1): -11, (0, 3): -10, (2, 3): -12},
                   {},
                   [1.0 for i in range(ndivx)],
                   [1.0 for i in range(ndivy)],
                   [1.0],
                   0, 1, 3, -1)
bks.append(ms.mesh_block())

bks[-1].set_coords(-2,
                   {2: (12.0, 4.0, 0.0), 3: (0.0, 4.0, 0.0), 4: (0.0, 8.0, 0.0), 5: (12.0, 8.0, 0.0)},
                   [(2, 3), (3, 4), (4, 5), (2, 5)],
                   {(4, 5): -13, (2, 3): -12, (2, 5): -10, (3, 4): -10},
                   {},
                   [1.0 for i in range(ndivx)],
                   [1.0 for i in range(ndivy)],
                   [1.0],
                   3, 2, 4, -1)
bks.append(ms.mesh_block())

bks[-1].set_coords(-3,
                   {4: (0.0, 8.0, 0.0), 5: (12.0, 8.0, 0.0), 6: (12.0, 12.0, 0.0), 7: (0.0, 12.0, 0.0)},
                   [(4, 5), (5, 6), (6, 7), (4, 7)],
                   {(5, 6): -10, (4, 5): -13, (4, 7): -10, (6, 7): -14},
                   {},
                   [1.0 for i in range(ndivx)],
                   [1.0 for i in range(ndivy)],
                   [1.0],
                   4, 5, 7, -1)

mesh = ms.mesh_structured (False)
mesh.set_blocks           (bks)
nels = mesh.generate      (True)

#////////////////////////////////////////////////////////////////////////////////////////// FEM /////

# Geometry
dat = ms.geom(2)

# Elements attributes
eatts = [[-1, "Quad4Biot", "LinElastic", "E=%f nu=%f gw=%f k=%f"%(E,nu,gw,k), "ZERO", "gam=20", True ],
         [-2, "Quad4Biot", "LinElastic", "E=%f nu=%f gw=%f k=%f"%(E,nu,gw,k), "ZERO", "gam=20", False],
         [-3, "Quad4Biot", "LinElastic", "E=%f nu=%f gw=%f k=%f"%(E,nu,gw,k), "ZERO", "gam=20", False]]

# Set geometry: nodes, elements, attributes, and boundaries
dat.set_nodes_elems (mesh, eatts)

# Solver
sol = ms.solver("ForwardEuler")
#sol.set_cte("dTini", 1.0)
sol.set_geom(dat)

# Open collection for output
out = ms.output()
out.open_collection("tembank02_py")

# Stage # -1 --------------------------------------------------------------
ebrys = [[-10, "ux",  0.0], [-11, "uy",  0.0], [-12, "pwp", 0.0]]
dat.set_brys          (mesh, [], ebrys, [])
dat.apply_body_forces   ()
sol.solve_with_info   (10, 1e+2, -1, "  Initial stress state due to self weight (zero displacements)\n")
dat.clear_displacements ()
out.vtu               (dat, sol.time())

# Stage # 0 ---------------------------------------------------------------
dat.activate            (-2)
ebrys = [[-10, "ux",  0.0], [-11, "uy",  0.0], [-13, "pwp", 0.0]]
dat.set_brys          (mesh, [], ebrys, [])
sol.solve_with_info   (10, 1e+2, 0, "  Construction of first layer\n")
out.vtu               (dat, sol.time())

# Stage # 1 ---------------------------------------------------------------
dat.activate            (-3)
ebrys = [[-10, "ux",  0.0], [-11, "uy",  0.0], [-14, "pwp", 0.0]]
dat.set_brys          (mesh, [], ebrys, [])
sol.solve_with_info   (10, 1e+2, 0, "  Construction of second layer\n")
out.vtu               (dat, sol.time())

# Close collection
out.close_collection()
