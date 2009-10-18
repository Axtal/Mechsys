from mechsys import *

########################################################################### Mesh ##

mesh = Generic(2)
mesh.SetSize   (6, 4)
mesh.SetVert   (0, -100, 0.0, 0.0, 0)
mesh.SetVert   (1, -100, 0.0, 2.0, 0)
mesh.SetVert   (2,    0, 2.0, 0.0, 0)
mesh.SetVert   (3,    0, 2.0, 1.5, 0)
mesh.SetVert   (4,    0, 4.0, 0.0, 0)
mesh.SetVert   (5,    0, 4.0, 1.0, 0)
mesh.SetCell   (0,   -1,  [0,2,3])
mesh.SetCell   (1,   -1,  [3,1,0])
mesh.SetCell   (2,   -1,  [2,4,5])
mesh.SetCell   (3,   -1,  [5,3,2])
mesh.SetBryTag (1, 0, -10)
mesh.SetBryTag (3, 0, -10)
mesh.WriteVTU  ("ex16_mesh_py")

############################################################################ FEM ##

# elements properties
prps = Dict()
prps.Set(-1, {"prob":PROB("Equilib"), "geom":GEOM("Tri3"), "h":0.25, "pse":True})

# models
mdls = Dict()
mdls.Set(-1, {"name":MODEL("LinElastic"), "E":1.0e+4, "nu":0.2, "pse":True})

# initial values
inis = Dict()
inis.Set(-1, {"sx":0.0, "sy":0.0, "sz":0.0, "sxy":0.0})

# domain
dom = FEM_Domain(mesh, prps, mdls, inis)
dom.SetOutNods ("ex16", [0,1,2,3,5])
dom.SetOutEles ("ex16", [3])

# solver
sol = FEM_Solver(dom)

# stage # 1 -----------------------------------------------------------
bcs = Dict()
bcs.Set( -10, {"qn":-20.0})
bcs.Set(-100, {"ux":0.0, "uy":0.0})
dom.SetBCs (bcs)
sol.Solve  (10)

########################################################################### Output ##

#print mesh
#print prps
#print mdls
#print inis
#print dom
dom.PrintResults ("%11.6g")
dom.WriteMPY     ("ex16_mesh_py")
dom.WriteVTU     ("ex16_py")

########################################################################### Check ##

# correct solution
nod_sol = Table()
nod_sol.Set([                    "ux",                    "uy"],
            [[ 0.000000000000000e+00,   0.000000000000000e+00],
             [ 0.000000000000000e+00,   0.000000000000000e+00],
             [-1.035527877607004e-02,  -2.552969847657423e-02],
             [ 4.727650463081949e-03,  -2.473565538172127e-02],
             [-1.313941349422282e-02,  -5.549310752960183e-02],
             [ 8.389015766816341e-05,  -5.556637423271112e-02]])

ele_sol = Table()
ele_sol.Set(["sx", "sy", "sz", "sxy",  "ex", "ey", "ez", "exy"],
            [[-5.283090599362460e+01, -5.272560566371797e+00, 0.000000000000000e+00, -1.128984616188524e+01, -5.177639388035024e-03,  5.293620632353122e-04,  1.162069331199928e-03, -2.709563078852457e-03/2.0],
             [ 2.462317949521848e+01,  4.924635899043697e+00, 0.000000000000000e+00, -5.153261537858599e+01,  2.363825231540974e-03,  0.000000000000000e+00, -5.909563078852436e-04, -1.236782769086064e-02/2.0],
             [-1.465334062185674e+01, -3.663335155464233e+00, 0.000000000000000e+00, -7.326670310928396e+00, -1.392067359076390e-03, -7.326670310928846e-05,  3.663335155464196e-04, -1.758400874622815e-03/2.0],
             [ 3.102227081237862e+00,  5.914066048600676e+00, 0.000000000000000e+00, -2.178221979271434e+01,  1.919413871517726e-04,  5.293620632353103e-04, -1.803258625967707e-04, -5.227732750251441e-03/2.0]])

# error tolerance
nod_tol = SDPair()
ele_tol = SDPair()
nod_tol.Set({"ux":1.0e-15, "uy":1.0e-15})
ele_tol.Set({"sx":1.0e-12, "sy":1.0e-12, "sz":1.0e-15, "sxy":1.0e-12, "ex":1.0e-15, "ey":1.0e-15, "ez":1.0e-15, "exy":1.0e-15})

# check error
dom.CheckError (nod_sol, ele_sol, nod_tol, ele_tol)
