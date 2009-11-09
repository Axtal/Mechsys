from numpy import pi, sin
from mechsys import DEM_Domain, pqt2L

# set the simulation domain ######################################
tag = -1   # tag of particles
R   = 0.1  # spheroradius
Lx  = 4.0  # length of cube with particles
Ly  = 4.0  # length of cube with particles
Lz  = 4.0  # length of cube with particles
nx  = 4    # number of particles per side
ny  = 4    # number of particles per side
nz  = 4    # number of particles per side
per = True # periodic ?
rho = 1.0  # density

# domain
d = DEM_Domain()
#d.SetCamPos ((0., 35., 0.)) # position of camera


# particles
Cf = 1.3 # walls length multiplier
#d.AddVoroPack (tag, R, Lx,Ly,Lz, nx,ny,nz, rho, per)
#d.AddRice     (-1, (0.0,0.0,0.0), 2.0, 0.1, 1.0, 0.0, (0.0,0.0,1.0))
#X = (0.0,0.0,0.0)
d.AddSphere   (-1, (0.0,0.0,0.0), 2.0, rho)
#d.GenSpheres  (-1,4,4,1.0)
d.GenBox      (-2, 6, 6, 6, R, True, Cf)
d.WriteBPY    ("test_triaxial01")

# stage 1: isotropic compresssion ###################################
tf     = 10.0
dt     = 0.001
dtOut  = 0.1
sigf   = (-0.1, -0.1, -0.1)    # final stress state
peps   = (False, False, False) # prescribed strain rates ?
depsdt = (0., 0., 0.)          # strain rate
d.SetTxTest (sigf, peps, depsdt)
d.Solve     (tf, dt, dtOut, "test_triaxial01a")

# stage 2: shearing wiht p-cte ####################################/
pf  = 0.1           # final pcam MPa
qf  = 0.15          # final qcam
thf = 30.0*pi/180.0 # final theta
ttf = sin(3.0*thf)  # final t = sin(3theta)

# calc principal values (lf)
sigf = pqt2L (pf,qf,ttf, 'cam')

# run
tf = 100.0
d.ResetEps  ()
d.SetTxTest (sigf, peps, depsdt)
d.Solve     (tf, dt, dtOut, "test_triaxial01b")
