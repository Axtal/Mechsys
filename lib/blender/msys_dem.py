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

from   multiprocessing import Process
import Blender
from   Blender import Mesh
from   bpy     import data
from   numpy   import pi
from   Blender.Mathutils import Vector
import msys_dict as di
import msys_mesh as me
import mechsys as ms
import subprocess, math


def gen_pkg(is_ttt, draw=True):
    edm = Blender.Window.EditMode()
    d   = di.load_dict()
    dom = ms.DEM_TTTDomain() if is_ttt else ms.DEM_Domain()
    if d['dem_pkg']==0: # Spheres
        dom.GenSpheres (-1,d['dem_Lx'],d['dem_Nx'],d['dem_rho'],'Normal',d['dem_seed'],d['dem_prob'])
    elif d['dem_pkg']==1: # Spheres HCP
        dom.GenSpheres (-1,d['dem_Lx'],d['dem_Nx'],d['dem_rho'],'HCP',d['dem_seed'],d['dem_prob'])
    elif d['dem_pkg']==2: # Voronoi
        dom.AddVoroPack (-1,d['dem_R'], d['dem_Lx'], d['dem_Ly'], d['dem_Lz'], d['dem_Nx'], d['dem_Ny'], d['dem_Nz'], d['dem_rho'], True, d['dem_seed'], d['dem_prob'])
    elif d['dem_pkg']==3: # Read Mesh
        mesh = me.gen_unstruct_mesh (False)
        dom.GenFromMesh (-1, mesh, d['dem_R'], d['dem_rho'])
    if draw:
        P = []
        dom.GetParticles (P)
        add_particles(P,d['dem_res'],d['dem_draw_verts'],d['dem_draw_edges'])
    if edm: Blender.Window.EditMode(1)
    Blender.Window.QRedrawAll()
    return dom


def gen_GSD(is_ttt):
    obj = di.get_obj()
    dom = gen_pkg(is_ttt,False)
    X, Y, D = [], [], []
    dom.GetGSD (X, Y, D)
    lin  = 'from numpy import sqrt, matrix, zeros, log, cos, pi, sin, tan\n'
    lin += 'from pylab import rc, subplot, plot, xlabel, ylabel, grid, axhline, axvline, axis, text, contour, show\n'
    lin += '\n'
    lin += 'X = '+X.__str__()+'\n'
    lin += 'Y = '+Y.__str__()+'\n'
    lin += 'D = '+D.__str__()+'\n'
    lin += '\n'
    lin += 'plot(X,Y, "ro")\n'
    lin += 'plot(X,Y, "r-", lw=2)\n'
    lin += 'xlabel("diameters")\n'
    lin += 'ylabel("% passing")\n'
    lin += 'grid()\n'
    lin += 'show()\n'
    f = open(obj.name+'_GSD.py', 'w')
    f.write (lin)
    f.close ()
    try: pid = subprocess.Popen(['python', obj.name+'_GSD.py']).pid
    except:
        Blender.Window.WaitCursor(0)
        raise Exception('Python (Matplotlib) command failed')


def load_hdf(filename):
    dom = ms.DEM_Domain()
    dom.Load (filename.replace('.hdf5',''))
    d = di.load_dict()
    P = []
    dom.GetParticles (P)
    add_particles(P,d['dem_res'],d['dem_draw_verts'],d['dem_draw_edges'])


def gen_script():
    Blender.Window.WaitCursor(1)
    d = di.load_dict()

    # C++ script
    if d['dem_cpp_script']:
        # create new text
        txt = Blender.Text.New('dem_ttt.cpp')

        # Header
        txt.write ('// MechSys\n')
        txt.write ('#include <mechsys/dem/domain.h>\n')
        txt.write ('#include <mechsys/util/fatal.h>\n')
        txt.write ('#include <mechsys/linalg/matvec.h>\n')
        txt.write ('\nint main(int argc, char **argv) try\n')
        txt.write ('{\n')

        # constants
        txt.write ('    // constants\n')
        txt.write ('    double Kn   = %g;\n' % d['dem_Kn'])
        txt.write ('    double Kt   = %g;\n' % d['dem_Kt'])
        txt.write ('    double Gn   = %g;\n' % d['dem_Gn'])
        txt.write ('    double Gt   = %g;\n' % d['dem_Gt'])
        txt.write ('    double mu   = %g;\n' % d['dem_mu'])
        txt.write ('    double beta = %g;\n' % d['dem_beta'])
        txt.write ('    double eta  = %g;\n' % d['dem_eta'])

        # input data
        txt.write ('\n    // input data\n')
        txt.write ('    double iso_pf     = %g;\n' % d['dem_iso_pf'])
        txt.write ('    double iso_timef  = %g;\n' % d['dem_iso_timef'])
        txt.write ('    double iso_dt     = %g;\n' % d['dem_iso_dt'])
        txt.write ('    double iso_dtout  = %g;\n' % d['dem_iso_dtout'])
        txt.write ('    bool   iso_render = %s;\n' % ('true' if d['dem_iso_render'] else 'false'))
        txt.write ('    bool   ttt_pcte   = %s;\n' % ('true' if d['dem_ttt_pcte'] else 'false'))
        txt.write ('    double ttt_pf     = %g;\n' % d['dem_ttt_pf'])
        txt.write ('    double ttt_qf     = %g;\n' % d['dem_ttt_qf'])
        txt.write ('    double ttt_thf    = %g;\n' % d['dem_ttt_thf'])
        txt.write ('    double ttt_pex    = %g;\n' % d['dem_ttt_pex'])
        txt.write ('    double ttt_pey    = %g;\n' % d['dem_ttt_pey'])
        txt.write ('    double ttt_pez    = %g;\n' % d['dem_ttt_pez'])
        txt.write ('    double ttt_exf    = %g;\n' % d['dem_ttt_exf'])
        txt.write ('    double ttt_eyf    = %g;\n' % d['dem_ttt_eyf'])
        txt.write ('    double ttt_ezf    = %g;\n' % d['dem_ttt_ezf'])
        txt.write ('    double ttt_timef  = %g;\n' % d['dem_ttt_timef'])
        txt.write ('    double ttt_dt     = %g;\n' % d['dem_ttt_dt'])
        txt.write ('    double ttt_dtout  = %g;\n' % d['dem_ttt_dtout'])
        txt.write ('    bool   ttt_render = %s;\n' % ('true' if d['dem_ttt_render'] else 'false'))

        # domain
        txt.write ('\n    // domain\n')
        txt.write ('    DEM::TriaxialDomain dom;\n')
        txt.write ('    Dict                prp;\n')
        txt.write ('    prp.Set (-1, "Kn Kt Gn Gt Mu Beta Eta", Kn,Kt,Gn,Gt, 0.0, beta,eta);\n')
        txt.write ('    prp.Set (-2, "Kn Kt Gn Gt Mu Beta Eta", Kn,Kt,Gn,Gt, 0.0, beta,eta);\n')
        txt.write ('    prp.Set (-3, "Kn Kt Gn Gt Mu Beta Eta", Kn,Kt,Gn,Gt, 0.0, beta,eta);\n')
        txt.write ('    prp.Set (-4, "Kn Kt Gn Gt Mu Beta Eta", Kn,Kt,Gn,Gt, 0.0, beta,eta);\n')
        txt.write ('    prp.Set (-5, "Kn Kt Gn Gt Mu Beta Eta", Kn,Kt,Gn,Gt, 0.0, beta,eta);\n')
        txt.write ('    prp.Set (-6, "Kn Kt Gn Gt Mu Beta Eta", Kn,Kt,Gn,Gt, 0.0, beta,eta);\n')
        txt.write ('    prp.Set (-7, "Kn Kt Gn Gt Mu Beta Eta", Kn,Kt,Gn,Gt, 0.0, beta,eta);\n')
        txt.write ('    dom.SetProps (prp);\n')

        # generate particles
        txt.write ('\n    // generate particles\n')
        if d['dem_pkg']==0: # Spheres
            txt.write ('    dom.GenSpheres (-1, %g, %d, %g, %s, %d, %g); // tag,Lx,Nx,rho,type,seed,prob\n' % (d['dem_Lx'],d['dem_Nx'],d['dem_rho'],'"Normal"',d['dem_seed'],d['dem_prob']))
        elif d['dem_pkg']==1: # Spheres HCP
            txt.write ('    dom.GenSpheres (-1, %g, %d, %g, %s, %d, %g); // tag,Lx,Nx,rho,type,seed,prob\n' % (d['dem_Lx'],d['dem_Nx'],d['dem_rho'],'"HCP"',   d['dem_seed'],d['dem_prob']))
        elif d['dem_pkg']==2: # Voronoi
            txt.write ('    dom.AddVoroPack (-1, %g, %g,%g,%g, %d,%d,%d, %g, True, %d, %g); // tag, R, Lx,Ly,Lz, Nx,Ny,Nz, rho,seed,prob\n' % (d['dem_R'],d['dem_Lx'],d['dem_Ly'],d['dem_Lz'],d['dem_Nx'],d['dem_Ny'],d['dem_Nz'],d['dem_rho'],d['dem_seed'],d['dem_prob']))
        elif d['dem_pkg']==3: # Read Mesh
            me.gen_unstruct_mesh (True, txt, True, False) # gen_script, txt, cpp, with_headers
            txt.write ('dom.GenFromMesh (-1, mesh, %g, %g);\n' % (d['dem_R'], d['dem_rho']))

        # generate bounding box
        txt.write ('\n    // generate bounding box\n')
        txt.write ('    dom.GenBoundingBox (-2, %g, %g); // ini_tag, R, scale_factor\n' % (d['dem_R'],1.2))

        # stage 1: isotropic compresssion
        txt.write('\n    // stage 1: isotropic compresssion\n')
        txt.write('    Vec3_t  iso_sigf   (-iso_pf, -iso_pf, -iso_pf); // final stress state\n')
        txt.write('    bVec3_t iso_peps   (false, false, false);       // prescribed strain rates ?\n')
        txt.write('    Vec3_t  iso_depsdt (0.0, 0.0, 0.0);             // strain rate\n')
        txt.write('    dom.ResetEps  ();\n')
        txt.write('    dom.SetTxTest (iso_sigf, iso_peps, iso_depsdt);\n')
        txt.write('    dom.Solve     (iso_timef/2.0, iso_dt, iso_dtout, "ttt_isocomp_a", iso_render);\n')
        txt.write('    dom.SetTxTest (iso_sigf, iso_peps, iso_depsdt);\n')
        txt.write('    dom.Solve     (iso_timef,     iso_dt, iso_dtout, "ttt_isocomp_b", iso_render);\n')

        # stage 2: triaxial test
        txt.write('\n    // stage 2: triaxial test\n')
        txt.write('    prp.Set (-1, "Kn Kt Gn Gt Mu Beta Eta", Kn,Kt,Gn,Gt, mu, beta,eta);\n')
        txt.write('    Vec3_t  ttt_sigf;  pqTh2L (ttt_pf, ttt_qf, ttt_thf, ttt_sigf, "cam");\n')
        txt.write('    bVec3_t ttt_peps   (ttt_pex, ttt_pey, ttt_pez);\n')
        txt.write('    Vec3_t  ttt_depsdt (ttt_exf/(ttt_timef-iso_timef), ttt_eyf/(ttt_timef-iso_timef), ttt_ezf/(ttt_timef-iso_timef));\n')
        if d['dem_ttt_pcte']:
            txt.write('    dom.Thf    = ttt_thf*M_PI/180.0;\n')
            txt.write('    dom.Pf     = iso_pf;\n')
            txt.write('    dom.IsPcte = true;\n')
        txt.write('    dom.SetProps  (prp);\n')
        txt.write('    dom.ResetEps  ();\n')
        txt.write('    dom.SetTxTest (ttt_sigf, ttt_peps, ttt_depsdt);\n')
        txt.write('    dom.Solve     (ttt_timef, ttt_dt, ttt_dtout, "ttt_shearing", ttt_render);\n')

        # footer
        txt.write ('\n    return 0;\n')
        txt.write ('}\n')
        txt.write ('MECHSYS_CATCH\n')

    # Python script
    else:
        # create new text
        txt = Blender.Text.New('dem_ttt.py')

        # load MechSys
        txt.write ('# load MechSys\n')
        txt.write ('from numpy   import pi\n')
        txt.write ('from mechsys import *\n')

        # constants
        txt.write ('\n# constants\n')
        txt.write ('Kn   = %g\n' % d['dem_Kn'])
        txt.write ('Kt   = %g\n' % d['dem_Kt'])
        txt.write ('Gn   = %g\n' % d['dem_Gn'])
        txt.write ('Gt   = %g\n' % d['dem_Gt'])
        txt.write ('mu   = %g\n' % d['dem_mu'])
        txt.write ('beta = %g\n' % d['dem_beta'])
        txt.write ('eta  = %g\n' % d['dem_eta'])

        # input data
        txt.write ('\n# input data\n')
        txt.write ('iso_pf     = %g\n' % d['dem_iso_pf'])
        txt.write ('iso_timef  = %g\n' % d['dem_iso_timef'])
        txt.write ('iso_dt     = %g\n' % d['dem_iso_dt'])
        txt.write ('iso_dtout  = %g\n' % d['dem_iso_dtout'])
        txt.write ('iso_render = %s\n' % ('True' if d['dem_iso_render'] else 'False'))
        txt.write ('ttt_pcte   = %s\n' % ('True' if d['dem_ttt_pcte'] else 'False'))
        txt.write ('ttt_pf     = %g\n' % d['dem_ttt_pf'])
        txt.write ('ttt_qf     = %g\n' % d['dem_ttt_qf'])
        txt.write ('ttt_thf    = %g\n' % d['dem_ttt_thf'])
        txt.write ('ttt_pex    = %g\n' % d['dem_ttt_pex'])
        txt.write ('ttt_pey    = %g\n' % d['dem_ttt_pey'])
        txt.write ('ttt_pez    = %g\n' % d['dem_ttt_pez'])
        txt.write ('ttt_exf    = %g\n' % d['dem_ttt_exf'])
        txt.write ('ttt_eyf    = %g\n' % d['dem_ttt_eyf'])
        txt.write ('ttt_ezf    = %g\n' % d['dem_ttt_ezf'])
        txt.write ('ttt_timef  = %g\n' % d['dem_ttt_timef'])
        txt.write ('ttt_dt     = %g\n' % d['dem_ttt_dt'])
        txt.write ('ttt_dtout  = %g\n' % d['dem_ttt_dtout'])
        txt.write ('ttt_render = %s\n' % ('True' if d['dem_ttt_render'] else 'False'))

        # domain
        txt.write ('\n# domain\n')
        txt.write ('dom = DEM_TTTDomain()\n')
        txt.write ('prp = Dict()\n')
        txt.write ('prp.Set (-1, {"Kn":Kn, "Kt":Kt, "Gn":Gn, "Gt":Gt, "Mu":0.0, "Beta":beta, "Eta":eta})\n')
        txt.write ('prp.Set (-2, {"Kn":Kn, "Kt":Kt, "Gn":Gn, "Gt":Gt, "Mu":0.0, "Beta":beta, "Eta":eta})\n')
        txt.write ('prp.Set (-3, {"Kn":Kn, "Kt":Kt, "Gn":Gn, "Gt":Gt, "Mu":0.0, "Beta":beta, "Eta":eta})\n')
        txt.write ('prp.Set (-4, {"Kn":Kn, "Kt":Kt, "Gn":Gn, "Gt":Gt, "Mu":0.0, "Beta":beta, "Eta":eta})\n')
        txt.write ('prp.Set (-5, {"Kn":Kn, "Kt":Kt, "Gn":Gn, "Gt":Gt, "Mu":0.0, "Beta":beta, "Eta":eta})\n')
        txt.write ('prp.Set (-6, {"Kn":Kn, "Kt":Kt, "Gn":Gn, "Gt":Gt, "Mu":0.0, "Beta":beta, "Eta":eta})\n')
        txt.write ('prp.Set (-7, {"Kn":Kn, "Kt":Kt, "Gn":Gn, "Gt":Gt, "Mu":0.0, "Beta":beta, "Eta":eta})\n')
        txt.write ('dom.SetProps (prp)\n')

        # generate particles
        txt.write ('\n# generate particles\n')
        if d['dem_pkg']==0: # Spheres
            txt.write ('dom.GenSpheres (-1, %g, %d, %g, %s, %d, %g) # tag,Lx,Nx,rho,type,seed,prob\n' % (d['dem_Lx'],d['dem_Nx'],d['dem_rho'],'"Normal"',d['dem_seed'],d['dem_prob']))
        elif d['dem_pkg']==1: # Spheres HCP
            txt.write ('dom.GenSpheres (-1, %g, %d, %g, %s, %d, %g) # tag,Lx,Nx,rho,type,seed,prob\n' % (d['dem_Lx'],d['dem_Nx'],d['dem_rho'],'"HCP"',   d['dem_seed'],d['dem_prob']))
        elif d['dem_pkg']==2: # Voronoi
            txt.write ('dom.AddVoroPack (-1, %g, %g,%g,%g, %d,%d,%d, %g, True, %d, %g) # tag, R, Lx,Ly,Lz, Nx,Ny,Nz, rho,seed,prob\n' % (d['dem_R'],d['dem_Lx'],d['dem_Ly'],d['dem_Lz'],d['dem_Nx'],d['dem_Ny'],d['dem_Nz'],d['dem_rho'],d['dem_seed'],d['dem_prob']))
        elif d['dem_pkg']==3: # Read Mesh
            me.gen_unstruct_mesh (True, txt, False, False) # gen_script, txt, cpp, with_headers
            txt.write ('dom.GenFromMesh (-1, mesh, %g, %g)\n' % (d['dem_R'], d['dem_rho']))

        # generate bounding box
        txt.write ('\n# generate bounding box\n')
        txt.write ('dom.GenBoundingBox (-2, %g, %g) # ini_tag, R, scale_factor\n' % (d['dem_R'],1.2))

        # stage 1: isotropic compresssion
        txt.write('\n# stage 1: isotropic compresssion\n')
        txt.write('iso_sigf   = (-iso_pf, -iso_pf, -iso_pf) # final stress state\n')
        txt.write('iso_peps   = (False, False, False)       # prescribed strain rates ?\n')
        txt.write('iso_depsdt = (0.0, 0.0, 0.0)             # strain rate\n')
        txt.write('dom.ResetEps  ()\n')
        txt.write('dom.SetTxTest (iso_sigf, iso_peps, iso_depsdt)\n')
        txt.write('dom.Solve     (iso_timef/2.0, iso_dt, iso_dtout, "ttt_isocomp_a", iso_render)\n')
        txt.write('dom.SetTxTest (iso_sigf, iso_peps, iso_depsdt)\n')
        txt.write('dom.Solve     (iso_timef,     iso_dt, iso_dtout, "ttt_isocomp_b", iso_render)\n')

        # stage 2: triaxial test
        txt.write('\n# stage 2: triaxial test\n')
        txt.write('prp.Set (-1, {"Kn":Kn, "Kt":Kt, "Gn":Gn, "Gt":Gt, "Mu":mu, "Beta":beta, "Eta":eta})\n')
        txt.write('ttt_sigf   = pqTh2L (ttt_pf, ttt_qf, ttt_thf, "cam")\n')
        txt.write('ttt_peps   = (ttt_pex, ttt_pey, ttt_pez)\n')
        txt.write('ttt_depsdt = (ttt_exf/(ttt_timef-iso_timef), ttt_eyf/(ttt_timef-iso_timef), ttt_ezf/(ttt_timef-iso_timef))\n')
        if d['dem_ttt_pcte']:
            txt.write('dom.Thf    = ttt_thf*pi/180.0\n')
            txt.write('dom.Pf     = iso_pf\n')
            txt.write('dom.IsPcte = True\n')
        txt.write('dom.SetProps  (prp)\n')
        txt.write('dom.ResetEps  ()\n')
        txt.write('dom.SetTxTest (ttt_sigf, ttt_peps, ttt_depsdt)\n')
        txt.write('dom.Solve     (ttt_timef, ttt_dt, ttt_dtout, "ttt_shearing", ttt_render)\n')

    Blender.Window.WaitCursor(0)


def add_particles(P,res,draw_verts,draw_edges):
    Blender.Window.WaitCursor(1)
    scn = data.scenes.active
    for p in P:
        # features
        R, V, E, F = p[0], p[1], p[2], p[3]
        objs = []

        # vertices
        if draw_verts:
            for v in V:
                m = Mesh.Primitives.UVsphere(res,res,R*2)
                objs.append (scn.objects.new (m,'vert'))
                objs[-1].setLocation(v[0],v[1],v[2])

        # edges
        if draw_edges:
            for e in E:
                x0  = Vector(V[e[0]])
                x1  = Vector(V[e[1]])
                dx  = x1-x0
                mid = (x0+x1)/2.0
                m   = Mesh.Primitives.Tube(res,R*2,dx.length)
                qua = dx.toTrackQuat ('z','x')
                objs.append (scn.objects.new (m,'edge'))
                objs[-1].setMatrix   (qua.toMatrix())
                objs[-1].setLocation (mid)

        # faces
        for f in F:
            ed0 = Vector(V[f[1]]) - Vector(V[f[0]])
            ed1 = Vector(V[f[2]]) - Vector(V[f[1]])
            n   = ed0.cross(ed1).normalize()

            # list of vertices
            vi = [Vector(V[v])-R*n for v in f]
            vo = [Vector(V[v])+R*n for v in f]

            # list of triangular faces
            fa = []
            nv = len(f) # number of vertices
            nf = nv-2   # number of faces
            for i in range(nf): fa.append ([   0,    i+1,    i+2])
            for i in range(nf): fa.append ([nv+0, nv+i+1, nv+i+2])

            # draw
            m = data.meshes.new ('face')
            m.verts.extend (vi+vo)
            m.faces.extend (fa)
            objs.append (scn.objects.new (m,'face'))

        # join objects
        msh = data.meshes.new ('particle')
        obj = scn.objects.new (msh,'particle')
        obj.join (objs)
        obj.getData(mesh=True).remDoubles (0.001)
        for o in objs: scn.objects.unlink(o)
    Blender.Window.WaitCursor(0)


def run_simulation(running, fatal):
    try:
        d = di.load_dict()

        # domain
        dom = ms.DEM_TTTDomain()
        prp = ms.Dict()
        prp.Set (-1, {"Kn":d['dem_Kn'], "Kt":d['dem_Kt'], "Gn":d['dem_Gn'], "Gt":d['dem_Gt'], "Mu":0.0, "Beta":d['dem_beta'], "Eta":d['dem_eta']})
        prp.Set (-2, {"Kn":d['dem_Kn'], "Kt":d['dem_Kt'], "Gn":d['dem_Gn'], "Gt":d['dem_Gt'], "Mu":0.0, "Beta":d['dem_beta'], "Eta":d['dem_eta']})
        prp.Set (-3, {"Kn":d['dem_Kn'], "Kt":d['dem_Kt'], "Gn":d['dem_Gn'], "Gt":d['dem_Gt'], "Mu":0.0, "Beta":d['dem_beta'], "Eta":d['dem_eta']})
        prp.Set (-4, {"Kn":d['dem_Kn'], "Kt":d['dem_Kt'], "Gn":d['dem_Gn'], "Gt":d['dem_Gt'], "Mu":0.0, "Beta":d['dem_beta'], "Eta":d['dem_eta']})
        prp.Set (-5, {"Kn":d['dem_Kn'], "Kt":d['dem_Kt'], "Gn":d['dem_Gn'], "Gt":d['dem_Gt'], "Mu":0.0, "Beta":d['dem_beta'], "Eta":d['dem_eta']})
        prp.Set (-6, {"Kn":d['dem_Kn'], "Kt":d['dem_Kt'], "Gn":d['dem_Gn'], "Gt":d['dem_Gt'], "Mu":0.0, "Beta":d['dem_beta'], "Eta":d['dem_eta']})
        prp.Set (-7, {"Kn":d['dem_Kn'], "Kt":d['dem_Kt'], "Gn":d['dem_Gn'], "Gt":d['dem_Gt'], "Mu":0.0, "Beta":d['dem_beta'], "Eta":d['dem_eta']})
        dom.SetProps (prp)

        # generate particles
        if d['dem_pkg']==0: # Spheres
            dom.GenSpheres (-1, d['dem_Lx'],d['dem_Nx'],d['dem_rho'], "Normal", d['dem_seed'],d['dem_prob'])
        elif d['dem_pkg']==1: # Spheres HCP
            dom.GenSpheres (-1, d['dem_Lx'],d['dem_Nx'],d['dem_rho'], "HCP", d['dem_seed'],d['dem_prob'])
        elif d['dem_pkg']==2: # Voronoi
            dom.AddVoroPack (-1, d['dem_R'],d['dem_Lx'],d['dem_Ly'],d['dem_Lz'],d['dem_Nx'],d['dem_Ny'],d['dem_Nz'],d['dem_rho'], True, d['dem_seed'],d['dem_prob'])
        elif d['dem_pkg']==3: # Read Mesh
            mesh = me.gen_unstruct_mesh (False)
            dom.GenFromMesh (-1, mesh, d['dem_R'], d['dem_rho'])

        # generate bounding box
        dom.GenBoundingBox (-2, d['dem_R'],1.2)

        # stage 1: isotropic compresssion
        iso_sigf   = (-d['dem_iso_pf'], -d['dem_iso_pf'], -d['dem_iso_pf']) # final stress state
        iso_peps   = (False, False, False)       # prescribed strain rates ?
        iso_depsdt = (0.0, 0.0, 0.0)             # strain rate
        dom.ResetEps  ()
        dom.SetTxTest (iso_sigf, iso_peps, iso_depsdt)
        dom.Solve     (d['dem_iso_timef']/2.0, d['dem_iso_dt'], d['dem_iso_dtout'], "ttt_isocomp_a", d['dem_iso_render'])
        dom.SetTxTest (iso_sigf, iso_peps, iso_depsdt)
        dom.Solve     (d['dem_iso_timef'],     d['dem_iso_dt'], d['dem_iso_dtout'], "ttt_isocomp_b", d['dem_iso_render'])

        # stage 2: triaxial test
        prp.Set (-1, {"Kn":d['dem_Kn'], "Kt":d['dem_Kt'], "Gn":d['dem_Gn'], "Gt":d['dem_Gt'], "Mu":d['dem_mu'], "Beta":d['dem_beta'], "Eta":d['dem_eta']})
        ttt_sigf   = ms.pqTh2L (d['dem_ttt_pf'], d['dem_ttt_qf'], d['dem_ttt_thf'], "cam")
        ttt_peps   = (d['dem_ttt_pex'], d['dem_ttt_pey'], d['dem_ttt_pez'])
        ttt_depsdt = (d['dem_ttt_exf']/(d['dem_ttt_timef']-d['dem_iso_timef']), d['dem_ttt_eyf']/(d['dem_ttt_timef']-d['dem_iso_timef']), d['dem_ttt_ezf']/(d['dem_ttt_timef']-d['dem_iso_timef']))
        if d['dem_ttt_pcte']:
            dom.Thf    = d['dem_ttt_thf']*pi/180.0
            dom.Pf     = d['dem_iso_pf']
            dom.IsPcte = True
        dom.SetProps  (prp)
        dom.ResetEps  ()
        dom.SetTxTest (ttt_sigf, ttt_peps, ttt_depsdt)
        dom.Solve     (d['dem_ttt_timef'], d['dem_ttt_dt'], d['dem_ttt_dtout'], "ttt_shearing", d['dem_ttt_render'])

        # notify parent that we have finished
        running.value = 0

    except Exception, inst:
        print '[1;34mMechSys[0m: Error: '+'[1;31m'+inst.args[0]+'[0m'
        print '>>>>>>>>>>> exception caught by msys_dem.run_simulation <<<<<<<<<<<'
        running.value = 0
        fatal  .value = 1


def run():
    d = di.load_dict()
    if d['dem_running'].value: raise Exception('Another DEM simulation is already running')
    d['dem_running'].value = 1
    d['dem_fatal']  .value = 0
    d['dem_process'] = Process(target=run_simulation, args=(d['dem_running'],d['dem_fatal']))
    print "\n#############################  DEM: running simulation  #################################\n"
    Blender.Window.QRedrawAll()
    d['dem_process'].start()


def stop():
    d = di.load_dict()
    if d['dem_running'].value:
        d['dem_process'].terminate()
        d['dem_running'].value = 0
        print "\n#############################  DEM: simulation stoped  ##################################\n"
        Blender.Window.QRedrawAll()
    else: raise Exception('There is no DEM simulation running')
