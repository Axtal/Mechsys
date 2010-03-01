/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////

    if (argc!=2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    String ptype;       // Particle type 
    bool   RenderVideo; // Decide is video should be render
    double fraction;    // Fraction of particles to be generated
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Beta;        // Rolling stiffness coefficient (only for spheres)
    double Eta;         // Plastic moment coefficient (only for spheres, 0 if rolling resistance is not used)
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt;          // Time step
    double dtOut;       // Time step for output
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double rho;         // rho
    double p0;          // Pressure for the isotropic compression
    double T0;          // Time span for the compression
    bool   isfailure;   // Flag for a failure stress path
    bool   pssrx;       // Prescribed strain rate in X ?
    bool   pssry;       // Prescribed strain rate in Y ?
    bool   pssrz;       // Prescribed strain rate in Z ?
    double srx;         // Final Strain x
    double sry;         // Final Strain y
    double srz;         // Final Strain z
    double pf;          // Final pressure p
    double qf;          // Final deviatoric stress q
    double thf;         // Angle of the stress path alpha
    double alpf;        // Angle of the p q plane
    double Tf;          // Final time for the test
    {
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> fraction;     infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu;           infile.ignore(200,'\n');
        infile >> Beta;         infile.ignore(200,'\n');
        infile >> Eta;          infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> seed;         infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> dtOut;        infile.ignore(200,'\n');
        infile >> Lx;           infile.ignore(200,'\n');
        infile >> Ly;           infile.ignore(200,'\n');
        infile >> Lz;           infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
        infile >> p0;           infile.ignore(200,'\n');
        infile >> T0;           infile.ignore(200,'\n');
        infile >> isfailure;    infile.ignore(200,'\n');
        infile >> pssrx;        infile.ignore(200,'\n');
        infile >> pssry;        infile.ignore(200,'\n');
        infile >> pssrz;        infile.ignore(200,'\n');
        infile >> srx;          infile.ignore(200,'\n');
        infile >> sry;          infile.ignore(200,'\n');
        infile >> srz;          infile.ignore(200,'\n');
        infile >> pf;           infile.ignore(200,'\n');
        infile >> qf;           infile.ignore(200,'\n');
        infile >> thf;          infile.ignore(200,'\n');
        infile >> alpf;         infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
    }

    // domain
    DEM::TriaxialDomain dom;
    dom.CamPos = Vec3_t(0.1*Lx, 0.7*(Lx+Ly+Lz), 0.15*Lz); // position of camera

    // particle
    if      (ptype=="sphere")  dom.GenSpheres  (-1, Lx, nx, rho, "HCP", seed, fraction);
    else if (ptype=="voronoi") dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, true, seed, fraction);
    else if (ptype=="tetra")
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/0.1*Lx*Ly*Lz,Lx,Ly,Lz);
        dom.GenFromMesh (-1,mesh,/*R*/R,/*rho*/rho);
    }
    else throw new Fatal("Packing for particle type not implemented yet");
    dom.GenBoundingBox (/*InitialTag*/-2, R, /*Cf*/1.3);

    // properties of particles prior the triaxial test
    Dict B;
    B.Set(-1,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,Gt,0.0,Beta,Eta);
    B.Set(-2,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,Gt,0.0,Beta,Eta);
    B.Set(-3,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,Gt,0.0,Beta,Eta);
    B.Set(-4,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,Gt,0.0,Beta,Eta);
    B.Set(-5,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,Gt,0.0,Beta,Eta);
    B.Set(-6,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,Gt,0.0,Beta,Eta);
    B.Set(-7,"Kn Kt Gn Gt Mu Beta Eta",Kn,Kt,Gn,Gt,0.0,Beta,Eta);
    dom.SetProps(B);

    // stage 1: isotropic compresssion  //////////////////////////////////////////////////////////////////////
    String fkey_a(filekey+"_a");
    String fkey_b(filekey+"_b");
    Vec3_t  sigf;                      // final stress state
    bVec3_t peps(false, false, false); // prescribed strain rates ?
    Vec3_t  depsdt(0.0,0.0,0.0);       // strain rate
    sigf =  Vec3_t(-p0,-p0,-p0);
    dom.ResetEps  ();
    dom.SetTxTest (sigf, peps, depsdt,false,0,0);
    dom.Solve     (/*tf*/T0/2.0, /*dt*/dt, /*dtOut*/dtOut, fkey_a.CStr(), RenderVideo);
    dom.SetTxTest (sigf, peps, depsdt,false,0,0);
    dom.Solve     (/*tf*/T0, /*dt*/dt, /*dtOut*/dtOut, fkey_b.CStr(), RenderVideo);

    // stage 2: The proper triaxial test /////////////////////////////////////////////////////////////////////////
    String fkey_c(filekey+"_c");
    Vec3_t lf;
    pqTh2L (pf, qf, thf, lf, "cam");
    sigf   = lf(0), lf(1), lf(2);
    peps   = bVec3_t(pssrx, pssry, pssrz);
    depsdt = Vec3_t(srx/(Tf-T0), sry/(Tf-T0), srz/(Tf-T0));

    // properties of particles at the start of  the triaxial test
    B.Set(-1,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu,Beta,Eta);
    dom.SetProps(B);
    
    // run
    dom.ResetEps  ();
    dom.SetTxTest (sigf, peps, depsdt, isfailure, thf*M_PI/180, alpf*M_PI/180);
    dom.ResetInteractons();
    dom.Solve     (/*tf*/Tf, /*dt*/dt, /*dtOut*/dtOut, fkey_c.CStr(), RenderVideo);

    return 0;
}
MECHSYS_CATCH
