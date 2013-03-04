/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
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

struct UserData
{
    double p0;
    Vec3_t Sig;
    std::ofstream      oss_ss;       ///< file for stress strain data
};

void Setup1 (DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    double Areax = (dom.GetParticle(-6)->x(2)-dom.GetParticle(-7)->x(2))*(dom.GetParticle(-4)->x(1)-dom.GetParticle(-5)->x(1));
    double Areay = (dom.GetParticle(-6)->x(2)-dom.GetParticle(-7)->x(2))*(dom.GetParticle(-2)->x(0)-dom.GetParticle(-3)->x(0));
    double Areaz = (dom.GetParticle(-4)->x(1)-dom.GetParticle(-5)->x(1))*(dom.GetParticle(-2)->x(0)-dom.GetParticle(-3)->x(0));
    
    dom.GetParticle(-2)->Ff = Vec3_t(-dat.p0*Areax,0.0          ,0.0          );
    dom.GetParticle(-3)->Ff = Vec3_t( dat.p0*Areax,0.0          ,0.0          );
    dom.GetParticle(-4)->Ff = Vec3_t(0.0          ,-dat.p0*Areay,0.0          );
    dom.GetParticle(-5)->Ff = Vec3_t(0.0          , dat.p0*Areay,0.0          );
    dom.GetParticle(-6)->Ff = Vec3_t(0.0          ,0.0          ,-dat.p0*Areaz);
    dom.GetParticle(-7)->Ff = Vec3_t(0.0          ,0.0          , dat.p0*Areaz);
}
void Setup2 (DEM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    Vec3_t trans = Vec3_t(dom.GetParticle(-6)->xb(0) - dom.GetParticle(-6)->x(0),0.0,0.0);
    dom.GetParticle(-6)->Translate(trans);
    trans = Vec3_t(dom.GetParticle(-7)->xb(0) - dom.GetParticle(-7)->x(0),0.0,0.0);
    dom.GetParticle(-7)->Translate(trans);
    double Area = (dom.Xmax-dom.Xmin)*(dom.GetParticle(-4)->x(1)-dom.GetParticle(-5)->x(1));
    dat.Sig = (dom.GetParticle(-6)->F-dom.GetParticle(-7)->F)/Area;
}

void Report (DEM::Domain & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sx" << Util::_8s << "sy" << Util::_8s << "sz" << Util::_8s << "Nc" << Util::_8s << "Nsc \n";
    }
    if (!dom.Finished) 
    {
        size_t Nc = 0;
        size_t Nsc = 0;
        for (size_t i=0; i<dom.CInteractons.Size(); i++)
        {
            Nc += dom.CInteractons[i]->Nc;
            Nsc += dom.CInteractons[i]->Nsc;
        }
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dat.Sig(0) << Util::_8s << dat.Sig(1) << Util::_8s << dat.Sig(2) << Util::_8s << Nc << Util::_8s << Nsc << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////

    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    //number of threads
    size_t Nproc = 1; 
    if (argc==3) Nproc=atoi(argv[2]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());
    
    double verlet;      // Verlet distance for optimization
    String ptype;       // Particle type 
    size_t  RenderVideo; // Decide is video should be render
    bool   Cohesion;    // Decide if coheison is going to be simulated
    double fraction;    // Fraction of particles to be generated
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Beta;        // Rolling stiffness coefficient (only for spheres)
    double Eta;         // Plastic moment coefficient (only for spheres, 0 if rolling resistance is not used)
    double Bn;          // Cohesion normal stiffness
    double Bt;          // Cohesion tangential stiffness
    double Bm;          // Cohesion torque stiffness
    double Eps;         // Threshold for breking bonds
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
    double str;         // Strain rate for shearing
    double T0;          // Time span for the compression
    double Tf;          // Final time for the test
    {
        infile >> verlet;       infile.ignore(200,'\n');
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> Cohesion;     infile.ignore(200,'\n');
        infile >> fraction;     infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu;           infile.ignore(200,'\n');
        infile >> Beta;         infile.ignore(200,'\n');
        infile >> Eta;          infile.ignore(200,'\n');
        infile >> Bn;           infile.ignore(200,'\n');
        infile >> Bt;           infile.ignore(200,'\n');
        infile >> Bm;           infile.ignore(200,'\n');
        infile >> Eps;          infile.ignore(200,'\n');
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
        infile >> str;          infile.ignore(200,'\n');
        infile >> T0;           infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
    }

    // domain and User data
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha=verlet;
    dom.CamPos = Vec3_t(0.0, 1.0*(Lx+Ly+Lz), 0.15*Lz); // position of camera
    dat.p0 = p0;

    // particle
    if      (ptype=="sphere")  dom.GenSpheres  (-1, Lx, nx, rho, "HCP", seed, fraction);
    else if (ptype=="voronoi") dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, Cohesion, !Cohesion, seed, fraction);
    else if (ptype=="tetra")
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        dom.GenFromMesh (mesh,/*R*/R,/*rho*/rho,Cohesion,false);
    }
    else if (ptype=="rice") dom.GenRice(-1,Lx,nx,R,rho,seed,fraction);
    else throw new Fatal("Packing for particle type not implemented yet");
    dom.GenBoundingBox (/*InitialTag*/-2, R, /*Cf*/1.3,Cohesion);
    dom.GetParticle(-2)->FixVeloc();
    dom.GetParticle(-3)->FixVeloc();
    dom.GetParticle(-4)->FixVeloc();
    dom.GetParticle(-5)->FixVeloc();
    dom.GetParticle(-6)->FixVeloc();
    dom.GetParticle(-7)->FixVeloc();

    dom.GetParticle(-2)->vxf = false;
    dom.GetParticle(-3)->vxf = false;
    dom.GetParticle(-4)->vyf = false;
    dom.GetParticle(-5)->vyf = false;
    dom.GetParticle(-6)->vzf = false;
    dom.GetParticle(-7)->vzf = false;

    // properties of particles prior the triaxial test
    Dict B1;
    B1.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt ,Mu     ,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    B1.Set(-2,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B1.Set(-3,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B1.Set(-4,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B1.Set(-5,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B1.Set(-6,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B1.Set(-7,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    dom.SetProps(B1);

    // stage 1: isotropic compresssion  //////////////////////////////////////////////////////////////////////
    String fkey_a  (filekey+"_a");
    String fkey_b  (filekey+"_b");
    String fkey_c  (filekey+"_c");
    String fkeybf_a(filekey+"bf_a");
    String fkeybf_b(filekey+"bf_b");
    String fkeybf_c(filekey+"bf_c");

    dom.Solve (/*tf*/0.5*T0, /*dt*/dt, /*dtOut*/dtOut, &Setup1, NULL, fkey_a.CStr(),RenderVideo,Nproc);
    dom.WriteXDMF(fkey_a.CStr());
    dom.WriteBF(fkeybf_a.CStr());


    dom.GetParticle(-4)->vyf = true;
    dom.GetParticle(-5)->vyf = true;
    dom.GetParticle(-4)->v = OrthoSys::O;
    dom.GetParticle(-5)->v = OrthoSys::O;

    dom.Xmax = dom.GetParticle(-2)->x(0) - dom.GetParticle(-2)->Props.R;
    dom.Xmin = dom.GetParticle(-3)->x(0) + dom.GetParticle(-3)->Props.R;

    dom.ClearInteractons();

    Array<int> idxs(2);
    idxs[0] = -2; idxs[1] = -3;
    dom.DelParticles(idxs);
    
    dom.Solve (/*tf*/   T0, /*dt*/dt, /*dtOut*/dtOut, NULL, NULL, fkey_b.CStr(),RenderVideo,Nproc);
    dom.WriteXDMF(fkey_b.CStr());
    dom.WriteBF(fkeybf_b.CStr());

    Dict B2;
    B2.Set(-1,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,Gt ,Mu     ,Beta,Eta,Bn,Bt ,Bm ,     Eps);
    B2.Set(-4,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B2.Set(-5,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,0.0    ,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B2.Set(-6,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,50.0*Mu,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    B2.Set(-7,"Kn Kt Gn Gt Mu Beta Eta Bn Bt Bm Eps",Kn,Kt,Gn,0.0,50.0*Mu,Beta,Eta,Bn,0.0,0.0,-0.1*Eps);
    dom.SetProps(B2);

    double vel = 0.5*str*(dom.GetParticle(-6)->x(2)-dom.GetParticle(-7)->x(2));
    dom.GetParticle(-6)->v = Vec3_t( vel,0.0,0.0);
    dom.GetParticle(-7)->v = Vec3_t(-vel,0.0,0.0);

    // stage 2: shearing stage         //////////////////////////////////////////////////////////////////////
    dom.Solve (/*tf*/T0+Tf, /*dt*/dt, /*dtOut*/dtOut, &Setup2, &Report, fkey_c.CStr(),RenderVideo,Nproc);
    dom.WriteXDMF(fkey_c.CStr());
    dom.WriteBF(fkeybf_c.CStr());

    return 0;
}
MECHSYS_CATCH
