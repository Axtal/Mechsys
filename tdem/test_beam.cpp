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

// Std lib
#include <math.h>

// GSL
#include <gsl/gsl_linalg.h>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/structured.h>

using std::cout;
using std::endl;

struct UserData
{
    Particle * p; // the array of particles at which the force is to be applied
};

void Setup (DEM::Domain & Dom, void * UD)
{
    // force at -3
    UserData & dat = (*static_cast<UserData *>(UD));
    dat.p->Ff=0.0,0.0,10.0*sin(0.3*Dom.Time);
}


int main(int argc, char **argv) try
{

    if (argc!=2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    double verlet;      // Verlet distance for optimization
    String ptype;       // Particle type 
    bool   RenderVideo; // Decide is video should be render
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Bn;          // Cohesion normal stiffness
    double Bt;          // Cohesion tangential stiffness
    double Bm;          // Cohesion torque stiffness
    double Eps;         // Threshold for breking bonds
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt;          // Time step
    double dtOut;       // Time step for output
    double Tf;          // Final time for the test
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double rho;         // rho
    {
        infile >> verlet;       infile.ignore(200,'\n');
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu;           infile.ignore(200,'\n');
        infile >> Bn;           infile.ignore(200,'\n');
        infile >> Bt;           infile.ignore(200,'\n');
        infile >> Bm;           infile.ignore(200,'\n');
        infile >> Eps;          infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> seed;         infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> dtOut;        infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
        infile >> Lx;           infile.ignore(200,'\n');
        infile >> Ly;           infile.ignore(200,'\n');
        infile >> Lz;           infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
    }



    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha = verlet;

    if (ptype=="voronoi") dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, true, false, seed, 1.0);
    else if (ptype=="cube")
    {
        Mesh::Structured mesh(3);
        mesh.GenBox (false, nx, ny, nz, Lx, Ly, Lz);
        dom.GenFromMesh (mesh, R, rho, true, false);
    }
    else throw new Fatal("Packing for particle type not implemented yet");

    //set the element properties
    dom.Center();
    dom.GenBoundingPlane(-2,R,1.0,true);
    Dict B;
    B.Set(-1,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn,Kt,Bn,Bt,Bm,Gn,Gt,Mu,Eps);
    B.Set(-2,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn,Kt,Bn,Bt,Bm,Gn,Gt,Mu,Eps);
    B.Set(-3,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn,Kt,Bn,Bt,Bm,Gn,Gt,Mu,Eps);
    dom.SetProps(B);
    // input
    double cam_x=2*Ly, cam_y=Ly/2.0, cam_z=0;
    dom.CamPos = cam_x, cam_y, cam_z;

     //connect
    dat.p = dom.GetParticle (-2);

    cout << dom.Particles[dom.Particles.Size()-1] << endl;
    // fix -2 particles at the left extreme of the beam
    Particle * p;
    p = dom.GetParticle (-3);
    p->FixVeloc();

    dom.Solve (Tf,dt,dtOut, &Setup, NULL, filekey.CStr(), RenderVideo);
}
MECHSYS_CATCH
