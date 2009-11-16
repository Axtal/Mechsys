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
#include "dem/domain.h"
#include "util/fatal.h"
#include "linalg/matvec.h"

using std::cout;
using std::endl;
using DEM::Domain;

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    
    ifstream infile("triaxialtestparameters.txt");
    double Kn;      //Normal stiffness
    double Kt;      //Tangential stiffness
    double Gn;      //Normal dissipative coefficient
    double Gt;      //Tangential dissipative coefficient
    double Mu;      //Microscopic friction coefficient
    double R;       //Spheroradius
    size_t seed;    //Seed of the ramdon generator
    double dt;      //Time step
    double dtOut;   //Time step for output
    double Lx;      //Lx
    double Ly;      //Ly
    double Lz;      //Lz
    size_t nx;      //nx
    size_t ny;      //ny
    size_t nz;      //nz
    double rho;     //rho
    double p0;      //Pressure for the isotropic compression
    double T0;      //Time span for the compression
    bool   pssrx;   //Prescribed strain rate in X ?
    bool   pssry;   //Prescribed strain rate in Y ?
    bool   pssrz;   //Prescribed strain rate in Z ?
    double srx;     //Final Strain x
    double sry;     //Final Strain y
    double srz;     //Final Strain z
    double pf;      //Final pressure p
    double qf;      //Final deviatoric stress q
    double thf;     //Angle of the stress path alpha
    double Tf;      //Final time for the test
    {
    infile >> Kn;       infile.ignore(200,'\n');
    infile >> Kt;       infile.ignore(200,'\n');
    infile >> Gn;       infile.ignore(200,'\n');
    infile >> Gt;       infile.ignore(200,'\n');
    infile >> Mu;       infile.ignore(200,'\n');
    infile >> R;        infile.ignore(200,'\n');
    infile >> seed;     infile.ignore(200,'\n');
    infile >> dt;       infile.ignore(200,'\n');
    infile >> dtOut;    infile.ignore(200,'\n');
    infile >> Lx;       infile.ignore(200,'\n');
    infile >> Ly;       infile.ignore(200,'\n');
    infile >> Lz;       infile.ignore(200,'\n');
    infile >> nx;       infile.ignore(200,'\n');
    infile >> ny;       infile.ignore(200,'\n');
    infile >> nz;       infile.ignore(200,'\n');
    infile >> rho;      infile.ignore(200,'\n');
    infile >> p0;       infile.ignore(200,'\n');
    infile >> T0;       infile.ignore(200,'\n');
    infile >> pssrx;    infile.ignore(200,'\n');
    infile >> pssry;    infile.ignore(200,'\n');
    infile >> pssrz;    infile.ignore(200,'\n');
    infile >> srx;      infile.ignore(200,'\n');
    infile >> sry;      infile.ignore(200,'\n');
    infile >> srz;      infile.ignore(200,'\n');
    infile >> pf;       infile.ignore(200,'\n');
    infile >> qf;       infile.ignore(200,'\n');
    infile >> thf;      infile.ignore(200,'\n');
    infile >> Tf;       infile.ignore(200,'\n');
    }



    // domain
    Domain d;
    d.CamPos = Vec3_t(0, 3*(Lx+Ly+Lz)/3.0, 0); // position of camera

    // particle
    //d.GenSpheres  (-1,4,10,1.0,"HCP",true);
    d.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho, true, seed);
    d.GenBoundingBox(/*InitialTag*/-2, R, /*Tx*/true, /*Cf*/1.3);
    d.WriteBPY    ("test_triaxial");

    // properties of particles
    Dict B;
    B.Set(-1,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
    B.Set(-2,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
    B.Set(-3,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
    B.Set(-4,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
    B.Set(-5,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
    B.Set(-6,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
    B.Set(-7,"Kn Kt Gn Gt Mu",Kn,Kt,Gn,Gt,Mu);
    d.SetProps(B);

    // stage 1: isotropic compresssion  //////////////////////////////////////////////////////////////////////
    Vec3_t  sigf;                      // final stress state
    bVec3_t peps(false, false, false); // prescribed strain rates ?
    Vec3_t  depsdt(0.0,0.0,0.0);       // strain rate
    sigf =  Vec3_t(-p0,-p0,-p0);
    d.SetTxTest (sigf, peps, depsdt);
    d.Solve     (/*tf*/T0/2.0, /*dt*/dt, /*dtOut*/dtOut, "test_triaxiala",true);
    d.SetTxTest (sigf, peps, depsdt);
    d.Solve     (/*tf*/T0, /*dt*/dt, /*dtOut*/dtOut, "test_triaxialb",true);


    // stage 2: The proper triaxial test /////////////////////////////////////////////////////////////////////////
    double tf  = sin(3.0*thf*M_PI/180.0);    // final t = sin(3theta)
    Vec3_t lf;
    pqt2L (pf, qf, tf, lf, "cam");
    sigf = lf(0), lf(1), lf(2);
    peps = bVec3_t(pssrx, pssry, pssrz);
    depsdt = Vec3_t(srx/(Tf-T0), sry/(Tf-T0), srz/(Tf-T0));

    // run
    d.ResetEps  ();
    d.SetTxTest (sigf, peps, depsdt);
    d.Solve     (/*tf*/Tf, /*dt*/dt, /*dtOut*/dtOut, "test_triaxialc",true);

    return 0;
}
MECHSYS_CATCH
