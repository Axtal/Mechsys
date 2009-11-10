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
    int    tag = -1;   // tag of particles
    double R   = 0.1;  // spheroradius
    double Lx  = 4.0;  // length of cube with particles
    double Ly  = 4.0;  // length of cube with particles
    double Lz  = 4.0;  // length of cube with particles
    double nx  = 4;    // number of particles per side
    double ny  = 4;    // number of particles per side
    double nz  = 4;    // number of particles per side
    bool   per = true; // periodic ?
    double rho = 1.0;  // density

    // domain
    Domain d;
    d.CamPos = 0, 35, 0; // position of camera

    // particles
    //d.AddVoroPack (tag, R, Lx,Ly,Lz, nx,ny,nz, rho, per);
    //d.AddRice     (-1,Vec3_t(0.0,0.0,0.0),2.0,0.1,1.0);
    //d.AddSphere   (-1, Vec3_t(0.0,0.0,0.0), /*R*/2.0, rho);
    d.GenSpheres  (-1,4,4,1.0,"HCP");
    d.GenBox      (/*InitialTag*/-2,/*Lx*/6,/*Ly*/6,/*Lz*/6, R, /*Tx*/true, /*Cf*/1.3);
    d.WriteBPY    ("test_triaxial01");

    //return 0;
    // stage 1: isotropic compresssion //////////////////////////////////////////////////////////////////////
    Vec3_t  sigf;                      // final stress state
    bVec3_t peps(false, false, false); // prescribed strain rates ?
    Vec3_t  depsdt(0.0,0.0,0.0);       // strain rate
    sigf =  -0.1,-0.1,-0.1;
    d.SetTxTest (sigf, peps, depsdt);
    d.Solve     (/*tf*/10, /*dt*/0.001, /*dtOut*/0.1, "test_triaxial01a");
    //return 0;

    // stage 2: shearing wiht p-cte /////////////////////////////////////////////////////////////////////////
    double pf  = 0.1;             // final p MPa
    double qf  = 0.15;            // final q
    double thf = 30.0*M_PI/180.0; // final theta
    double tf  = sin(3.0*thf);    // final t = sin(3theta)
    
    // calc principal values (lf)
    Vec3_t lf;
    pqt2L (pf,qf,tf, lf, "cam");
    sigf = lf(0), lf(1), lf(2);
    //cout << sigf << endl;

    // run
    d.ResetEps  ();
    d.SetTxTest (sigf, peps, depsdt);
    d.Solve     (/*tf*/100, /*dt*/0.001, /*dtOut*/0.1, "test_triaxial01b");

    return 0;
}
MECHSYS_CATCH
