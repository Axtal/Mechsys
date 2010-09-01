/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

// Std Lib
#include <iostream>

// Google
#include <google/dense_hash_map>
#include <google/dense_hash_set>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/maps.h>
#include <mechsys/dem/domain.h>

#define MACH_EPS 1.0e-16

using std::cout;
using std::endl;

typedef google::dense_hash_set<int>                 BoxSet_t;
typedef google::dense_hash_set<Particle*>           PartSet_t;
typedef google::dense_hash_map<int,PartSet_t>       Box2Part_t;
typedef google::dense_hash_map<Particle*,Particle*> Neighbours_t;

int main(int argc, char **argv) try
{
    MECHSYS_CATCH_PARALLEL = true;
    MECHSYS_MPI_INIT
    int my_id  = MPI::COMM_WORLD.Get_rank();
    int nprocs = MPI::COMM_WORLD.Get_size();

    // number:  nx ny nz
    Array<int> N(11, 11, 2);
    //Array<int> N(2, 2, 2);
    double dt = 0.001;
    if (argc>1) N[0]  = atoi(argv[1]);
    if (argc>2) N[1]  = atoi(argv[2]);
    if (argc>3) N[2]  = atoi(argv[3]);
    if (argc>4) dt    = atof(argv[4]);
    if (N[0]<2) throw new Fatal("nx must be greater than 1");
    if (N[1]<2) throw new Fatal("ny must be greater than 1");
    if (N[2]<2) throw new Fatal("nz must be greater than 1");

    // limits
    Array<double> L(6);
    //     0   1     2   3     4   5
    //   xmi xma   ymi yma   zmi zma
    //L =    0,  1.01,    0,  1.01,    0,0.101;
    L =  -2,  2,    -2,  2,    0,0.1;

    // set particle with limits
    //Particle::SetLimits (N,L);

    // read data
    Table tab;
    tab.Read ("parts1.dat");
    Array<double> const & X = tab("Xc");
    Array<double> const & Y = tab("Yc");
    Array<double> const & Z = tab("Zc");
    Array<double> const & R = tab("R");

    DEM::Domain dom;

    //double Rmax = 0.0; // largest halo radius
    double Ekin = 0.0;
    for (size_t i=0; i<X.Size(); ++i)
    {
        double vx = static_cast<double>(rand())/static_cast<double>(RAND_MAX)-0.5;
        //double vy = static_cast<double>(rand())/static_cast<double>(RAND_MAX)-0.5;
        //double vx = 0.0;
        double vy = 0.0;
        double vz = 0.0;
        Vec3_t v(vx,vy,vz);
        Vec3_t x(X[i],Y[i],Z[i]);

        dom.AddSphere(-1, x, R[i], 1.0);
        dom.Particles.Last()->v = v;

        dom.Particles[i]->Initialize();
        Ekin += 0.5*dom.Particles.Last()->Props.m*dot(dom.Particles.Last()->v,dom.Particles.Last()->v);
    }

    Dict prps;
    prps.Set (-1,"Kn Kt Gn Gt Mu Beta Eta", 1000.0, 0., 0., 0., 0. ,0.,0.);
    dom.SetProps (prps);
    printf("\nEkin (before) = %16.8e\n", Ekin);


    double Tf      = 1.0;
    double dtout   = 0.1;
    double tout    = 0.1;
    int    stp_out = 0;

    dom.Solve (Tf, dt, dtout);


    // energy
    Ekin = 0.0;
    for (size_t i=0; i<dom.Particles.Size(); ++i)
    {
        Ekin += 0.5*dom.Particles[i]->Props.m*dot(dom.Particles[i]->v,dom.Particles[i]->v);
    }
    printf("Ekin (after ) = %16.8e\n\n", Ekin);

    // end
    MPI::Finalize();
    return 0;
}
MECHSYS_CATCH
