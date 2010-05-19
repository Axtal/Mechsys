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
#include <iostream>
#include <math.h>

// GSL
#include <gsl/gsl_linalg.h>

// MechSys
#include <mechsys/util/fatal.h>
#include <mechsys/dem/particle.h>
#include <mechsys/dem/domain.h>

using std::cout;
using std::endl;

#define TAG_PARTICLE 10000

int main(int argc, char **argv) try
{
#ifdef USE_MPI
    // init
    MPI::Init (argc, argv);
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors

    // create particle data type
    BuildParticleDataType (MPI_Particle_Type);

    // check
    cout << "hi, I'm processor # " << my_id << endl;

    // destinations
    if (nprocs<2) throw new Fatal("The number of processors needs to be greater than 1");
    Array<int> destinations;
    for (int i=0; i<nprocs; ++i)
    {
        if (my_id!=i) destinations.Push(i);
    }

    // gen particle
    DEM::Domain dom;
    dom.AddRice (-(my_id+1), /*x*/Vec3_t(0.5,0.5,0.5), /*R*/0.1, /*L*/1.0, /*rho*/1.0);
    Particle * p_snd = dom.GetParticle (-(my_id+1));

    // dummy particle
    Particle p_rcv;

    for (int k=0; k<nprocs-1; ++k)
    {
        MPI::Request snd_part = MPI::COMM_WORLD.Isend (p_snd,  /*number*/1, MPI_Particle_Type, destinations[k], TAG_PARTICLE);
        MPI::Request rcv_part = MPI::COMM_WORLD.Irecv (&p_rcv, /*number*/1, MPI_Particle_Type, MPI::ANY_SOURCE, TAG_PARTICLE);

        // Wait for all requests to finish
        MPI::Status st;
        snd_part.Wait();
        rcv_part.Wait(st);

        cout << "processor # " << my_id << " sent particle "     << p_snd->Tag << " to proc # "   << destinations[k] << endl;
        cout << "processor # " << my_id << " received particle " << p_rcv.Tag  << endl;
    }

    cout << "processor # " << my_id << " has finished" << endl;

    // end
    MPI::Finalize();
    return 0;
#else
    throw new Fatal("This program needs MPI");
#endif
}
MECHSYS_CATCH
