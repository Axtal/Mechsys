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
#define TAG_VERTSIZE 10001
#define TAG_VERTICES 10002

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

    size_t verts_size = 0;

    if (my_id==0)
    {
        // gen particle
        DEM::Domain dom;
        dom.AddCube (-1, /*x*/Vec3_t(0.5,0.5,0.5), /*R*/0.1, /*L*/1.0, /*rho*/1.0);
        Particle * p = dom.GetParticle (-1);
        verts_size = p->Verts.Size();

        // send data
        for (int k=1; k<nprocs; ++k)
        {
            MPI::COMM_WORLD.Send (p,  /*number of part.*/1, MPI_Particle_Type,  /*destination*/k, TAG_PARTICLE);
            MPI::COMM_WORLD.Send (&verts_size, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/k, TAG_VERTSIZE);
            for (size_t i=0; i<verts_size; ++i)
                MPI::COMM_WORLD.Send (p->Verts[i]->data(), /*number*/3, MPI::DOUBLE, /*destination*/k, TAG_VERTICES);
        }
    }
    else
    {
        // dummy particle
        Particle p;

        // send data
        MPI::COMM_WORLD.Recv (&p, /*number of part.*/1, MPI_Particle_Type,  /*destination*/0, TAG_PARTICLE);
        MPI::COMM_WORLD.Recv (&verts_size, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/0, TAG_VERTSIZE);

        for (size_t i=0; i<verts_size; ++i)
        {
            p.Verts.Push (new Vec3_t(0,0,0));
            MPI::COMM_WORLD.Recv (p.Verts[i]->data(), /*number*/3, MPI::DOUBLE, /*destination*/0, TAG_VERTICES);
        }

        cout << "processor # " << my_id << " has got the following particle: \n" << p;
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
