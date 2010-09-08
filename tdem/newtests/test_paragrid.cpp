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

// MechSys
#include <mechsys/util/fatal.h>
#include <mechsys/util/maps.h>
#include <mechsys/mesh/paragrid3d.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // initialize MPI
    MECHSYS_CATCH_PARALLEL = true;
    MECHSYS_MPI_INIT
    int my_id  = MPI::COMM_WORLD.Get_rank();

    Array<int> N(20, 30, 1); // number of cells/boxes along each side of grid
    if (argc>1) N[0]  = atoi(argv[1]);
    if (argc>2) N[1]  = atoi(argv[2]);
    if (argc>3) N[2]  = atoi(argv[3]);
    if (N[0]<1) throw new Fatal("nx must be greater than 0");
    if (N[1]<1) throw new Fatal("ny must be greater than 0");
    if (N[2]<1) throw new Fatal("nz must be greater than 0");

    // limits of grid
    Array<double> L(6);
    //     0    1      2    3      4    5
    //   xmi  xma    ymi  yma    zmi  zma
    L =  -2.,  2.,   -2.,  2.,    0., 1.0;

    // grid
    String buf;
    buf.Printf("test_paragrid_proc_%d",my_id);
    ParaGrid3D grid(N, L, buf.CStr());
    cout << "File <" << buf << ".vtu> written" << endl;

    // neighbour partitions
    cout << "Proc # " << my_id << ": neighbours = ";
    for (size_t i=0; i<grid.NeighProcs.Size(); ++i) cout << grid.NeighProcs[i] << " ";
    cout << endl;

    // end
    MPI::Finalize();
    return 0;
}
MECHSYS_CATCH
