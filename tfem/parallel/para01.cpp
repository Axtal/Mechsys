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

// STL
#include <iostream>
#include <cmath>

// MechSys
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    // input
    bool parallel  = true;
    int  nx        = 3;
    int  ny        = 3;
    bool part_full = false;
    if (argc>1) parallel  = atoi(argv[1]);
    if (argc>2) nx        = atoi(argv[2]);
    if (argc>3) ny        = atoi(argv[3]);
    if (argc>4) part_full = atoi(argv[4]);

    // mpi
    int my_id  = -1;
    int nprocs = 1;
    if (parallel)
    {
#ifdef USE_MPI
        MPI::Init (argc, argv);
        my_id  = MPI::COMM_WORLD.Get_rank();
        nprocs = MPI::COMM_WORLD.Get_size();
        cout << "\n========================= parallel =========================" << endl;
#else
        throw new Fatal("main.cpp: this code wasn't compiled with USE_MPI ==> parallel version is not available");
#endif
    }
    else cout << "\n========================= serial ===========================" << endl;

    // mesh
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -1.0,  0.0, 0.0,
                 -2.0,  1.0, 0.0,
                 -3.0,  1.0, 1.0,
                 -4.0,  0.0, 1.0,  -10.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (nx);
    blks[0].SetNy (ny);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/false);
    if (parallel) mesh.PartDomain (nprocs, part_full);
    mesh.WriteVTU ("para01_mesh");

    // domain
    Dict prps, mdls, inis;
    prps.Set (-1, "prob geom psa", PROB("Equilib"), GEOM("Quad4"), 1.);
    mdls.Set (-1, "name E nu psa", MODEL("LinElastic"), 1000.0, 0.2, 1.);
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.Parallel = parallel;

    //cout << dom << endl;

    // output
    String buf;
    buf.Printf ("para01_%d",my_id);

    FEM::Solver sol(dom);
    sol.SetScheme("FE");
    //sol.Initialize ();
    //sol.AssembleKA ();
    //sol.A11.WriteSMAT (buf.CStr());

    Dict bcs;
    bcs.Set (-1, "ux uy", 0.0, 0.0);
    bcs.Set (-2, "ux uy", 0.0, 0.0);
    bcs.Set (-3, "fy",    -10.0);
    bcs.Set (-4, "fy",    -10.0);
    dom.SetBCs (bcs);
    sol.Solve  ();

    dom.WriteVTU (buf.CStr());

    // end
#ifdef USE_MPI
    if (parallel) MPI::Finalize();
#endif
    return 0;
}
MECHSYS_CATCH
