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
    // mpi
    MPI::Init (argc, argv);
    int my_id  = MPI::COMM_WORLD.Get_rank();
    int nprocs = MPI::COMM_WORLD.Get_size();

    // input
    int  nx     = 3;
    int  ny     = 3;
    bool full   = false;
    if (argc>1) nx   = atoi(argv[1]);
    if (argc>2) ny   = atoi(argv[2]);
    if (argc>3) full = atoi(argv[3]);

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
    mesh.Generate   (blks,/*O2*/false);
    mesh.PartDomain (nprocs, full);
    mesh.WriteVTU   ("para01_mesh");

    // domain
    Dict prps, mdls, inis;
    prps.Set (-1, "prob geom psa", PROB("Equilib"), GEOM("Quad4"), 1.);
    mdls.Set (-1, "name E nu psa", MODEL("LinElastic"), 1000.0, 0.2, 1.);
    FEM::Domain dom(mesh, prps, mdls, inis);

    //cout << dom << endl;

    // output
    String buf;
    buf.Printf ("para01_%d",my_id);
    dom.WriteVTU (buf.CStr());

    FEM::Solver sol(dom);
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

    // end
    MPI::Finalize();
    return 0;
}
MECHSYS_CATCH
