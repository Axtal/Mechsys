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
    bool nonlin    = false;
    int  nx        = 3;
    int  ny        = 3;
    bool FE        = true;
    bool NR        = false;
    int  nincs     = 1;
    bool part_full = false;
    if (argc>1) parallel  = atoi(argv[1]);
    if (argc>2) nonlin    = atoi(argv[2]);
    if (argc>3) nx        = atoi(argv[3]);
    if (argc>4) ny        = atoi(argv[4]);
    if (argc>5) FE        = atoi(argv[5]);
    if (argc>6) NR        = atoi(argv[6]);
    if (argc>7) nincs     = atoi(argv[7]);
    if (argc>8) part_full = atoi(argv[8]);
    MECHSYS_CATCH_PARALLEL = parallel;

    // mpi
    int  my_id  = -1;
    int  nprocs = 1;
    bool root   = true;
    if (parallel)
    {
#ifdef HAS_MPI
        MECHSYS_MPI_INIT
        my_id  = MPI::COMM_WORLD.Get_rank();
        nprocs = MPI::COMM_WORLD.Get_size();
        if (my_id!=0) root = false;
        if (root) printf("\n%s===================================== Parallel =====================================%s\n",TERM_YELLOW_BLUE,TERM_RST);
#else
        throw new Fatal("main.cpp: this code wasn't compiled with HAS_MPI ==> parallel version is not available");
#endif
    }
    else printf("\n%s====================================== Serial ======================================%s\n",TERM_BLACK_WHITE,TERM_RST);

    // fkey
    String fkey, buf;
    fkey.Printf ("para01_%d", my_id);

#ifdef USE_MTL4
    if (root) printf("\n%s--------------------------------------- MTL4 ---------------------------------------%s\n",TERM_BLACK_WHITE,TERM_RST);
    fkey.append("_MTL4");
#else
    if (root) printf("\n%s---------------------------------- Raul's LaExpr -----------------------------------%s\n",TERM_BLACK_WHITE,TERM_RST);
#endif

    // mesh
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -100.,  0.0, 0.0,
                 -200.,  1.0, 0.0,
                 -300.,  1.0, 1.0,
                 -400.,  0.0, 1.0,  -10.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (nx);
    blks[0].SetNy (ny);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.WithInfo = root;
    if (parallel) mesh.OnlyRoot = true;
    mesh.Generate (blks,/*O2*/false);
    if (parallel) mesh.PartDomain (nprocs, part_full);
    buf = fkey + "_mesh";
    mesh.WriteVTU (buf.CStr());

    // domain
    FEM::Domain::PARA = parallel;
    Array<int> out_verts(-300,true);
    Dict prps, mdls, inis;
    prps.Set (-1, "prob geom psa", PROB("Equilib"), GEOM("Quad4"), 1.);
    if (nonlin) mdls.Set (-1, "name K0 G0 alp bet psa", MODEL("NLElastic"), 4000.0, 4000.0, 0.4, 0.4, 1.);
    else        mdls.Set (-1, "name E nu psa", MODEL("LinElastic"), 1000.0, 0.2, 1.);
    FEM::Domain dom(mesh, prps, mdls, inis, fkey.CStr(), &out_verts);

    /*
    buf = fkey + "_dom_before.txt";
    std::ofstream of1(buf.CStr(), std::ios::out);
    of1 << dom << endl;
    of1.close();
    */

    // solver
    FEM::Solver sol(dom);
    if (FE) sol.SetScheme("FE");
    if (NR) sol.SetScheme("NR");
    //sol.CorR = false;
    //sol.Initialize ();
    //sol.AssembleKA ();
    //sol.A11.WriteSMAT (fkey.CStr());

    // boundary conditions for stage # 1
    Dict bcs;
    /*
    bcs.Set (-100, "ux uy", 0.0, 0.0);
    //bcs.Set (-200, "uy",    0.0);
    bcs.Set (-10, "uy",    0.0);
    bcs.Set (-300, "fy",    -10.0);
    bcs.Set (-400, "fy",    -20.0);
    */
    bcs.Set (-300, "fy",    -30.0);
    bcs.Set (-400, "ux uy", 0.0, 0.0);
    bcs.Set (-10,  "uy", -0.1);
    dom.SetBCs (bcs);

    // output domain
    /*
    buf = fkey + "_dom_after.txt";
    of1.open(buf.CStr(), std::ios::out);
    of1 << dom << endl;
    of1.close();
    */

    // solve
    sol.Solve        (nincs);
    dom.WriteVTU     (fkey.CStr());
    //dom.PrintResults ("%12.6g", /*with_elems*/false);

    // end
#ifdef HAS_MPI
    if (parallel) MPI::Finalize();
#endif
    return 0;
}
MECHSYS_CATCH
