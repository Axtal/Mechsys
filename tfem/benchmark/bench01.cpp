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
    bool parallel  = false;
    bool mixed     = false;
    bool usigeps   = true;
    bool nonlin    = false;
    int  nxyz      = 3;
    bool FE        = true;
    bool NR        = false;
    int  nincs     = 1;
    bool part_full = false; // use full neighbours during PartDom ?
    if (argc>1) parallel  = atoi(argv[1]);
    if (argc>2) mixed     = atoi(argv[2]);
    if (argc>3) usigeps   = atoi(argv[3]);
    if (argc>4) nonlin    = atoi(argv[4]);
    if (argc>5) nxyz      = atoi(argv[5]);
    if (argc>6) FE        = atoi(argv[6]);
    if (argc>7) NR        = atoi(argv[7]);
    if (argc>8) nincs     = atoi(argv[8]);
    if (argc>9) part_full = atoi(argv[9]);

    // mpi
    int my_id  = -1;
    int nprocs = 1;
    if (parallel)
    {
#ifdef HAS_MPI
        MPI::Init (argc, argv);
        my_id  = MPI::COMM_WORLD.Get_rank();
        nprocs = MPI::COMM_WORLD.Get_size();
        printf("\n%s===================================== Parallel =====================================%s\n",TERM_YELLOW_BLUE,TERM_RST);
#else
        throw new Fatal("main.cpp: this code wasn't compiled with HAS_MPI ==> parallel version is not available");
#endif
    }
    else printf("\n%s====================================== Serial ======================================%s\n",TERM_BLACK_WHITE,TERM_RST);

    // fkey
    String fkey, buf;
    fkey.Printf ("bench01_%d", my_id);

#ifdef USE_MTL4
    printf("\n%s--------------------------------------- MTL4 ---------------------------------------%s\n",TERM_BLACK_WHITE,TERM_RST);
    fkey.append("MTL4");
#else
    printf("\n%s---------------------------------- Raul's LaExpr -----------------------------------%s\n",TERM_BLACK_WHITE,TERM_RST);
#endif

    // mesh
    /*                4----------------7
                    ,'|              ,'|
                  ,'  |            ,'  |
                ,'    | -60   -10,'    |
              ,'      |        ,'      |
            5'===============6'        |
            |         |      |    -40  |
            |    -30  |      |         |
            |         0- - - | -  - - -3
            |       ,'       |       ,'
            |     ,' -20     |     ,'
            |   ,'        -50|   ,'
            | ,'             | ,'
            1----------------2'                   */
    Mesh::Structured mesh(/*NDim*/3);
    mesh.GenBox (/*O2*/true, nxyz,nxyz,nxyz, 1.0,1.0,1.0);
    if (parallel) mesh.PartDomain (nprocs, part_full);
    buf = fkey + "_mesh";
    mesh.WriteVTU (buf.CStr());

    // domain
    FEM::Domain::PARA = parallel;
    Dict prps, mdls, inis;
    double prob = (mixed ? (usigeps ? PROB("USigEps") : PROB("USig")) : PROB("Equilib"));
    prps.Set (-1, "prob geom d3d", prob, GEOM("Hex20"), 1.);
    if (nonlin) mdls.Set (-1, "name K0 G0 alp bet d3d", MODEL("NLElastic"), 4000.0, 4000.0, 0.4, 0.4, 1.);
    else        mdls.Set (-1, "name E nu d3d", MODEL("LinElastic"), 1000.0, 0.2, 1.);
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    FEM::Solver sol(dom);
    if (FE) sol.SetScheme("FE");
    if (NR) sol.SetScheme("NR");

    // boundary conditions for stage # 1
    Dict bcs;
    bcs.Set (-10, "ux",  0.0);
    bcs.Set (-30, "uy",  0.0);
    bcs.Set (-50, "uz",  0.0);
    bcs.Set (-60, "uz", -0.2);
    dom.SetBCs (bcs);

    // output domain
    /*
    buf = fkey + "_dom.txt";
    std::ofstream of1(buf.CStr(), std::ios::out);
    of1 << dom << endl;
    of1.close();
    */

    // solve
    sol.Solve    (nincs);
    dom.WriteVTU (fkey.CStr());

    // end
#ifdef HAS_MPI
    if (parallel) MPI::Finalize();
#endif
    return 0;
}
MECHSYS_MPI_CATCH
