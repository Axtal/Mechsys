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
#include <cstdlib> // for std::rand

// MechSys
#define PARALLEL_DEBUG
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    // time mpirun -np 4 ./bench01 1 0 1 0 10
    // time ./bench01 0 0 1 0 10
    // time ./bench01 0 0 1 0 17/18            => UMFPACK out of memory
    // time mpirun -np 4 ./bench01 1 0 1 0 19  => Error with Array (own implementation gives wrong residual. Sometimes UMFPACK fails)
    // time mpirun -np 4 ./bench01 1 0 1 0 20  => Error with Array (own implementation fails eventually, std::vector gives wrong residual)
    // time mpirun -np 4 ./bench01 1 0 1 0 15  => Segfault (with std::vector)
    // time mpirun -np 4 ./bench01 1 0 1 0 13  => wrong results

    // mpirun -np 7 ./bench01 1 1 123  0 0 0  10 => segfault
    // mpirun -np 4 ./bench01 1 1 123  0 0 0  6 => wrong results
    // mpirun -np 3 ./bench01 1 1 123  0 0 0  6 => wrong results

    // input
    bool parallel  = false;
    bool part_rnd  = true;
    bool seed      = 123;
    bool mixed     = false;
    bool usigeps   = true;
    bool nonlin    = false;
    int  nxyz      = 3;
    bool FE        = true;
    bool NR        = false;
    int  nincs     = 1;
    bool part_full = false; // use full neighbours during PartDom ?
    if (argc> 1) parallel  = atoi(argv[ 1]);
    if (argc> 2) part_rnd  = atoi(argv[ 2]);
    if (argc> 3) seed      = atoi(argv[ 3]);
    if (argc> 4) mixed     = atoi(argv[ 4]);
    if (argc> 5) usigeps   = atoi(argv[ 5]);
    if (argc> 6) nonlin    = atoi(argv[ 6]);
    if (argc> 7) nxyz      = atoi(argv[ 7]);
    if (argc> 8) FE        = atoi(argv[ 8]);
    if (argc> 9) NR        = atoi(argv[ 9]);
    if (argc>10) nincs     = atoi(argv[10]);
    if (argc>11) part_full = atoi(argv[11]);
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
    fkey.Printf ("bench01_%d", my_id);

#ifdef USE_MTL4
    if (root) printf("\n%s--------------------------------------- MTL4 ---------------------------------------%s\n",TERM_BLACK_WHITE,TERM_RST);
    fkey.append("_MTL4");
#else
    if (root) printf("\n%s---------------------------------- Raul's LaExpr -----------------------------------%s\n",TERM_BLACK_WHITE,TERM_RST);
#endif

    // mesh
    /*                4----------------7
                    ,'|              ,'|
                  ,'  |            ,'  |
                ,'    | -60   -10,'    |
              ,'      |        ,'      |
            5'===============6(-7)     |
            |         |      |    -40  |
            |    -30  |      |         |
            |         0- - - | -  - - -3
            |       ,'       |       ,'
            |     ,' -20     |     ,'
            |   ,'        -50|   ,'
            | ,'             | ,'
            1----------------2'                   */
    Mesh::Structured mesh(/*NDim*/3);
    mesh.WithInfo = root;
    if (parallel) mesh.OnlyRoot = true;
    mesh.GenBox (/*O2*/true, nxyz,nxyz,nxyz, 1.0,1.0,1.0);
    /*
    Array<int> part(27);
    part = 0,1,2, 3,2,1, 0,1,2,
           1,2,3, 2,1,0, 1,2,3,
           2,3,2, 1,0,1, 2,3,0;
           */
    //if (parallel) mesh.PartDomain (nprocs, part_full, part.GetPtr());
    if (parallel)
    {
        if (part_rnd)
        {
            std::srand(seed);
            Array<int> part(mesh.Cells.Size());
            for (size_t i=0; i<part.Size(); ++i)
            {
                part[i] = (std::rand() % nprocs);
                if (part[i]<0 || part[i]>nprocs-1) throw new Fatal("part[i]=%d is wrong",part[i]);
            }
            mesh.PartDomain (nprocs, part_full, part.GetPtr());
        }
        else mesh.PartDomain (nprocs, part_full);
    }
    buf = fkey + "_mesh";
    mesh.WriteVTU (buf.CStr());

    // domain
    FEM::Domain::PARA = parallel;
    Array<int> out_verts(-7,true);
    Dict prps, mdls, inis;
    double prob = (mixed ? (usigeps ? PROB("USigEps") : PROB("USig")) : PROB("Equilib"));
    prps.Set (-1, "prob geom d3d", prob, GEOM("Hex20"), 1.);
    if (nonlin) mdls.Set (-1, "name K0 G0 alp bet d3d", MODEL("NLElastic"), 4000.0, 4000.0, 0.4, 0.4, 1.);
    else        mdls.Set (-1, "name E nu d3d", MODEL("LinElastic"), 1000.0, 0.2, 1.);
    FEM::Domain dom(mesh, prps, mdls, inis, fkey.CStr(), &out_verts);

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
    printf("Norm(dU) = %.8f\n",Norm(sol.dU));

    // end
#ifdef HAS_MPI
    if (parallel) MPI::Finalize();
#endif
    return 0;
}
MECHSYS_CATCH