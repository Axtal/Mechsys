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
    // to debug: mpirun -np 2 xterm -e gdb --args ./activation01 1

    // input
    bool parallel = false;
    bool nonlin   = false;
    int  nx       = 4;
    int  ny       = 2;
    bool FE       = true;
    bool NR       = false;
    int  nincs    = 1;
    if (argc>1) parallel  = atoi(argv[1]);
    if (argc>2) nonlin    = atoi(argv[2]);
    if (argc>3) nx        = atoi(argv[3]);
    if (argc>4) ny        = atoi(argv[4]);
    if (argc>5) FE        = atoi(argv[5]);
    if (argc>6) NR        = atoi(argv[6]);
    if (argc>7) nincs     = atoi(argv[7]);
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
    fkey.Printf ("activation01_%d", my_id);

    // mesh
    Array<Mesh::Block> blks(2);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -100.,  0.0, 0.0,
                 -200.,  1.0, 0.0,
                 -300.,  1.0, 0.5,
                 -400.,  0.0, 0.5,  -10.0,-20.0,0.0,-40.0);
    blks[1].Set (/*NDim*/2, /*Tag*/-2, /*NVert*/4,
                 -500.,  0.0, 0.5,
                 -600.,  1.0, 0.5,
                 -700.,  1.0, 1.0,
                 -800.,  0.0, 1.0,  0.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (nx);
    blks[0].SetNy (ny);
    blks[1].SetNx (nx);
    blks[1].SetNy (ny);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.WithInfo = root;
    mesh.Generate (blks,/*O2*/false);
    if (parallel) mesh.PartDomain (nprocs);
    buf = fkey + "_mesh";
    mesh.WriteMPY (buf.CStr());
    mesh.WriteVTU (buf.CStr());

    // domain
    FEM::Domain::PARA = parallel;
    Array<int> out_verts(-300,/*justone*/true);
    Dict prps, mdls, inis;
    prps.Set (-1, "prob geom psa active rho", PROB("Equilib"), GEOM("Quad4"), 1., 1., 1.0);
    prps.Set (-2, "prob geom psa active rho", PROB("Equilib"), GEOM("Quad4"), 1., 0., 1.0);
    if (nonlin)
    {
        mdls.Set (-1, "name K0 G0 alp bet psa", MODEL("NLElastic"), 4000.0, 4000.0, 0.4, 0.4, 1.);
        mdls.Set (-2, "name K0 G0 alp bet psa", MODEL("NLElastic"), 4000.0, 4000.0, 0.4, 0.4, 1.);
    }
    else
    {
        mdls.Set (-1, "name E nu psa", MODEL("LinElastic"), 1000.0, 0.2, 1.);
        mdls.Set (-2, "name E nu psa", MODEL("LinElastic"), 1000.0, 0.2, 1.);
    }
    FEM::Domain dom(mesh, prps, mdls, inis, fkey.CStr(), &out_verts);

    // solver
    FEM::Solver sol(dom);
    if (FE) sol.SetScheme("FE");
    if (NR) sol.SetScheme("NR");

    // stage # 1 -------------------------------------
    buf = fkey + "_stage_1";
    Dict bcs;
    bcs.Set      (-10, "uy", 0.0);
    bcs.Set      (-40, "ux", 0.0);
    bcs.Set      (-1, "gravity", 10.0);
    //bcs.Set      (-30, "uy", -0.2);
    dom.SetBCs   (bcs);
    cout << dom << endl;
    sol.Solve    (nincs);
    dom.WriteVTU (buf.CStr());
    dom.PrintResults ("%15.6e", /*onlysummary*/false, /*withelems*/false);

    // stage # 2 -------------------------------------
    buf = fkey + "_stage_2";
    bcs.clear();
    bcs.Set      (-10, "uy", 0.0);
    bcs.Set      (-40, "ux", 0.0);
    bcs.Set      (-2, "activate gravity", 1., 10.0);
    //bcs.Set      (-30, "uy", -0.2);
    dom.SetBCs   (bcs);
    cout << dom << endl;
    sol.Solve    (nincs);
    dom.WriteVTU (buf.CStr());
    dom.PrintResults ("%15.6e", /*onlysummary*/false, /*withelems*/false);

    // end
#ifdef HAS_MPI
    if (parallel) MPI::Finalize();
#endif
    return 0;
}
MECHSYS_CATCH
