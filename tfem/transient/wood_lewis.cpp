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

/*  Wood & Lewis (1975) A comparison of time marching schemes for
 *  the transient heat conduction equation. IJNME. 9, 679-689
 *  ============================================================= */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/flowelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/models/linflow.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    double theta = 0.5;
    double dt    = 2.0;
    if (argc>1) theta = atof(argv[1]);
    if (argc>2) dt    = atof(argv[2]);

    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    double L  = 4.0;
    size_t nx = 10;
    double H  = L/nx;
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 0.,  0.0, 0.0,
                 0.,    L, 0.0,
                 0.,    L,   H,
                 0.,  0.0,   H,   -10.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (nx);
    blks[0].SetNy (1);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/true);
    mesh.WriteMPY ("wood_lewis");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom m", PROB("Flow"), GEOM("Quad8"), 1.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name k", MODEL("LinFlow"), 1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "vx vy", 0.0,0.0);

    // domain
    Array<int> out_verts(0,10,24,34);
    FEM::Domain dom(mesh, prps, mdls, inis, "wood_lewis", &out_verts);

    // solver
    FEM::Solver sol(dom);
    sol.Theta = theta;
    
    // stage # 1 -----------------------------------------------------------
    Dict   bcs;
    bcs.Set(-30, "flux", 0.0);
    bcs.Set(-40, "H",    1.0);
    dom.SetBCs (bcs);
    sol.TransSolve (/*tf*/30.0, /*dt*/dt, /*dtOut*/2.0);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    //nod_sol.Set("H", dom.Nods.Size(),

    Table ele_sol;
    //ele_sol.Set("gx gy", dom.Eles.Size(),

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("H", 1.0e-12);
    ele_tol.Set("gx gy", 1.0e-11, 1.0e-11);

    // return error flag
    return dom.CheckError (nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
