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

/*  Bhatti (2005): Example 8.1 p552 *
 *  =============================== */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "fem/elems/tri3.h"
#include "fem/flowelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linflow.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (4/*verts*/, 2/*cells*/);
    mesh.SetVert   ( 0,    0,  0.0/100,  0.0/100);
    mesh.SetVert   ( 1, -100,  0.0/100,  2.0/100);
    mesh.SetVert   ( 2,    0,  2.0/100,  0.0/100);
    mesh.SetVert   ( 3, -100,  2.0/100,  4.0/100);
    mesh.SetCell   ( 0,   -1, /*NVerts*/3, 0,2,3);
    mesh.SetCell   ( 1,   -1, /*NVerts*/3, 0,3,1);
    mesh.SetBryTag (0, 0, -10);
    mesh.SetBryTag (0, 1, -20);
    mesh.SetBryTag (1, 2, -20);
    mesh.WriteMPY  ("ex81", /*OnlyMesh*/false);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom m", PROB("Flow"), GEOM("Tri3"), 1280000.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name k", MODEL("LinFlow"), 3.0);

    // initial values
    Dict inis;
    inis.Set(-1, "vx vy", 0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.SetOutNods ("ex81",/*NNod*/2,/*IDs*/0,2);

    // solver
    FEM::Solver sol(dom);
    sol.CteTg  = true;
    //sol.Scheme = FEM::Solver::FE_t;
    sol.nSS    = 100;
    
    // stage # 1 -----------------------------------------------------------
    SDPair uvs;
    Dict   bcs;
    uvs.Set("H", 50.0);
    bcs.Set(-100, "H",           300.0)
       .Set( -10, "conv h Tinf", 1.0, 200.0, 50.0)
       .Set( -20, "flux",        0.0);
    dom.SetBCs     (bcs);
    dom.SetUVals   (uvs);
    sol.TransSolve (/*tf*/300.0, /*dt*/1.0, /*dtOut*/10.0);

    // check matrices
    {
        double tol   = 1.0e-11;
        double error = 0.0;
        Mat_t M0c(3,3),M1c(3,3),K0c(3,3),K1c(3,3);
        Mat_t K0,K1,M0,M1;
        dom.Eles[0]->CalcK(K0);  dom.Eles[1]->CalcK(K1);
        dom.Eles[0]->CalcM(M0);  dom.Eles[1]->CalcM(M1);
        M0c =
          8.5333333333333329e+01,  4.2666666666666664e+01,  4.2666666666666664e+01,
          4.2666666666666664e+01,  8.5333333333333329e+01,  4.2666666666666664e+01,
          4.2666666666666664e+01,  4.2666666666666664e+01,  8.5333333333333329e+01;
        M1c =
          4.2666666666666664e+01,  2.1333333333333332e+01,  2.1333333333333332e+01,
          2.1333333333333332e+01,  4.2666666666666664e+01,  2.1333333333333332e+01,
          2.1333333333333332e+01,  2.1333333333333332e+01,  4.2666666666666664e+01;
        K0c =
          4.3333333333333330e+00, -2.3333333333333335e+00,  0.0000000000000000e+00,
         -2.3333333333333335e+00,  5.0833333333333330e+00, -7.5000000000000000e-01,
          0.0000000000000000e+00, -7.5000000000000000e-01,  7.5000000000000000e-01;
        K1c =
          3.0000000000000000e+00,  1.5000000000000000e+00, -4.5000000000000000e+00,
          1.5000000000000000e+00,  1.5000000000000000e+00, -3.0000000000000000e+00,
         -4.5000000000000000e+00, -3.0000000000000000e+00,  7.5000000000000000e+00;
        error += CompareMatrices (M0,M0c);
        error += CompareMatrices (M1,M1c);
        error += CompareMatrices (K0,K0c);
        error += CompareMatrices (K1,K1c);
        cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    }

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults (cout, Util::_12_6);

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
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
