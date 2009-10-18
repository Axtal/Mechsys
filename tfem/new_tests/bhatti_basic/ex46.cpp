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

/*  Bhatti (2005): Example 4.6, p253  *
 *  ================================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "fem/beam.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    double L = 300.0;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (4/*nodes*/, 3/*cells*/);
    mesh.SetVert  (0, -100,  0.0, 0.0);
    mesh.SetVert  (1, -200,    L, 0.0);
    mesh.SetVert  (2, -300, 2.*L, 0.0);
    mesh.SetVert  (3, -400, 3.*L, 0.0);
    mesh.SetCell  (0,   -1, Array<int>(0, 1));
    mesh.SetCell  (1,   -1, Array<int>(1, 2));
    mesh.SetCell  (2,   -2, Array<int>(2, 3));
    mesh.AddPin   (-300);
    mesh.WriteMPY ("ex46",/*OnlyMesh*/false);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob E A Izz fra", PROB("Beam"), 1.0, 3.0, 180000.0, TRUE);
    prps.Set(-2, "prob E A Izz fra", PROB("Beam"), 2.0, 3.0, 180000.0, TRUE);

    // domain
    FEM::Domain dom(mesh, prps, /*mdls*/Dict(), /*inis*/Dict());

    // check matrices
    {
        //double tol   = 1.0e-9;
        //double error = 0.0;
        Mat_t K0c(6,6),K1c(6,6),K2c(6,6);
        Mat_t K0,K1,K2;
        dom.Eles[0]->CalcK(K0);
        dom.Eles[1]->CalcK(K1);
        dom.Eles[2]->CalcK(K2);
        cout << "K0 = \n" << PrintMatrix(K0);
        cout << "K1 = \n" << PrintMatrix(K1);
        cout << "K2 = \n" << PrintMatrix(K2);
        //error += CompareMatrices (K0,K0c);
        //error += CompareMatrices (K1,K1c);
        //error += CompareMatrices (K2,K2c);
        //cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        //cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    }

    // solver
    FEM::Solver sol(dom);
    //sol.CteTg  = true;
    sol.Scheme = FEM::Solver::FE_t;
    //sol.Scheme = FEM::Solver::NR_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy wz", 0.0,   0.0, 0.0);
    bcs.Set(-200, "ux uy",    0.0,   0.0);
    bcs.Set(-400, "ux uy wz", 0.0, -10.0, 0.0);
    dom.SetBCs (bcs);
    sol.Solve  (1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%11.6g");
    dom.WriteMPY     ("ex46_res", /*sf*/0.01);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    //nod_sol.Set("                  ux                      uy", /*NRows*/4,

    Table ele_sol;
    //ele_sol.Set("                   ea                      sa                      fa", /*NRows*/5,

    // return error flag
    SDPair nod_tol, ele_tol;
    return dom.CheckError (nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
