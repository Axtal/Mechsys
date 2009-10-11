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

double BendingMoment (FEM::Cells const & C, double L, double x)
{
    double 
}

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
    mesh.SetCell  (0,   -1, /*NVerts*/2, 0, 1);
    mesh.SetCell  (1,   -1, /*NVerts*/2, 1, 2);
    mesh.SetCell  (2,   -2, /*NVerts*/2, 2, 3);
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
        double tol   = 1.0e-9;
        double error = 0.0;
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
    //sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy wz", 0.0,   0.0, 0.0);
    bcs.Set(-200, "ux uy",    0.0,   0.0);
    bcs.Set(-400, "ux uy wz", 0.0, -10.0, 0.0);
    dom.SetBCs (bcs);
    sol.Solve  (1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults (cout);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    //nod_sol.Set("                  ux                      uy", /*NRows*/4,

    Table ele_sol;
    //ele_sol.Set("                   ea                      sa                      fa", /*NRows*/5,

    // bending moment
    Mat_t M(7,2);
    M = 0.0, -3.1642857142857132e+01,
        1.0, -1.1035714285714279e+01,
        2.0,  9.5714285714285747e+00,
        2.0,  9.5714285714285641e+00,
        3.0,  1.2178571428571418e+01,
        4.0,  1.4785714285714272e+01,
        4.0,  1.4785714285714279e+01,
        5.0,  1.2392857142857141e+01,
        6.0,  0.0000000000000000e+00;
    for (size_t i=0; i<7; ++i)
    {
        int    idx;
        double r;
        if      (x<L)    { idx = 0;  r = x/L;        }
        else if (x<2.*L) { idx = 1;  r = (x-L)/L;    }
        else             { idx = 3;  r = (x-2.*L)/L; }
        cout << "r = " << r << "  cell = " << idx;
    }

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy",    1.0e-15,1.0e-15);
    ele_tol.Set("ea sa fa", 1.0e-15,1.0e-13,1.0e-10);

    // return error flag
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
