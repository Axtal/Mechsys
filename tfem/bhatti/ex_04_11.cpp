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

/*  Bhatti (2005): Example 4.11 p274  *
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

    double sq2 = sqrt(2.0);
    double L   = 15.*12.;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (3/*nodes*/, 2/*cells*/);
    mesh.SetVert  (0, -100,  0.0,     0.0);
    mesh.SetVert  (1,    0,  L/sq2,   L/sq2);
    mesh.SetVert  (2, -100,  L/sq2+L, L/sq2);
    mesh.SetCell  (0,   -1, Array<int>(0, 1));
    mesh.SetCell  (1,   -2, Array<int>(1, 2));
    mesh.WriteMPY ("ex411");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob  E A Izz  fra", PROB("Beam"),  3.0e+4, 100.0, 1000.0,  TRUE);
    prps.Set(-2, "prob  E A Izz  fra", PROB("Beam"),  3.0e+4, 100.0, 1000.0,  TRUE);

    // domain
    FEM::Domain dom(mesh, prps, /*mdls*/Dict(), /*inis*/Dict());

    // check matrices
    {
        double tol   = 1.0e-9;
        double error = 0.0;
        Mat_t K0c(6,6),K1c(6,6);
        K0c =
          8.3641975308641959e+03,  8.3024691358024684e+03, -3.9283710065919304e+03, -8.3641975308641959e+03, -8.3024691358024684e+03, -3.9283710065919304e+03,
          8.3024691358024684e+03,  8.3641975308641959e+03,  3.9283710065919304e+03, -8.3024691358024684e+03, -8.3641975308641959e+03,  3.9283710065919304e+03,
         -3.9283710065919304e+03,  3.9283710065919304e+03,  6.6666666666666663e+05,  3.9283710065919304e+03, -3.9283710065919304e+03,  3.3333333333333331e+05,
         -8.3641975308641959e+03, -8.3024691358024684e+03,  3.9283710065919304e+03,  8.3641975308641959e+03,  8.3024691358024684e+03,  3.9283710065919304e+03,
         -8.3024691358024684e+03, -8.3641975308641959e+03, -3.9283710065919304e+03,  8.3024691358024684e+03,  8.3641975308641959e+03, -3.9283710065919304e+03,
         -3.9283710065919304e+03,  3.9283710065919304e+03,  3.3333333333333331e+05,  3.9283710065919304e+03, -3.9283710065919304e+03,  6.6666666666666663e+05;
        K1c =
          1.6666666666666668e+04,  0.0000000000000000e+00,  0.0000000000000000e+00, -1.6666666666666668e+04,  0.0000000000000000e+00,  0.0000000000000000e+00,
          0.0000000000000000e+00,  6.1728395061728392e+01,  5.5555555555555557e+03,  0.0000000000000000e+00, -6.1728395061728392e+01,  5.5555555555555557e+03,
          0.0000000000000000e+00,  5.5555555555555557e+03,  6.6666666666666663e+05,  0.0000000000000000e+00, -5.5555555555555557e+03,  3.3333333333333331e+05,
         -1.6666666666666668e+04,  0.0000000000000000e+00,  0.0000000000000000e+00,  1.6666666666666668e+04,  0.0000000000000000e+00,  0.0000000000000000e+00,
          0.0000000000000000e+00, -6.1728395061728392e+01, -5.5555555555555557e+03,  0.0000000000000000e+00,  6.1728395061728392e+01, -5.5555555555555557e+03,
          0.0000000000000000e+00,  5.5555555555555557e+03,  3.3333333333333331e+05,  0.0000000000000000e+00, -5.5555555555555557e+03,  6.6666666666666663e+05;
        Mat_t K0,K1;
        dom.Eles[0]->CalcK(K0);
        dom.Eles[1]->CalcK(K1);
        error += CompareMatrices (K0,K0c);
        error += CompareMatrices (K1,K1c);
        cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    }

    // solver
    FEM::Solver sol(dom);

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy wz", 0.0,0.0,0.0);
    bcs.Set(  -1, "qn", -1./12.);
    dom.SetBCs (bcs);
    sol.Solve  (1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%11.6g");

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("                   ux                      uy                      wz", /*NRows*/dom.Nods.Size(),
                 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
                 6.016073731432885e-04, -1.254737159610633e-03,  1.685087639678406e-04,
                 0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00);

    Table ele_sol;
    ele_sol.Set("                    N                       V                        M", /*NRows*/dom.Eles.Size(),
                -7.697208350099685e+00,  1.017186578051514e+00,   1.405847939946401e+02,
                -1.002678955238814e+01,  8.587068887342608e-01,  -2.808479399464011e+01);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy wz", 1.0e-15,1.0e-15,1.0e-15);
    ele_tol.Set("N  V  M",  1.0e-14,1.0e-14,1.0e-15);

    // return error flag
    return dom.CheckError (nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
