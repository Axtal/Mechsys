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

/*  Bhatti (2005): Example 1.5, p28  *
 *  ===============================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/elems/tri3.h>
#include <mechsys/fem/flowelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solvers/stdsolver.h>
#include <mechsys/models/linflow.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (5/*verts*/, 4/*cells*/);
    mesh.SetVert   (0, -100, 0.0, 0.0);
    mesh.SetVert   (1,    0, 0.2, 0.0);
    mesh.SetVert   (2,    0, 0.2, 0.3);
    mesh.SetVert   (3, -100, 0.0, 0.1);
    mesh.SetVert   (4,    0, 0.1, 0.1);
    mesh.SetCell   (0,   -1, Array<int>(0,1,4));
    mesh.SetCell   (1,   -1, Array<int>(1,2,4));
    mesh.SetCell   (2,   -1, Array<int>(2,3,4));
    mesh.SetCell   (3,   -1, Array<int>(0,4,3));
    mesh.SetBryTag (1, 0, -10);
    mesh.WriteMPY  ("ex15");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom", PROB("Flow"), GEOM("Tri3"));

    // models
    Dict mdls;
    mdls.Set(-1, "name k", MODEL("LinFlow"), 1.4);

    // initial values
    Dict inis;
    inis.Set(-1, "vx vy", 0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    SDPair flags;
    FEM::STDSolver sol(dom, flags);

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "conv h Tinf", 1.0, 27.0, 20.0);
    bcs.Set(-100, "hh", 300.0);
    dom.SetBCs (bcs);
    sol.Solve  ();

    // check matrices
    if (false)
    {
        double tol   = 1.0e-10;
        double error = 0.0;
        Mat_t K0c(3,3),K1c(3,3),K2c(3,3),K3c(3,3);
        K0c =
          0.7,  0.0, -0.7,
          0.0,  0.7, -0.7,
         -0.7, -0.7,  1.4;
        K1c =
          3.8666666666666663e+00,  1.5833333333333333e+00, -1.4,
          1.5833333333333333e+00,  3.1666666666666665e+00, -0.7,
         -1.4,                    -0.7,                     2.1;
        K2c =
          0.35,  0.35, -0.7,
          0.35,  1.75, -2.1,
         -0.7,  -2.1,   2.8;
        K3c =
          0.7,  0.0, -0.7,
          0.0,  0.7, -0.7,
         -0.7, -0.7,  1.4;
        Mat_t K0,K1,K2,K3;
        dom.Eles[0]->CalcK(K0);
        dom.Eles[1]->CalcK(K1);
        dom.Eles[2]->CalcK(K2);
        dom.Eles[3]->CalcK(K3);
        error += CompareMatrices (K0,K0c);
        error += CompareMatrices (K1,K1c);
        error += CompareMatrices (K2,K2c);
        error += CompareMatrices (K3,K3c);
        cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    }

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%11.6g");
    //dom.WriteVTU     ("ex15");

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("                  hh                    Rhh", dom.Nods.Size(),
                3.000000000000000e+02, 8.201709351675694e+01,
                9.354661202985511e+01, 0.0,
                2.384369969266794e+01, 0.0,
                3.000000000000000e+02, 2.314136689594616e+02,
                1.828327235474901e+02, 0.0);

    Table ele_sol;
    ele_sol.Set("gx  gy", dom.Eles.Size(),
                -1.032266939850724e+03,  -1.394058246743742e+02,
                -1.125204156300308e+03,  -2.323430411239573e+02,
                -1.171672764525099e+03,  -2.091087370115617e+02,
                -1.171672764525099e+03,   0.000000000000000e+00);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("hh Rhh", 1.0e-13,1.0e-12);
    ele_tol.Set("gx gy",  1.0e-12,1.0e-12);

    // return error flag
    bool err1 = dom.CheckErrorNods(nod_sol, nod_tol);
    bool err2 = dom.CheckErrorEles(ele_sol, ele_tol);
    return (err1 || err2);
}
MECHSYS_CATCH
