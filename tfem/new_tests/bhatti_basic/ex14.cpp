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

/*  Bhatti (2005): Example 1.4, p25  *
 *  ===============================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "fem/rod.h"
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

    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (4/*nodes*/, 5/*cells*/);
    mesh.SetVert  (0, -100,  0.0,  0.0);
    mesh.SetVert  (1, -200, 1500, 3500);
    mesh.SetVert  (2,    0,  0.0, 5000);
    mesh.SetVert  (3, -100, 5000, 5000);
    mesh.SetCell  (0,   -1, /*NVerts*/2, 0, 1);
    mesh.SetCell  (1,   -1, /*NVerts*/2, 1, 3);
    mesh.SetCell  (2,   -2, /*NVerts*/2, 0, 2);
    mesh.SetCell  (3,   -2, /*NVerts*/2, 2, 3);
    mesh.SetCell  (4,   -3, /*NVerts*/2, 1, 2);
    //mesh.WriteVTU ("ex14");
    //cout << mesh << endl;

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob active E A fra", PROB("Rod"), TRUE, 200000.0, 4000.0, TRUE)
        .Set(-2, "prob active E A fra", PROB("Rod"), TRUE, 200000.0, 3000.0, TRUE)
        .Set(-3, "prob active E A fra", PROB("Rod"), TRUE,  70000.0, 2000.0, TRUE);

    // domain
    FEM::Domain dom(mesh, prps, /*mdls*/Dict(), /*inis*/Dict());

    // check matrices
    {
        double tol   = 1.0e-9;
        double error = 0.0;
        Mat_t K0c(4,4),K1c(4,4),K2c(4,4),K3c(4,4),K4c(4,4);
        K0c =
          3.2600217813448358e+04,  7.6067174898046171e+04, -3.2600217813448358e+04, -7.6067174898046171e+04,
          7.6067174898046171e+04,  1.7749007476210775e+05, -7.6067174898046171e+04, -1.7749007476210775e+05,
         -3.2600217813448358e+04, -7.6067174898046171e+04,  3.2600217813448358e+04,  7.6067174898046171e+04,
         -7.6067174898046171e+04, -1.7749007476210775e+05,  7.6067174898046171e+04,  1.7749007476210775e+05;
        K1c =
          1.7749007476210775e+05,  7.6067174898046171e+04, -1.7749007476210775e+05, -7.6067174898046171e+04,
          7.6067174898046171e+04,  3.2600217813448358e+04, -7.6067174898046171e+04, -3.2600217813448358e+04,
         -1.7749007476210775e+05, -7.6067174898046171e+04,  1.7749007476210775e+05,  7.6067174898046171e+04,
         -7.6067174898046171e+04, -3.2600217813448358e+04,  7.6067174898046171e+04,  3.2600217813448358e+04;
        K2c =
          0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
          0.0000000000000000e+00,  1.2000000000000000e+05,  0.0000000000000000e+00, -1.2000000000000000e+05,
          0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
          0.0000000000000000e+00, -1.2000000000000000e+05,  0.0000000000000000e+00,  1.2000000000000000e+05;
        K3c =
          1.2000000000000000e+05,  0.0000000000000000e+00, -1.2000000000000000e+05,  0.0000000000000000e+00,
          0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,
         -1.2000000000000000e+05,  0.0000000000000000e+00,  1.2000000000000000e+05,  0.0000000000000000e+00,
          0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00;
        K4c =
          3.2998316455372231e+04, -3.2998316455372231e+04, -3.2998316455372231e+04,  3.2998316455372231e+04,
         -3.2998316455372231e+04,  3.2998316455372231e+04,  3.2998316455372231e+04, -3.2998316455372231e+04,
         -3.2998316455372231e+04,  3.2998316455372231e+04,  3.2998316455372231e+04, -3.2998316455372231e+04,
          3.2998316455372231e+04, -3.2998316455372231e+04, -3.2998316455372231e+04,  3.2998316455372231e+04;
        Mat_t K0,K1,K2,K3,K4;
        dom.Eles[0]->CalcK(K0);
        dom.Eles[1]->CalcK(K1);
        dom.Eles[2]->CalcK(K2);
        dom.Eles[3]->CalcK(K3);
        dom.Eles[4]->CalcK(K4);
        error += CompareMatrices (K0,K0c);
        error += CompareMatrices (K1,K1c);
        error += CompareMatrices (K2,K2c);
        error += CompareMatrices (K3,K3c);
        error += CompareMatrices (K4,K4c);
        cout << "\n[1;37m--- Matrices: Error ----------------------------------------------------------[0m\n";
        cout << "error (K) = " << (error>tol ? "[1;31m" : "[1;32m") << error << "[0m" << endl;
    }

    // solver
    FEM::Solver sol(dom);
    //sol.CteTg  = true;
    //sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy", 0.0,0.0)
       .Set(-200, "fy", -150000.0);
    dom.SetBCs (bcs);
    sol.Solve  (1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults (cout);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("                  ux                      uy", /*NRows*/4,
                0.000000000000000e+00,  0.000000000000000e+00,
                5.389536380057676e-01, -9.530613006371175e-01,
                2.647036149579491e-01, -2.647036149579490e-01,
                0.000000000000000e+00,  0.000000000000000e+00);

    Table ele_sol;
    ele_sol.Set("                   ea                      sa                      fa", /*NRows*/5,
                -1.742954548428455e-04, -3.485909096856910e+01, -1.394363638742764e+05,
                -3.149970910789726e-05, -6.299941821579453e+00, -2.519976728631781e+04,
                -5.294072299158981e-05, -1.058814459831796e+01, -3.176443379495388e+04,
                -5.294072299158982e-05, -1.058814459831796e+01, -3.176443379495389e+04,
                 3.208692362423290e-04,  2.246084653696303e+01,  4.492169307392606e+04);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy",    1.0e-15,1.0e-15);
    ele_tol.Set("ea sa fa", 1.0e-15,1.0e-13,1.0e-10);

    // return error flag
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
