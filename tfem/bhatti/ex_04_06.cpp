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
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/beam.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

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

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob E A Izz fra", PROB("Beam"), 1.0, 100.0, 180.0, 1.0);
    prps.Set(-2, "prob E A Izz fra", PROB("Beam"), 2.0, 100.0, 180.0, 1.0);

    // domain and solver
    FEM::Domain dom(mesh, prps, /*mdls*/Dict(), /*inis*/Dict());
    FEM::Solver sol(dom);
    //cout << dom << endl;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy wz", 0.0,   0.0, 0.0);
    bcs.Set(-200, "ux uy",    0.0,   0.0);
    bcs.Set(-400, "ux uy wz", 0.0, -10.0, 0.0);
    dom.SetBCs (bcs);
    sol.Solve  (1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%11.6g");
    //dom.WriteMPY     ("ex_04_06_res");

    int i = sol.NEq-sol.NLag;
    cout << "\n[1;31mlambda0 = " << sol.U(i) << ",    lambda1 = " << sol.U(i+1) << endl << "[0m\n";

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set(" ux    uy         wz         M", /*NRows*/5,
                 0.0,  0.0,       0.0,        1.0/75.,
                 0.0,  0.0,      -1.0/90.0,  -2.0/75.,
                 0.0, -70.0/9.0, -1.0/30.0,   0.0,
                 0.0, -10.0,      0.0,        2.0/75.,
                 0.0, -70.0/9.0, -1.0/90.0,   0.0);
    SDPair nod_tol;
    nod_tol.Set("ux uy wz M",1.0e-15,1.0e-14,1.0e-15,1.0e-15);

    // return error flag
    return dom.CheckErrorNods (nod_sol, nod_tol);
}
MECHSYS_CATCH
