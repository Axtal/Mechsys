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

/*  Bhatti (2005): Example 1.17, p76  *
 *  ================================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/rod.h>
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

    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (4/*nodes*/, 5/*cells*/);
    mesh.SetVert  (0, -100,     0.0,     0.0);
    mesh.SetVert  (1, -200,     0.0,  3000.0);
    mesh.SetVert  (2, -300,  5000.0, -3000.0);
    mesh.SetVert  (3,    0,  5000.0,  3000.0);
    mesh.SetCell  (0,   -1, Array<int>(0, 2));
    mesh.SetCell  (1,   -1, Array<int>(0, 3));
    mesh.SetCell  (2,   -1, Array<int>(0, 1));
    mesh.SetCell  (3,   -1, Array<int>(1, 3));
    mesh.SetCell  (4,   -1, Array<int>(2, 3));
    //mesh.WriteMPY ("ex117");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob E A fra", PROB("Rod"), 70000.0, 1000.0, 1.0);

    // domain and solver
    FEM::Domain dom(mesh, prps, /*mdls*/Dict(), /*inis*/Dict());
    FEM::Solver sol(dom);

    // stage # 1 -----------------------------------------------------------
    double alpha = 150.0*Util::PI/180.0;
    Dict bcs;
    bcs.Set(-100, "inclsupport alpha", 1.0, alpha);
    bcs.Set(-200, "ux uy", 0.0, 0.0);
    bcs.Set(-300, "fx", 20000.0);
    dom.SetBCs (bcs);
    sol.Solve  (1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%15.6g");
    //dom.WriteVTU     ("ex117");

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("ux uy", /*NRows*/4,
             5.14286,  -2.96923,
             0.0,       0.0,
            16.8629,   12.788,
            -1.42857,  11.7594);
    SDPair nod_tol;
    nod_tol.Set("ux uy",1.0e-4,1.0e-4);

    Table ele_sol;
    ele_sol.Set("N", 5, 23323.8, 23323.8, 69282.0, -20000.0, -12000.0);
    SDPair ele_tol;
    ele_tol.Set("N",1.0e-1);

    // return error flag
    bool err1 = dom.CheckErrorNods(nod_sol, nod_tol);
    bool err2 = dom.CheckErrorEles(ele_sol, ele_tol);
    return (err1 || err2);
}
MECHSYS_CATCH
