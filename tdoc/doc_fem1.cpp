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
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/rod.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solvers/stdsolver.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try {
    // mesh
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize  (4/*nodes*/, 5/*cells*/);
    mesh.SetVert  (0, -100,  0.0,  0.0);
    mesh.SetVert  (1, -200, 1500, 3500);
    mesh.SetVert  (2,    0,  0.0, 5000);
    mesh.SetVert  (3, -100, 5000, 5000);
    mesh.SetCell  (0,   -1, Array<int>(0, 1));
    mesh.SetCell  (1,   -1, Array<int>(1, 3));
    mesh.SetCell  (2,   -2, Array<int>(0, 2));
    mesh.SetCell  (3,   -2, Array<int>(2, 3));
    mesh.SetCell  (4,   -3, Array<int>(1, 2));
    mesh.WriteMPY ("doc_fem1");

    // elements properties
    Dict prps;
    prps.Set(-1, "prob active E A fra", PROB("Rod"), 1.0, 200000.0, 4000.0, 1.0);
    prps.Set(-2, "prob active E A fra", PROB("Rod"), 1.0, 200000.0, 3000.0, 1.0);
    prps.Set(-3, "prob active E A fra", PROB("Rod"), 1.0,  70000.0, 2000.0, 1.0);

    // domain
    FEM::Domain dom(mesh, prps, /*mdls*/Dict(), /*inis*/Dict());

    // solver
    SDPair flags;
    FEM::STDSolver sol(dom, flags);

    // set boundary conditions and solve
    Dict bcs;
    bcs.Set(-100, "ux uy", 0.0,0.0);
    bcs.Set(-200, "fy", -150000.0);
    dom.SetBCs (bcs);
    sol.Solve  (1);

    // output
    dom.PrintResults ("%11.6g");
    dom.WriteVTU     ("doc_fem1"); // for Paraview
} MECHSYS_CATCH