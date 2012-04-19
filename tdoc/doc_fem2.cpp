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

/*  Bhatti (2005): Example 1.6, p32  *
 *  ===============================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/elems/tri3.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solvers/stdsolver.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try {
    // mesh
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (6/*verts*/, 4/*cells*/);
    mesh.SetVert   (0, -100, 0.0, 0.0, 0);
    mesh.SetVert   (1, -101, 0.0, 2.0, 0);
    mesh.SetVert   (2,    0, 2.0, 0.0, 0);
    mesh.SetVert   (3,    0, 2.0, 1.5, 0);
    mesh.SetVert   (4,    0, 4.0, 0.0, 0);
    mesh.SetVert   (5,    0, 4.0, 1.0, 0);
    mesh.SetCell   (0,   -1, Array<int>(0,2,3));
    mesh.SetCell   (1,   -1, Array<int>(3,1,0));
    mesh.SetCell   (2,   -1, Array<int>(2,4,5));
    mesh.SetCell   (3,   -1, Array<int>(5,3,2));
    mesh.SetBryTag (1, 0, -10);
    mesh.SetBryTag (3, 0, -10);
    mesh.WriteMPY  ("doc_fem2");

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active  h pse", PROB("Equilib"), GEOM("Tri3"), 1.0, 0.25, 1.0);

    // material models
    Dict mdls;
    mdls.Set(-1, "name E nu pse rho", MODEL("LinElastic"), 1.0e+4, 0.2, 1.0, 1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    SDPair flags;
    FEM::STDSolver sol(dom, flags);
    sol.WithInfo = false;

    // solve
    Dict bcs;
    bcs.Set( -10, "qn",   -20.0);
    bcs.Set(-100, "ux uy", 0.0,0.1);
    bcs.Set(-101, "ux uy", 0.33,0.0);
    dom.SetBCs (bcs);
    sol.Solve  (/*NDiv*/1);

    // output
    dom.PrintResults ("%11.6g");
} MECHSYS_CATCH
