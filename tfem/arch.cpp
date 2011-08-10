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

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/unstructured.h>
#include <mechsys/fem/elems/tri3.h>
#include <mechsys/fem/elems/tri6.h>
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

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    double A  = 0.5;
    bool   o2 = false;
    #include "arch.h"
    mesh.WriteMPY ("arch", false);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    if (o2) prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Tri6"), 1.0);
    else    prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Tri3"), 1.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 10.0, 0.2, 1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    Array<int> out_verts(-5,-6,-7); // 0,192
    FEM::Domain dom(mesh, prps, mdls, inis, "arch", &out_verts);

    // solver
    SDPair flags;
    FEM::STDSolver sol(dom, flags);

    ////////////////////////////////////////////////////////////////////////////////////////// Run /////
    
    Dict bcs;
    bcs.Set(-50, "uy",  0.0);
    bcs.Set(-51, "ux",  0.0);
    bcs.Set( -5, "fy", -10.0);
    dom.SetBCs (bcs);
    sol.Solve  (/*NDiv*/1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.WriteVTU ("arch");

    return 0;
}
MECHSYS_CATCH
