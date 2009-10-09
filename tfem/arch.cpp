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
#include "mesh/unstructured.h"
#include "fem/elems/tri3.h"
#include "fem/elems/tri6.h"
#include "fem/equilibelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "util/maps.h"
#include "util/fatal.h"

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

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    if (o2) prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Tri6"), TRUE);
    else    prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Tri3"), TRUE);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 10.0, 0.2, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    //dom.SetOutNods ("arch", /*NNods*/2, /*WithTags*/false, /*Ids*/0,192);
    dom.SetOutNods ("arch", /*NNods*/3, /*WithTags*/true, /*Tags*/-5,-6,-7);

    // solver
    FEM::Solver sol(dom);
    sol.Scheme = FEM::Solver::FE_t;
    //sol.Scheme = FEM::Solver::NR_t;

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
