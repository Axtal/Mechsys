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
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad4.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/beam.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
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

    bool o2 = false;
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -1.0,   0.0, 0.0,
                 -2.0,  10.0, 0.0,
                 -3.0,  10.0, 1.0,
                 -4.0, 0.0, 1.0,  -10.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (10);
    blks[0].SetNy (1);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate    (blks,/*O2*/o2);
    mesh.AddLinCells (Array<int>(-10, /*JustOne*/true));
    mesh.WriteMPY    ("beam02_mesh");

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    if (o2) prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Quad8"), 1.0);
    else    prps.Set(-1, "prob geom psa", PROB("Equilib"), GEOM("Quad4"), 1.0);
    prps.Set(-10, "prob fra rho E A Izz", PROB("Beam"), 1.0, 1.0, 1.0, 1.0, 1.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 1.0, 0.2, 1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    //dom.SetOutNods ("arch", /*NNods*/2, /*WithTags*/false, /*Ids*/0,192);
    //dom.SetOutNods ("beam", /*NNods*/3, /*WithTags*/true, /*Tags*/-5,-6,-7);

    // solver
    FEM::Solver sol(dom);
    //sol.Scheme = FEM::Solver::FE_t;
    sol.Scheme = FEM::Solver::NR_t;

    ////////////////////////////////////////////////////////////////////////////////////////// Run /////
    
    Dict bcs;
    bcs.Set(-1,  "ux uy", 0.0);
    bcs.Set(-2,  "uy",    0.0);
    bcs.Set(-30, "qn",   -1.0);
    dom.SetBCs (bcs);
    sol.Solve  (/*NDiv*/1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.WriteVTU ("beam02_res");
    dom.WriteMPY ("beam02_res", /*sf*/0.03);

    return 0;
}
MECHSYS_CATCH
