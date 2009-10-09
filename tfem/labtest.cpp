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
#include "mesh/structured.h"
#include "fem/elems/hex8.h"
#include "fem/elems/hex20.h"
#include "fem/equilibelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "models/elastoplastic.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    double qx = -10.0;
    double qy = -10.0;
    double uz = -0.02;

    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    bool   o2 = false;
    size_t nd = 1;
    Mesh::Structured mesh(/*NDim*/3);
    mesh.GenBox (/*O2*/o2,/*Nx*/nd,/*Ny*/nd,/*Nz*/nd);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    if (o2) prps.Set(-1, "prob geom", PROB("Equilib"), GEOM("Hex20"));
    else    prps.Set(-1, "prob geom", PROB("Equilib"), GEOM("Hex8"));

    // models
    Dict mdls;
    //mdls.Set(-1, "name E nu", MODEL("LinElastic"), 10.0, 0.2);
    mdls.Set(-1, "name E nu fc sY", MODEL("ElastoPlastic"), 1.0, 0.3, FAILCRIT("VM"), 2.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy syz szx", -10.0,-10.0,-10.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.SetOutNods ("labtest", /*NNod*/1, /*WithTags*/false, /*ID*/(nd==3 ?  6 : 6));
    dom.SetOutEles ("labtest", /*NEle*/1, /*ID*/(nd==3 ? 13 : 0));

    // solver
    FEM::Solver sol(dom);
    sol.Scheme = FEM::Solver::FE_t;
    //sol.Scheme = FEM::Solver::NR_t;

    ////////////////////////////////////////////////////////////////////////////////////////// Run /////
    
    Dict bcs;
    bcs.Set(-10, "ux", 0.0);
    bcs.Set(-30, "uy", 0.0);
    bcs.Set(-50, "uz", 0.0);
    bcs.Set(-20, "qn", qx);
    bcs.Set(-40, "qn", qy);
    bcs.Set(-60, "uz", uz);
    dom.SetBCs (bcs);
    sol.Solve  (/*NDiv*/10);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.WriteVTU ("labtest");

    return 0;
}
MECHSYS_CATCH
