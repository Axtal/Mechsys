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

/*  Smith & Griffiths (2004): Figure 11.1 p470 
 *  ==========================================  */

// STL
#include <iostream>
#include <cmath>    // for cos

// MechSys
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "fem/beam.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "util/maps.h"
#include "util/util.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

double Multiplier (double t)
{
    if (t>1.0) return 0.0;
    else       return 3.194*sin(Util::PI*t);
}

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize    (2/*verts*/, 1/*cells*/);
    mesh.SetVert    ( 0, -100,  0.0,  0.0, 0);
    mesh.SetVert    ( 1, -200,  1.0,  0.0, 0);
    mesh.SetCell    ( 0,   -1, Array<int>(0,1));
    mesh.WriteMPY   ("fig_11_01");
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob fra rho E A Izz", PROB("Beam"), TRUE, 1.0, 3.194, 1.0, 1.0);

    // domain
    FEM::Domain dom(mesh, prps, Dict(), Dict());
    dom.SetOutNods ("fig_11_01", Array<int>(1,/*JustOne*/true));
    dom.MFuncs[-200] = &Multiplier; // set database of callbacks

    // solver
    FEM::Solver sol(dom);
    sol.DScheme = FEM::Solver::GN22_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "ux uy wz", 0.0,0.0,0.0);
    bcs.Set(-200, "fy mfunc", 1.0, 0.0);
    dom.SetBCs (bcs);
    sol.DynSolve (/*tf*/1.8, /*dt*/0.05, /*dtOut*/0.05, "fig_11_01");

    return 0.0;
}
MECHSYS_CATCH
