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
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/beam.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solvers/stdsolver.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
const double TRUE = 1.0;

class BCF : public FEM::BCFuncs
{
public:
    double fm (double t) { return (t<1.0 ? 3.194*sin(Util::PI*t) : 0.0); }
};

int main(int argc, char **argv) try
{
    // input
    bool rksolver = true;
    if (argc>1) rksolver = atoi(argv[1]);

    // filekey
    String fnkey;
    if (rksolver) fnkey = "fig_11_01";
    else          fnkey = "fig_11_01_GN22";

    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize    (2/*verts*/, 1/*cells*/);
    mesh.SetVert    ( 0, -100,  0.0,  0.0, 0);
    mesh.SetVert    ( 1, -200,  1.0,  0.0, 0);
    mesh.SetCell    ( 0,   -1, Array<int>(0,1));
    mesh.WriteMPY   (fnkey.CStr());
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob fra rho E A Izz", PROB("Beam"), TRUE, 1.0, 3.194, 1.0, 1.0);

    // domain
    BCF bcf;
    Array<int> out_nods(1, /*justone*/true);
    FEM::Domain dom(mesh, prps, Dict(), Dict(), fnkey.CStr(), &out_nods);
    dom.MFuncs["-200_fy"] = &bcf; // set database of callbacks

    // boundary conditions
    Dict bcs;
    bcs.Set(-100, "ux uy wz", 0.0,0.0,0.0);
    bcs.Set(-200, "FY", 1.0);
    dom.SetBCs (bcs);

    // solver
    SDPair flags;
    FEM::STDSolver sol(dom, flags);
    sol.DynSolve (/*tf*/5, /*dt*/0.05, /*dtOut*/0.05, fnkey.CStr());

    return 0.0;
}
MECHSYS_CATCH
