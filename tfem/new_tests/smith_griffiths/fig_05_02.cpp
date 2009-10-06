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

/*  Smith & Griffiths (2004): Figure 5.2 p171 
 *  =========================================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "fem/elems/tri3.h"
#include "fem/equilibelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "util/maps.h"
#include "util/util.h"
#include "util/fatal.h"
#include "draw.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (9/*verts*/, 8/*cells*/);
    mesh.SetVert   (0, -100, 0.0, 0.0, 0);
    mesh.SetVert   (1, -200, 0.5, 0.0, 0);
    mesh.SetVert   (2, -100, 1.0, 0.0, 0);
    mesh.SetVert   (3,    0, 0.0,-0.5, 0);
    mesh.SetVert   (4,    0, 0.5,-0.5, 0);
    mesh.SetVert   (5,    0, 1.0,-0.5, 0);
    mesh.SetVert   (6,    0, 0.0,-1.0, 0);
    mesh.SetVert   (7,    0, 0.5,-1.0, 0);
    mesh.SetVert   (8,    0, 1.0,-1.0, 0);
    mesh.SetCell   (0,   -1, /*NVerts*/3, 1,0,3);
    mesh.SetCell   (1,   -1, /*NVerts*/3, 1,3,4);
    mesh.SetCell   (2,   -1, /*NVerts*/3, 2,1,4);
    mesh.SetCell   (3,   -1, /*NVerts*/3, 2,4,5);
    mesh.SetCell   (4,   -1, /*NVerts*/3, 4,3,6);
    mesh.SetCell   (5,   -1, /*NVerts*/3, 4,6,7);
    mesh.SetCell   (6,   -1, /*NVerts*/3, 5,4,7);
    mesh.SetCell   (7,   -1, /*NVerts*/3, 5,7,8);
    mesh.SetBryTag (0, 1, -10);
    mesh.SetBryTag (4, 1, -10);
    mesh.SetBryTag (5, 1, -20);
    mesh.SetBryTag (7, 1, -20);
    mesh.WriteMPY  ("fig_05_02",/*OnlyMesh*/false);
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa nip", PROB("Equilib"), GEOM("Tri3"), TRUE, TRUE, 3.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"),  1.0e+6, 0.3,  TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    /*
    Mat_t M;
    dom.Eles[0]->CalcM(M);
    M *= 12./0.125;
    cout << PrintMatrix(M);
    */

    // solver
    FEM::Solver sol(dom);
    //sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "ux",  0.0)
       .Set( -20, "uy",  0.0)
       .Set(-100, "fy", -0.25)
       .Set(-200, "fy", -0.5);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.Solve ();

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults (cout, Util::_6s);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////
    
    Table nod_sol;
    nod_sol.Set("ux uy", dom.Nods.Size(),
                 0.000E+00, -9.100E-07,
                 1.950E-07, -9.100E-07,
                 3.900E-07, -9.100E-07,
                 0.000E+00, -4.550E-07,
                 1.950E-07, -4.550E-07,
                 3.900E-07, -4.550E-07,
                 0.000E+00,  0.000E+00,
                 1.950E-07,  0.000E+00,
                 3.900E-07,  0.000E+00);

    Table ele_sol;
    ele_sol.Set("sx sy sxy  ex ey exy", dom.Eles.Size(),
                 0.0, -1.0, 0.0,  3.90E-07, -9.10E-07, 0.0,
                 0.0, -1.0, 0.0,  3.90E-07, -9.10E-07, 0.0,
                 0.0, -1.0, 0.0,  3.90E-07, -9.10E-07, 0.0,
                 0.0, -1.0, 0.0,  3.90E-07, -9.10E-07, 0.0,
                 0.0, -1.0, 0.0,  3.90E-07, -9.10E-07, 0.0,
                 0.0, -1.0, 0.0,  3.90E-07, -9.10E-07, 0.0,
                 0.0, -1.0, 0.0,  3.90E-07, -9.10E-07, 0.0,
                 0.0, -1.0, 0.0,  3.90E-07, -9.10E-07, 0.0);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy", 1.0e-15, 1.0e-15);
    ele_tol.Set("sx sy sz sxy  ex ey ez exy", 1.0e-15,1.0e-14,1.0e-15,1.0e-15, 1.0e-15,1.0e-15,1.0e-15,1.0e-15);

    // return error flag
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
