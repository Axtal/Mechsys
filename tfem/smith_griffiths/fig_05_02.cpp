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
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/tri3.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>
#include <mechsys/util/fatal.h>
#include <mechsys/draw.h>

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
    mesh.SetCell   (0,   -1, Array<int>(1,0,3));
    mesh.SetCell   (1,   -1, Array<int>(1,3,4));
    mesh.SetCell   (2,   -1, Array<int>(2,1,4));
    mesh.SetCell   (3,   -1, Array<int>(2,4,5));
    mesh.SetCell   (4,   -1, Array<int>(4,3,6));
    mesh.SetCell   (5,   -1, Array<int>(4,6,7));
    mesh.SetCell   (6,   -1, Array<int>(5,4,7));
    mesh.SetCell   (7,   -1, Array<int>(5,7,8));
    mesh.SetBryTag (0, 1, -10);
    mesh.SetBryTag (4, 1, -10);
    mesh.SetBryTag (5, 1, -20);
    mesh.SetBryTag (7, 1, -20);
    //mesh.WriteMPY  ("fig_05_02");
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa nip", PROB("Equilib"), GEOM("Tri3"), 1.0, 1.0, 3.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa rho", MODEL("LinElastic"),  1.0e+6, 0.3,  1.0, 1.0);

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
    sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "ux",  0.0);
    bcs.Set( -20, "uy",  0.0);
    bcs.Set(-100, "fy", -0.25);
    bcs.Set(-200, "fy", -0.5);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.Solve ();

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
    ele_sol.Set("sx sy sxy", dom.Eles.Size()*3/*nip*/,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,

        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,

        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,

        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,

        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,

        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,

        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,

        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0,
        0.0, -1.0, 0.0);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy",      1.0e-15, 1.0e-15);
    ele_tol.Set("sx sy sxy ", 1.0e-15, 1.0e-14, 1.0e-15);

    // return error flag
    bool err_nods = dom.CheckErrorNods (nod_sol, nod_tol);
    bool err_eles = dom.CheckErrorIPs  (ele_sol, ele_tol);
    return (err_nods || err_eles);
}
MECHSYS_CATCH
