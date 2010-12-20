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

/*  Smith & Griffiths (2004): Figure 11.4 p475 
 *  ==========================================  */

// STL
#include <iostream>
#include <cmath>    // for cos

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/rksolver.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/util.h>
#include <mechsys/util/fatal.h>
#include <mechsys/draw.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
const double TRUE = 1.0;

class BCF : public FEM::BCFuncs
{
public:
    double fm (double t) { return cos(0.3*t); }
};

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    double H = 1.0;  // height
    double h = 0.5;  // half height
    double L = 4.0;  // length
    double l = L/3.; // element length
    double m = l/2.; // half length (for mid nodes)
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize    (18/*verts*/, 3/*cells*/);
    mesh.SetVert    ( 0,    0,  0.0,    H, 0);
    mesh.SetVert    ( 1,    0,  0.0,    h, 0);
    mesh.SetVert    ( 2,    0,  0.0,  0.0, 0);
    mesh.SetVert    ( 3,    0,    m,    H, 0);
    mesh.SetVert    ( 4,    0,    m,  0.0, 0);
    mesh.SetVert    ( 5,    0,    l,    H, 0);
    mesh.SetVert    ( 6,    0,    l,    h, 0);
    mesh.SetVert    ( 7,    0,    l,  0.0, 0);
    mesh.SetVert    ( 8,    0,  l+m,    H, 0);
    mesh.SetVert    ( 9,    0,  l+m,  0.0, 0);
    mesh.SetVert    (10,    0,  l+l,    H, 0);
    mesh.SetVert    (11,    0,  l+l,    h, 0);
    mesh.SetVert    (12,    0,  l+l,  0.0, 0);
    mesh.SetVert    (13,    0,  L-m,    H, 0);
    mesh.SetVert    (14,    0,  L-m,  0.0, 0);
    mesh.SetVert    (15,    0,    L,    H, 0);
    mesh.SetVert    (16,    0,    L,    h, 0);
    mesh.SetVert    (17, -100,    L,  0.0, 0);
    mesh.SetCell    ( 0,   -1, Array<int>( 2, 7, 5, 0,  4, 6, 3, 1));
    mesh.SetCell    ( 1,   -1, Array<int>( 7,12,10, 5,  9,11, 8, 6));
    mesh.SetCell    ( 2,   -1, Array<int>(12,17,15,10, 14,16,13,11));
    mesh.SetBryTag  ( 0, 3, -10);
    mesh.WriteMPY   ("fig_11_04");
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa rho", PROB("Equilib"), GEOM("Quad8"), 1.0, 1.0, 1.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 1.0, 0.3, 1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    BCF bcf;
    Array<int> out_nods(17, /*justone*/true);
    FEM::Domain dom(mesh, prps, mdls, inis, "fig_11_04", &out_nods);
    dom.MFuncs[-100] = &bcf; // set database of callbacks

    // solver
    FEM::RKSolver sol(dom);
    sol.DampAm = 0.005;
    sol.DampAk = 0.272;
    sol.DampTy = FEM::RKSolver::Rayleigh_t;
    sol.Scheme = "RK4I";
    sol.STOL   = 1.0e-3;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "ux uy",  0.0, 0.0);
    bcs.Set(-100, "fy bcf", 1.0, TRUE);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.DynSolve (/*tf*/100, /*dt*/1.0, /*dtOut*/1.0, "fig_11_04");

    return 0.0;
}
MECHSYS_CATCH
