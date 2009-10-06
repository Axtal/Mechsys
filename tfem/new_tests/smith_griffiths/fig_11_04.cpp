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
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "fem/elems/quad8.h"
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

void calc_F(double t, double & fx, double & fy, double & fz)
{
    double omega = 0.3;
    fx = 0.0;
    fy = cos(omega*t);
    fz = 0.0;
}

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
    mesh.SetCell    ( 0,   -1, /*NVerts*/8, 2,7,5,0, 4,6,3,1);
    mesh.SetCell    ( 1,   -1, /*NVerts*/8, 7,12,10,5, 9,11,8,6);
    mesh.SetCell    ( 2,   -1, /*NVerts*/8, 12,17,15,10, 14,16,13,11);
    mesh.SetBryTag  ( 0, 3, -10);
    mesh.WriteMPY   ("fig_11_04",/*OnlyMesh*/false);
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa rho", PROB("Equilib"), GEOM("Quad8"), TRUE, TRUE, 1.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 1.0, 0.3, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.SetOutNods ("fig_11_04",/*NNod*/1, /*IDs*/17);
    dom.FFuncs[-100] = &calc_F; // set database of callbacks

    // solver
    FEM::Solver sol(dom);
    //sol.DScheme = FEM::Solver::SS22_t;
    sol.DScheme = FEM::Solver::GN22_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "ux uy", 0.0)
       .Set(-100, "ffunc", 0.0);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.DynSolve (/*tf*/100.0, /*dt*/1.0, /*dtOut*/1.0);

    return 0.0;
}
MECHSYS_CATCH
