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

/*  Bhatti (2005): Example 5.3, p337 *
 *  ================================ */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "fem/elems/quad4.h"
#include "fem/flowelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linflow.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (13/*verts*/, 6/*cells*/);
    mesh.SetVert   ( 0,    0,  0.0,   0.03);
    mesh.SetVert   ( 1,    0,  0.015, 0.03);
    mesh.SetVert   ( 2,    0,  0.03,  0.03);
    mesh.SetVert   ( 3,    0,  0.0,   0.015);
    mesh.SetVert   ( 4,    0,  0.015, 0.015);
    mesh.SetVert   ( 5,    0,  0.03,  0.015);
    mesh.SetVert   ( 6,    0,  0.045, 0.015);
    mesh.SetVert   ( 7,    0,  0.06,  0.015);
    mesh.SetVert   ( 8, -100,  0.0,   0.0);
    mesh.SetVert   ( 9, -100,  0.015, 0.0);
    mesh.SetVert   (10, -100,  0.03,  0.0);
    mesh.SetVert   (11, -100,  0.045, 0.0);
    mesh.SetVert   (12, -100,  0.06,  0.0);
    mesh.SetCell   ( 0,   -1, Array<int>( 8, 9, 4, 3));
    mesh.SetCell   ( 1,   -1, Array<int>( 3, 4, 1, 0));
    mesh.SetCell   ( 2,   -1, Array<int>( 9,10, 5, 4));
    mesh.SetCell   ( 3,   -1, Array<int>( 4, 5, 2, 1));
    mesh.SetCell   ( 4,   -1, Array<int>(10,11, 6, 5));
    mesh.SetCell   ( 5,   -1, Array<int>(11,12, 7, 6));
    mesh.SetBryTag (0, 3, -40);
    mesh.SetBryTag (1, 2, -30);  mesh.SetBryTag (1, 3, -40);
    mesh.SetBryTag (3, 1, -30);  mesh.SetBryTag (3, 2, -30);
    mesh.SetBryTag (4, 2, -30);
    mesh.SetBryTag (5, 1, -20);  mesh.SetBryTag (5, 2, -30);
    mesh.WriteMPY  ("ex53", /*OnlyMesh*/false);

    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom", PROB("Flow"), GEOM("Quad4"));

    // models
    Dict mdls;
    mdls.Set(-1, "name k", MODEL("LinFlow"), 45.0);

    // initial values
    Dict inis;
    inis.Set(-1, "vx vy", 0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    FEM::Solver sol(dom);

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-100, "H",           110.0);
    bcs.Set( -20, "flux",        0.0);
    bcs.Set( -30, "conv h Tinf", 1.0, 55.0, 20.0);
    bcs.Set( -40, "flux",        8000.0);
    bcs.Set(  -1, "s",           5.0e+6);
    dom.SetBCs (bcs);
    sol.Solve  ();

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%11.6g");

    //////////////////////////////////////////////////////////////////////////////////////// Check /////

    // correct solution
    Table nod_sol;
    nod_sol.Set("H", dom.Nods.Size(),
                 1.549619804450071e+02,
                 1.512282609195495e+02,
                 1.486731378550244e+02,
                 1.454325133938252e+02,
                 1.425208024332478e+02,
                 1.348705268525763e+02,
                 1.224358952463445e+02,
                 1.210878374484298e+02,
                 1.100000000000000e+02,
                 1.100000000000000e+02,
                 1.100000000000000e+02,
                 1.100000000000000e+02,
                 1.100000000000000e+02);

    Table ele_sol;
    ele_sol.Set("gx  gy", dom.Eles.Size(),
                 -9.705703201924695e+01,  2.265110527569101e+03,
                 -2.215143495345008e+02,  6.078975179161185e+02,
                 -2.550091860223839e+02,  1.913044309527470e+03,
                 -3.401799548398884e+02,  7.503356496249911e+02,
                 -4.144877202077272e+02,  1.243547403297359e+03,
                 -4.493525993048971e+01,  7.841244231591417e+02);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("H", 1.0e-12);
    ele_tol.Set("gx gy", 1.0e-11, 1.0e-11);

    // return error flag
    return dom.CheckError (nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
