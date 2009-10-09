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

/*  Smith & Griffiths (2004): Figure 5.11 p178 
 *  ==========================================  */

// STL
#include <iostream>

// MechSys
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "fem/elems/quad4.h"
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
    mesh.SetSize   (12/*verts*/, 6/*cells*/);
    mesh.SetVert   ( 0,  0,  0.0,   0.0, 0);
    mesh.SetVert   ( 1,  0,  0.0,  -5.0, 0);
    mesh.SetVert   ( 2,  0,  0.0, -10.0, 0);
    mesh.SetVert   ( 3,  0, 10.0,   0.0, 0);
    mesh.SetVert   ( 4,  0, 10.0,  -5.0, 0);
    mesh.SetVert   ( 5,  0, 10.0, -10.0, 0);
    mesh.SetVert   ( 6,  0, 20.0,   0.0, 0);
    mesh.SetVert   ( 7,  0, 20.0,  -5.0, 0);
    mesh.SetVert   ( 8,  0, 20.0, -10.0, 0);
    mesh.SetVert   ( 9,  0, 30.0,   0.0, 0);
    mesh.SetVert   (10,  0, 30.0,  -5.0, 0);
    mesh.SetVert   (11,  0, 30.0, -10.0, 0);
    mesh.SetCell   ( 0, -1, /*NVerts*/4, 0,1,4,3);
    mesh.SetCell   ( 1, -1, /*NVerts*/4, 1,2,5,4);
    mesh.SetCell   ( 2, -1, /*NVerts*/4, 3,4,7,6);
    mesh.SetCell   ( 3, -1, /*NVerts*/4, 4,5,8,7);
    mesh.SetCell   ( 4, -1, /*NVerts*/4, 6,7,10,9);
    mesh.SetCell   ( 5, -1, /*NVerts*/4, 7,8,11,10);
    mesh.SetBryTag ( 0, 0, -10);
    mesh.SetBryTag ( 0, 3, -30);
    mesh.SetBryTag ( 1, 0, -10);
    mesh.SetBryTag ( 1, 1, -20);
    mesh.SetBryTag ( 3, 1, -20);
    mesh.SetBryTag ( 4, 2, -10);
    mesh.SetBryTag ( 5, 1, -20);
    mesh.SetBryTag ( 5, 2, -10);
    mesh.WriteMPY  ("fig_05_11",/*OnlyMesh*/false);
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa nip", PROB("Equilib"), GEOM("Quad4"), TRUE, TRUE, 4.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"),  1.0e+6, 0.3,  TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    FEM::Solver sol(dom);
    //sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "ux",     0.0);
    bcs.Set( -20, "ux uy",  0.0);
    bcs.Set( -30, "uy",    -1.0e-5);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.Solve ();

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults (cout, Util::_6s);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////
    
    Table nod_sol;
    nod_sol.Set("ux uy", dom.Nods.Size(),
       0.000000000000000E+00, -1.000000000000000E-05,
       0.000000000000000E+00, -5.152429576289283E-06,
       0.000000000000000E+00,  0.000000000000000E+00,
       8.101475206125646E-08, -1.000000000000000E-05,
       1.582094832251061E-06, -4.593608955706439E-06,
       0.000000000000000E+00,  0.000000000000000E+00,
       1.240701017391965E-07,  1.257758501770675E-06,
       1.472140401416949E-06,  1.953453413631710E-07,
       0.000000000000000E+00,  0.000000000000000E+00,
       0.000000000000000E+00,  2.815484477779523E-07,
       0.000000000000000E+00,  3.474895306354692E-07,
       0.000000000000000E+00,  0.000000000000000E+00);


    Table ele_sol;
    ele_sol.Set("sx sy sxy  ex ey exy", dom.Eles.Size(),
       -4.796346390494155E-01, -1.332366671768454E+00, -4.698729954241076E-02,  8.315548045472907E-08, -1.025396162080021E-06, -1.221669788102680E-07/2.0,
       -4.557843083899548E-01, -1.266329393400422E+00,  7.159635270311812E-02,  7.910474279130555E-08, -9.746038677223013E-07,  1.861505170281071E-07/2.0,
       -2.551169352573997E-01, -5.866960435296958E-01,  1.990079267943947E-01, -3.344954107652315E-09, -4.343977948616375E-07,  5.174206096654264E-07/2.0,
       -2.611467606404129E-01, -5.952457193824300E-01,  2.095658639767290E-01, -5.497721623628055E-09, -4.398263679882503E-07,  5.448712463394955E-07/2.0,
       -4.994847244422382E-02,  8.809593954572581E-02, -6.769627916405944E-02, -7.981052634707675E-08,  9.964720923985780E-08, -1.760103258265546E-07/2.0,
       -6.776897772067103E-02,  3.060833703042781E-02,  5.954663535089467E-02, -7.360702116767749E-08,  5.428348800875101E-08,  1.548212519123262E-07/2.0);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy", 1.0e-14, 1.0e-13);
    ele_tol.Set("sx sy sz sxy  ex ey ez exy", 1.0e-8,1.0e-7,1.0e-8,1.0e-8, 1.0e-13,1.0e-13,1.0e-13,1.0e-13);

    // return error flag
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
