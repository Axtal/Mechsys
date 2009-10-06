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

/*  Smith & Griffiths (2004): Figure 5.15 p181 
 *  ==========================================  */

// STL
#include <iostream>

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

int main(int argc, char **argv) try
{
    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize    (12/*verts*/, 6/*cells*/);
    mesh.SetVert    ( 0,  0,  0.0,   0.0, 0);
    mesh.SetVert    ( 1,  0,  0.0,  -3.0, 0);
    mesh.SetVert    ( 2,  0,  0.0,  -6.0, 0);
    mesh.SetVert    ( 3,  0,  0.0,  -9.0, 0);
    mesh.SetVert    ( 4,  0,  3.0,   0.0, 0);
    mesh.SetVert    ( 5,  0,  3.0,  -3.0, 0);
    mesh.SetVert    ( 6,  0,  3.0,  -6.0, 0);
    mesh.SetVert    ( 7,  0,  3.0,  -9.0, 0);
    mesh.SetVert    ( 8,  0,  6.0,   0.0, 0);
    mesh.SetVert    ( 9,  0,  6.0,  -3.0, 0);
    mesh.SetVert    (10,  0,  6.0,  -6.0, 0);
    mesh.SetVert    (11,  0,  6.0,  -9.0, 0);
    mesh.SetCell    ( 0, -1, /*NVerts*/4, 0,1,5,4);
    mesh.SetCell    ( 1, -1, /*NVerts*/4, 1,2,6,5);
    mesh.SetCell    ( 2, -1, /*NVerts*/4, 2,3,7,6);
    mesh.SetCell    ( 3, -1, /*NVerts*/4, 4,5,9,8);
    mesh.SetCell    ( 4, -1, /*NVerts*/4, 5,6,10,9);
    mesh.SetCell    ( 5, -1, /*NVerts*/4, 6,7,11,10);
    mesh.SetBryTag  ( 0, 0, -10);
    mesh.SetBryTag  ( 0, 3, -30);
    mesh.SetBryTag  ( 1, 0, -10);
    mesh.SetBryTag  ( 2, 0, -10);
    mesh.SetBryTag  ( 2, 1, -20);
    mesh.SetBryTag  ( 3, 2, -10);
    mesh.SetBryTag  ( 4, 2, -10);
    mesh.SetBryTag  ( 5, 1, -20);
    mesh.SetBryTag  ( 5, 2, -10);
    mesh.GenO2Verts ();
    mesh.WriteMPY   ("fig_05_15",/*OnlyMesh*/false);
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa nip", PROB("Equilib"), GEOM("Quad8"), TRUE, TRUE, 4.0);

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
    bcs.Set( -10, "ux",     0.0)
       .Set( -20, "ux uy",  0.0)
       .Set( -30, "qn",    -1.0);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.Solve ();

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults (cout, Util::_6s);

    //////////////////////////////////////////////////////////////////////////////////////// Check /////
    
    Table nod_sol;
    nod_sol.Set("ux uy", dom.Nods.Size(),
          0.000000000000000E+00, -5.310671749739340E-06,
          0.000000000000000E+00, -3.243122913926113E-06,
          0.000000000000000E+00, -1.378773528167206E-06,
          0.000000000000000E+00,  0.000000000000000E+00,
         -7.221882229919595E-07, -3.342856970059000E-06,
          3.669995660268174E-07, -2.228571333618901E-06,
          1.912334843273510E-07, -1.114285662470972E-06,
          0.000000000000000E+00,  0.000000000000000E+00,
          0.000000000000000E+00, -1.375042190378655E-06,
          0.000000000000000E+00, -1.214019753311685E-06,
          0.000000000000000E+00, -8.497977967747343E-07,
          0.000000000000000E+00,  0.000000000000000E+00,
          0.000000000000000E+00, -4.288483503711751E-06,
         -4.211153193865630E-07, -5.041231308584792E-06,
          0.000000000000000E+00, -2.217378283206796E-06,
          2.708147610339125E-07, -2.873165860694200E-06,
          0.000000000000000E+00, -6.453528854650004E-07,
          1.370485290061350E-07, -1.298706937331023E-06,
          0.000000000000000E+00,  0.000000000000000E+00,
          3.774202249726477E-07, -2.785714138358060E-06,
         -4.211153193865630E-07, -1.644482664256252E-06,
          2.996313570447692E-07, -1.671428489956403E-06,
          2.708147610339123E-07, -1.583976767620262E-06,
          1.121994376475524E-07, -5.571428285393078E-07,
          1.370485290061355E-07, -9.298643811646859E-07,
          0.000000000000000E+00,  0.000000000000000E+00,
          0.000000000000000E+00, -1.282944773004366E-06,
          0.000000000000000E+00, -1.125478696706005E-06,
          0.000000000000000E+00, -4.689327716136146E-07);

    Table ele_sol;
    ele_sol.Set("sx sy sxy  ex ey exy", dom.Eles.Size(),
        -2.475804412799619E-01, -9.002690947990966E-01,  1.039537575766807E-01,  1.258067454068823E-07, -7.226885041679928E-07,  2.702797696993699E-07/2.0,
        -1.683305992971648E-01, -6.488665829536232E-01,  8.714308276093419E-02,  9.987712199149316E-08, -5.248196567619029E-07,  2.265720151784289E-07/2.0,
        -1.994054385363477E-01, -5.611763146970815E-01,  2.887930674936624E-02,  3.739981366378542E-08, -4.329023253451686E-07,  7.508619754835223E-08/2.0,
        -1.809909839014195E-01, -9.973089729079350E-02,  1.039537575766811E-01, -1.258067454068823E-07, -2.016863281306848E-08,  2.702797696993709E-07/2.0,
        -2.602408192096339E-01, -3.511333935622404E-01,  8.714308276093476E-02, -9.987712199149316E-08, -2.180374686498816E-07,  2.265720151784303E-07/2.0,
        -2.291659816390965E-01, -4.388236657122882E-01,  2.887930674936637E-02, -3.739981366378542E-08, -3.099548029589347E-07,  7.508619754835257E-08/2.0);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy", 1.0e-13, 1.0e-12);
    ele_tol.Set("sx sy sz sxy  ex ey ez exy", 1.0e-8,1.0e-7,1.0e-7,1.0e-8, 1.0e-14,1.0e-13,1.0e-15,1.0e-14);

    // return error flag
    return dom.CheckError (cout, nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
