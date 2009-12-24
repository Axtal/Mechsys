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

/*  Smith & Griffiths (2004): Figure 5.17 p182 
 *  ==========================================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad4.h>
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
    mesh.ReadMesh ("fig_05_17");
    mesh.WriteMPY ("fig_05_17");
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active axs nip", PROB("Equilib"), GEOM("Quad4"), TRUE, TRUE, 9.0);
    prps.Set(-2, "prob geom active axs nip", PROB("Equilib"), GEOM("Quad4"), TRUE, TRUE, 9.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu axs", MODEL("LinElastic"),  100.0, 0.3,  TRUE);
    mdls.Set(-2, "name E nu axs", MODEL("LinElastic"), 1000.0, 0.45, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);
    inis.Set(-2, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    FEM::Solver sol(dom);
    sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set(-10, "qn",   -1.0);
    bcs.Set(-11, "ux",    0.0);
    bcs.Set(-12, "ux uy", 0.0,0.0);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.Solve ();

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.PrintResults ("%11.6g");
    dom.WriteVTU ("fig_05_17");

    //////////////////////////////////////////////////////////////////////////////////////// Check /////
    
    Table nod_sol;
    nod_sol.Set("                   ux                      uy", /*NRows*/dom.Nods.Size(),
                 0.000000000000000E+00, -3.176060471584131E-02,
                 0.000000000000000E+00, -3.231272276712264E-03,
                 0.000000000000000E+00,  0.000000000000000E+00,
                 1.394996413258417E-03, -3.990499561989564E-02,
                 1.164812317196636E-03, -2.497995580677533E-03,
                 0.000000000000000E+00,  0.000000000000000E+00,
                 1.703556484195066E-03, -6.045921606453851E-03,
                 1.330158368942836E-03, -4.421423576906878E-04,
                 0.000000000000000E+00,  0.000000000000000E+00,
                 0.000000000000000E+00,  2.587551628166839E-03,
                 0.000000000000000E+00,  3.090608328409015E-04,
                 0.000000000000000E+00,  0.000000000000000E+00);

    Table ele_sol;
    ele_sol.Set("sx sy sz sxy  ex ey ez exy", /*NRows*/dom.Eles.Size(),
                -4.139685339675263E-01, -1.072585276359855E+00, -4.139685339675263E-01, -3.452370246133565E-02,  3.199760913068816E-04, -8.242041559793394E-03,  3.199760913068816E-04, -8.976162639947269E-04/2.0,
                -4.775586993482320E-01, -9.072418317580158E-01, -4.775586993482319E-01,  6.507837344720818E-02,  1.456015396495795E-04, -4.774390023446070E-04,  1.456015396495795E-04,  1.887272829969037E-04/2.0,
                -2.933325170823173E-01, -7.099355623949982E-01, -2.810035539017913E-01,  1.180137872556022E-01,  3.949217806719598E-05, -5.376347410997659E-03,  1.997686994140341E-04,  3.068358468645657E-03/2.0,
                -4.315798060253994E-01, -6.100559027302879E-01, -3.796299730065065E-01,  1.307708058836621E-01,  1.377883805615801E-05, -2.450115021659299E-04,  8.910609593355258E-05,  3.792353370626199E-04/2.0,
                -3.200228788954709E-02, -5.814214821424391E-02, -2.325118726471492E-02,  1.081891071274616E-02, -7.584287245859439E-05, -4.156610566796532E-04,  3.792143566422378E-05,  2.812916785314001E-04/2.0,
                -1.089523748330700E-01, -9.366697323461043E-02, -7.455172701867076E-02,  4.469883244546728E-02, -3.325395971909350E-05, -1.109012740132706E-05,  1.662697961178545E-05,  1.296266140918551E-04/2.0);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy", 1.0e-8, 1.0e-6);
    ele_tol.Set("sx sy sz sxy  ex ey ez exy", 1.0e-5,1.0e-5,1.0e-5,1.0e-5, 1.0e-8,1.0e-7,1.0e-8,1.0e-7);

    // return error flag
    return dom.CheckError (nod_sol, ele_sol, nod_tol, ele_tol);
}
MECHSYS_CATCH
