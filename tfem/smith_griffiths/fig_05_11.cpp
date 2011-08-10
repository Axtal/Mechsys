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
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad4.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solvers/stdsolver.h>
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
    mesh.SetCell   ( 0, -1, Array<int>(0,1, 4, 3));
    mesh.SetCell   ( 1, -1, Array<int>(1,2, 5, 4));
    mesh.SetCell   ( 2, -1, Array<int>(3,4, 7, 6));
    mesh.SetCell   ( 3, -1, Array<int>(4,5, 8, 7));
    mesh.SetCell   ( 4, -1, Array<int>(6,7,10, 9));
    mesh.SetCell   ( 5, -1, Array<int>(7,8,11,10));
    mesh.SetBryTag ( 0, 0, -10);
    mesh.SetBryTag ( 0, 3, -30);
    mesh.SetBryTag ( 1, 0, -10);
    mesh.SetBryTag ( 1, 1, -20);
    mesh.SetBryTag ( 3, 1, -20);
    mesh.SetBryTag ( 4, 2, -10);
    mesh.SetBryTag ( 5, 1, -20);
    mesh.SetBryTag ( 5, 2, -10);
    //mesh.WriteMPY  ("fig_05_11");
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa nip", PROB("Equilib"), GEOM("Quad4"), 1.0, 1.0, 4.0);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa rho", MODEL("LinElastic"),  1.0e+6, 0.3,  1.0, 1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);

    // solver
    SDPair flags;
    FEM::STDSolver sol(dom, flags);

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "ux",     0.0);
    bcs.Set( -20, "ux uy",  0.0);
    bcs.Set( -30, "uy",    -1.0e-5);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.Solve ();

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

    /*
    for (size_t i=0; i<dom.Eles.Size(); ++i)
    {
        for (size_t j=0; j<4; ++j)
        {
            Vec_t X;
            SDPair res;
            dom.Eles[i]->GetState (res, j);
            dom.Eles[i]->CoordsOfIP (j, X);
            cout << Util::_4<<i << Util::_8s<<X(0) << Util::_8s<<X(1) << Util::_8s<<res("sx") << Util::_8s<<res("sy") << Util::_8s<<res("sxy") << endl;
        }
        cout << endl;
    }
    */

    Table ele_sol;
    ele_sol.Set("sx sy sxy", dom.Eles.Size()*4/*nip*/,
        -5.193531824099671E-01, -1.313934421944146E+00, -1.985916864383694E-02,
        -4.026889013713688E-01, -1.263935444356176E+00, -7.450121739119428E-03,
        -5.565803231241198E-01, -1.400797750277169E+00, -8.652447209446459E-02,
        -4.399160420855213E-01, -1.350798772689198E+00, -7.411542518974709E-02,

        -4.129174703940255E-01, -1.283412251283522E+00,  4.266922479239042E-02,
        -5.358782361621740E-01, -1.336109722327015E+00,  3.026017788767373E-02,
        -3.756903296798754E-01, -1.196548922950506E+00,  1.129325195170467E-01,
        -4.986510954480239E-01, -1.249246393993998E+00,  1.005234726123300E-01,

        -4.646381934406058E-01, -1.086904650735379E+00,  2.674326309495864E-01,
        -4.765301479544879E-01, -1.092001202669900E+00,  1.237877978189764E-01,
        -3.370369404877578E-02, -8.139081882110877E-02,  2.742280335289477E-01,
        -4.559564856265794E-02, -8.648737075562969E-02,  1.305832003983376E-01,

        -4.249332749103382E-01, -9.692755076905859E-01,  2.651787051092106E-01,
        -4.163875918454830E-01, -9.656130720913623E-01,  1.588362468890717E-01,
        -1.059059002499217E-01, -2.248783001496140E-01,  2.602954576435790E-01,
        -9.736021718506645E-02, -2.212158645503904E-01,  1.539529994234401E-01,

         4.002167015041122E-02,  1.982430927124814E-01, -1.101592179972651E-01,
        -6.475065637430502E-02,  1.533406670590315E-01, -8.510323363647247E-02,
        -3.514628293196669E-02,  2.285120218693292E-02, -5.028931712599868E-02,
        -1.399186094566829E-01, -2.205122346651690E-02, -2.523333276520604E-02,

        -1.300442355289804E-01, -5.733862970360981E-03,  9.392589573416293E-02,
        -1.562915282568706E-02,  4.330117247390757E-02,  9.054741557179707E-02,
        -1.199087950418828E-01,  1.791549816620003E-02,  2.854584847513820E-02,
        -5.493712338589483E-03,  6.695053361046858E-02,  2.516736831277234E-02);

    // error tolerance
    SDPair nod_tol, ele_tol;
    nod_tol.Set("ux uy",     1.0e-14, 1.0e-13);
    ele_tol.Set("sx sy sxy", 1.0e-7,  1.0e-7,  1.0e-7);

    // return error flag
    bool err_nods = dom.CheckErrorNods (nod_sol, nod_tol);
    bool err_eles = dom.CheckErrorIPs  (ele_sol, ele_tol);
    return (err_nods || err_eles);
}
MECHSYS_CATCH
