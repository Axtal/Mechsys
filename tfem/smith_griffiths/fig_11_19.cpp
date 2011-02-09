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

/*  Smith & Griffiths (2004): Figure 11.19 p500 
 *  ===========================================  */

// STL
#include <iostream>
#include <cmath>    // for cos

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/equilibelem.h>
#include <mechsys/fem/domain.h>
#include <mechsys/fem/solver.h>
#include <mechsys/fem/rksolver.h>
#include <mechsys/models/linelastic.h>
#include <mechsys/models/elastoplastic.h>
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
    // input
    bool rksolver = true;
    if (argc>1) rksolver = atoi(argv[1]);

    // filekey
    String fnkey;
    if (rksolver) fnkey = "fig_11_19";
    else          fnkey = "fig_11_19_GN22";

    ///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
    
    size_t nc = 6;             // number of cells
    size_t nv = 3*(nc+1)+2*nc; // number of vertices
    double H  = 1.0;           // height
    double h  = 0.5;           // half height
    double L  = 15.0;          // length
    double l  = L/nc;          // element length
    double m  = l/2.;          // half length (for mid nodes)
    double x  = 0.0;
    size_t k  = 0;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize (nv/*verts*/, nc/*cells*/);
    mesh.SetVert (k,    0, x, H  );  k++;
    mesh.SetVert (k,    0, x, h  );  k++;
    mesh.SetVert (k, -100, x, 0.0);  k++;
    for (size_t i=0; i<nc; ++i)
    {
        int tag = (i==nc-1 ? -200 : 0);
        mesh.SetVert (k,    0, x+m, H  );  k++;
        mesh.SetVert (k, -200, x+m, 0.0);  k++;
        mesh.SetVert (k,  tag, x+l, H  );  k++;
        mesh.SetVert (k,  tag, x+l, h  );  k++;
        mesh.SetVert (k, -200, x+l, 0.0);  k++;
        mesh.SetCell (i, -1, Array<int>(i*5, i*5+2, i*5+7, i*5+5,  i*5+1, i*5+4, i*5+6, i*5+3));
        mesh.SetBryTag (i, 3, -10);
        x += l;
    }
    mesh.WriteMPY (fnkey.CStr());
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom psa rho", PROB("Equilib"), GEOM("Quad8"), 1.0, 7.33e-4);

    // models
    Dict mdls;
    //mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 3.0e7, 0.3, 1.0);
    mdls.Set(-1, "name E nu sY psa", MODEL("ElastoPlastic"), 3.0e7, 0.3, 5.0e4, 1.0);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    Array<int> out_nods(30, /*justone*/true);
    FEM::Domain dom(mesh, prps, mdls, inis, fnkey.CStr(), &out_nods);

    // boundary conditions
    Dict bcs;
    bcs.Set( -10, "qn",   -180.0);
    bcs.Set(-100, "ux uy", 0.0,0.0);
    bcs.Set(-200, "ux",    0.0);
    dom.SetBCs (bcs);
    //cout << dom << endl;

    // solver
    if (rksolver)
    {
        FEM::RKSolver sol(dom);
        sol.DynSolve (/*tf*/0.01, /*dt*/1.0e-4, /*dtOut*/1.0e-4, fnkey.CStr());
    }
    else
    {
        FEM::Solver sol(dom);
        sol.DScheme = FEM::Solver::GN22_t;
        sol.DynSolve (/*tf*/0.01, /*dt*/1.0e-4, /*dtOut*/1.0e-4, fnkey.CStr());
    }

    return 0.0;
}
MECHSYS_CATCH
