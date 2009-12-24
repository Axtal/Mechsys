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

/*  Smith & Griffiths (2004): Figure 9.8 p421 
 *  =========================================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/mesh/mesh.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/fem/elems/quad8.h>
#include <mechsys/fem/hydromechelem.h>
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
    
    size_t nc = 4;             // number of cells
    size_t nv = 3*(nc+1)+2*nc; // number of vertices
    double H  = 1.0;           // height
    double W  = 0.25;          // width
    double h  = H/nc;          // element height
    double m  = h/2.0;         // half element height
    double w  = W/2.0;         // half element width
    double y  = H;
    size_t k  = 0;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize (nv/*verts*/, nc/*cells*/);
    mesh.SetVert (k, 0, 0.0, y);  k++;
    mesh.SetVert (k, 0,   w, y);  k++;
    mesh.SetVert (k, 0,   W, y);  k++;
    for (size_t i=0; i<nc; ++i)
    {
        int a = i*5;
        mesh.SetVert (k, 0, 0.0, y-m);  k++;
        mesh.SetVert (k, 0,   W, y-m);  k++;
        mesh.SetVert (k, 0, 0.0, y-h);  k++;
        mesh.SetVert (k, 0,   w, y-h);  k++;
        mesh.SetVert (k, 0,   W, y-h);  k++;
        mesh.SetCell (i, -1, Array<int>(a+5, a+7, a+2, a,  a+6, a+4, a+1, a+3));
        if (i==0)    mesh.SetBryTag (i, 2, -10);
        if (i==nc-1) mesh.SetBryTag (i, 0, -20);
        mesh.SetBryTag (i, 1, -30);
        mesh.SetBryTag (i, 3, -30);
        y -= h;
    }
    mesh.WriteMPY ("fig_09_08");
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom active psa nip", PROB("HMEquilib"), GEOM("Quad8"), TRUE, TRUE, 4.0);

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
    sol.DampTy = FEM::Solver::Rayleigh_t;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "qn pw", -1.0, 0.0);
    bcs.Set( -20, "ux uy",  0.0);
    bcs.Set( -30, "ux",     0.0);
    dom.SetBCs (bcs);
    //cout << dom << endl;
    sol.DynSolve (3.0, 0.01, 0.1);

    //////////////////////////////////////////////////////////////////////////////////////// Output ////

    dom.WriteVTU ("fig_09_08");

    //////////////////////////////////////////////////////////////////////////////////////// Check /////
    
    return 0;
}
MECHSYS_CATCH
