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
#include "mesh/mesh.h"
#include "mesh/structured.h"
#include "fem/elems/quad8.h"
#include "fem/equilibelem.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "models/elastoplastic.h"
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
    mesh.WriteMPY ("fig_11_19",/*OnlyMesh*/false);
    
    ////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // elements properties
    Dict prps;
    prps.Set(-1, "prob geom psa rho", PROB("Equilib"), GEOM("Quad8"), TRUE, 7.33e-4);

    // models
    Dict mdls;
    mdls.Set(-1, "name E nu psa", MODEL("LinElastic"), 3.0e7, 0.3, TRUE);
    //mdls.Set(-1, "name E nu sY psa", MODEL("ElastoPlastic"), 3.0e7, 0.3, 5.0e4, TRUE);

    // initial values
    Dict inis;
    inis.Set(-1, "sx sy sz sxy", 0.0,0.0,0.0,0.0);

    // domain
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.SetOutNods ("fig_11_19", Array<int>(30,/*JustOne*/true));

    // solver
    FEM::Solver sol(dom);
    //sol.DScheme = FEM::Solver::SS22_t;
    //sol.DScheme = FEM::Solver::GN22_t;
    sol.DScheme = FEM::Solver::SG113_t;    sol.DampTy=FEM::Solver::Rayleigh_t;    sol.DampAm=0.0;    sol.DampAk=0.0;

    // stage # 1 -----------------------------------------------------------
    Dict bcs;
    bcs.Set( -10, "qn", -180.0);
    bcs.Set(-100, "ux uy", 0.0,0.0);
    bcs.Set(-200, "ux",    0.0);
    dom.SetBCs (bcs);
    cout << dom << endl;
    sol.DynSolve (/*tf*/0.004, /*dt*/1.0e-4, /*dtOut*/1.0e-4);
    //sol.DynSolve (/*tf*/0.01, /*dt*/1.0e-6, /*dtOut*/1.0e-4);

    return 0.0;
}
MECHSYS_CATCH
