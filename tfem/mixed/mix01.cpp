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

// STL
#include <iostream>

// MechSys
#include <mechsys/fem/fem.h>
#include <mechsys/fem/usigelem.h>
#include <mechsys/models/nlelastic.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using FEM::Domain;
using FEM::Solver;

int main(int argc, char **argv) try
{
    // input
    bool mixed  = true;
    bool nonlin = true;
    int  ninc   = 10;
    bool FE     = false;
    bool NR     = false;
    if (argc>1) mixed  = atoi(argv[1]);
    if (argc>2) nonlin = atoi(argv[2]);
    if (argc>3) ninc   = atoi(argv[3]);
    if (argc>4) FE     = atoi(argv[4]);
    if (argc>5) NR     = atoi(argv[5]);

    // mesh
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -1.0,  0.0, 0.0,
                 -2.0,  2.0, 0.0,
                 -3.0,  2.0, 1.0,
                 -4.0,  0.0, 1.0,  -10.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (2);
    blks[0].SetNy (1);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/true);
    mesh.WriteMPY ("mix01");

    // elements properties
    Dict prps;
    if (mixed) prps.Set (-1, "prob geom psa nip", PROB("USig"),    GEOM("Quad8"), 1., 4.0);
    else       prps.Set (-1, "prob geom psa nip", PROB("Equilib"), GEOM("Quad8"), 1., 4.0);

    // model parameters
    Dict mdls;
    double E  = 1000.0;
    double nu = 0.2;
    double K  = E/(3.0*(1.0-2.0*nu));
    double G  = E/(2.0*(1.0+nu));
    if (nonlin) mdls.Set (-1, "name K0 G0 alp bet psa", MODEL("NLElastic"),  K, G, 1.6, 1.0, 1.);
    else        mdls.Set (-1, "name E nu psa",          MODEL("LinElastic"), E, nu, 1.);

    // initial values
    Dict inis;
    inis.Set (-1, "sx sy", 0.0, 0.0);

    // domain
    Domain dom(mesh, prps, mdls, inis);
    dom.SetOutNods ("mix01", Array<int>(2,5,12));
    dom.SetOutEles ("mix01", Array<int>(0,1));

    //Mat_t A,Q;
    //FEM::USigElem * ele = static_cast<FEM::USigElem*>(dom.Eles[0]);
    //ele->Matrices(A,Q);
    //cout << PrintMatrix(A);
    //return 0;

    // solver
    Solver sol(dom);
    sol.STOL = 1.0e-7;
    sol.TolR = 1.0e-10;
    if (FE) sol.SetScheme ("FE");
    if (NR) sol.SetScheme ("NR");

    sol.Initialize ();
    sol.AssembleKA ();
    sol.A11.WriteSMAT ("mix01");

    // solve stage # 1
    Dict bcs;
    bcs.Set      (-1,  "ux uy", 0.0, 0.0);
    bcs.Set      (-40, "ux",    0.0);
    bcs.Set      (-20, "qn",    10.0);
    dom.SetBCs   (bcs);
    sol.Solve    (ninc);
    //dom.WriteVTU ("mix01");

    cout << endl;
    cout << "Node # 5: U = " << dom.Nods[5]->U << endl;
    cout << endl;
    
    dom.PrintResults();

    return 0;
}
MECHSYS_CATCH
