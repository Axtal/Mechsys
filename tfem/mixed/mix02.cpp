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
#include <mechsys/fem/usigepselem.h>
#include <mechsys/models/nlelastic.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using FEM::Domain;
using FEM::Solver;
using Util::_8s;

int main(int argc, char **argv) try
{
    // input
    bool   mixed    = true;
    bool   nonlin   = true;
    int    ninc     = 10;
    bool   FE       = false;
    bool   NR       = false;
    bool   usigeps  = true;
    String geom_sig = "Quad8";
    double a        = 1.0;
    if (argc>1) mixed    = atoi(argv[1]);
    if (argc>2) nonlin   = atoi(argv[2]);
    if (argc>3) ninc     = atoi(argv[3]);
    if (argc>4) FE       = atoi(argv[4]);
    if (argc>5) NR       = atoi(argv[5]);
    if (argc>6) usigeps  = atoi(argv[6]);
    if (argc>7) geom_sig =      argv[7];
    if (argc>8) a        = atof(argv[8]);

    // mesh
    double L = 10.0;
    double l = 2.0;
    double m = L/2.0+a/2.0;
    double n = L-m;
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (13/*verts*/, 2/*cells*/);
    mesh.SetVert   ( 0, -100, 0.0,    0.0,  0.0);
    mesh.SetVert   ( 1,    0, m/2.,   0.0,  0.0);
    mesh.SetVert   ( 2,    0, m,      0.0,  0.0);
    mesh.SetVert   ( 3,    0, m+n/2., 0.0,  0.0);
    mesh.SetVert   ( 4, -300, L,      0.0,  0.0);
    mesh.SetVert   ( 5, -200, 0.0,    l/2., 0.0);
    mesh.SetVert   ( 6,    0, L/2.,   l/2., 0.0);
    mesh.SetVert   ( 7,    0, L,      l/2., 0.0);
    mesh.SetVert   ( 8, -200, 0.0,    l,    0.0);
    mesh.SetVert   ( 9,    0, n/2.,   l,    0.0);
    mesh.SetVert   (10,    0, n,      l,    0.0);
    mesh.SetVert   (11,    0, n+m/2., l,    0.0);
    mesh.SetVert   (12, -400, L,      l,    0.0);
    mesh.SetCell   (0, -1, Array<int>(0,2,10,8,  1,6,9,5));
    mesh.SetCell   (1, -1, Array<int>(2,4,12,10, 3,7,11,6));
    mesh.SetBryTag (0, 3, -40);
    mesh.SetBryTag (1, 1, -20);
    mesh.WriteMPY ("mix02");

    // elements properties
    double prob = (mixed ? (usigeps ? PROB("USigEps") : PROB("USig")) : PROB("Equilib"));
    Dict prps;
    prps.Set (-1, "prob geom psa nip geom_sig", prob, GEOM("Quad8"), 1., 4.0, GEOM(geom_sig));

    // model parameters
    Dict mdls;
    double E  = 75.0;
    double nu = 0.0;
    double K  = E/(3.0*(1.0-2.0*nu));
    double G  = E/(2.0*(1.0+nu));
    if (nonlin) mdls.Set (-1, "name K0 G0 alp bet psa", MODEL("NLElastic"),  K, G, 2.0, 1.0, 1.);
    else        mdls.Set (-1, "name E nu psa",          MODEL("LinElastic"), E, nu, 1.);

    // initial values
    Dict inis;
    //inis.Set (-1, "sx sy sz", -10.0, -10.0, -10.0);

    // domain
    Array<int> out_verts(4,7,12);
    Domain dom(mesh, prps, mdls, inis, "mix02", &out_verts);

    // solver
    SDPair flags;
    flags.Set ("stol tolr", 1.0e-7, 1.0e-10);
    if (FE) flags.Set ("fe", 1.);
    if (NR) flags.Set ("nr", 1.);
    FEM::STDSolver sol(dom, flags);

    /*
    sol.Initialize ();
    sol.AssembleKA ();
    sol.A11.WriteSMAT ("mix02");
    cout << "DetA11 = " << UMFPACK::Det(sol.A11) << endl;
    */

    // solve stage # 1
    Dict bcs;
    bcs.Set      (-100, "ux uy", 0.0, 0.0);
    bcs.Set      (-200, "ux",    0.0);
    bcs.Set      (-300, "ux",    0.5);
    bcs.Set      (-400, "ux",   -0.5);
    dom.SetBCs   (bcs);
    sol.Solve    (ninc);
    dom.WriteVTU ("mix02");

    //std::ofstream of("mix02.res", std::ios::out);
    //of << _8s << "a" << _8s<< "v" << endl;

    // end
    return 0;
}
MECHSYS_CATCH
