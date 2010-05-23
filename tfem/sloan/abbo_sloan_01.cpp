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

/*  Abbo & Sloan (1996): Figure 1 p1747 
 *  ===================================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    // input
    bool quad8 = false;
    if (argc>1) quad8 = atoi(argv[1]);
    String fkey("abbo_sloan_01");
    if (quad8) 
    {
        cout << "\n[1;33m=================================== Quad8 ======================================[0m\n";
        fkey.append("_quad8");
    }
    else
    {
        cout << "\n[1;33m=================================== Tri15 ======================================[0m\n";
        fkey.append("_tri15");
    }
    //fkey.append("_ref");

    // constants
    double a  = 1;
    double b  = a*2;
    int    nx = 5;
    double t  = (b-a)/5;

    // mesh
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 0.,  a, 0.0,
                 0.,  b, 0.0,
                 0.,  b,   t,
                 0.,  a,   t,  -10., 0., -20., -30.);
    blks[0].SetNx (nx);
    blks[0].SetNy (1);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/true);
    if (!quad8)
    {
        mesh.Quad8ToTri6 ();
        mesh.Tri6ToTri15 ();
    }
    mesh.Check    ();
    mesh.WriteMPY (fkey.CStr(), /*tags*/true, /*ids*/true, /*shares*/true);

    // parameters
    double E   = 100000.;
    double nu  = 0.3;
    double sY  = 100.;
    double Hp  = 0.0;
    double c   = 10.0;
    double phi = 30.0;
    double del = 0.3*a*1.0e-3;

    // props and domain
    double geom = (quad8 ? GEOM("Quad8") : GEOM("Tri15"));
    double nip  = (quad8 ? 4.0           : 16.0);
    Dict prps, mdls;
    prps.Set(-1, "prob geom axs nip", PROB("Equilib"), geom, 1.0, nip);
    //mdls.Set(-1, "name E nu fc sY Hp axs", MODEL("ElastoPlastic"), E, nu, FAILCRIT("VM"), sY, Hp, 1.0);
    mdls.Set(-1, "name E nu fc c phi axs", MODEL("ElastoPlastic"), E, nu, FAILCRIT("MC"), c, phi, 1.0);
    //mdls.Set(-1, "name E nu axs", MODEL("LinElastic"), E, nu, 1.0);
    FEM::Domain dom(mesh, prps, mdls, /*inis*/Dict());

    // solver
    FEM::Solver sol(dom);
    //sol.SetScheme("NR");
    sol.STOL = 1.0e-4;

    // solve
    if (quad8) dom.SetOutNods (fkey.CStr(), Array<int>(0,6,22));
    else       dom.SetOutNods (fkey.CStr(), Array<int>(0,6,22,44,45));
    Dict bcs;
    bcs.Set      (-10, "uy", 0.0);
    bcs.Set      (-20, "uy", 0.0);
    bcs.Set      (-30, "ux", del);
    dom.SetBCs   (bcs);
    sol.Solve    (10);
    dom.WriteVTU (fkey.CStr());

    //dom.PrintResults();

    // end
    return 0;
}
MECHSYS_CATCH
