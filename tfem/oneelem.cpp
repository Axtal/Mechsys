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

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
using Util::_6_3;
using Util::_8s;

int main(int argc, char **argv) try
{
    double a   = 0.5;
    double E   = 100000.;
    double nu  = 0.3;
    double sY  = 100.;
    double Hp  = 0.0;
    double c   = 10.0;
    double phi = 30.0;
    double del = 0.3*a*1.0e-2;

    Mesh::Generic mesh(2);
    mesh.SetSize  (4, 1);
    mesh.SetVert  (0, -100, 0.0, 0.0);
    mesh.SetVert  (1, -200,   a, 0.0);
    mesh.SetVert  (2, -300,   a,   a);
    mesh.SetVert  (3, -400, 0.0,   a);
    mesh.SetCell  (0, -1, Array<int>(0,1,2,3));
    mesh.WriteMPY ("oneelem");

    // props and domain
    Array<int> out_verts(0,1,2,3);
    Dict prps, mdls;
    prps.Set(-1, "prob geom psa nip", PROB("Equilib"), GEOM("Quad4"), 1.0, 4.0);
    mdls.Set(-1, "name E nu fc sY Hp psa", MODEL("ElastoPlastic"), E, nu, FAILCRIT("VM"), sY, Hp, 1.0);
    //mdls.Set(-1, "name E nu fc c phi pse", MODEL("ElastoPlastic"), E, nu, FAILCRIT("MC"), c, phi, 1.0);
    //mdls.Set(-1, "name E nu pse", MODEL("LinElastic"), E, nu, 1.0);
    FEM::Domain dom(mesh, prps, mdls, /*inis*/Dict(), "oneelem", &out_verts);

    // solver
    FEM::Solver sol(dom);
    //sol.SetScheme("FE");
    sol.SetScheme("NR");
    //sol.STOL = 1.0e-4;

    // solve
    Dict bcs;
    bcs.Set          (-100, "ux uy", 0.0, 0.0);
    bcs.Set          (-200, "ux",    del);
    bcs.Set          (-300, "ux",    del);
    bcs.Set          (-400, "ux",    0.0);
    dom.SetBCs       (bcs);
    sol.Solve        (100);
    dom.WriteVTU     ("oneelem");
    dom.PrintResults ();
}
MECHSYS_CATCH
