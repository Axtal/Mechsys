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
#include <cmath>

// MechSys
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
const double TRUE = 1.0;

int main(int argc, char **argv) try
{
    /////////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    // generate mesh
    Mesh::Generic mesh(/*NDim*/2);
    mesh.GenGroundSG (/*Nx*/5, /*Ny*/5);

    // layer 1
    mesh.Cells[ 8]->Tag = -2;
    mesh.Cells[12]->Tag = -2;

    // layer 2
    mesh.Cells[ 9]->Tag = -3;
    mesh.Cells[13]->Tag = -3;

    // output
    mesh.WriteMPY ("fig_06_38");
    
    /////////////////////////////////////////////////////////////////////////////////////////// FEM /////

    // data
    double gra = 9.81;    // m/s2
    double gam = 20.0;    // kN/m3
    double K0  = 1.0;

    // elements properties
    Dict prps;
    prps.Set (-1, "prob geom psa  geosta K0 surf", PROB("Equilib"),GEOM("Quad8"),TRUE, TRUE,K0,0.0);
    prps.Set (-2, "prob geom psa  geosta K0 surf", PROB("Equilib"),GEOM("Quad8"),TRUE, TRUE,K0,0.0);
    prps.Set (-3, "prob geom psa  geosta K0 surf", PROB("Equilib"),GEOM("Quad8"),TRUE, TRUE,K0,0.0);

    // parameters
    Dict prms;
    prms.Set (-1, "name psa  E nu  MC c phi", MODEL("LinElastic"),TRUE,  1.0e+5, 0.49, TRUE, 9.0, 0.0);
    prms.Set (-2, "name psa  E nu  MC c phi", MODEL("LinElastic"),TRUE,  1.0e+5, 0.49, TRUE, 9.0, 0.0);
    prms.Set (-3, "name psa  E nu  MC c phi", MODEL("LinElastic"),TRUE,  1.0e+5, 0.49, TRUE, 9.0, 0.0);

    // constants
    prms.Set (-1, "gamW grav gamNat", 9.81, 9.81, gam);
    prms.Set (-2, "gamW grav gamNat", 9.81, 9.81, gam);
    prms.Set (-3, "gamW grav gamNat", 9.81, 9.81, gam);

    // domain
    FEM::Domain dom(mesh, prps, prms, Dict());
    dom.WriteVTU   ("fig_06_38_stg_0");

    // solver
    FEM::Solver sol(dom);
    sol.SetScheme ("FE");

    // stage # 1 ==============================================================
    Dict bcs;
    bcs.Set        (-10, "ux",     0.0);
    bcs.Set        (-20, "ux",     0.0);
    bcs.Set        (-30, "ux uy",  0.0, 0.0);
    bcs.Set        (-2,  "deactivate", gra);
    dom.SetBCs     (bcs);
    sol.Solve      (1);
    dom.WriteVTU   ("fig_06_38_stg_1");

    // stage # 2 ==============================================================
    bcs.clear();
    bcs.Set        (-10, "ux",     0.0);
    bcs.Set        (-20, "ux",     0.0);
    bcs.Set        (-30, "ux uy",  0.0, 0.0);
    bcs.Set        (-3,  "deactivate", gra);
    dom.SetBCs     (bcs);
    sol.Solve      (1);
    dom.WriteVTU   ("fig_06_38_stg_2");
}
MECHSYS_CATCH
