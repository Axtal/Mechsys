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

    // elements properties
    double gam = 20.0;     // kN/m3
    double rho = gam/9.81; // Mg/m3
    double K0  = 1.0;
    Dict prps;
    prps.Set(-1, "prob geom psa  rho", PROB("Equilib"),GEOM("Quad8"),1.0, rho);
    prps.Set(-2, "prob geom psa  rho", PROB("Equilib"),GEOM("Quad8"),1.0, rho);
    prps.Set(-3, "prob geom psa  rho", PROB("Equilib"),GEOM("Quad8"),1.0, rho);

    // models
    Dict mdls;
    mdls.Set(-1, "name psa  E nu", MODEL("LinElastic"),1.0,  1.0e+5, 0.49);
    mdls.Set(-2, "name psa  E nu", MODEL("LinElastic"),1.0,  1.0e+5, 0.49);
    mdls.Set(-3, "name psa  E nu", MODEL("LinElastic"),1.0,  1.0e+5, 0.49);

    // domain
    Dict inis;
    inis.Set(-1, "geostatic K0 gam y_surf", 1.0, K0, gam, 0.0);
    inis.Set(-2, "geostatic K0 gam y_surf", 1.0, K0, gam, 0.0);
    inis.Set(-3, "geostatic K0 gam y_surf", 1.0, K0, gam, 0.0);
    FEM::Domain dom(mesh, prps, mdls, inis);
    dom.WriteVTU  ("fig_06_38_stg_0");


    // solver
    FEM::Solver sol(dom);
    //sol.Scheme = FEM::Solver::FE_t;

    // stage # 1 ==============================================================

    // boundary conditions
    Dict bcs;
    bcs.Set(-10, "ux",     0.0);
    bcs.Set(-20, "ux",     0.0);
    bcs.Set(-30, "ux uy",  0.0, 0.0);
    dom.SetBCs (bcs);
    dom.Deactivate (-2);

    // solve
    sol.Solve (1);

    // output
    dom.WriteVTU ("fig_06_38_stg_1");
}
MECHSYS_CATCH
