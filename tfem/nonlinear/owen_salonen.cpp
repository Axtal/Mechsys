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

/*  Owen & Salonen (1975): Fig 5 pag 216                      *
 *  Three-dimensional elasto-plastic finite element analysis  *
 *  Int J Numerical Methods Engineering, Vol 9, 209-218       *
 *  ========================================================  */

// STL
#include <iostream>

// MechSys
#include <mechsys/fem/fem.h>

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;
const double TRUE = 1.0;

int main(int argc, char **argv) try
{
    bool is2d    = true;
    bool unitest = false; // uniaxial test ?
    bool NR      = true;

    cout << "\nUsage:\n\t" << argv[0] << "  [0,1](is2d)  [0,1](unitest)\n" << endl;
    if (argc>1) is2d    = atoi(argv[1]);
    if (argc>2) unitest = atoi(argv[2]);
    if (argc>3) NR      = atoi(argv[3]);
    if (unitest) is2d = false;

    double r  = 2.5;
    double R  = 3.5;
    double t  = 0.5;

    // mesh
	bool o2 = true;
    Mesh::Structured mesh((is2d?2:3));
    if (unitest)
    {
        mesh.GenBox   (o2, 1,1,1);
        mesh.WriteVTU ("owen_salonen_uni_mesh", /*shell*/1);
    }
    else
    {
        if (is2d)
        {
            int nx = 2;
            int ny = 9;
            mesh.GenQRing (o2, nx, ny, r, R, /*Nb*/1, /*Ax*/0.0);
            mesh.WriteMPY ("owen_salonen_mesh");
        }
        else
        {
            int nx = 2;
            int ny = 2;
            int nz = 9;
            mesh.GenQRing (o2, nx, ny, nz, r, R, t);
            mesh.WriteVTU ("owen_salonen_mesh", /*shell*/1);
        }
    }

    // parameters
    double nu = 0.3237;
    
    // calibrated
    //double E  = 12462596.911;
    //double sY = 27373.6759576;
    //double Hp = 408698.716781;

    // from reference
    double E  = 12.458e+6;
    double sY = 27.4e+3;
    double Hp = 4.58e+5;

    // props, domain, and solver
    String fkey;
    Array<int> out_verts;
    Dict prps, mdls;
    if (unitest)
    {
        out_verts.Resize(1);
        fkey      = "owen_salonen_uni";
        out_verts = 0;
        prps.Set(-1, "prob geom", PROB("Equilib"), (o2 ? GEOM("Hex20") : GEOM("Hex8")));
        mdls.Set(-1, "name E nu VM sY Hp", MODEL("ElastoPlastic"), E, nu, TRUE, sY, Hp);
        //mdls.Set(-1, "name E nu", MODEL("LinElastic"), E, nu);
    }
    else
    {
        if (is2d)
        {
            out_verts.Resize(6);
            fkey      = "owen_salonen_2d";
            out_verts = 0,1,2, 27,28,29;
            prps.Set(-1, "prob geom pse h", PROB("Equilib"), (o2 ? GEOM("Quad8") : GEOM("Quad4")), TRUE, t);
            mdls.Set(-1, "name E nu VM sY Hp pse", MODEL("ElastoPlastic"), E, nu, TRUE, sY, Hp, TRUE);
            //mdls.Set(-1, "name E nu pse", MODEL("LinElastic"), E, nu, 1.0);
        }
        else
        {
            out_verts.Resize(8);
            fkey      = "owen_salonen_3d";
            out_verts = -1,-2,-3,-4,-5,-6,-7,-8;
            prps.Set(-1, "prob geom", PROB("Equilib"), (o2 ? GEOM("Hex20") : GEOM("Hex8")));
            mdls.Set(-1, "name E nu VM sY Hp", MODEL("ElastoPlastic"), E, nu, TRUE, sY, Hp);
            //mdls.Set(-1, "name E nu", MODEL("LinElastic"), E, nu);
        }
    }
    FEM::Domain dom(mesh, prps, mdls, /*inis*/Dict(), fkey.CStr(), &out_verts);
    FEM::Solver sol(dom);
    if (NR) sol.SetScheme ("NR");
    sol.SSOut = true;

    // solve
    if (unitest)
    {
        Dict bcs;
        bcs.Set      (-10, "ux", 0.0);
        bcs.Set      (-20, "uy", 0.0);
        bcs.Set      (-30, "uz", 0.0);
        bcs.Set      (-31, "uz", 0.012);
        dom.SetBCs   (bcs);
        sol.Solve    (10);
        dom.WriteVTU ("owen_salonen_uni");
    }
    else
    {
        double DelP = -18000.0;
        if (is2d)
        {
            Dict bcs;
            bcs.Set      (-10, "uy", 0.0);
            bcs.Set      (-30, "ux", 0.0);
            bcs.Set      (-40, "qn", DelP);
            dom.SetBCs   (bcs);
            sol.Solve    (10);
            dom.WriteVTU ("owen_salonen_2d");
        }
        else
        {
            Dict bcs;
            bcs.Set      (-10, "ux", 0.0);
            bcs.Set      (-50, "uz", 0.0);
            bcs.Set      (-60, "uy", 0.0);
            bcs.Set      (-30, "qn", DelP);
            dom.SetBCs   (bcs);
            sol.Solve    (10);
            dom.WriteVTU ("owen_salonen_3d");
        }
    }

    // end
	return 0;
}
MECHSYS_CATCH
