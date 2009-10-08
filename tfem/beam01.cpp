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
#include "mesh/mesh.h"
#include "fem/equilibelem.h"
#include "fem/beam.h"
#include "fem/rod.h"
#include "fem/domain.h"
#include "fem/solver.h"
#include "models/linelastic.h"
#include "util/maps.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using FEM::PROB;
using FEM::GEOM;

int main(int argc, char **argv) try
{
    // data
    int    tst = 1;
    double L   = 1.0;
    double sf  = 1.0;
    if (argc>1) tst = atoi(argv[1]);
    if (argc>2) L   = atof(argv[2]);

    // elements properties
    Dict prps;
    prps.Set (-1, "prob fra rho E A Izz", PROB("Beam"), TRUE, 1.0, 1.0, 1.0, 1.0);

    // boundary conditions
    Dict bcs;

    Mesh::Generic mesh(/*NDim*/2);
    switch (tst)
    {
        case 1:
        {
            mesh.SetSize (2, 1);
            mesh.SetVert (0, -100, 0.0, 0.0);
            mesh.SetVert (1, -200,   L, 0.0);
            mesh.SetCell (0,    -1, 2, 0,1);

            bcs.Set(-100, "ux uy", 0.0);
            bcs.Set(-200, "uy",    0.0);
            bcs.Set(  -1, "qn",   -1.0);
            break;
        }
        case 2:
        {
            mesh.SetSize (2, 1);
            mesh.SetVert (0, -100, 0.0, 0.0);
            mesh.SetVert (1, -200,   L, 0.0);
            mesh.SetCell (0,    -1, 2, 0,1);

            bcs.Set(-100, "ux uy wz", 0.0);
            bcs.Set(  -1, "qn",      -1.0);
            break;
        }
        case 3:
        {
            mesh.SetSize (3, 2);
            mesh.SetVert (0, -100,  0.0, 0.0);
            mesh.SetVert (1,    0, L/2., 0.0);
            mesh.SetVert (2, -200,    L, 0.0);
            mesh.SetCell (0,   -1, 2, 0,1);
            mesh.SetCell (1,   -1, 2, 1,2);

            bcs.Set(-100, "ux uy wz", 0.0);
            bcs.Set(  -1, "qn",      -1.0);
            break;
        }
        case 4:
        {
            mesh.SetSize (6, 5);
            mesh.SetVert (0, -100,  0.0, 2.5);
            mesh.SetVert (1, -100,  2.5, 0.0);
            mesh.SetVert (2, -100, 10.0, 5.0);
            mesh.SetVert (3,    0,  7.5, 2.5);
            mesh.SetVert (4,    0,  2.5, 2.5);
            mesh.SetVert (5, -200,  5.0, 2.5);
            mesh.SetCell (0,   -3, 2, 0,4);
            mesh.SetCell (1,   -1, 2, 1,4);
            mesh.SetCell (2,   -1, 2, 3,2);
            mesh.SetCell (3,   -2, 2, 3,5);
            mesh.SetCell (4,   -2, 2, 5,4);

            bcs.Set (-100, "ux uy", 0.0);
            bcs.Set (-200, "mz",   12.5);
            bcs.Set (  -2, "qn",   -9.2);

            prps.Set (-2, "prob fra rho E A Izz", PROB("Beam"), TRUE, 1.0, 1.0, 1.0, 1.0);
            prps.Set (-3, "prob fra E A",         PROB("Rod"),  TRUE, 1.0, 1.0);

            sf = 0.005;
            break;
        }
        case 5:
        {
            mesh.SetSize  (6, 5);
            mesh.SetVert  (0, -100, 0.0, 0.0);
            mesh.SetVert  (1,    0, 0.0,   L);
            mesh.SetVert  (2,    0,   L, 2*L);
            mesh.SetVert  (3,    0, 2*L, 2*L);
            mesh.SetVert  (4,    0, 3*L,   L);
            mesh.SetVert  (5, -100, 3*L, 0.0);
            mesh.SetCell  (0,   -1, 2, 0,1);
            mesh.SetCell  (1,   -1, 2, 1,2);
            mesh.SetCell  (2,   -2, 2, 2,3);
            mesh.SetCell  (3,   -1, 2, 3,4);
            mesh.SetCell  (4,   -1, 2, 4,5);
            //mesh.SetCell  (0,   -1, 2, 1,0);
            //mesh.SetCell  (1,   -1, 2, 2,1);
            //mesh.SetCell  (2,   -2, 2, 3,2);
            //mesh.SetCell  (3,   -1, 2, 3,4);
            //mesh.SetCell  (4,   -1, 2, 5,4);
            
            bcs.Set (-100, "ux uy wz", 0.0);
            bcs.Set (  -2, "qn",      -1.0);

            prps.Set (-2, "prob fra rho E A Izz", PROB("Beam"), TRUE, 1.0, 1.0, 1.0, 1.0);

            sf = 0.5;
            break;
        }
        default: throw new Fatal("main: Test = %d is not available",tst);
    }
    mesh.WriteMPY ("beam01_mesh",/*OnlyMesh*/false);

    // domain
    FEM::Domain dom(mesh, prps, Dict(), Dict());

    // solver
    FEM::Solver sol(dom);
    //sol.Scheme = FEM::Solver::FE_t;
    sol.Scheme = FEM::Solver::NR_t;

    // run
    dom.SetBCs (bcs);
    sol.Solve  (/*NDiv*/1);

    // output
    dom.WriteMPY ("beam01_res", sf);

    return 0;
}
MECHSYS_CATCH
