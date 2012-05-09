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
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    //////////////////////////////////////// 2D: structured //////////////////////////////////////////////////
    
    // 1) mesh03_quad8_to_tri6
    {
        Array<Mesh::Block> blks(1);
        blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                     -1.0,  0.0, 0.0,
                     -2.0,  3.0, 0.0,
                     -3.0,  3.0, 2.0,
                     -4.0,  0.0, 2.0,  -10.0,-20.0,-30.0,-40.0);
        blks[0].SetNx (3);
        blks[0].SetNy (2);
        Mesh::Structured mesh(/*NDim*/2);
        mesh.Generate    (blks,/*O2*/true);
        mesh.Quad8ToTri6 ();
        mesh.FindNeigh   ();
        mesh.Check       ();
        mesh.WriteMPY ("mesh03_quad8_to_tri6", /*tags*/true, /*ids*/true, /*shares*/true);
        cout << " File <mesh03_quad8_to_tri6.mpy> generated\n";
    }

    // 2) mesh03_tri6_to_tri15
    {
        Array<Mesh::Block> blks(1);
        blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                     -1.0,  0.0, 0.0,
                     -2.0,  3.0, 0.0,
                     -3.0,  3.0, 2.0,
                     -4.0,  0.0, 2.0,  -10.0,-20.0,-30.0,-40.0);
        blks[0].SetNx (3);
        blks[0].SetNy (2);
        Mesh::Structured mesh(/*NDim*/2);
        mesh.Generate    (blks,/*O2*/true);
        mesh.Quad8ToTri6 ();
        mesh.Tri6ToTri15 ();
        mesh.Check       ();
        mesh.WriteMPY ("mesh03_tri6_to_tri15", /*tags*/true, /*ids*/true, /*shares*/true);
        cout << " File <mesh03_tri6_to_tri15.mpy> generated\n";
    }

    return 0;
}
MECHSYS_CATCH
