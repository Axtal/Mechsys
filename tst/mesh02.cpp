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
#include "util/fatal.h"
#include "mesh/structured.h"
#include "mesh/unstructured.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // 2D: structured
    {
        Array<Mesh::Block> blks(1);
        blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                     -1.0,   0.0, 0.0,
                     -2.0,  10.0, 0.0,
                     -3.0,  10.0, 1.0,
                     -4.0,   0.0, 1.0,  -10.0,-20.0,-30.0,-40.0);
        blks[0].SetNx (2);
        blks[0].SetNy (1);
        Mesh::Structured mesh(/*NDim*/2);
        mesh.Generate    (blks,/*O2*/true);
        mesh.AddLinCells (/*NTagsOrPairs*/1, /*WTags*/true, /*Tags*/-10);
        mesh.WriteMPY    ("mesh02",/*OnlyMesh*/false);
        mesh.WriteVTU    ("mesh02");
    }

    return 0;
}
MECHSYS_CATCH
