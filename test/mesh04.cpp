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

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // input
    int  nx     = 3;
    int  ny     = 3;
    int  nparts = 2;
    bool full   = false;
    if (argc>1) nx     = atoi(argv[1]);
    if (argc>2) ny     = atoi(argv[2]);
    if (argc>3) nparts = atoi(argv[3]);
    if (argc>4) full   = atoi(argv[4]);

    // mesh
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -1.0,  0.0, 0.0,
                 -2.0,  1.0, 0.0,
                 -3.0,  1.0, 1.0,
                 -4.0,  0.0, 1.0,  -10.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (nx);
    blks[0].SetNy (ny);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/false);

    // partitions
    mesh.PartDomain (nparts, full);
    //for (size_t i=0; i<mesh.Cells.Size(); ++i) cout << mesh.Cells[i]->DomID << " ";
    //cout << endl;

    // output
    mesh.WriteMPY ("mesh04");
    mesh.WriteVTU ("mesh04");

    return 0;
}
MECHSYS_CATCH
