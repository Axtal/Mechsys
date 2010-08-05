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

// ParMETIS
#ifdef USE_PMETIS
extern "C" {
  #include <metis.h>
}
#endif

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    int nx = 3;
    int ny = 2;
    if (argc>1) nx = atoi(argv[1]);
    if (argc>2) ny = atoi(argv[2]);

    Array<int> Xadj, Adjncy;
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                 -1.0,  0.0, 0.0,
                 -2.0,  1.0, 0.0,
                 -3.0,  1.0, 1.0,
                 -4.0,  0.0, 1.0,  -10.0,-20.0,-30.0,-40.0);
    blks[0].SetNx (nx);
    blks[0].SetNy (ny);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate  (blks,/*O2*/false);
    mesh.WriteMPY  ("mesh04");
    mesh.Adjacency (Xadj, Adjncy);
    cout << "Xadj   = " << Xadj   << endl;
    cout << "Adjncy = " << Adjncy << endl;
    cout << " File <mesh04.mpy> generated\n";

#ifdef USE_PMETIS
    int n          = mesh.Cells.Size();
    int wgtflag    = 0; // No weights
    int numflag    = 0; // zero numbering
    int nparts     = 2; // num of partitions
    int options[5] = {0,0,0,0,0};
    int edgecut;
    int part[n];
    METIS_PartGraphKway (&n, Xadj.GetPtr(), Adjncy.GetPtr(), NULL, NULL, &wgtflag, &numflag, &nparts, options, &edgecut, part);
    cout << "edgecut = " << edgecut << endl;
    cout << "part = ";
    for (int i=0; i<n; ++i) cout << part[i] << ", ";
    cout << endl;
#endif

    return 0;
}
MECHSYS_CATCH
