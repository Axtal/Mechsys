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

// Std Lib
#include <iostream>

// MechSys
#include <mechsys/util/tree.h>

using std::cout;
using std::endl;

int main(int argc, char **argv)
{
    // dot cluster01_before.dot -Tpng -o cluster01_before.png
    // dot cluster01_after.dot  -Tpng -o cluster01_after.png

    Array<int> edges(18);
    edges = 0,1, 0,5, 3,5, 1,2, 2,5, 3,4, 2,3, 1,3, 8,3;
    Util::Tree tree(edges);

    tree.WriteDOT ("cluster01_before");
    cout << "before: Number of edges = " << tree.nEdges() << endl;
    cout << tree << endl;

    tree.DelEdge(3,5);
    tree.DelEdge(3,2);
    tree.DelEdge(3,1);

    tree.WriteDOT ("cluster01_after");
    cout << endl;
    cout << "after: Number of edges = " << tree.nEdges() << endl;
    cout << tree << endl;

    Array< Array<int> > C;
    tree.GetClusters(C);

    cout << endl;
    for (size_t i=0; i<C.Size(); ++i)
    for (size_t j=0; j<C[i].Size(); ++j)
        cout << "Cluster # " << i << " contains vertex => "  << C[i][j] << endl;

    return 0;
}
