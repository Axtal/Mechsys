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
    /////////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    // footing: mesh input
    Array<double> X(9);
    Array<double> Y(5);
    X = 0.0, 1.0, 2.0, 3.0, 4.0, 5.5, 7.0, 9.0, 12.0;
    Y = 0.0, -1.25, -2.5, -3.75, -5.0;
    double L = 2.0; // length of footing

    // constants
    size_t nx = X.Size();
    size_t ny = Y.Size();
    size_t nc = (nx-1)*(ny-1);
    size_t nv = nx*ny + nx*(ny-1) + (nx-1)*ny;

    // allocate mesh
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize (nv/*verts*/, nc/*cells*/);

    // generate mesh
    size_t idx_vert = 0;
    size_t idx_cell = 0;
    for (size_t i=0; i<nx-1; ++i)
    {
        // vertices: first column
        if (i==0)
        {
            for (size_t j=0; j<ny; ++j)
            {
                mesh.SetVert (idx_vert, 0, X[i], Y[j]);
                idx_vert++;
                if (j!=ny-1) // intermediate nodes
                {
                    double dy = Y[j+1] - Y[j];
                    mesh.SetVert (idx_vert, 0, X[i], Y[j]+dy/2.0);
                    idx_vert++;
                }
            }
        }

        // vertices: second column
        double dx = X[i+1] - X[i];
        for (size_t j=0; j<ny; ++j)
        {
            mesh.SetVert (idx_vert, 0, X[i]+dx/2.0, Y[j]);
            idx_vert++;
        }

        // vertices: third column
        for (size_t j=0; j<ny; ++j)
        {
            mesh.SetVert (idx_vert, 0, X[i]+dx, Y[j]);
            idx_vert++;
            if (j!=ny-1)
            {
                double dy = Y[j+1] - Y[j];
                mesh.SetVert (idx_vert, 0, X[i]+dx, Y[j]+dy/2.0);
                idx_vert++;
            }
        }

        // set cells
        for (size_t j=0; j<ny-1; ++j)
        {
            int a = (3*ny-1)*i + 2*j;
            int b = (3*ny-1)*i + (2*ny-1) + j;
            int c = a + 3*ny - 1;
            mesh.SetCell (idx_cell, -1, Array<int>(a+2, c+2, c, a,  b+1, c+1, b, a+1));
            if (i==0)    mesh.SetBryTag (idx_cell, 3, -10);
            if (i==nx-2) mesh.SetBryTag (idx_cell, 1, -20);
            if (j==ny-2) mesh.SetBryTag (idx_cell, 0, -30);
            if (j==0)
            {
                if (X[i]+0.95*dx<=L)
                    mesh.SetBryTag (idx_cell, 2, -40);
            }
            idx_cell++;
        }
    }

    // write mesh
    mesh.WriteMPY ("fig_06_09");
    
    /////////////////////////////////////////////////////////////////////////////////////////// FEM /////

}
MECHSYS_CATCH
