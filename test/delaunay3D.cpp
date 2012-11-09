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
    Mesh::Unstructured mesh(3); // 2D
    Array<double> X(6), Y(6), Z(6);
    X = 0.0, 1.0, 1.0, 0.0, 0.5, 2.0;
    Y = 0.0, 0.0, 1.0, 1.0, 0.5, 2.0;
    Z = 1.0, 1.0, 1.0, 1.0, 0.5, 2.0;
    mesh.Delaunay    (X, Y, Z,/*Tag*/-1);
    mesh.WriteVTU    ("delaunay3D");
    cout << " File <delaunay3D.vtu> generated\n";
    return 0;
}
MECHSYS_CATCH
