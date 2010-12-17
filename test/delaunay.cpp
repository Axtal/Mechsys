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
    Mesh::Unstructured mesh(2); // 2D
    Array<double> X(5), Y(5);
    X = 0.0, 1.0, 1.0, 0.0, 0.5;
    Y = 0.0, 0.0, 1.0, 1.0, 0.5;
    mesh.Delaunay    (X, Y, /*Tag*/-1);
    mesh.TagSegsLine (-33, /*x0*/0.0, /*y0*/0.0, /*AlpDeg*/ 45.0);
    mesh.TagSegsLine (-44, /*x0*/0.0, /*y0*/1.0, /*AlpDeg*/-45.0);
    mesh.TagHorzSegs (-55, /*y*/1.0, /*xMin*/0.0, /*xMax*/1.0);
    mesh.GroundTags  (/*L*/-10, /*R*/-20, /*B*/-30);
    mesh.WriteMPY    ("delaunay01");
    cout << " File <delaunay01.mpy> generated\n";
    return 0;
}
MECHSYS_CATCH
