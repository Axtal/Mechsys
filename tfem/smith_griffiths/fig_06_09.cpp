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
    /////////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    // coordinates
    Array<double> X(9);
    Array<double> Y(5);
    X = 0.0, 1.0, 2.0, 3.0, 4.0, 5.5, 7.0, 9.0, 12.0;
    Y = 0.0, -1.25, -2.5, -3.75, -5.0;
    double Lx = 2.0; // length of footing

    // generate mesh
    Mesh::Generic mesh(/*NDim*/2);
    mesh.GenGroundSG (X, Y, Lx);

    // write mesh
    mesh.WriteMPY ("fig_06_09");
    
    /////////////////////////////////////////////////////////////////////////////////////////// FEM /////

}
MECHSYS_CATCH
