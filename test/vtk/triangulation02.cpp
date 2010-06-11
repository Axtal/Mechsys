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
#include <mechsys/vtk/triangulate.h>
#include <mechsys/util/meshgrid.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    // Create a meshgrid
    MeshGrid mg(0.0,1.0,6,  // X
                0.0,1.0,6,  // Y
                0.0,1.0,6); // Z

    VTK::Triangulate (mg.X, mg.Y, mg.Z, mg.Length(), "triangulation02.vtk");

    cout << "Triangulate done. File <triangulation02.vtk> saved\n";

    return 0;
}
MECHSYS_CATCH
