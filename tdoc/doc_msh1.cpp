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
#include <mechsys/mesh/mesh.h>
#include <mechsys/util/fatal.h>

int main(int argc, char **argv) try {
    Mesh::Generic mesh(/*NDim*/2);
    mesh.SetSize   (6/*verts*/, 4/*cells*/);
    mesh.SetVert   (0, -100, 0.0, 0.0, 0);
    mesh.SetVert   (1, -101, 0.0, 2.0, 0);
    mesh.SetVert   (2,    0, 2.0, 0.0, 0);
    mesh.SetVert   (3,    0, 2.0, 1.5, 0);
    mesh.SetVert   (4,    0, 4.0, 0.0, 0);
    mesh.SetVert   (5,    0, 4.0, 1.0, 0);
    mesh.SetCell   (0,   -1, Array<int>(0,2,3));
    mesh.SetCell   (1,   -1, Array<int>(3,1,0));
    mesh.SetCell   (2,   -1, Array<int>(2,4,5));
    mesh.SetCell   (3,   -1, Array<int>(5,3,2));
    mesh.SetBryTag (1, 0, -10);
    mesh.SetBryTag (3, 0, -10);
    mesh.WriteVTU  ("doc_msh1");
    mesh.WriteMPY  ("doc_msh1");
} MECHSYS_CATCH
