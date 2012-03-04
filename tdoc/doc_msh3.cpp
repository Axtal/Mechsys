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
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try {
    /*           -20
          -4@-----------@-3
            | -1        |
            |   @---@   |
         -30|   | h |   |-20
            |   @---@   |
            |           |
          -1@-----------@-2
                 -10
    */
    Mesh::Unstructured mesh(/*NDim*/2);
    mesh.Set    (8, 8, 1, 1);             // 8 points, 8 segmts, 1 reg, 1 hole
    mesh.SetReg (0, -1, 0.005, 0.2, 0.8); // id, tag, max{area}, x, y <<< regs
    mesh.SetHol (0, 0.7, 0.7);            // id, x, y        <<<<<<<<<<< holes
    mesh.SetPnt (0, -1, 0.0, 0.0);        // id, vtag, x, y  <<<<<<<<<< points
    mesh.SetPnt (1, -2, 1.5, 0.0);        // id, vtag, x, y
    mesh.SetPnt (2, -3, 1.5, 1.5);        // id, vtag, x, y
    mesh.SetPnt (3, -4, 0.0, 1.5);        // id, vtag, x, y
    mesh.SetPnt (4,  0, 0.5, 0.5);        // id, vtag, x, y
    mesh.SetPnt (5,  0, 1.0, 0.5);        // id, vtag, x, y
    mesh.SetPnt (6,  0, 1.0, 1.0);        // id, vtag, x, y
    mesh.SetPnt (7,  0, 0.5, 1.0);        // id, vtag, x, y
    mesh.SetSeg (0, -10,  0, 1);          // id, etag, L, R  <<<<<<<< segments
    mesh.SetSeg (1, -20,  1, 2);          // id, etag, L, R
    mesh.SetSeg (2, -30,  2, 3);          // id, etag, L, R
    mesh.SetSeg (3, -40,  3, 0);          // id, etag, L, R
    mesh.SetSeg (4,   0,  4, 5);          // id, etag, L, R
    mesh.SetSeg (5,   0,  5, 6);          // id, etag, L, R
    mesh.SetSeg (6,   0,  6, 7);          // id, etag, L, R
    mesh.SetSeg (7,   0,  7, 4);          // id, etag, L, R
    mesh.Generate ();
    mesh.WriteMPY ("doc_msh3", false, false, false);
} MECHSYS_CATCH
