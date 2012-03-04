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

// MechSys
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/structured.h>

using std::cout;
using std::endl;
using Util::PI;
using Util::SQ2;

int main(int argc, char **argv) try {
    Array<Mesh::Block> blks(2);
    blks[0].Set(/*NDim*/2, /*Tag*/-1, /*NVert*/8,
                -1.,  1.0,                 0.0,
                -2.,  2.0,                 0.0,
                -3.,  2.0,                 2.0,
                -4.,  cos(PI/4.),          sin(PI/4.),
                -5.,  1.5,                 0.0,
                -6.,  2.0,                 1.0,
                -7.,  (SQ2+.5)*cos(PI/4.), (SQ2+.5)*sin(PI/4.),
                -8.,  cos(PI/8.),          sin(PI/8.),
                -10., -11., 0., -14.);
    blks[1].Set(/*NDim*/2, /*Tag*/-1, /*NVert*/8,
                -4.,  cos(PI/4.),          sin(PI/4.),
                -3.,  2.0,                 2.0,
                -9.,  0.0,                 2.0,
                -10., 0.0,                 1.0,
                -7.,  (SQ2+.5)*cos(PI/4.), (SQ2+.5)*sin(PI/4.),
                -11., 1.0,                 2.0,
                -12., 0.0,                 1.5,
                -13., cos(3.*PI/8.),       sin(3.*PI/8.),
                0., -12., -13., -14.);
    blks[0].SetNx (10);
    blks[0].SetNy (10);
    blks[1].SetNx (10);
    blks[1].SetNy (10);
    Mesh::Structured mesh(/*NDim*/2);
    mesh.Generate (blks,/*O2*/false);
    mesh.WriteMPY ("doc_msh2", false, false, false);
} MECHSYS_CATCH
