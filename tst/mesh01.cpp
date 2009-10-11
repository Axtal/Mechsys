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
    // 2D: structured
    {
        Array<Mesh::Block> blks(2);
        blks[0].Set (/*NDim*/2, /*Tag*/-1, /*NVert*/4,
                     -1.0,  0.0, 0.0,
                     -2.0,  1.0, 0.0,
                     -3.0,  1.0, 1.0,
                     -4.0,  0.0, 1.0,  -10.0,-20.0,-30.0,-40.0);
        blks[1].Set (/*NDim*/2, /*Tag*/-2, /*NVert*/4,
                     -5.0,  1.0, 0.0,
                     -6.0,  2.0, 0.0,
                     -7.0,  2.0, 1.0,
                     -8.0,  1.0, 1.0,  -11.0,-22.0,-33.0,-44.0);
        blks[0].SetNx (2);
        blks[0].SetNy (3);
        blks[1].SetNx (4);
        blks[1].SetNy (3);
        Mesh::Structured mesh(/*NDim*/2);
        mesh.Generate (blks,/*O2*/true);
        mesh.WriteVTU ("mesh01_quad");
        cout << " File <mesh01_quad.vtu> generated\n";
    }

    // 2D: structured
    {
        Mesh::Structured mesh(/*NDim*/2);
        mesh.GenQRing (/*O2*/true,/*Nx*/4,/*Ny*/1,/*r*/100.,/*R*/200.,/*Nb*/6);
        mesh.WriteMPY ("mesh01_quad_ring", /*OnlyMesh*/false);
        mesh.WriteVTU ("mesh01_quad_ring", /*VolSurfOrBoth*/0);
        cout << " File <mesh01_quad_ring.vtu> generated\n";
    }

    // 3D: structured
    {
        Mesh::Structured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/true);
        mesh.WriteVTU ("mesh01_hex_box", /*VolSurfOrBoth*/0);
        cout << " File <mesh01_hex_box.vtu> generated\n";
    }

    // 2D: unstructured
    {
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
        
        Mesh::Unstructured mesh(2);         // 2D
        mesh.Set (8, 8, 1, 1,               // 8 points, 8 segments, 1 region, 1 hole
                  0., -1.0, 0.0, 0.0,       // id, vtag, x, y, <<<<<< points
                  1., -2.0, 1.5, 0.0,       // id, vtag, x, y,
                  2., -3.0, 1.5, 1.5,       // id, vtag, x, y,
                  3., -4.0, 0.0, 1.5,       // id, vtag, x, y,
                  4.,  0.0, 0.5, 0.5,       // id, vtag, x, y,
                  5.,  0.0, 1.0, 0.5,       // id, vtag, x, y,
                  6.,  0.0, 1.0, 1.0,       // id, vtag, x, y,
                  7.,  0.0, 0.5, 1.0,       // id, vtag, x, y,
                      -1.0, 0.2, 0.8, -1.0, //      tag, x, y, max{area} <<<<<<< regions
                            0.7, 0.7);      //           x, y, <<<<<<< holes
        mesh.SetSeg (0, -10,  0, 1);        // id, etag, L, R <<<<<<<<<<<< segments
        mesh.SetSeg (1, -20,  1, 2);        // id, etag, L, R
        mesh.SetSeg (2, -30,  2, 3);        // id, etag, L, R
        mesh.SetSeg (3, -40,  3, 0);        // id, etag, L, R
        mesh.SetSeg (4,   0,  4, 5);        // id, etag, L, R
        mesh.SetSeg (5,   0,  5, 6);        // id, etag, L, R
        mesh.SetSeg (6,   0,  6, 7);        // id, etag, L, R
        mesh.SetSeg (7,   0,  7, 4);        // id, etag, L, R
        mesh.Generate ();
        mesh.WriteVTU ("mesh01_tri", /*VolSurfOrBoth*/0);
        cout << " File <mesh01_tri.vtu> generated\n";
    }

    // 3D: unstructured
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.Set (4,4,1,0,
                  0., -1., 0.0, 0.0, 0.0,
                  1., -2., 1.0, 0.0, 0.0,
                  2., -3., 0.0, 1.0, 0.0,
                  3., -4., 0.0, 0.0, 1.0,
                      -1., 0.1, 0.1, 0.1, -1.0);
        mesh.SetFac (0, -1, 1,  3., 0.,2.,3.);
        mesh.SetFac (1, -2, 1,  3., 0.,3.,1.);
        mesh.SetFac (2, -3, 1,  3., 0.,1.,2.);
        mesh.SetFac (3, -4, 1,  3., 1.,2.,3.);
        mesh.Generate (/*O2*/true);
        mesh.WriteVTU ("mesh01_1tet");
        cout << " File <mesh01_1tet.vtu> generated\n";
    }

    // 3D: unstructured
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/true,/*V*/0.1);
        mesh.WriteVTU ("mesh01_tet_box", /*VolSurfOrBoth*/0);
        cout << " File <mesh01_tet_box.vtu> generated\n";
    }

    // 3D: unstructured
    {
        Mesh::Unstructured mesh(3);                      // 3D
        mesh.Set (16, 12, 1, 1,                          // 18 points, 12 facets, 1 region, 1 hole
                   0., -1.,  0.0, 0.0, 0.0,              // id, vtag, x, y, z, <<<<<< points
                   1., -2.,  1.5, 0.0, 0.0,              // id, vtag, x, y, z,
                   2., -3.,  1.5, 1.5, 0.0,              // id, vtag, x, y, z,
                   3., -4.,  0.0, 1.5, 0.0,              // id, vtag, x, y, z,
                   4.,  0.,  0.0, 0.0, 1.5,              // id, vtag, x, y, z, <<<<<< points
                   5.,  0.,  1.5, 0.0, 1.5,              // id, vtag, x, y, z,
                   6.,  0.,  1.5, 1.5, 1.5,              // id, vtag, x, y, z,
                   7.,  0.,  0.0, 1.5, 1.5,              // id, vtag, x, y, z,
                   8.,  0.,  0.5, 0.5, 0.5,              // id, vtag, x, y, z,
                   9.,  0.,  1.0, 0.5, 0.5,              // id, vtag, x, y, z,
                  10.,  0.,  1.0, 1.0, 0.5,              // id, vtag, x, y, z,
                  11.,  0.,  0.5, 1.0, 0.5,              // id, vtag, x, y, z,
                  12.,  0.,  0.5, 0.5, 1.0,              // id, vtag, x, y, z,
                  13.,  0.,  1.0, 0.5, 1.0,              // id, vtag, x, y, z,
                  14.,  0.,  1.0, 1.0, 1.0,              // id, vtag, x, y, z,
                  15.,  0.,  0.5, 1.0, 1.0,              // id, vtag, x, y, z,
                       -1.,  0.2, 0.2, 0.2, -1.0,        //      tag, x, y, z, max{volume} <<<<<<< regions
                             0.7, 0.7, 0.7);             //           x, y, z, <<<<<<< holes
        mesh.SetFac ( 0, -1, 1,  4.,  0.,3.,7.,4.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 1, -2, 1,  4.,  1.,2.,6.,5.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 2, -3, 1,  4.,  0.,1.,5.,4.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 3, -4, 1,  4.,  2.,3.,7.,6.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 4, -5, 1,  4.,  0.,1.,2.,3.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 5, -6, 1,  4.,  4.,5.,6.,7.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 6,  0, 1,  4.,  8.,11.,15.,12.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 7,  0, 1,  4.,  9.,10.,14.,13.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 8,  0, 1,  4.,  8., 9.,13.,12.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 9,  0, 1,  4., 10.,11.,15.,14.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac (10,  0, 1,  4.,  8., 9.,10.,11.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac (11,  0, 1,  4., 12.,13.,14.,15.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.Generate ();
        mesh.WritePLY ("mesh01_tet_hole");
        mesh.WriteVTU ("mesh01_tet_hole");
        cout << " File <mesh01_tet_hole.vtu> generated\n";
    }

    return 0;
}
MECHSYS_CATCH
