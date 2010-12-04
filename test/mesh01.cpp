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
    //////////////////////////////////////// 2D: structured //////////////////////////////////////////////////
    
    // 1) mesh01_quad
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

    // 2) mesh01_quad_ring
    {
        Mesh::Structured mesh(/*NDim*/2);
        mesh.GenQRing (/*O2*/true,/*Nx*/4,/*Ny*/1,/*r*/100.,/*R*/200.,/*Nb*/6);
        mesh.WriteMPY ("mesh01_quad_ring", /*OnlyMesh*/false);
        mesh.WriteVTU ("mesh01_quad_ring", /*VolSurfOrBoth*/0);
        cout << " File <mesh01_quad_ring.vtu> generated\n";
    }

    // 3) mesh01_qring2d
    {
        bool   o2 = true;
        double r  = 2.5;
        double R  = 3.5;
        int    nx = 2;
        int    ny = 20;
        int    nb = 1;
        double ax = 0.0;
        bool   nl = false;
        Mesh::Structured mesh(2);
        mesh.GenQRing (o2,nx,ny,r,R,nb,ax,nl);
        mesh.WriteMPY ("mesh01_qring2d");
        cout << " File <mesh01_qring2d.mpy> generated\n";
    }

    // 4) mesh01_qdisk2d
    {
        bool   o2  = true;
        double r   = 2.5;
        double R   = 3.5;
        int    nx1 = 10;
        int    nx2 = 5;
        int    ny  = 10;
        Mesh::Structured mesh(2);
        mesh.GenQDisk (o2, nx1, nx2, ny, r, R);
        mesh.WriteMPY ("mesh01_qdisk2d", false, false);
        cout << " File <mesh01_qdisk2d.mpy> generated\n";
    }

    // 5) mesh01_qplatehole
    {
        String extra("\
from msys_fig import *\n\
Arc (0., 0., 2.5, 0., pi/2., zorder=10)\n\
Arc (0., 0., 3.5, 0., pi/2., zorder=10)\n");
        bool   o2  = false;
        double r   = 2.5;
        double R   = 3.5;
        double L   = 5.0;
        int    nx1 = 3;
        int    nx2 = 5;
        int    ny1 = 2;
        int    ny2 = 1;
        Mesh::Structured mesh(2);
        mesh.GenQPlateHole (r,R,L,nx1,nx2,ny1,ny2,o2);
        mesh.WriteMPY ("mesh01_qplatehole", /*tags*/true, /*ids*/false, /*shares*/false, extra.CStr());
        cout << " File <mesh01_qplatehole.mpy> generated\n";
    }

    //////////////////////////////////////// 3D: structured //////////////////////////////////////////////////
    
    // 6) mesh01_blk_box
    {
        /*
                          4----------------7
                        ,'|              ,'|
                      ,'  |            ,'  |
                    ,'    | -60   -10,'    |
                  ,'      |        ,'      |
                5'===============6'        |
                |         |      |    -40  |
                |    -30  |      |         |
                |         0- - - | -  - - -3
                |       ,'       |       ,'
                |     ,' -20     |     ,'
                |   ,'        -50|   ,'
                | ,'             | ,'
                1----------------2'
        */
        double Lx=1. , Ly=1. , Lz=2.;
        size_t Nx=2  , Ny=2  , Nz=3;
        Array<Mesh::Block> blks(1);
        blks[0].Set (/*NDim*/3, /*Tag*/-1, /*NVert*/8,
                     -1.,  0.0, 0.0, 0.0,  // tag, x, y, z
                     -2.,   Lx, 0.0, 0.0, 
                     -3.,   Lx,  Ly, 0.0, 
                     -4.,  0.0,  Ly, 0.0,
                     -5.,  0.0, 0.0,  Lz,  // tag, x, y, z
                     -6.,   Lx, 0.0,  Lz, 
                     -7.,   Lx,  Ly,  Lz, 
                     -8.,  0.0,  Ly,  Lz,
                     -10.,-20.,-30.,-40.,-50.,-60.); // face tags
        blks[0].SetNx (Nx);
        blks[0].SetNy (Ny);
        blks[0].SetNz (Nz);
        Mesh::Structured mesh(/*NDim*/3);
        mesh.Generate (blks,true);
        mesh.WriteVTU ("mesh01_blk_box");
        cout << " File <mesh01_blk_box.vtu> generated\n";
    }

    // 7) mesh01_hex_box
    {
        Mesh::Structured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/true);
        mesh.WriteVTU ("mesh01_hex_box", /*VolSurfOrBoth*/0);
        cout << " File <mesh01_hex_box.vtu> generated\n";
    }

    // 8) mesh01_qring3d
    {
        bool   o2 = true;
        double r  = 2.5;
        double R  = 3.5;
        double t  = 0.5;
        int    nx = 2;
        int    ny = 2;
        int    nz = 9;
        Mesh::Structured mesh(3);
        mesh.GenQRing (o2, nx, ny, nz, r, R, t);
        mesh.WriteVTU ("mesh01_qring3d");
        cout << " File <mesh01_qring3d.vtu> generated\n";
    }

    //////////////////////////////////////// 2D: unstructured ////////////////////////////////////////////////
    
    // 9) mesh01_tri
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
        
        Mesh::Unstructured mesh(/*NDim*/2);
        mesh.Set    (8, 8, 1, 1);            // 8 points, 8 segments, 1 region, 1 hole
        mesh.SetReg (0, -1, -1.0, 0.2, 0.8); // id, tag, max{area}, x, y <<<<<<< regions
        mesh.SetHol (0, 0.7, 0.7);           // id, x, y <<<<<<< holes
        mesh.SetPnt (0, -1, 0.0, 0.0);       // id, vtag, x, y <<<<<< points
        mesh.SetPnt (1, -2, 1.5, 0.0);       // id, vtag, x, y
        mesh.SetPnt (2, -3, 1.5, 1.5);       // id, vtag, x, y
        mesh.SetPnt (3, -4, 0.0, 1.5);       // id, vtag, x, y
        mesh.SetPnt (4,  0, 0.5, 0.5);       // id, vtag, x, y
        mesh.SetPnt (5,  0, 1.0, 0.5);       // id, vtag, x, y
        mesh.SetPnt (6,  0, 1.0, 1.0);       // id, vtag, x, y
        mesh.SetPnt (7,  0, 0.5, 1.0);       // id, vtag, x, y
        mesh.SetSeg (0, -10,  0, 1);         // id, etag, L, R <<<<<<<<<<<< segments
        mesh.SetSeg (1, -20,  1, 2);         // id, etag, L, R
        mesh.SetSeg (2, -30,  2, 3);         // id, etag, L, R
        mesh.SetSeg (3, -40,  3, 0);         // id, etag, L, R
        mesh.SetSeg (4,   0,  4, 5);         // id, etag, L, R
        mesh.SetSeg (5,   0,  5, 6);         // id, etag, L, R
        mesh.SetSeg (6,   0,  6, 7);         // id, etag, L, R
        mesh.SetSeg (7,   0,  7, 4);         // id, etag, L, R
        mesh.Generate ();
        mesh.FindNeigh ();
        mesh.WriteVTU ("mesh01_tri", /*VolSurfOrBoth*/0);
        cout << " File <mesh01_tri.vtu> generated\n";
    }

    //////////////////////////////////////// 3D: unstructured ////////////////////////////////////////////////
    
    // 10) mesh01_1tet
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.Set    (4,4,1,0);
        mesh.SetReg (0, -1, -1.0, 0.1, 0.1, 0.1);
        mesh.SetPnt (0, -1, 0.0, 0.0, 0.0);
        mesh.SetPnt (1, -2, 1.0, 0.0, 0.0);
        mesh.SetPnt (2, -3, 0.0, 1.0, 0.0);
        mesh.SetPnt (3, -4, 0.0, 0.0, 1.0);
        mesh.SetFac (0, -1, Array<int>(0,2,3));
        mesh.SetFac (1, -2, Array<int>(0,3,1));
        mesh.SetFac (2, -3, Array<int>(0,1,2));
        mesh.SetFac (3, -4, Array<int>(1,2,3));
        mesh.Generate (/*O2*/true);
        mesh.WriteVTU ("mesh01_1tet");
        cout << " File <mesh01_1tet.vtu> generated\n";
    }

    // 11) mesh01_tet_box
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/true,/*V*/0.1);
        mesh.WriteVTU ("mesh01_tet_box", /*VolSurfOrBoth*/0);
        cout << " File <mesh01_tet_box.vtu> generated\n";
    }

    // 12) mesh01_tet_hole
    {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.Set    (16, 12, 1, 1);                      // 16 points, 12 facets, 1 region, 1 hole
        mesh.SetReg ( 0, -1, -1.0, 0.2, 0.2, 0.2);       // id, tag, max{volume}, x, y, z <<<<<<< regions
        mesh.SetHol ( 0, 0.7, 0.7, 0.7);                 // id, x, y, z <<<<<<< holes
        mesh.SetPnt ( 0, -1,  0.0, 0.0, 0.0);            // id, vtag, x, y, z <<<<<< points
        mesh.SetPnt ( 1, -2,  1.5, 0.0, 0.0);            // id, vtag, x, y, z
        mesh.SetPnt ( 2, -3,  1.5, 1.5, 0.0);            // id, vtag, x, y, z
        mesh.SetPnt ( 3, -4,  0.0, 1.5, 0.0);            // id, vtag, x, y, z
        mesh.SetPnt ( 4,  0,  0.0, 0.0, 1.5);            // id, vtag, x, y, z <<<<<< points
        mesh.SetPnt ( 5,  0,  1.5, 0.0, 1.5);            // id, vtag, x, y, z
        mesh.SetPnt ( 6,  0,  1.5, 1.5, 1.5);            // id, vtag, x, y, z
        mesh.SetPnt ( 7,  0,  0.0, 1.5, 1.5);            // id, vtag, x, y, z
        mesh.SetPnt ( 8,  0,  0.5, 0.5, 0.5);            // id, vtag, x, y, z
        mesh.SetPnt ( 9,  0,  1.0, 0.5, 0.5);            // id, vtag, x, y, z
        mesh.SetPnt (10,  0,  1.0, 1.0, 0.5);            // id, vtag, x, y, z
        mesh.SetPnt (11,  0,  0.5, 1.0, 0.5);            // id, vtag, x, y, z
        mesh.SetPnt (12,  0,  0.5, 0.5, 1.0);            // id, vtag, x, y, z
        mesh.SetPnt (13,  0,  1.0, 0.5, 1.0);            // id, vtag, x, y, z
        mesh.SetPnt (14,  0,  1.0, 1.0, 1.0);            // id, vtag, x, y, z
        mesh.SetPnt (15,  0,  0.5, 1.0, 1.0);            // id, vtag, x, y, z
        mesh.SetFac ( 0, -1, Array<int>( 0, 3, 7, 4));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 1, -2, Array<int>( 1, 2, 6, 5));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 2, -3, Array<int>( 0, 1, 5, 4));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 3, -4, Array<int>( 2, 3, 7, 6));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 4, -5, Array<int>( 0, 1, 2, 3));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 5, -6, Array<int>( 4, 5, 6, 7));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 6,  0, Array<int>( 8,11,15,12));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 7,  0, Array<int>( 9,10,14,13));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 8,  0, Array<int>( 8, 9,13,12));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 9,  0, Array<int>(10,11,15,14));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac (10,  0, Array<int>( 8, 9,10,11));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac (11,  0, Array<int>(12,13,14,15));   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.Generate ();
        mesh.WritePLY ("mesh01_tet_hole");
        mesh.WriteVTU ("mesh01_tet_hole");
        cout << " File <mesh01_tet_hole.vtu> generated\n";
    }

    return 0;
}
MECHSYS_CATCH
