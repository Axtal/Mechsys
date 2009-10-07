/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

// Std lib
#include <math.h>

// GSL
#include <gsl/gsl_linalg.h>

// MechSys
#include "dem/domain.h"
#include "dem/distance.h"
#include "util/fatal.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using std::map;

int main(int argc, char **argv) try
{
    
    /////////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    double L = 1.0;
    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/3, /*Tag*/-1, /*NVert*/8,
                 0.0,  0.0, 0.0, 0.0,  // tag, x, y, z
                 0.0,    L, 0.0, 0.0, 
                 0.0,    L,   L, 0.0, 
                 0.0,  0.0,   L, 0.0,
                 0.0,  0.0, 0.0,   L,  // tag, x, y, z
                 0.0,    L, 0.0,   L, 
                 0.0,    L,   L,   L, 
                 0.0,  0.0,   L,   L,
                 0.0,0.0,0.0,0.0,0.0,0.0); // face tags
    blks[0].SetNx (3);
    blks[0].SetNy (3);
    blks[0].SetNz (3);
    Mesh::Structured mesh(/*NDim*/3);
    mesh.Generate (blks);
    mesh.WriteVTU ("test_27cubes");

    /////////////////////////////////////////////////////////////////////////////////////////// Domain /////
    
    Domain d;
    d.GenFromMesh (-1,mesh,/*R*/0.08);

    //////////////////////////////////////////////////////////////////////////////////// First timestep /////
    
    Vec3_t cam_pos(4.0,3.0,3.0);
    d.WritePOV ("init_test_27cubes",cam_pos);
    d.WriteBPY ("init_test_27cubes");

    /////////////////////////////////////////////////////////////////////////////////////////////// Solve /////

    d.Particles[0]->w = Vec3_t(0.01,0.02,0.03);
    //d.Solve(/*tf*/30.0, /*dt*/0.001, /*dtOut*/0.1, "test_27cubes", cam_pos);

    return 0;    
}
MECHSYS_CATCH
