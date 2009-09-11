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

    Array<Mesh::Block> blks(1);
    blks[0].Set (/*NDim*/3, /*Tag*/-1, /*NVert*/8,
                 0.0,  0.0, 0.0, 0.0,  // tag, x, y, z
                 0.0,  1.0, 0.0, 0.0, 
                 0.0,  1.0, 1.0, 0.0, 
                 0.0,  0.0, 1.0, 0.0,
                 0.0,  0.0, 0.0, 1.0,  // tag, x, y, z
                 0.0,  1.0, 0.0, 1.0, 
                 0.0,  1.0, 1.0, 1.0, 
                 0.0,  0.0, 1.0, 1.0,
                 0.0,0.0,0.0,0.0,0.0,0.0); // face tags
    blks[0].SetNx (3);
    blks[0].SetNy (3);
    blks[0].SetNz (3);
    Mesh::Structured mesh(/*NDim*/3);
    mesh.Generate (blks);
    mesh.WriteVTU ("test_27cubes");

    /////////////////////////////////////////////////////////////////////////////////////////// Domain /////
    
    Domain d;
    d.GenFromMesh (mesh);

    //////////////////////////////////////////////////////////////////////////////////// First timestep /////
    
    std::ofstream of("test_27cubes.pov",std::ios::out);
    PovHeader     (of);
    PovSetCam     (of,Vec3_t(2,1.5,1.5),OrthoSys::O);
    d.WritePov    (of,"Blue");
    of.close      ();

    /////////////////////////////////////////////////////////////////////////////////////////////// Solve /////

    double dt = 0.001;
    d.Particles[13]->w = Vec3_t(0,.0,.0);
    d.Initialize(dt);
    d.Solve(0,30,dt,1.,"test_27cubes");

    return 0;    
}
MECHSYS_CATCH
