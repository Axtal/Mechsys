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
#include "mesh/unstructured.h"

using std::cout;
using std::endl;
using std::map;
using DEM::Domain;

int main(int argc, char **argv) try
{
    /////////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Unstructured mesh(3);
    mesh.GenBox   (/*O2*/false,/*A*/-1,/*L*/1.0);
    mesh.WriteVTU ("mesh01_cube", /*VolSurfOrBoth*/1);

    ///////////////////////////////////////////////////////////////////////////////////////// Domain /////
    
    Domain d;
    d.GenFromMesh (-1, mesh,0.05);

    ///////////////////////////////////////////////////////////////////////////////// First timestep /////
    
    d.WriteBPY ("test_tet");

    ////////////////////////////////////////////////////////////////////////////////////////// Solve /////

    double dt = 0.001;
    d.Particles[13]->w = Vec3_t(0,1.,1.);
    d.Solve(/*tf*/1, dt, /*dtOut*/.1, "test_27cubes", /*CamPos*/Vec3_t(0,10,0));

    return 0;    
}
MECHSYS_CATCH
