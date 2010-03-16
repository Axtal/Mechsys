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
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;
using std::map;
using DEM::Domain;

int main(int argc, char **argv) try
{
    
    /////////////////////////////////////////////////////////////////////////////////////////// Mesh /////

    Mesh::Generic mesh(/*NDim*/3);
    mesh.SetSize(5,2);
    mesh.SetVert(0,-1,-7.0, 0.0, 3.0);
    mesh.SetVert(1,-1, 4.0, 0.0, 0.0);
    mesh.SetVert(2,-1, 0.0, 0.0, 6.0);
    mesh.SetVert(3,-1, 0.0, 4.0, 0.0);
    mesh.SetVert(4,-1, 0.0,-4.0, 0.0);
    mesh.SetCell(0,-1,Array<int>(3,4,1,2));
    mesh.SetCell(1,-1,Array<int>(3,0,4,2));

    /////////////////////////////////////////////////////////////////////////////////////////// Domain /////
    
    Domain d;
    d.GenFromMesh (-1,mesh,/*R*/0.3,/*rho*/1.0,true,false);
    d.Center(Vec3_t(0.0,0.0,3.5));
    d.AddPlane(-2,OrthoSys::O,0.3,100,100,1.0);
    d.Initialize();
    d.Save("test_tets");
    //////////////////////////////////////////////////////////////////////////////////// First timestep /////
    
    //Fix the plane 
    Dict B;
    B.Set(-2,"vx vy vz",0.0,0.0,0.0);
    d.SetBC(B);

    d.CamPos= 4.0,3.0,3.0;
    d.WritePOV ("test_tets");
    d.WriteBPY ("test_tets");

    // Initialize the gravity on the particles
    for (size_t i=0;i<d.FreeParticles.Size();i++)
    {
        d.FreeParticles[i]->Ff = d.FreeParticles[i]->m*Vec3_t(0.0,0.0,-9.8);
    }

    d.CamPos = Vec3_t(0.0, 20.0, 2.5); // position of camera
    d.ResetInteractons();
    d.Solve     (/*tf*/5.0, /*dt*/0.0005, /*dtOut*/0.1, "test_tets", true);
    return 0;    
}
MECHSYS_CATCH
