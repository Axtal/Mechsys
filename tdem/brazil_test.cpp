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

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/mesh/unstructured.h>

using std::cout;
using std::endl;
using std::ofstream;
using DEM::Domain;

int main(int argc, char **argv) try
{
    if (argc!=2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".hdf5");
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    
    Domain d;
    d.CamPos = Vec3_t(10.0, 30.0,5.0); // position of camera
    Mesh::Unstructured mesh(3);                  // 3D

    size_t n_divisions = 30;
    double thickness = 5.0;
    double radius = 10.0;
    mesh.Set    (2*n_divisions, 2+n_divisions, 1, 0);                  // 18 points, 12 facets, 1 region, 1 hole
    mesh.SetReg (0,  -1,  -1,  0.0, 0.0, 0.0);  // id, tag, max{volume}, x, y, z <<<<<<< regions
    Array<int> Front;
    Array<int> Back;
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetPnt(i  , 0, radius*cos(2*i*M_PI/n_divisions),-thickness/2.0,radius*sin(2*i*M_PI/n_divisions));
        mesh.SetPnt(n_divisions+i, 0, radius*cos(2*i*M_PI/n_divisions), thickness/2.0,radius*sin(2*i*M_PI/n_divisions));
    }
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetFac(i, 0, Array<int>(i,(i+1)%n_divisions,(i+1)%n_divisions+n_divisions,i+n_divisions));
        Front.Push(i);
        Back.Push(i+n_divisions);
    }
    mesh.SetFac(n_divisions, 0,Front);
    mesh.SetFac(n_divisions+1, 0,Back);

    mesh.Generate();
    d.GenFromMesh(mesh,0.1,3.0,true,false);
    d.Center();
    Vec3_t Xmin,Xmax;
    d.BoundingBox(Xmin,Xmax);
    d.AddPlane(-2, Vec3_t(0.0,0.0,Xmin(2)-0.1), 0.1, 2.4*radius, 1.2*thickness, 1.0);
    d.AddPlane(-3, Vec3_t(0.0,0.0,Xmax(2)+0.1), 0.1, 2.4*radius, 1.2*thickness, 1.0);

    Vec3_t velocity(0.0,0.0,0.1);

    Particle * p1 = d.GetParticle(-2);
    Particle * p2 = d.GetParticle(-3);
    p1->Fixvelocities();
    p1->v =  velocity;
    p2->Fixvelocities();
    p2->v = -velocity;

    d.WriteBPY(filekey.CStr());

    d.Solve(10.0, 0.001, 0.1, NULL, NULL, filekey.CStr());


    return 0;
}
MECHSYS_CATCH
