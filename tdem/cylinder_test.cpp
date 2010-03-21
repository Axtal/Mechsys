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
    d.CamPos = Vec3_t(0.0, 20.0, 2.5); // position of camera
    Mesh::Unstructured mesh(3);                  // 3D
    mesh.Set    (16, 10, 1, 1);                  // 18 points, 12 facets, 1 region, 1 hole
    mesh.SetReg (0,  0.,  3.0,  0.2, 0.2, 0.2);  // id, tag, max{volume}, x, y, z <<<<<<< regions
    mesh.SetHol (0,  2.5, 1.5, 2.5);             // id, x, y, z, <<<<<<< holes
    mesh.SetPnt ( 0,  0,  0.0, 0.0, 0.0);        // id, vtag, x, y, z, <<<<<< points
    mesh.SetPnt ( 1,  0,  5.0, 0.0, 0.0);        // id, vtag, x, y, z,
    mesh.SetPnt ( 2,  0,  5.0, 0.0, 5.0);        // id, vtag, x, y, z,
    mesh.SetPnt ( 3,  0,  0.0, 0.0, 5.0);        // id, vtag, x, y, z,
    mesh.SetPnt ( 4,  0,  0.0, 3.0, 0.0);        // id, vtag, x, y, z, <<<<<< points
    mesh.SetPnt ( 5,  0,  5.0, 3.0, 0.0);        // id, vtag, x, y, z,
    mesh.SetPnt ( 6,  0,  5.0, 3.0, 5.0);        // id, vtag, x, y, z,
    mesh.SetPnt ( 7,  0,  0.0, 3.0, 5.0);        // id, vtag, x, y, z,
    mesh.SetPnt ( 8,  0,  2.0, 0.0, 2.0);        // id, vtag, x, y, z,
    mesh.SetPnt ( 9,  0,  3.0, 0.0, 2.0);        // id, vtag, x, y, z,
    mesh.SetPnt (10,  0,  3.0, 0.0, 3.0);        // id, vtag, x, y, z,
    mesh.SetPnt (11,  0,  2.0, 0.0, 3.0);        // id, vtag, x, y, z,
    mesh.SetPnt (12,  0,  2.0, 3.0, 2.0);        // id, vtag, x, y, z,
    mesh.SetPnt (13,  0,  3.0, 3.0, 2.0);        // id, vtag, x, y, z,
    mesh.SetPnt (14,  0,  3.0, 3.0, 3.0);        // id, vtag, x, y, z,
    mesh.SetPnt (15,  0,  2.0, 3.0, 3.0);        // id, vtag, x, y, z,
    mesh.SetFac ( 0,  0, Array<int>( 0, 1, 2, 3), Array<int>( 8, 9,10,11));
    mesh.SetFac ( 1,  0, Array<int>( 4, 5, 6, 7), Array<int>(12,13,14,15));
    mesh.SetFac ( 2,  0, Array<int>( 0, 3, 7, 4));
    mesh.SetFac ( 3,  0, Array<int>( 1, 2, 6, 5));
    mesh.SetFac ( 4,  0, Array<int>( 0, 1, 5, 4));
    mesh.SetFac ( 5,  0, Array<int>( 2, 3, 7, 6));
    mesh.SetFac ( 6,  0, Array<int>( 8,11,15,12));
    mesh.SetFac ( 7,  0, Array<int>( 9,10,14,13));
    mesh.SetFac ( 8,  0, Array<int>( 8, 9,13,12));
    mesh.SetFac ( 9,  0, Array<int>(10,11,15,14));
    mesh.Generate();

    d.GenFromMesh(-1,mesh,0.1,1.0,true,false);
    d.Center(Vec3_t(0.0,1.5,8.0));

    Mesh::Unstructured mesh2(3);
    mesh2.Set    (16, 10, 1, 1);
    mesh2.SetReg (0, 0, 3.0,  0.2, 0.2, 0.2);
    mesh2.SetHol (0, 2.5, 1.5, 2.5);
    mesh2.SetPnt ( 0,  0,  0.0, 0.0, 0.0);
    mesh2.SetPnt ( 1,  0,  5.0, 0.0, 0.0);
    mesh2.SetPnt ( 2,  0,  5.0, 0.0, 5.0);
    mesh2.SetPnt ( 3,  0,  0.0, 0.0, 5.0);
    mesh2.SetPnt ( 4,  0,  0.0, 3.0, 0.0);
    mesh2.SetPnt ( 5,  0,  5.0, 3.0, 0.0);
    mesh2.SetPnt ( 6,  0,  5.0, 3.0, 5.0);
    mesh2.SetPnt ( 7,  0,  0.0, 3.0, 5.0);
    mesh2.SetPnt ( 8,  0,  2.0, 0.0, 2.0);
    mesh2.SetPnt ( 9,  0,  3.0, 0.0, 2.0);
    mesh2.SetPnt (10,  0,  3.0, 0.0, 3.0);
    mesh2.SetPnt (11,  0,  2.0, 0.0, 3.0);
    mesh2.SetPnt (12,  0,  2.0, 3.0, 2.0);
    mesh2.SetPnt (13,  0,  3.0, 3.0, 2.0);
    mesh2.SetPnt (14,  0,  3.0, 3.0, 3.0);
    mesh2.SetPnt (15,  0,  2.0, 3.0, 3.0);
    mesh2.SetFac ( 0,  0, Array<int>( 0, 1, 2, 3), Array<int>( 8, 9,10,11));
    mesh2.SetFac ( 1,  0, Array<int>( 4, 5, 6, 7), Array<int>(12,13,14,15));
    mesh2.SetFac ( 2,  0, Array<int>( 0, 3, 7, 4));
    mesh2.SetFac ( 3,  0, Array<int>( 1, 2, 6, 5));
    mesh2.SetFac ( 4,  0, Array<int>( 0, 1, 5, 4));
    mesh2.SetFac ( 5,  0, Array<int>( 2, 3, 7, 6));
    mesh2.SetFac ( 6,  0, Array<int>( 8,11,15,12));
    mesh2.SetFac ( 7,  0, Array<int>( 9,10,14,13));
    mesh2.SetFac ( 8,  0, Array<int>( 8, 9,13,12));
    mesh2.SetFac ( 9,  0, Array<int>(10,11,15,14));
    mesh2.Generate();
    d.GenFromMesh(-1,mesh2,0.1,1.0,true,false);


    d.AddPlane(-2,Vec3_t(0.0,0.0,-0.2),0.2,100,100,1.0);
    d.WriteBPY(filekey.CStr());
    d.Initialize();
    d.Save(filekey.CStr());
    d.WriteBPY(filekey.CStr());
    Array<double> X,Y,D;
    d.GetGSD(X,Y,D);

    ofstream fg("granulometry.txt");
    for (size_t i = 0;i < X.Size() ; i++ )
    {
        fg << Util::_10_6 << X[i] << Util::_8s << Y[i] << std::endl;
    }
    fg.close();


    //Fix the plane 
    Dict B;
    B.Set(-2,"vx vy vz",0.0,0.0,0.0);
    d.SetBC(B);

    // Initialize the gravity on the particles
    for (size_t i=0;i<d.FreeParticles.Size();i++)
    {
        d.FreeParticles[i]->Ff = d.FreeParticles[i]->m*Vec3_t(0.0,0.0,-9.8);
    }

    d.Solve     (/*tf*/10.0, /*dt*/0.00005, /*dtOut*/0.1, filekey.CStr(), true);

    return 0;
}
MECHSYS_CATCH
