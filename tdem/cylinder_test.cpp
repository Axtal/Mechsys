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
    if (!Util::FileExists(filename))
    {
        Mesh::Unstructured mesh(3);                      // 3D
        mesh.Set (16, 10, 1, 1,                          // 18 points, 12 facets, 1 region, 1 hole
                   0.,  0.,  0.0, 0.0, 0.0,               // id, vtag, x, y, z, <<<<<< points
                   1.,  0.,  5.0, 0.0, 0.0,               // id, vtag, x, y, z,
                   2.,  0.,  5.0, 0.0, 5.0,               // id, vtag, x, y, z,
                   3.,  0.,  0.0, 0.0, 5.0,               // id, vtag, x, y, z,
                   4.,  0.,  0.0, 3.0, 0.0,               // id, vtag, x, y, z, <<<<<< points
                   5.,  0.,  5.0, 3.0, 0.0,               // id, vtag, x, y, z,
                   6.,  0.,  5.0, 3.0, 5.0,               // id, vtag, x, y, z,
                   7.,  0.,  0.0, 3.0, 5.0,               // id, vtag, x, y, z,
                   8.,  0.,  2.0, 0.0, 2.0,               // id, vtag, x, y, z,
                   9.,  0.,  3.0, 0.0, 2.0,               // id, vtag, x, y, z,
                  10.,  0.,  3.0, 0.0, 3.0,               // id, vtag, x, y, z,
                  11.,  0.,  2.0, 0.0, 3.0,               // id, vtag, x, y, z,
                  12.,  0.,  2.0, 3.0, 2.0,               // id, vtag, x, y, z,
                  13.,  0.,  3.0, 3.0, 2.0,               // id, vtag, x, y, z,
                  14.,  0.,  3.0, 3.0, 3.0,               // id, vtag, x, y, z,
                  15.,  0.,  2.0, 3.0, 3.0,               // id, vtag, x, y, z,
                        0.,  0.2, 0.2, 0.2, 3.0,        //      tag, x, y, z, max{volume} <<<<<<< regions
                             2.5, 1.5, 2.5);             //           x, y, z, <<<<<<< holes
        mesh.SetFac ( 0,  0, 2,  4.,  0.,1.,2.,3.,4., 8.,9.,10.,11.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 1,  0, 2,  4.,  4.,5.,6.,7.,4., 12.,13.,14.,15.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 2,  0, 1,  4.,  0.,3.,7.,4.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 3,  0, 1,  4.,  1.,2.,6.,5.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 4,  0, 1,  4.,  0.,1.,5.,4.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 5,  0, 1,  4.,  2.,3.,7.,6.);      // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 6,  0, 1,  4.,  8.,11.,15.,12.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 7,  0, 1,  4.,  9.,10.,14.,13.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 8,  0, 1,  4.,  8., 9.,13.,12.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.SetFac ( 9,  0, 1,  4., 10.,11.,15.,14.);   // id, ftag, npolygons,  npoints, point0,point1,point2,point3
        mesh.Generate();
        //mesh.FindNeigh();
        //cout<<mesh<<endl;

        d.GenFromMesh(-1,mesh,0.1,1.0,true);
        d.Center(Vec3_t(0.0,0.0,5.0));
        d.AddPlane(-2,OrthoSys::O,0.1,100,100,1.0);
        d.WriteBPY(filekey.CStr());
        d.Initialize();
        //d.Save(filekey.CStr());
    }
    else 
    {
        //d.Load(filekey.CStr());
    }
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

    d.Solve     (/*tf*/5.0, /*dt*/0.00005, /*dtOut*/0.05, filekey.CStr(), true);

    return 0;
}
MECHSYS_CATCH
