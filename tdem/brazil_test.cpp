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

struct UserData
{
    Particle         * p1;  // Upper plane
    Particle         * p2;  // Lower plane
    Vec3_t          force;  // Force on planes
    std::ofstream  oss_ss;  // File to store the forces
};

void Setup (DEM::Domain const & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    dat.force = 0.5*(dat.p2->F-dat.p1->F);
}

void Report (DEM::Domain const & dom, void *UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_walls.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "fx" << Util::_8s << "fy" << Util::_8s << "fz \n";
    }
    else 
    {
        if (!dom.Finished) 
        {
            dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << fabs(dat.force(0)) << Util::_8s << fabs(dat.force(1)) << Util::_8s << fabs(dat.force(2)) << std::endl;
        }
        else
        {
            dat.oss_ss.close();
        }
    }
}

int main(int argc, char **argv) try
{
    if (argc!=2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    ifstream infile(filename.CStr());
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    UserData dat; 
    Domain d(&dat);
    Mesh::Unstructured mesh(2);
    size_t n_divisions = 30;
    double thickness   = 5.0;
    double radius      = 20.0;
    double Bn          = 1.0e5;
    double Bt          = 3.3e4;
    double Bm          = 3.3e4;
    double Gn          = 16.0;
    double Gt          = 8.0;
    double eps         = 0.01;
    double Amax        = 10.0;

    infile >> n_divisions;        infile.ignore(200,'\n');
    infile >> thickness;          infile.ignore(200,'\n');
    infile >> radius;             infile.ignore(200,'\n');
    infile >> Bn;                 infile.ignore(200,'\n');
    infile >> Bt;                 infile.ignore(200,'\n');
    infile >> Bm;                 infile.ignore(200,'\n');
    infile >> Gn;                 infile.ignore(200,'\n');
    infile >> Gt;                 infile.ignore(200,'\n');
    infile >> eps;                infile.ignore(200,'\n');
    infile >> Amax;               infile.ignore(200,'\n');


    d.CamPos = Vec3_t(0.0, 0.0, 3*radius); // position of camera
    mesh.Set    (n_divisions, n_divisions, 1, 0);  // division, edges, regions and holes
    mesh.SetReg (0,  -1,  Amax,  0.0, 0.0);  // id, tag, max{volume}, x, y, z <<<<<<< regions
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetPnt(i  , 0, radius*cos(2*i*M_PI/n_divisions),radius*sin(2*i*M_PI/n_divisions));
    }
    for(size_t i=0; i<n_divisions; i++)
    {
        mesh.SetSeg(i, 0, i, (i+1)%n_divisions);
    }

    mesh.Generate();
    d.GenFromMesh(mesh,0.1*sqrt(Amax/10),3.0,true,false,thickness);
    d.Center();
    Vec3_t Xmin,Xmax;
    d.BoundingBox(Xmin,Xmax);
    d.AddPlane(-2, Vec3_t(0.0,Xmin(1)-0.5,0.0), 0.5, 0.5*radius, 1.2*thickness, 1.0, M_PI/2.0, &OrthoSys::e0);
    d.AddPlane(-3, Vec3_t(0.0,Xmax(1)+0.5,0.0), 0.5, 0.5*radius, 1.2*thickness, 1.0, M_PI/2.0, &OrthoSys::e0);
    
    // properties of particles prior the brazilian test
    Dict B;
    B.Set(-1,"Bn Bt Bm Gn Gt eps",Bn,Bt,Bm,Gn,Gt,eps);
    d.SetProps(B);

    Vec3_t velocity(0.0,0.1*radius/10.0,0.0);

    Particle * p1 = d.GetParticle(-2);
    Particle * p2 = d.GetParticle(-3);
    p1->FixVeloc();
    p1->v =  velocity;
    p2->FixVeloc();
    p2->v = -velocity;
    dat.p1=p1;
    dat.p2=p2;

    d.WriteBPY(filekey.CStr());

    d.Solve(10.0, 0.0005, 0.1, &Setup, &Report, filekey.CStr());


    return 0;
}
MECHSYS_CATCH
