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
#include <mechsys/mesh/structured.h>

using std::cout;
using std::endl;

struct UserData
{
    Particle * p; // the particle at which the force is to be applied
};

void Setup (DEM::Domain & Dom, void * UD)
{
    // force at -3
    //UserData & dat = (*static_cast<UserData *>(UD));
    //dat.p->Ff=0.0,0.0,3.0*cos(0.3*Dom.Time);
}

//void Report (DEM::Domain const & Dom, void * UD)
//{
    //UserData & dat = (*static_cast<UserData *>(UD));
    //cout << "hello\n";
//}

int main(int argc, char **argv) try
{

    // input
    double cam_x=20, cam_y=0, cam_z=0;
    if (argc>1) cam_x = atof(argv[1]);
    if (argc>2) cam_y = atof(argv[2]);
    if (argc>3) cam_z = atof(argv[3]);

    // gen beam
    size_t nx = 1;
    size_t ny = 2;
    size_t nz = 1;
    double lx = 1.;
    double ly = 2.;
    double lz = 1.;
    Mesh::Structured mesh(3);
    mesh.GenBox (false, nx, ny, nz, lx, ly, lz);

    // change tags
    size_t nc = mesh.Cells.Size();
    mesh.Cells[0]->Tag    = -2;
    mesh.Cells[nc-1]->Tag = -3;

    // write
    mesh.WriteVTU("test_beam");

    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha = 0.001;

    // gen particles
    bool cohesion   = true;
    bool montecarlo = false;
    dom.GenFromMesh (mesh, 0.1, 1., cohesion, montecarlo);
    dom.Center();
    dom.Save ("test_beam");
    Dict B;
    B.Set(-2,"Gn Gt Mu eps",0.0,0.0,0.0,0.01);
    B.Set(-3,"Gn Gt Mu eps",0.0,0.0,0.0,0.01);
    dom.SetProps(B);
    dom.CamPos = cam_x, cam_y, cam_z;

    // connect
    dat.p = dom.GetParticle (-3);
    dat.p->w=0.0,0.0,0.0;

    // fix -2
    Particle * p = dom.GetParticle (-2);
    p->FixVeloc();

    double tf        = 20.;
    double dt        = 0.0001;
    double dtOut     = 0.5;
    //char   filekey[] = "test_beam";
    //bool   render    = false;
    dom.Solve (tf,dt,dtOut, &Setup, NULL, "test_beam",true);
}
MECHSYS_CATCH
