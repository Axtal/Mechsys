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
    Array <Particle *> p; // the array of particles at which the force is to be applied
};

void Setup (DEM::Domain & Dom, void * UD)
{
    // force at -3
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dat.p.Size();i++) dat.p[i]->Ff=0.0,0.0,10.0*sin(0.3*Dom.Time);
}


int main(int argc, char **argv) try
{


    // gen beam
    size_t nx = 2;
    size_t ny = 20;
    size_t nz = 2;
    double lx = 2.;
    double ly = 20.;
    double lz = 2.;
    
    // input
    double cam_x=2*ly, cam_y=ly/2.0, cam_z=0;
    Mesh::Structured mesh(3);

    mesh.GenBox (false, nx, ny, nz, lx, ly, lz);


    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha = 0.05;

    // gen particles
    bool cohesion   = true;
    bool montecarlo = false;
    double Kn       = 1.0e4;
    double Kt       = 5.0e3;
    dom.GenFromMesh (mesh, 0.1, 1., cohesion, montecarlo);

    // change tags and identify the particles at the beam's end
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        if      (dom.Particles[i]->x(1)>ly*(1-1.0/ny)) 
        {
            mesh.Cells[i]->Tag    = -3;
            dom.Particles[i]->Tag = -3;
        }
        else if (dom.Particles[i]->x(1)<ly/ny) 
        {
            mesh.Cells[i]->Tag    = -2;
            dom.Particles[i]->Tag = -2;
        }
        else
        {
            mesh.Cells[i]->Tag    = -1;
            dom.Particles[i]->Tag = -1;
        }
        
    }

    // write
    mesh.WriteVTU("test_beam");

    //set the element properties
    dom.Center();
    Dict B;
    B.Set(-1,"Kn Kt Gn Gt Mu Eps",Kn,Kt,8.0,4.0,0.0,-0.01);
    B.Set(-2,"Kn Kt Gn Gt Mu Eps",Kn,Kt,8.0,4.0,0.0,-0.01);
    B.Set(-3,"Kn Kt Gn Gt Mu Eps",Kn,Kt,8.0,4.0,0.0,-0.01);
    dom.SetProps(B);
    dom.CamPos = cam_x, cam_y, cam_z;

     //connect
    dom.GetParticles (-3,dat.p);

    // fix -2 particles at the left extreme of the beam
    Array <Particle *> p;
    dom.GetParticles (-2,p);
    for (size_t i=0;i<p.Size();i++) p[i]->FixVeloc();

    double tf        = 100.;
    double dt        = 1.0e-4;
    double dtOut     = 0.5;
    dom.Solve (tf,dt,dtOut, &Setup, NULL, "test_beam",true);
}
MECHSYS_CATCH
