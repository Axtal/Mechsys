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

/*
class TTT
{
public:
    Particle * LidZ, LidX,LidY;
    void SayHi() { cout << "TTT says hi\n"; }
    double Lx0, Ly0, Lz0;
    double Lx,  Ly,  Lz;
private:
};
*/


/*
void Setup (DEM::Domain const & Dom, void * Data)
{
    //TTT const & ttt = (*static_cast<TTT const *>(Data));
    TTT & ttt = (*static_cast<TTT *>(Data));
    ttt.LidX->x - ttt.LidX->R;
    ttt.SayHi();
}
*/

struct UserData
{
    Particle * p; // the particle at which the force is to be applied
};

void Setup (DEM::Domain const & Dom, void * UD)
{
    // force at -3
    UserData & dat = (*static_cast<UserData *>(UD));
    dat.p->Ff=0.0,0.0,cos(0.3*Dom.Time);
}

//void Report (DEM::Domain const & Dom, void * UD)
//{
    //UserData & dat = (*static_cast<UserData *>(UD));
    //cout << "hello\n";
//}

int main(int argc, char **argv) try
{
    // input
    double cam_x=4, cam_y=5, cam_z=3;
    if (argc>1) cam_x = atof(argv[1]);
    if (argc>2) cam_y = atof(argv[2]);
    if (argc>3) cam_z = atof(argv[3]);

    // gen beam
    size_t nx = 1;
    size_t ny = 4;
    size_t nz = 1;
    double lx = 1.;
    double ly = 4.;
    double lz = 1.;
    Mesh::Structured mesh(3);
    mesh.GenBox (false, nx, ny, nz, lx, ly, lz);

    // change tags
    size_t nc = mesh.Cells.Size();
    mesh.Cells[0]->Tag   = -2;
    mesh.Cells[nc-1]->Tag = -3;

    // write
    mesh.WriteVTU("test_beam");

    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);

    // gen particles
    bool cohesion   = true;
    bool montecarlo = false;
    dom.GenFromMesh (mesh, 0.1, 1., cohesion, montecarlo);
    dom.Save ("test_beam");
    dom.CamPos = cam_x, cam_y, cam_z;

    // connect
    dat.p = dom.GetParticle (-3);

    // fix -2
    Particle * p = dom.GetParticle (-2);
    p->v=0.0,0.0,0.0;  p->vxf=true;  p->vyf=true;  p->vzf=true;

    double tf        = 100.;
    double dt        = 0.001;
    double dtOut     = 0.5;
    //char   filekey[] = "test_beam";
    //bool   render    = false;
    dom.Solve (tf,dt,dtOut, &Setup, NULL, "test_beam");
}
MECHSYS_CATCH
