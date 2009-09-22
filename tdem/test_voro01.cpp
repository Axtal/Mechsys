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
#include "dem/domain.h"
#include "util/fatal.h"

using std::cout;
using std::endl;

// Random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main(int argc, char **argv) try
{
    // generate voro cell
    double x,y,z,rsq,r;
    voronoicell v,u;
    v.init(-1,1,-1,1,-1,1);
    for(int i=0;i<20;i++) 
    {
        x=2*rnd()-1;
        y=2*rnd()-1;
        z=2*rnd()-1;
        rsq=x*x+y*y+z*z;
        if(rsq>0.01&&rsq<1) 
        {
            r=1/sqrt(rsq);x*=r;y*=r;z*=r;
            v.plane(x,y,z,1);
        }
    }

    // generate voro cell
    u.init(-1,1,-1,1,-1,1);
    for(int i=0;i<20;i++) 
    {
        x=2*rnd()-1;
        y=2*rnd()-1;
        z=2*rnd()-1;
        rsq=x*x+y*y+z*z;
        if(rsq>0.01&&rsq<1) 
        {
            r=1/sqrt(rsq);x*=r;y*=r;z*=r;
            u.plane(x,y,z,1);
        }
    }

    // domain
    Domain d;
    d.AddVoroCell (u,0.05);
    d.AddVoroCell (v,0.05);

    // tranlate particles
    Vec3_t trans(2.5,0,0);
    d.Particles[0]->Translate (trans);
    trans = -trans;
    d.Particles[1]->Translate (trans);
    d.Particles[0]->v = Vec3_t(-0.5,0,0);

    d.WriteBPY ("test_voro01");
    // solve
    double dt = 0.001;
    //d.Solve (/*tf*/30, dt, /*dtOut*/0.1, "test_voro01", /*CamPos*/Vec3_t(0,10,0));
}
MECHSYS_CATCH
