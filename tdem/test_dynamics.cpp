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

using std::cout;
using std::endl;
using std::map;

int main(int argc, char **argv) try
{
	//This test the Dynamic engine.
    Quaternion_t q;
	Domain D;
	Vec3_t r(-10,0,0),p(0,50,0);
	D.AddCube(r,0.3,3.,1.);
    Vec3_t ome(0,M_PI/50,0),vel(1.,0,0);
    D.Particles[0]->v = vel;
    D.Particles[0]->w = ome;
    D.Particles[0]->CalcMassProperties();
    r = 10 , 0 , 0;
    ome = 0 , 0 , 0;
    vel = -1. , 0 ,0;
    D.AddTetra(r,0.8,8.,1.);
    D.Particles[1]->v = vel;
    D.Particles[1]->w = ome;
    D.Particles[1]->CalcMassProperties();
    double dt = 0.001;
    for(size_t i = 0 ;i < 14; ++i)
    {
        r = Vec3_t((-10.*(i+2)), 0 ,0);
	    D.AddCube(r,0.3,3.,1.);
        D.Particles[2*i+2]->x = r;
        D.Particles[2*i+2]->I = D.Particles[0]->I;
        D.Particles[2*i+2]->Q = 1,0,0,0;
        D.Particles[2*i+2]->v = 1.,0,0;
        D.Particles[2*i+2]->w = D.Particles[0]->w;
        D.Particles[2*i+2]->m = D.Particles[0]->m;
        D.Particles[2*i+2]->Dmax = D.Particles[0]->Dmax;


        r = Vec3_t(10.*(i+2), 0 ,0);
	    D.AddTetra(r,0.44,4.4,1.);
        D.Particles[2*i+3]->x = r;
        D.Particles[2*i+3]->I = D.Particles[1]->I;
        D.Particles[2*i+3]->Q = 1,0,0,0;
        D.Particles[2*i+3]->v = -1.,0,0;
        D.Particles[2*i+3]->w = D.Particles[1]->w;
        D.Particles[2*i+3]->m = D.Particles[1]->m;
        D.Particles[2*i+3]->Dmax = D.Particles[1]->Dmax;
    }
    D.Initialize(dt);
    r=0,0,0; 
    Vec3_t L0(0,0,0),P0(0,0,0);
    double E0 = D.TotalEnergy();
    D.LinearMomentum(P0);
    D.AngularMomentum(L0);
    size_t i =0;
    for (double t = 0;t < 140;t+=dt,i++)
    {
        D.OneStep(dt);
        if(i%100==0)
        {
	        String         fn;  
            fn.Printf("test_dynamics_%04d.pov",i/100);
            std::ofstream ws(fn.CStr());
            PovHeader(ws);
            PovSetCam(ws,p,OrthoSys::O);
            D.WritePov(ws,"Blue");
            ws.close();
        }
    }
    Vec3_t L1(0,0,0),P1(0,0,0);
    double E1 = D.TotalEnergy();
    D.LinearMomentum(P1);
    D.AngularMomentum(L1);
    double error = fabs(E1-E0)/E0+norm(P1-P0)/norm(P0)+norm(L1-L0)/norm(L0),tol = 0.1;
    if (error>tol) return 1;
    else return 0;
}
MECHSYS_CATCH
