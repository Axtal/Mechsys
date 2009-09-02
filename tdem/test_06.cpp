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
#include "dem/graph.h"
#include "dem/featuredistance.h"
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
	D.AddCube(r,0.7,7.,1.);
    Vec3_t ome(0,M_PI/5,0),vel(1,0,0);
    D.Particles(0)->v() = vel;
    D.Particles(0)->w() = ome;
    r=10,0,0;
    ome=0,0,0;
    vel=-1,0,0;
    D.AddTetra(r,1.,10.,1.);
    D.Particles(1)->v() = vel;
    D.Particles(1)->w() = ome;
    double dt = 0.001;
    D.InitializeSimulation(dt);
    size_t i= 0;
    r=0,0,0; 
    Graph gp("drawing",false);
    gp.DrawEntireDomain(D,"Blue");
    gp.Close();
    Vec3_t L0(0,0,0),P0(0,0,0);
    double E0 = D.TotalEnergy();
    D.LinearMomentum(P0);
    D.AngularMomentum(L0);
    for (double t = 0;t < 30;t+=dt,i++)
    {
        D.OneStep(dt);
        if (i%100==0)
        {
            String         fn;  fn.Printf("mu%.4d",i/100);
	        Graph gp(fn.CStr(),true);
            gp.SetCamera(p,r);
	        gp.DrawEntireDomain(D,"Blue");
	        gp.Close();
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
