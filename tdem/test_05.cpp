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
//#include <gsl/gsl_linalg.h>

// MechSys
#include "dem/graph.h"
#include "dem/featuredistance.h"
#include "util/fatal.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	//This test the varius Domain and Particle constructors.
	Domain D;

	//D.GenerateSpheres(1000,0,10,0,10,0,10,1,0.5);
	Vec3_t r(0,0,0),p(0,0,0);
	D.AddRice(r,1.,10.,1.);
	D.Particles(0)->CalcMassProperties(50000);
    double error1 = fabs(D.Particles(0)->Volume()-(4./3.)*M_PI-M_PI*10),tol1=0.1;
    double error2 = 0,tol2=1.;
    Vec3_t Ireal((1./3.)*M_PI*100+(1./12.)*M_PI*1000+0.75*M_PI*10+(8./15.)*M_PI,(1./3.)*M_PI*100+(1./12.)*M_PI*1000+0.75*M_PI*10+(8./15.)*M_PI,0.5*M_PI*10+(8./15.)*M_PI);
    for (size_t i = 0;i<3;i++)
    {
        error2+=fabs(D.Particles(0)->I()(i)-Ireal(i));
    }

	cout << "Volume " << D.Particles(0)->Volume() << endl;
	cout << "Center of mass " << D.Particles(0)->r() <<endl;
	cout << "Moment of inertia " << D.Particles(0)->I() <<endl;
	cout << "Quaternion " << D.Particles(0)->Q() << endl;



	Graph gp("drawing",false);
	gp.DrawEntireDomain(D,"Blue");
	gp.Close();

    if ((error1>tol1)||(error2>tol2)) return 1;
    else return 0;
}
MECHSYS_CATCH
