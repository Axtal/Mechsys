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

int main(int argc, char **argv) try
{
	//This test the varius Domain and Particle constructors.
	Domain D;
	//D.GenerateSpheres(1000,0,10,0,10,0,10,1,0.5);
	Vec3_t r(0,0,0),p(0,0,0);
	D.AddTetra(r,1.,10,1.);
	cout << "inside \? " << D.Particles(0)->IsInside(p) <<endl;
	Graph gp("drawing",false);
	gp.DrawEntireDomain(D,"Blue");
	gp.Close();
}
MECHSYS_CATCH
