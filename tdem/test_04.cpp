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

// MechSys
//#include "dem/graph.h"
#include "dem/featuredistance.h"
#include "util/fatal.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	Vec3_t a(1,0,0),b(0,0,0),c(0,1,0),d(1,0.5,0),I(0,0,0),F(0,0,1);
	Vec3_t v[4];
	v[0]=a;
	v[1]=b;
	v[2]=c;
	v[3]=d;
	Face Fa(v,4);
	Quaternion_t q;
	NormalizeRotation(M_PI/3,F,q);
	Fa.Rotate(q,I);

    /*
	Graph g("drawing",false);
	g.DrawFace(Fa,0.1,"Blue");
	g.Close();
    */
}
MECHSYS_CATCH
