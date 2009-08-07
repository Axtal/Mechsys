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
#include "dem/graph.h"
#include "dem/featuredistance.h"
#include "util/fatal.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	//This tests the distance between two edges in 3D, the edges are given by the points a,b and c,d. I and F are auxiliary vectors to define the distance edge and draw it.
	Vec3_t a(1,1,1),b(0,0,0),c(0.5,1,0),d(2,0,0),I,F;
	Edge E1(a,b),E2(c,d);
	Distance(E1,E2,I,F);
	Edge D(I,F);


	Graph g("drawing",false);
	g.DrawEdge(E1,0.2,"Blue");
	g.DrawEdge(E2,0.2,"Blue");
	g.DrawEdge(D,0.1,"Blue");
	g.Close();
}
MECHSYS_CATCH
