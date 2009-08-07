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
	//This tests the distance between a point and an face in 3D, the face is given by the points a,b,c and the point is d. I and F are auxiliary vectors to define the distance edge and draw it. the edges of the face have been draw to for comprehension.
	Vec3_t a(1,0,0),b(0,0,0),c(0,1,0),d(1,1,1),I,F;
	Vec3_t v[3];
	v[0]=a;
	v[1]=b;
	v[2]=c;
	Face Fa(v,3);
	Distance(d,Fa,I,F);
	Edge D(I,F);
	Graph g("drawing",false);
	I=10,0,0;
	F=0,0,0;
	g.SetCamera(I,F);
	g.DrawPoint(d,0.2,"Blue");
	g.DrawFace(Fa,0.2,"Blue");
	g.DrawEdge(*Fa.Edges(0),0.2,"Blue");
	g.DrawEdge(*Fa.Edges(1),0.2,"Blue");
	g.DrawEdge(*Fa.Edges(2),0.2,"Blue");
	g.DrawEdge(D,0.05,"Blue");
	g.Close();
}
MECHSYS_CATCH
