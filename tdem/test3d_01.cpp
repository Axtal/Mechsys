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
#include "util/exception.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	Vec3_t a(0,1,0),b(1,0,0),c(1,1,1),In(5,0,0),Vi(0,0,0),v[3];
	v[0]=Vec3_t(0,0,-1);
	v[1]=Vec3_t(1,0,-1);
	v[2]=Vec3_t(0,1,-1);
	Edge3D E(a,b);
	Face3D F(v,3);
	Distance(c,E,a,b);
	Edge3D E1(a,b);


	Graph g("drawing");
	g.SetCamera(In,Vi);
	g.DrawPoint(c,0.1,"Blue");
	g.DrawEdge3D(E,0.1,"Red");
	g.DrawEdge3D(E1,0.05,"Black");
	g.DrawFace3D(F,0.1,"Green");
	g.DrawPolygon(v,3,"Green");
	g.Close();
	return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
