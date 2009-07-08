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
 * MERCHANTABILITY or finNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef DEM_GRAPH_H
#define DEM_GRAPH_H

// Std lib
#include <math.h>
#include <fstream>
#include <sstream>
#include <iostream>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/face3d.h"

using namespace std;

class Graph
{
public:
	// Constructor
	Graph(char *f);     ///< Filename

	// Set methods
	void SetCamera (const Vec3_t & x,const Vec3_t & v);  ///< Put the Camera at position x looking at point v
	void DrawPoint (const Vec3_t & r,double R,char *c);  ///< Draw a single point as a sphere at position r with radius R and color c
	void DrawEdge3D ( Edge3D & E,double R,char *c);      ///< Draw a edge as a cylinder with radius R and color c
	void DrawFace3D ( Face3D & F,double R,char *c);      ///< Draw a face as a polygonal face with radius R and color c
	void DrawPolygon (const Vec3_t *v,size_t N,char *c); ///< Draw a polygon with N sides and color c whose vertices are stored in v
	void Close ();                                       ///< Flushes the working string into the final file




protected:
	String _fn;              ///< The exporting File Name
	std::ostringstream _oss; ///< Working String

};



/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Graph::Graph(char *f)
{
	_fn.Printf("%s.pov",f);
	_oss << "#include \"colors.inc\" \n";
	_oss << "background {color White} \n";
	_oss << "light_source{<3,0,0> color White shadowless}  \n";
	
}	

inline void Graph::SetCamera (const Vec3_t & x,const Vec3_t & v)
{
	_oss << "camera { location <"<<x(0)<<","<<x(1)<<","<<x(2)<<"> sky <0,0,1> look_at <"<<v(0)<<","<<v(1)<<","<<v(2)<<"> }\n";
}

inline void Graph::DrawPoint (const Vec3_t & r,double R,char *c)
{
	_oss << "sphere  { <"<<r(0)<<","<<r(1)<<","<<r(2)<<">,"<<R<<"\n pigment { color "<<c<<" } }\n";
}

inline void Graph::DrawEdge3D (Edge3D & E,double R,char *c)
{
	Vec3_t ini = E.ri();
	Vec3_t fin = E.rf();
	_oss << "cylinder { <"<<ini(0)<<","<<ini(1)<<","<<ini(2)<<">,<"<<fin(0)<<","<<fin(1)<<","<<fin(2)<<">,"<<R<<"\n pigment { color "<<c<<" } }\n";
}

inline void Graph::DrawPolygon (const Vec3_t *v,size_t N,char *c)
{
	size_t i=0;
	_oss << "polygon {"<<N<<", \n";
	_oss << "<"<<v[i](0)<<","<<v[i](1)<<","<<v[i](2)<<">";
	for(i=1;i<N;i++) {
		_oss << ",<"<<v[i](0)<<","<<v[i](1)<<","<<v[i](2)<<">";
	}
	_oss <<"\n pigment { color "<<c<<" } }\n";
}

inline void Graph::DrawFace3D (Face3D & F,double R,char *c)
{
	Vec3_t * vs;
	Vec3_t * vi;
	Vec3_t n = cross(F.Edge(0)->dr(),F.Edge(1)->dr());
	n=n/norm(n);
	size_t ns = F.NumberofSides();
	vs = new Vec3_t [ns];
	vi = new Vec3_t [ns];
	for(size_t i=0;i<ns;i++) {
		vi[i]=F.Edge(i%ns)->ri()-R*n;
        	vs[i]=F.Edge(i%ns)->ri()+R*n;
		//cout << vs [i] << vi[i] << endl;
	}
	DrawPolygon(vi,ns,c);
    	DrawPolygon(vs,ns,c);
}



inline void Graph::Close ()
{
	std::ofstream file(_fn.CStr());
	file << _oss.str();
	file.close();
}

#endif //DEM_GRAPH_H

