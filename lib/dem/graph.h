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
#include "dem/particle.h"

using namespace std;

class Graph
{
public:
	// Constructor
	Graph(char const *f,           ///< Filename
	      bool IsPovray=true);     ///< Flag for the user to choose between Povray or Blender for visualization

	// Methods
	void SetCamera (const Vec3_t & x,const Vec3_t & v);        ///< Put the Camera at position x looking at point v
	void DrawPoint (const Vec3_t & r,double R,char const *c);  ///< Draw a single point as a sphere at position r with radius R and color c
	void DrawEdge ( Edge & E,double R,char const *c);          ///< Draw a edge as a cylinder with radius R and color c
	void DrawFace ( Face & F,double R,char const *c);          ///< Draw a face as a polygonal face with radius R and color c
	void DrawPolygon (const Vec3_t *v,size_t N,char const *c); ///< Draw a polygon with N sides and color c whose vertices are stored in v
	void Close ();                                             ///< Flushes the working string into the final file




protected:
	String _fn;              ///< The exporting File Name
	std::ostringstream _oss; ///< Working String
	bool _IsPovray;		 ///< Povray-Blender Flag
};



/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Graph::Graph(char const *f,bool IsPovray)
{
	_IsPovray=IsPovray;
	if(IsPovray)
	{
		_fn.Printf("%s.pov",f);
		_oss << "#include \"colors.inc\" \n";
		_oss << "background {color White} \n";
		_oss << "light_source{<3,0,0> color White shadowless}  \n";
	}
	else
	{
		_fn.Printf("%s.py",f);
		_oss << "from Blender import * \n";
		_oss << "from bpy import * \n";
		_oss << "from Blender import Mathutils \n";
		_oss << "from Blender.Mathutils import * \n";
		_oss << "s = data.scenes.active \n";
	}
}	

inline void Graph::SetCamera (const Vec3_t & x,const Vec3_t & v)
{
	if(_IsPovray)
	{
		_oss << "camera { location <"<<x(0)<<","<<x(1)<<","<<x(2)<<"> sky <0,0,1> look_at <"<<v(0)<<","<<v(1)<<","<<v(2)<<"> }\n";
	}
}

inline void Graph::DrawPoint (const Vec3_t & r,double R,char const *c)
{
	if(_IsPovray)
	{
		_oss << "sphere  { <"<<r(0)<<","<<r(1)<<","<<r(2)<<">,"<<R<<"\n pigment { color "<<c<<" } }\n";
	}
	else
	{
		_oss << "m = Mesh.Primitives.UVsphere(32,32,"<<R*2.0<<") \no = s.objects.new(m,'Sphere') \no.setLocation("<<r(0)<<","<<r(1)<<","<<r(2)<<") \n";
	}
}

inline void Graph::DrawEdge (Edge & E,double R,char const *c)
{
	Vec3_t ini = E.ri();
	Vec3_t fin = E.rf();
	if(_IsPovray)
	{
		_oss << "cylinder { <"<<ini(0)<<","<<ini(1)<<","<<ini(2)<<">,<"<<fin(0)<<","<<fin(1)<<","<<fin(2)<<">,"<<R<<"\n pigment { color "<<c<<" } }\n";
	}
	else
	{
		Vec3_t middle=(ini+fin)/2;
		double L = norm(fin-ini);
		Vec3_t Bv(0,0,L/2);
		Vec3_t axis;
		axis = fin-middle;
		axis = cross(Bv,axis);
		if(norm(axis)==0) axis=1,0,0;
		double angle = acos(4*dot(Bv,fin-middle)/(L*L));
		_oss << "m = Mesh.Primitives.Cylinder(32,"<<R*2.0<<","<<L<<") \no = s.objects.new(m,'Cylinder') \n";
		_oss << "axis = ["<<axis(0)<<","<<axis(1)<<","<<axis(2)<<"] \nquat = Quaternion(axis,"<<angle*180/M_PI<<") \nquat.normalize() \no.setMatrix(quat.toMatrix()) \no.setLocation("<<middle(0)<<","<<middle(1)<<","<<middle(2)<<") \n";
	}
}

inline void Graph::DrawPolygon (const Vec3_t *v,size_t N,char const *c)
{
	if (_IsPovray)
	{
		Vec3_t middle(0,0,0);
		for(size_t i=0;i<N;i++) 
		{
			middle += v[i];
		}
		middle /= N;
		for(size_t i=0;i<N;i++) {
			_oss << "polygon {"<<3<<", \n";
			_oss << "<"<<v[i](0)<<","<<v[i](1)<<","<<v[i](2)<<">";
			_oss << ",<"<<v[(i+1)%N](0)<<","<<v[(i+1)%N](1)<<","<<v[(i+1)%N](2)<<">";
			_oss << ",<"<<middle(0)<<","<<middle(1)<<","<<middle(2)<<">";
		}
		_oss <<"\n pigment { color "<<c<<" } }\n";
	}
	else
	{
		Vec3_t middle(0,0,0);
		for(size_t i=0;i<N;i++) 
		{
			middle += v[i];
		}
		middle /= N;
		for(size_t i=0;i<N;i++) 
		{	
			_oss << "m = data.meshes.new('Face') \no = s.objects.new(m,'Face') \nv = [["<<v[i](0)<<","<<v[i](1)<<","<<v[i](2)<<"],["<<v[(i+1)%N](0)<<","<<v[(i+1)%N](1)<<","<<v[(i+1)%N](2)<<"],["<<middle(0)<<","<<middle(1)<<","<<middle(2)<<"]] \nf = [[0,1,2]] \nm.verts.extend(v) \nm.faces.extend(f) \n";
		}
		
	}

}

inline void Graph::DrawFace (Face & F,double R,char const *c)
{
	Vec3_t * vs;
	Vec3_t * vi;
	Vec3_t n = cross(F.Edges(0)->dr(),F.Edges(1)->dr());
	n=n/norm(n);
	size_t ns = F.NumberofSides();
	vs = new Vec3_t [ns];
	vi = new Vec3_t [ns];
	for(size_t i=0;i<ns;i++) {
		vi[i]=F.Edges(i%ns)->ri()-R*n;
        	vs[i]=F.Edges(i%ns)->ri()+R*n;
		//cout << vs [i] << vi[i] << endl;
	}
	DrawPolygon(vi,ns,c);
    	DrawPolygon(vs,ns,c);
	delete [] vs;
	delete [] vi;
}



inline void Graph::Close ()
{
	std::ofstream file(_fn.CStr());
	file << _oss.str();
	file.close();
}

#endif //DEM_GRAPH_H

