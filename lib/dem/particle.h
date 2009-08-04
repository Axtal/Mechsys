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

#ifndef DEM_PARTICLE_H
#define DEM_PARTICLE_H

// Std lib
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/face.h"
#include "util/array.h"

class Particle
{
public:
	// Constructor
	Particle(const Array<Vec3_t *> &,          ///< The list of vertices
	         const Array<Array <int> > &,      ///< The list of edges with connectivity
	         const Array<Array <int> > &,      ///< The list of faces with connectivity
	         const double R,                   ///< The spheroradius
	         const double rho0,                ///< The density of the material
                 const Vec3_t & v0,                ///< Initial velocity
                 const Vec3_t & w0);               ///< Initial angular velocity

	

	// Methods
	void StartForce(Vec3_t F);          ///< Start the force of the particle a value F, use for external forces like gravity
	void Start(double dt);              ///< Initialize the particle for the Verlet algorithm
	void DynamicRotation(double dt);    ///< Apply rotation on the particle once the total torque is found
	void DynamicTranslation(double dt); ///< Apply translation once the total force is found

	// Access Methods
        size_t NumberVertices ( )         { return _vertex.Size();} ///< Return the number of vertices.
        size_t NumberEdges    ( )         { return _edges.Size();}  ///< Return the number of edges.
        size_t NumberFaces    ( )         { return _faces.Size();}  ///< Return the number of faces.
        double Radius         ( )         { return _R;}             ///< Return the spheroradius
        Vec3_t * Vertex       ( size_t i) { return _vertex[i];}     ///< Return pointer to the i-th vertex
        Edge * Edges          ( size_t i) { return _edges[i];}      ///< Return pointer to the i-th Edge
        Face * Faces          ( size_t i) { return _faces[i];}      ///< Return pointer to the i-th vertex


protected:

	bool            _IsLarge; ///< Flag to see if it is a large particle
	double          _m;       ///< Mass of the particle
	double          _R;       ///< Spheroradius of the particle
	Vec3_t          _r;       ///< Position of the particle
	Vec3_t          _rb;      ///< Former position for the Verlet algorithm
	Vec3_t          _v;       ///< Velocity of the particle
	Vec3_t          _f;       ///< Force over the particle
	Vec3_t          _T;       ///< Torque over the particle
	Vec3_t          _w;       ///< Angular velocity
	Vec3_t          _wb;      ///< Former angular velocity for the leap frog algorithm
	Vec3_t          _I;       ///< Vector containing the principal components of the inertia tensor
	Quaternion_t    _Q;       ///< The quaternion representing the rotation
	Array<Vec3_t *> _vertex;  ///< The set of vertices defining the geometry of the particle
	Array<Edge *>   _edges;   ///< The set of edges defining the geometry of the particle
	Array<Face *>   _faces;   ///< The set of faces defining the geometry of the particle
	
	
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Particle::Particle(const Array<Vec3_t *> & V,const Array<Array <int> > & E,const Array<Array <int> > & F,const double R,const double rho0,const Vec3_t & v0,const Vec3_t & w0)
{
	if (V.Size()==1&&E.Size()==0&&F.Size()==0)
	{
		_R                 = R;
		_m                 = (4./3.)*rho0*_R*_R*_R;
		_r                 = *V[0];
		_w                 = 
		_I                 = Vec3_t(0.4*_m*_R*_R,0.4*_m*_R*_R,0.4*_m*_R*_R);
		_Q                 = 1,0,0,0;
		_vertex.Resize(1);
		_edges.Resize(0);
		_faces.Resize(0);
		_vertex[0]         = new Vec3_t(0,0,0);
		*_vertex[0]        = _r;
	}
}




#endif // DEM_PARTICLE_H
