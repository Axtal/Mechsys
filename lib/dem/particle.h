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
	void CalcMassProperties();          ///< Calculate the mass, center of mass and moment of inertia
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
	double          _V;       ///< Volume of the particle
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
	for (size_t i = 0;i<V.Size();i++)
	{
		_vertex.Push(new Vec3_t(0,0,0));
		*_vertex[i] = *V[i];
	}
	for (size_t i = 0;i<E.Size();i++)
	{
		size_t n = E[i][0],m = E[i][1];
		_edges.Push(new Edge(*_vertex[n],*_vertex[m]));
	}
	for (size_t i = 0;i<F.Size();i++)
	{
		Vec3_t *v;
		v = new Vec3_t [F[i].Size()];
		for (size_t j = 0;j<F[i].Size();j++) 
		{
			v[j] = *_vertex[F[i][j]];
		}
		_faces.Push(new Face(v,F[i].Size()));
		delete [] v;
	}
	_R = R;
	_v = v0;
	_w = w0;
	CalcMassProperties();
	_m = rho0*_V;
	_I = rho0*_I;
	
}

inline void Particle::CalcMassProperties()
{
	if (_vertex.Size()==1&&_edges.Size()==0&&_faces.Size()==0)
	{
		_V = (4./3.)*M_PI*_R*_R*_R;
		_I = Vec3_t((8./15.)*M_PI*pow(_R,5.),(8./15.)*M_PI*pow(_R,5.),(8./15.)*M_PI*pow(_R,5.));
		_r = *_vertex[0];
		_Q = 1,0,0,0;
	}
	else if (_vertex.Size()==4&&_edges.Size()==6&&_faces.Size()==4) {
		double l = norm(*_vertex[0]-*_vertex[1]);
		_V = sqrt(2)*pow(l,3.)/12.;
		_I = Vec3_t((16./(135.*sqrt(3)))*pow(l,5.),(16./(135.*sqrt(3)))*pow(l,5.),(16./(135.*sqrt(3)))*pow(l,5.));
		_r = 0.25*(*_vertex[0] + *_vertex[1] + *_vertex[2] + *_vertex[3]);
		_Q = 1,0,0,0;
	}
	else 
	{
		_V = 1;
		_I = (1,1,1);
		_r = (0,0,0);
		_Q = 1,0,0,0;
	}
 /* :TODO:08/05/2009 04:32:49 PM:: The rest of the mass properties with Monte Carlo integration */

}



#endif // DEM_PARTICLE_H
