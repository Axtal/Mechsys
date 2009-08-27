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

// GSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>


// MechSys
#include "dem/face.h"
#include "dem/featuredistance.h"
#include "util/array.h"
#include "numerical/montecarlo.h"

using Numerical::MonteCarlo;

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
    void CalcMassProperties(size_t NCALCS = 5000);          ///< Calculate the mass, center of mass and moment of inertia
    void StartForce(Vec3_t F);                              ///< Start the force of the particle a value F, use for external forces like gravity
    void Start(double dt);                                  ///< Initialize the particle for the Verlet algorithm
    void DynamicRotation(double dt);                        ///< Apply rotation on the particle once the total torque is found
    void DynamicTranslation(double dt);                     ///< Apply translation once the total force is found
    bool IsInside(Vec3_t & v);                              ///< Enquire if the point v is inside the particle.
    double IsInside(double *v);                             ///< Wrapper to the previous function
    double MaxX();                                          ///< Find Maximun X coordinate
    double MaxY();                                          ///< Find Maximun Y coordinate
    double MaxZ();                                          ///< Find Maximun Y coordinate
    double MinX();                                          ///< Find Minimun X coordinate
    double MinY();                                          ///< Find Minimun Y coordinate
    double MinZ();                                          ///< Find Minimun Y coordinate

    //Integrants for mass properties
    double V(double *r);                ///< Calculate the volume of the sample
    double Xc(double *r);               ///< Calculate the coordinates of the center of mass
    double Yc(double *r);               ///< Calculate the coordinates of the center of mass
    double Zc(double *r);               ///< Calculate the coordinates of the center of mass
    double Ixx(double *r);              ///< Calculate the inertia tensor
    double Iyy(double *r);              ///< Calculate the inertia tensor
    double Izz(double *r);              ///< Calculate the inertia tensor
    double Ixy(double *r);              ///< Calculate the inertia tensor
    double Ixz(double *r);              ///< Calculate the inertia tensor
    double Iyz(double *r);              ///< Calculate the inertia tensor

    // Access Methods
        size_t NumberVertices ( )         { return _vertex.Size();} ///< Return the number of vertices.
        size_t NumberEdges    ( )         { return _edges.Size();}  ///< Return the number of edges.
        size_t NumberFaces    ( )         { return _faces.Size();}  ///< Return the number of faces.
        double Radius         ( )         { return _R;}             ///< Return the spheroradius
        double Volume         ( )         { return _V;}             ///< Return the Volume
        Vec3_t * Vertex       ( size_t i) { return _vertex[i];}     ///< Return pointer to the i-th vertex
    Vec3_t & r            ( )         { return _r;}             ///< Return center of mass
    Vec3_t & I            ( )         { return _I;}             ///< Return the principal moments of inertia
    Quaternion_t & Q      ( )         { return _Q;}             ///< Return the quaternion of the particle
        Edge * Edges          ( size_t i) { return _edges[i];}      ///< Return pointer to the i-th Edge
        Face * Faces          ( size_t i) { return _faces[i];}      ///< Return pointer to the i-th vertex


protected:

    bool            _IsLarge; ///< Flag to see if it is a large particle
    double          _m;       ///< Mass of the particle
    double          _rho;     ///< Density of the particle
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
        _rho = rho0;
    
}

inline void Particle::CalcMassProperties(size_t NCALCS)
{
    if (_vertex.Size()==1&&_edges.Size()==0&&_faces.Size()==0)
    {
        _V = (4./3.)*M_PI*_R*_R*_R;
        _I = Vec3_t((8./15.)*M_PI*pow(_R,5.),(8./15.)*M_PI*pow(_R,5.),(8./15.)*M_PI*pow(_R,5.));
        _r = *_vertex[0];
        _Q = 1,0,0,0;
    }
/*  else if (_vertex.Size()==4&&_edges.Size()==6&&_faces.Size()==4) {
 *      double l = norm(*_vertex[0]-*_vertex[1]);
 *      _V = sqrt(2)*pow(l,3.)/12.;
 *      _I = Vec3_t((16./(135.*sqrt(3)))*pow(l,5.),(16./(135.*sqrt(3)))*pow(l,5.),(16./(135.*sqrt(3)))*pow(l,5.));
 *      _r = 0.25*(*_vertex[0] + *_vertex[1] + *_vertex[2] + *_vertex[3]);
 *      _Q = 1,0,0,0;
 *  }
 */
    else 
    {
        double ri[3] = { MinX() , MinY() , MinZ() };
        double rs[3] = { MaxX() , MaxY() , MaxZ() };
        Mat3_t I;
        MonteCarlo<Particle> *MC;
        MC = new MonteCarlo<Particle> (this,&Particle::V,Numerical::VEGAS,NCALCS);
        _V = MC->Integrate(ri,rs);
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Xc,Numerical::VEGAS,NCALCS);
        _r(0) = MC->Integrate(ri,rs)/_V;
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Yc,Numerical::VEGAS,NCALCS);
        _r(1) = MC->Integrate(ri,rs)/_V;
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Zc,Numerical::VEGAS,NCALCS);
        _r(2) = MC->Integrate(ri,rs)/_V;
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Ixx,Numerical::VEGAS,NCALCS);
        I(0,0) = MC->Integrate(ri,rs);
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Iyy,Numerical::VEGAS,NCALCS);
        I(1,1) = MC->Integrate(ri,rs);
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Izz,Numerical::VEGAS,NCALCS);
        I(2,2) = MC->Integrate(ri,rs);
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Ixy,Numerical::VEGAS,NCALCS);
        I(0,1) = I(1,0) = MC->Integrate(ri,rs);
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Ixz,Numerical::VEGAS,NCALCS);
        I(0,2) = I(2,0) = MC->Integrate(ri,rs);
        delete MC;
        MC = new MonteCarlo<Particle> (this,&Particle::Iyz,Numerical::VEGAS,NCALCS);
        I(1,2) = I(2,1) = MC->Integrate(ri,rs);
        delete MC;
        Vec3_t xp,yp,zp,In(1,0,0),Jn(0,1,0),Kn(0,0,1);
        Eig(I,_I,xp,yp,zp);
        Vec3_t axis = cross(Kn,zp);
        double angle = acos(dot(Jn,zp));
        NormalizeRotation(angle,axis,_Q);
/*      _Q(0) = 0.5*sqrt(1+gsl_matrix_get(evec,0,0)+gsl_matrix_get(evec,1,1)+gsl_matrix_get(evec,2,2));
 *      _Q(1) = (gsl_matrix_get(evec,2,1)-gsl_matrix_get(evec,1,2))/(4*_Q(0));
 *      _Q(2) = (gsl_matrix_get(evec,0,2)-gsl_matrix_get(evec,2,0))/(4*_Q(0));
 *      _Q(3) = (gsl_matrix_get(evec,1,0)-gsl_matrix_get(evec,0,1))/(4*_Q(0));
 */
        
    }

}

inline bool Particle::IsInside (Vec3_t & v)
{
    bool inside = false,insideface = false;
    size_t nv = NumberVertices(),ne = NumberEdges(),nf = NumberFaces();
    for (size_t i = 0; i < nv; i++)
    {
        if (Distance(v,*Vertex(i)) < _R) {
            inside = true;
            return inside;
        }
    }

    for (size_t i = 0; i < ne; i++)
    {
        if (Distance(v,*Edges(i)) < _R) {
            inside = true;
            return inside;
        }
    }
    int numberofintercepts = 0;
    for (size_t i = 0; i < nf; i++)
    {
        if (Distance(v,*Faces(i)) < _R) {
            inside = true;
            return inside;
        }
        Vec3_t D(0,0,1),nor,x,b;
        Mat3_t m;
        for (size_t j = 0;j < 3;j++)
        {
            m(j,0) = Faces(i)->Edges(0)->dr()(j);
            m(j,1) = Faces(i)->Edges(1)->dr()(j);
            m(j,2) = -D(j);
            b(j) = v(j)-Faces(i)->Edges(0)->ri()(j);
        }
        Sol(m,b,x);
        D = v + x(2)*D;
        nor = cross(Faces(i)->Edges(0)->dr(),Faces(i)->Edges(1)->dr());
        insideface = true;
        for(size_t j=0;j < Faces(i)->NumberofSides();j++) 
        {
            Vec3_t tmp = D-Faces(i)->Edges(j)->ri();
            if (dot(cross(Faces(i)->Edges(j)->dr(),tmp),nor)<0) insideface = false;
        }
        if ((insideface)&&(x(2)>0)) 
        {
            numberofintercepts++;
        }
    }
    if (numberofintercepts%2!=0) inside = true;
    return inside;
}

inline double Particle::IsInside (double *r)
{
    Vec3_t p;
    p(0)=r[0];
    p(1)=r[1];
    p(2)=r[2];
    return double(IsInside(p));
}

inline double Particle::MaxX()
{
    double result = (*Vertex(0))(0)+_R;
    for (size_t i = 1; i < NumberVertices(); i++)
    {
        if ((*Vertex(i))(0)+_R > result) result = (*Vertex(i))(0)+_R; 
    }
    return result;
}

inline double Particle::MaxY()
{
    double result = (*Vertex(0))(1)+_R;
    for (size_t i = 1; i < NumberVertices(); i++)
    {
        if ((*Vertex(i))(1)+_R > result) result = (*Vertex(i))(1)+_R; 
    }
    return result;
}

inline double Particle::MaxZ()
{
    double result = (*Vertex(0))(2)+_R;
    for (size_t i = 1; i < NumberVertices(); i++)
    {
        if ((*Vertex(i))(2)+_R > result) result = (*Vertex(i))(2)+_R; 
    }
    return result;
}

inline double Particle::MinX()
{
    double result = (*Vertex(0))(0)-_R;
    for (size_t i = 1; i < NumberVertices(); i++)
    {
        if ((*Vertex(i))(0)-_R < result) result = (*Vertex(i))(0)-_R; 
    }
    return result;
}

inline double Particle::MinY()
{
    double result = (*Vertex(0))(1)-_R;
    for (size_t i = 1; i < NumberVertices(); i++)
    {
        if ((*Vertex(i))(1)-_R < result) result = (*Vertex(i))(1)-_R; 
    }
    return result;
}

inline double Particle::MinZ()
{
    double result = (*Vertex(0))(2)-_R;
    for (size_t i = 1; i < NumberVertices(); i++)
    {
        if ((*Vertex(i))(2)-_R < result) result = (*Vertex(i))(2)-_R; 
    }
    return result;
}

inline double Particle::V (double *r)
{
    return IsInside(r);
}

inline double Particle::Xc (double *r)
{
    return r[0]*IsInside(r);
}

inline double Particle::Yc (double *r)
{
    return r[1]*IsInside(r);
}

inline double Particle::Zc (double *r)
{
    return r[2]*IsInside(r);
}

inline double Particle::Ixx (double *r)
{
    //std::cout << (r[1]-_r(1)) << " " << (r[2]-_r(2)) << std::endl;
    return ((r[1]-_r(1))*(r[1]-_r(1))+(r[2]-_r(2))*(r[2]-_r(2)))*IsInside(r);
}

inline double Particle::Iyy (double *r)
{
    return ((r[0]-_r(0))*(r[0]-_r(0))+(r[2]-_r(2))*(r[2]-_r(2)))*IsInside(r);
}

inline double Particle::Izz (double *r)
{
    return ((r[0]-_r(0))*(r[0]-_r(0))+(r[1]-_r(1))*(r[1]-_r(1)))*IsInside(r);
}

inline double Particle::Ixy (double *r)
{
    return -(r[0]-_r(0))*(r[1]-_r(1))*IsInside(r);
}

inline double Particle::Ixz (double *r)
{
    return -(r[0]-_r(0))*(r[2]-_r(2))*IsInside(r);
}

inline double Particle::Iyz (double *r)
{
    return -(r[1]-_r(1))*(r[2]-_r(2))*IsInside(r);
}

#endif // DEM_PARTICLE_H
