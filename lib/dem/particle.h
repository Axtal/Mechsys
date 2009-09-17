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

#ifndef MECHSYS_DEM_PARTICLE_H
#define MECHSYS_DEM_PARTICLE_H

// Std lib
#include <iostream>

// MechSys
#include "dem/face.h"
#include "dem/distance.h"
#include "util/array.h"
#include "numerical/montecarlo.h"

class Particle
{
public:
    // Constructor
    Particle(Array<Vec3_t>       const & V,            ///< The list of vertices
             Array<Array <int> > const & E,            ///< The list of edges with connectivity
             Array<Array <int> > const & F,            ///< The list of faces with connectivity
             Vec3_t              const & v0,           ///< Initial velocity
             Vec3_t              const & w0,           ///< Initial angular velocity
             double                      R,            ///< The spheroradius
             double                      rho=1.0,      ///< The density of the material
             bool                        Open=false);  ///< Ask if the particle should opened or just dilated

    // Destructor
    ~Particle ();


    // Data
    Vec3_t         x;     ///< Position of the center of mass
    Vec3_t         xb;    ///< Former position for the Verlet algorithm
    Vec3_t         v;     ///< Velocity
    Vec3_t         w;     ///< Angular velocity
    Vec3_t         wb;    ///< Former angular velocity for the leap frog algorithm
    Vec3_t         F;     ///< Force over the particle
    Vec3_t         T;     ///< Torque over the particle
    Vec3_t         I;     ///< Vector containing the principal components of the inertia tensor
    Quaternion_t   Q;     ///< The quaternion representing the rotation
    double         R;     ///< Spheroradius
    double         rho;   ///< Density
    double         V;     ///< Volume
    double         m;     ///< Mass
    double         Erot;  ///< Rotational energy of the particle
    double         Ekin;  ///< Kinetical energy of the particle
    double         Dmax;  ///< Maximal distance from the center of mass to the surface of the body
    Array<Vec3_t*> Verts; ///< Vertices
    Array<Edge*>   Edges; ///< Edges
    Array<Face*>   Faces; ///< Faces
    
    // Methods
    void Draw (std::ostream & os, char const * Color="Blue", bool Blender=false); ///< Draw the particle
    void CalcMassProperties(size_t NCALCS = 5000);                                ///< Calculate the mass, center of mass and moment of inertia
    void StartForce(Vec3_t F) { F =  F;  T = 0,0,0;}                              ///< Start the force of the particle a value F, use for external forces like gravity
    void StartForce() { F = 0,0,0;  T = 0,0,0;}                                   ///< Start the force of the particle a value F, use for external forces like gravity
    void Start(double dt) { xb =  x -  v*dt;  wb =  w;}                           ///< Initialize the particle for the Verlet algorithm
    void DynamicRotation(double dt);                                              ///< Apply rotation on the particle once the total torque is found
    void DynamicTranslation(double dt);                                           ///< Apply translation once the total force is found
    void QuaternionRotation(Quaternion_t & Q,Vec3_t & V);                         ///< Apply rotation given by Quaternion Q at point v
    void Translation(Vec3_t & t);                                                 ///< Apply translation by vector t
    bool IsInside(Vec3_t & V);                                                    ///< Enquire if the point v is inside the particle.
    double IsInside(double *V);                                                   ///< Wrapper to the previous function
    double MaxX();                                                                ///< Find Maximun X coordinate
    double MaxY();                                                                ///< Find Maximun Y coordinate
    double MaxZ();                                                                ///< Find Maximun Y coordinate
    double MinX();                                                                ///< Find Minimun X coordinate
    double MinY();                                                                ///< Find Minimun Y coordinate
    double MinZ();                                                                ///< Find Minimun Y coordinate

    // Integrants for mass properties
    double Vol(double *r);                                      ///< Calculate the volume of the sample
    double Xc(double *r);                                       ///< Calculate the coordinates of the center of mass
    double Yc(double *r);                                       ///< Calculate the coordinates of the center of mass
    double Zc(double *r);                                       ///< Calculate the coordinates of the center of mass
    double Ixx(double *r);                                      ///< Calculate the inertia tensor
    double Iyy(double *r);                                      ///< Calculate the inertia tensor
    double Izz(double *r);                                      ///< Calculate the inertia tensor
    double Ixy(double *r);                                      ///< Calculate the inertia tensor
    double Ixz(double *r);                                      ///< Calculate the inertia tensor
    double Iyz(double *r);                                      ///< Calculate the inertia tensor

    
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Particle::Particle (Array<Vec3_t> const & V, Array<Array <int> > const & E, Array<Array <int> > const & F, Vec3_t const & v0, Vec3_t const & w0, double TheR, double TheRho, bool Open)
    : v(v0), w(w0), R(TheR), rho(TheRho)
{
    for (size_t i=0; i<V.Size(); i++) Verts.Push (new Vec3_t(V[i]));
    for (size_t i = 0;i<F.Size();i++)
    {
        Array<Vec3_t> verts(F[i].Size());
        for (size_t j=0; j<F[i].Size(); ++j) verts[j] = (*Verts[F[i][j]]);
        Faces.Push (new Face(verts));
    }
    if (Open)
    {
        Array<int> Faux(3);
        Array<Vec3_t> Normal(3);
        for (size_t i=0; i<V.Size(); i++)
        {
            size_t k=0;
            for (size_t j = 0;j<F.Size();j++)
            {
                if (F[j].Find(i)!=-1)
                {
                    Faux[k] = j;
                    Normal[k] = cross(Faces[j]->Edges[0]->dL,Faces[j]->Edges[1]->dL);
                    Normal[k] = Normal[k]/norm(Normal[k]);
                    Normal[k] = Faces[j]->Edges[0]->X0 - R*Normal[k];
                    k++;
                }
            }
            Mat_t M(9,9);
            Vec_t X(9),B(9);
            for (size_t j = 0;j<9;j++)
            {
                for (size_t k = 0;k<9;k++)
                {
                    M(j,k) = 0;
                }
                X(j) = 0;
                B(j) = 0;
            }
            M(0,0) = Faces[Faux[0]]->Edges[0]->dL(0);
            M(0,1) = Faces[Faux[0]]->Edges[1]->dL(0);
            M(0,6) = -1;
            M(1,0) = Faces[Faux[0]]->Edges[0]->dL(1);
            M(1,1) = Faces[Faux[0]]->Edges[1]->dL(1);
            M(1,7) = -1;
            M(2,0) = Faces[Faux[0]]->Edges[0]->dL(2);
            M(2,1) = Faces[Faux[0]]->Edges[1]->dL(2);
            M(2,8) = -1;
            M(3,2) = Faces[Faux[1]]->Edges[0]->dL(0);
            M(3,3) = Faces[Faux[1]]->Edges[1]->dL(0);
            M(3,6) = -1;
            M(4,2) = Faces[Faux[1]]->Edges[0]->dL(1);
            M(4,3) = Faces[Faux[1]]->Edges[1]->dL(1);
            M(4,7) = -1;
            M(5,2) = Faces[Faux[1]]->Edges[0]->dL(2);
            M(5,3) = Faces[Faux[1]]->Edges[1]->dL(2);
            M(5,8) = -1;
            M(6,4) = Faces[Faux[2]]->Edges[0]->dL(0);
            M(6,5) = Faces[Faux[2]]->Edges[1]->dL(0);
            M(6,6) = -1;
            M(7,4) = Faces[Faux[2]]->Edges[0]->dL(1);
            M(7,5) = Faces[Faux[2]]->Edges[1]->dL(1);
            M(7,7) = -1;
            M(8,4) = Faces[Faux[2]]->Edges[0]->dL(2);
            M(8,5) = Faces[Faux[2]]->Edges[1]->dL(2);
            M(8,8) = -1;

            B(0) = -Normal[0](0);
            B(1) = -Normal[0](1);
            B(2) = -Normal[0](2);
            B(3) = -Normal[1](0);
            B(4) = -Normal[1](1);
            B(5) = -Normal[1](2);
            B(6) = -Normal[2](0);
            B(7) = -Normal[2](1);
            B(8) = -Normal[2](2);

            Sol(M,B,X);
            Vec3_t Inter;
            Inter(0) = X(6);
            Inter(1) = X(7);
            Inter(2) = X(8);

            (*Verts[i]) = (Inter);
        }
        for (size_t i = 0;i<F.Size();i++) delete Faces[i];
        Faces.Clear();
        for (size_t i = 0;i<F.Size();i++)
        {
            Array<Vec3_t> verts(F[i].Size());
            for (size_t j=0; j<F[i].Size(); ++j) verts[j] = (*Verts[F[i][j]]);
            Faces.Push (new Face(verts));
        }
    }
    for (size_t i=0; i<E.Size(); i++) Edges.Push (new Edge((*Verts[E[i][0]]), (*Verts[E[i][1]])));
}

inline Particle::~Particle()
{
    for (size_t i=0; i<Verts.Size(); ++i) delete Verts[i];
    for (size_t i=0; i<Edges.Size(); ++i) delete Edges[i];
    for (size_t i=0; i<Faces.Size(); ++i) delete Faces[i];
}

inline void Particle::Draw (std::ostream & os, char const * Color, bool Blender)
{
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        if (Blender) BlenderDrawVert((*Verts[i]),os,R);
        else PovDrawVert((*Verts[i]),os,R,Color); 
    }
    for (size_t i=0; i<Edges.Size(); ++i)
    {
        Edges[i]->Draw(os,R,Color,Blender);
    }
    for (size_t i=0; i<Faces.Size(); ++i)
    {
        Faces[i]->Draw(os,R,Color,Blender);
    }
}

inline void Particle::CalcMassProperties(size_t NCALCS)
{
    if (Verts.Size()==1&&Edges.Size()==0&&Faces.Size()==0)
    {
        V = (4./3.)*M_PI*R*R*R;
        I = Vec3_t((8./15.)*M_PI*pow(R,5.),(8./15.)*M_PI*pow(R,5.),(8./15.)*M_PI*pow(R,5.));
        x = *Verts[0];
        Q = 1,0,0,0;
        m = rho*V;
    }
    else 
    {
        double Xi[3] = { MinX() , MinY() , MinZ() };
        double Xs[3] = { MaxX() , MaxY() , MaxZ() };
        Mat3_t It;
        Numerical::MonteCarlo<Particle> MC;
        V = MC.Integrate(this,&Particle::Vol,Xi,Xs,Numerical::VEGAS,NCALCS);
        x(0) = MC.Integrate(this,&Particle::Xc,Xi,Xs,Numerical::VEGAS,NCALCS)/V;
        x(1) = MC.Integrate(this,&Particle::Yc,Xi,Xs,Numerical::VEGAS,NCALCS)/V;
        x(2) = MC.Integrate(this,&Particle::Zc,Xi,Xs,Numerical::VEGAS,NCALCS)/V;
        It(0,0) = MC.Integrate(this,&Particle::Ixx,Xi,Xs,Numerical::VEGAS,NCALCS);
        It(1,1) = MC.Integrate(this,&Particle::Iyy,Xi,Xs,Numerical::VEGAS,NCALCS);
        It(2,2) = MC.Integrate(this,&Particle::Izz,Xi,Xs,Numerical::VEGAS,NCALCS);
        It(0,1) = It(1,0) = MC.Integrate(this,&Particle::Ixy,Xi,Xs,Numerical::VEGAS,NCALCS);
        It(0,2) = It(2,0) = MC.Integrate(this,&Particle::Ixz,Xi,Xs,Numerical::VEGAS,NCALCS);
        It(1,2) = It(2,1) = MC.Integrate(this,&Particle::Iyz,Xi,Xs,Numerical::VEGAS,NCALCS);
        Vec3_t xp,yp,zp;
        Eig(It,I,xp,yp,zp);
        Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
        Q(1) = (yp(2)-zp(1))/(4*Q(0));
        Q(2) = (zp(0)-xp(2))/(4*Q(0));
        Q(3) = (xp(1)-yp(0))/(4*Q(0));
        Q = Q/norm(Q);
        Rotation(w,Q,wb);
        w = wb;
        m = rho*V;
        Ekin = 0.5*m*dot(v,v);
        Erot = 0.5*(I(0)*w(0)*w(0)+I(1)*w(1)*w(1)+I(2)*w(2)*w(2));
        Dmax = Distance(x,(*Verts[0]))+R;
        for (size_t i=1; i<Verts.Size(); ++i)
        {
            if (Distance(x,(*Verts[i]))+R > Dmax) Dmax = Distance(x,(*Verts[i]))+R;
        }

    }
}  

inline void Particle::DynamicRotation(double dt)
{
    double q0,q1,q2,q3,wx,wy,wz;
    q0 = 0.5*Q(0);
    q1 = 0.5*Q(1);
    q2 = 0.5*Q(2);
    q3 = 0.5*Q(3);

    Vec3_t Td; 
    Td(0)=(T(0)+(I(1)-I(2))*wb(1)*wb(2))/I(0);
    Td(1)=(T(1)+(I(2)-I(0))*wb(0)*wb(2))/I(1);
    Td(2)=(T(2)+(I(0)-I(1))*wb(1)*wb(0))/I(2);
    w = wb+0.5*dt*Td;
    wx = w(0);
    wy = w(1);
    wz = w(2);
    Quaternion_t dq(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz),qm;

    wb  = wb+Td*dt;
    qm  = Q+dq*(0.5*dt);
    q0  = 0.5*qm(0);
    q1  = 0.5*qm(1);
    q2  = 0.5*qm(2);
    q3  = 0.5*qm(3);
    wx  = wb(0);
    wy  = wb(1);
    wz  = wb(2);
    dq  = Quaternion_t(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz);
    Quaternion_t Qd = (qm+dq*0.5*dt),temp;
    Conjugate(Q,temp);
    QuaternionRotation(temp,x);
    Q  = Qd/norm(Qd);
    QuaternionRotation(Q,x);    
    Erot=0.5*(I(0)*wx*wx+I(1)*wy*wy+I(2)*wz*wz);
}

inline void Particle::DynamicTranslation(double dt)
{
    Vec3_t temp,xa;
    xa    = 2*x - xb + F*(dt*dt/m);
    temp  = xa - x;
    v    = 0.5*(xa - xb)/dt;
    xb   = x;
    x    = xa;
    Ekin = 0.5*m*dot(v,v);
    Translation(temp);
}

inline void Particle::QuaternionRotation(Quaternion_t & Q,Vec3_t & V)
{
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size();
    for (size_t i = 0; i < nv; i++)
    {
        Vec3_t xt = *Verts[i]-V;
        Rotation(xt,Q,*Verts[i]);
        *Verts[i] += V;
    }

    for (size_t i = 0; i < ne; i++)
    {
        Edges[i]->Rotate(Q,V);
    }

    for (size_t i = 0; i < nf; i++)
    {
        Faces[i]->Rotate(Q,V);
    }
}

inline void Particle::Translation(Vec3_t & V)
{
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size();
    for (size_t i = 0; i < nv; i++)
    {
        *Verts[i] += V;
    }

    for (size_t i = 0; i < ne; i++)
    {
        Edges[i]->Translate(V);
    }

    for (size_t i = 0; i < nf; i++)
    {
        Faces[i]->Translate(V);
    }
}

inline bool Particle::IsInside (Vec3_t & V)
{
    bool inside = false;
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size();
    for (size_t i = 0; i < nv; i++)
    {
        if (Distance(V,*Verts[i]) < R) {
            inside = true;
            return inside;
        }
    }

    for (size_t i = 0; i < ne; i++)
    {
        if (Distance(V,*Edges[i]) < R) {
            inside = true;
            return inside;
        }
    }
    for (size_t i = 0; i < nf; i++)
    {
        if (Distance(V,*Faces[i]) < R) {
            inside = true;
            return inside;
        }
    }
    if (nf>2)
    {
        size_t k = 0;
        double Mindistance = Distance(V,*Faces[k]);
        for (size_t i = 1; i < nf; i++)
        {
            if (Distance(V,*Faces[i])<Mindistance) 
            {
                k = i;
                Mindistance = Distance(V,*Faces[k]);
            }
        }
        Vec3_t ct(0,0,0);
        for (size_t i = 0; i < Faces[k]->Edges.Size(); i++)
        {
            ct += Faces[k]->Edges[i]->X0;
        }
        ct = ct/double(Faces[k]->Edges.Size());
        Vec3_t pro = V - ct;
        Vec3_t nor = cross(Faces[k]->Edges[0]->dL,Faces[k]->Edges[1]->dL);
        nor = nor/norm(nor);
        if (dot(pro,nor)<0) inside =true;
    }


    return inside;
}

inline double Particle::MaxX()
{
    double result = (*Verts[0])(0)+R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(0)+R > result) result = (*Verts[i])(0)+R; 
    }
    return result;
}

inline double Particle::MaxY()
{
    double result = (*Verts[0])(1)+R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(1)+R > result) result = (*Verts[i])(1)+R; 
    }
    return result;
}

inline double Particle::MaxZ()
{
    double result = (*Verts[0])(2)+R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(2)+R > result) result = (*Verts[i])(2)+R; 
    }
    return result;
}

inline double Particle::MinX()
{
    double result = (*Verts[0])(0)-R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(0)-R < result) result = (*Verts[i])(0)-R; 
    }
    return result;
}

inline double Particle::MinY()
{
    double result = (*Verts[0])(1)-R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(1)-R < result) result = (*Verts[i])(1)-R; 
    }
    return result;
}

inline double Particle::MinZ()
{
    double result = (*Verts[0])(2)-R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(2)-R < result) result = (*Verts[i])(2)-R; 
    }
    return result;
}

inline double Particle::IsInside (double *r)
{
    Vec3_t p;
    p(0)=r[0];
    p(1)=r[1];
    p(2)=r[2];
    return double(IsInside(p));
}

inline double Particle::Vol (double *Thex)
{
    return IsInside(Thex);
}

inline double Particle::Xc (double *Thex)
{
    return Thex[0]*IsInside(Thex);
}

inline double Particle::Yc (double *Thex)
{
    return Thex[1]*IsInside(Thex);
}

inline double Particle::Zc (double *Thex)
{
    return Thex[2]*IsInside(Thex);
}

inline double Particle::Ixx (double *Thex)
{
    return ((Thex[1]-x(1))*(Thex[1]-x(1))+(Thex[2]-x(2))*(Thex[2]-x(2)))*IsInside(Thex);
}

inline double Particle::Iyy (double *Thex)
{
    return ((Thex[0]-x(0))*(Thex[0]-x(0))+(Thex[2]-x(2))*(Thex[2]-x(2)))*IsInside(Thex);
}

inline double Particle::Izz (double *Thex)
{
    return ((Thex[0]-x(0))*(Thex[0]-x(0))+(Thex[1]-x(1))*(Thex[1]-x(1)))*IsInside(Thex);
}

inline double Particle::Ixy (double *Thex)
{
    return -(Thex[0]-x(0))*(Thex[1]-x(1))*IsInside(Thex);
}

inline double Particle::Ixz (double *Thex)
{
    return -(Thex[0]-x(0))*(Thex[2]-x(2))*IsInside(Thex);
}

inline double Particle::Iyz (double *Thex)
{
    return -(Thex[1]-x(1))*(Thex[2]-x(2))*IsInside(Thex);
}

#endif // MECHSYS_DEM_PARTICLE_H
