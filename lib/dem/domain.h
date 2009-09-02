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

#ifndef MECHSYS_DEM_DOMAIN_H
#define MECHSYS_DEM_DOMAIN_H

// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI

// MechSys
#include "dem/particle.h"
//#include "dem/interacton.h"
#include "util/array.h"
#include "util/util.h"

class Domain
{
public:
    // Destructor
    ~Domain();

    // Methods to generate particles
    void GenerateSpheres (size_t N,     ///< Number of spheres
                          double Xmin,  ///< Left boundary
                          double Xmax,  ///< Right boundary
                          double Ymin,  ///< Back boundary
                          double Ymax,  ///< Front boundary
                          double Zmin,  ///< Bottom boundary
                          double Zmax,  ///< Top boundary
                          double rho,   ///< Density of the material
                          double Rmin); ///< Minimun radius in units of the maximun radius

    void GenerateBox (double Xmin,       ///< Left boundary
                      double Xmax,       ///< Right boundary
                      double Ymin,       ///< Back boundary
                      double Ymax,       ///< Front boundary
                      double Zmin,       ///< Bottom boundary
                      double Zmax,       ///< Top boundary
                      double Thickness); ///< Thickness of the wall, cannot be zero

    void AddTetra (Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a tetrahedron at position X with spheroradius R, side of length L and density rho
    void AddRice  (Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a rice at position X with spheroradius R, side of length L and density rho
    void AddCube  (Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho

    void WritePov     (char const * Filename);
    void WriteBlender (char const * Filename);

    // Methods
    //void CopyParticle (const Particle & P);  ///< Create a new particle as a copy of particle P, it should be translated to another position.
    //void Initialize   (double dt);           ///< Set the particles to a initial state and asign the possible insteractions
    //void OneStep      (double dt);           ///< One simualtion step

    /*
    void   LinearMomentum  (Vec3_t & L); ///< Return total momentum of the system
    void   AngularMomentum (Vec3_t & L); ///< Return total angular momentum of the system
    double TotalEnergy     ();           ///< Return total energy of the system
    */

    // Data
    Array<Particle*>   Particles;   ///< All particles in domain
    //Array<Interacton*> Interactons; ///< All interactons in domain
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Domain::~Domain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    //for (size_t i=0; i<Interactons.Size(); ++i) if (Interactons[i]!=NULL) delete Interactons[i];
}

inline void Domain::GenerateSpheres (size_t N, double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, double rho, double Rmin)
{
    double Lx   = Xmax-Xmin;
    double Ly   = Ymax-Ymin;
    double Lz   = Zmax-Zmin;
    double Rmax = pow(Lx*Ly*Lz/(8*N),1./3.);
    size_t nx   = Lx/(2*Rmax);
    size_t ny   = Ly/(2*Rmax);
    Array <Array <int> > Empty;
    for (size_t i=0; i<N; i++)
    {
        Array<Vec3_t> V(1);
        V[0] = Xmin + Rmax+2*Rmax*(i%nx), Ymin + Rmax+2*Rmax*((i/nx)%ny), Zmin + Rmax+2*Rmax*(i/(nx*ny));
        double R = Rmax*Rmin + (1.*rand())/RAND_MAX*Rmax*(1-Rmin);
        Particles.Push (new Particle(V,Empty,Empty, OrthoSys::O,OrthoSys::O, R,rho));
    }
}

inline void Domain::GenerateBox (double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, double Thickness)
{
}

inline void Domain::AddTetra (Vec3_t const & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    double sq8 = sqrt(8.0);
    Array<Vec3_t> V(4);
    V[0] =  L/sq8,  L/sq8, L/sq8;
    V[1] = -L/sq8, -L/sq8, L/sq8;
    V[2] = -L/sq8,  L/sq8,-L/sq8;
    V[3] =  L/sq8, -L/sq8,-L/sq8;

    // edges
    Array<Array <int> > E(6);
    for (size_t i=0; i<6; ++i) E[i].Resize(2);
    E[0] = 0, 1;
    E[1] = 1, 2;
    E[2] = 2, 0;
    E[3] = 0, 3;
    E[4] = 1, 3;
    E[5] = 2, 3;

    // face
    Array<Array <int> > F;
    F.Resize(4);
    for (size_t i=0; i<4; ++i) F[i].Resize(2);
    F[0] = 0, 3, 2;
    F[1] = 0, 1, 3;
    F[2] = 0, 2, 1;
    F[3] = 1, 2, 3;

    // random orientation
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        (*Axis) = (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX;
    }

    // calculate the rotation
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
}

inline void Domain::AddRice (const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(2);
    V[0] = 0.0, 0.0,  L/2;
    V[1] = 0.0, 0.0, -L/2;

    // edges
    Array<Array <int> > E(1);
    E[0].Resize(2);
    E[0] = 0, 1;

    // faces
    Array<Array <int> > F; // no faces

    // random orientation
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        (*Axis) = (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX;
    }

    // calculate the rotation
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
}

inline void Domain::AddCube (const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(8);
    double l = L/2.0;
    V[0] = -l, -l, -l;
    V[1] =  l, -l, -l;
    V[2] =  l,  l, -l;
    V[3] = -l,  l, -l;
    V[4] = -l, -l,  l;
    V[5] =  l, -l,  l;
    V[6] =  l,  l,  l;
    V[7] = -l,  l,  l;

    // edges
    Array<Array <int> > E(12);
    for (size_t i=0; i<12; ++i) E[i].Resize(2);
    E[ 0] = 0, 1;
    E[ 1] = 1, 2;
    E[ 2] = 2, 3;
    E[ 3] = 3, 0;
    E[ 4] = 4, 5;
    E[ 5] = 5, 6;
    E[ 6] = 6, 7;
    E[ 7] = 7, 4;
    E[ 8] = 0, 4;
    E[ 9] = 1, 5;
    E[10] = 2, 6;
    E[11] = 3, 7;

    // faces
    Array<Array <int> > F(6);
    for (size_t i=0; i<6; i++) F[i].Resize(4);
    F[0] = 0, 3, 7, 4;
    F[1] = 1, 2, 6, 5;
    F[2] = 0, 1, 5, 4;
    F[3] = 2, 3, 7, 6;
    F[4] = 0, 3, 2, 1;
    F[5] = 4, 5, 6, 7;

    // calculate the rotation
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
}

/*
inline void Domain::CopyParticle(const Particle & P)
{

}

inline void Domain::InitializeSimulation(double dt)
{
    for(size_t i = 0;i < NumberParticles();i++)
    {
        _Particles[i]->CalcMassProperties(5000);
        _Particles[i]->Start(dt);
    }

    for(size_t i = 0;i < NumberParticles()-1;i++)
    {
        for(size_t j = i+1;j < NumberParticles();j++)
        {
            _Interactons.Push(new Interacton(_Particles[i],_Particles[j]));
        }
    }
}

inline void Domain::OneStep(double dt)
{
    for(size_t i = 0;i < NumberParticles();i++)
    {
        _Particles[i]->StartForce();
    }
    for(size_t i = 0;i < _Interactons.Size();i++)
    {
        _Interactons[i]->CalcForce(dt);
    }
    for(size_t i = 0;i < NumberParticles();i++)
    {
        _Particles[i]->DynamicRotation(dt);
        _Particles[i]->DynamicTranslation(dt);
    }
}

inline void Domain::LinearMomentum(Vec3_t & L)
{
    L = 0.,0.,0.;
    for(size_t i = 0;i < NumberParticles();i++)
    {
        L+=_Particles[i]->Mass()*_Particles[i]->v();
    }
}

inline void Domain::AngularMomentum(Vec3_t & L)
{
    L = 0.,0.,0.;
    for(size_t i = 0;i < NumberParticles();i++)
    {
        Vec3_t t1,t2;
        t1 = _Particles[i]->I()(0)*_Particles[i]->w()(0),_Particles[i]->I()(1)*_Particles[i]->w()(1),_Particles[i]->I()(2)*_Particles[i]->w()(2);
        Rotation(t1,_Particles[i]->Q(),t2);
        L += _Particles[i]->Mass()*cross(_Particles[i]->r(),_Particles[i]->v())+t2;
    }
    
}

inline double Domain::TotalEnergy()
{
    double E = 0;
    for(size_t i = 0;i < NumberParticles();i++)
    {
        E += _Particles[i]->KinEnergy();
    }

    for(size_t i = 0;i < _Interactons.Size();i++)
    {
        E += _Interactons[i]->_Epot;
    }
    return E;
}
*/

#endif // MECHSYS_DEM_DOMAIN_H
