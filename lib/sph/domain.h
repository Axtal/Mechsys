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

#ifndef MECHSYS_SPH_DOMAIN_H
#define MECHSYS_SPH_DOMAIN_H

// Mechsys
#include <mechsys/sph/particle.h>
#include <mechsys/sph/special_functions.h>
#include <mechsys/dem/special_functions.h>
#include <mechsys/dem/graph.h>


class SPHDomain
{
public:

    // Constructor
    SPHDomain();

    // Methods
    void AddBox(Vec3_t const & x, size_t nx, size_t ny, size_t nz, double R, double rho0, bool Fixed); ///< Add a box of SPHparticles
    void StartAcceleration (Vec3_t const & a = Vec3_t(0.0,0.0,0.0));                                   ///< Add a fixed acceleration
    void ComputeAcceleration ();                                                                       ///< Compute the accleration due to the other particles
    void Move                (double dt);                                                              ///< Compute the accleration due to the other particles
    void WriteBPY (char const * FileKey);                                                              ///< Draw the entire domain in a POV file
    void WritePOV (char const * FileKey);                                                              ///< Draw the entire domain in a blender file

    // Data
    Vec3_t CamPos;                    ///< Camera position
    Array <SPHParticle*> Particles;   ///< Array of SPH particles

};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// Constructor
inline SPHDomain::SPHDomain ()
{
    CamPos = 1.0,2.0,3.0;
}


// Methods
inline void SPHDomain::AddBox(Vec3_t const & V, size_t nx, size_t ny, size_t nz, double R, double rho0, bool Fixed)
{
    Vec3_t C(V);
    C -= Vec3_t((nx-1)*R,(ny-1)*R,(nz-1)*R);
    for (size_t i = 0;i<nx;i++)
    for (size_t j = 0;j<ny;j++)
    for (size_t k = 0;k<nz;k++)
    {
        Vec3_t x(2*i*R,2*j*R,2*k*R);
        x -=C;
        Particles.Push(new SPHParticle(x,OrthoSys::O,rho0,R,Fixed));
    }
}

inline void SPHDomain::StartAcceleration (Vec3_t const & a)
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->a = a;
        Particles[i]->dDensity = 0.0;
    }
}

inline void SPHDomain::ComputeAcceleration ()
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
        for (size_t j=0; j<Particles.Size(); j++)
        {
            if (Particles[i]->IsFree)
            {
                double di = Particles[i]->Density;
                double dj = Particles[j]->Density;
                double d0 = Particles[j]->Density0;
                double h  = 2*ReducedValue(Particles[i]->h,Particles[j]->h);
                Vec3_t vij = Particles[j]->v - Particles[i]->v;
                Vec3_t rij = Particles[j]->x - Particles[i]->x;

                //Calculate the fictional viscosity
                double alphacij = 0.1;
                double beta = 0.1;
                double muij = h*dot(vij,rij)/(dot(rij,rij)+0.01*h*h);
                double piij;
                if (dot(vij,rij)>0) piij = 2*(-alphacij*muij+beta*muij)/(di+dj);
                else                piij = 0.0;
                Particles[i]->a -= d0*(Pressure(di)/(di*di)+Pressure(dj)/(dj*dj)+piij)*rij*GradSPHKernel(norm(rij),h)/norm(rij);
                Particles[i]-> dDensity -= d0*dot(vij,rij)*GradSPHKernel(norm(rij),h)/norm(rij);
            }
        }
    }
}

inline void SPHDomain::Move (double dt)
{
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Move(dt);
}

inline void SPHDomain::WritePOV (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".pov");
    std::ofstream of(fn.CStr(), std::ios::out);
    POVHeader (of);
    POVSetCam (of, CamPos, OrthoSys::O);
    for (size_t i=0; i<Particles.Size(); i++) 
    {
        if (Particles[i]->IsFree) POVDrawVert(Particles[i]->x,of,Particles[i]->h,"Blue");
        else                      POVDrawVert(Particles[i]->x,of,Particles[i]->h,"Col_Glass_Ruby");
    }
    of.close();
}

inline void SPHDomain::WriteBPY (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".bpy");
    std::ofstream of(fn.CStr(), std::ios::out);
    BPYHeader(of);
    for (size_t i=0; i<Particles.Size(); i++) BPYDrawVert(Particles[i]->x,of,Particles[i]->h);
    of.close();
}


#endif // MECHSYS_SPH_DOMAIN_H
