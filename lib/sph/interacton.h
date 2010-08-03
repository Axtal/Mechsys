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

#ifndef MECHSYS_SPH_INTERACTON_H
#define MECHSYS_SPH_INTERACTON_H

// Mechsys
#include <mechsys/sph/particle.h>
#include <mechsys/sph/special_functions.h>
#include <mechsys/dem/special_functions.h>

class SPHInteracton
{   
public:
    // Constructor and destructor
    SPHInteracton (SPHParticle * Pt1, SPHParticle * Pt2); ///< Default constructor

    // Methods
    bool UpdateContacts (double alpha);          ///< Update contacts by verlet algorithm
    void CalcForce      (double dt = 0.0);       ///< Calculates the contact force between particles
    

    // Data
    SPHParticle * P1;                            ///< Pointer to first particle
    SPHParticle * P2;                            ///< Pointer to second particle
    double alpha;                                ///< Coefficient of bulk viscosity
    double beta;                                 ///< Coefficient of Neumann - Richtmyer viscosity
    double h;                                    ///< Smoothing length
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline SPHInteracton::SPHInteracton (SPHParticle * Pt1, SPHParticle * Pt2)
{
    P1 = Pt1;
    P2 = Pt2;
    h  = 2*ReducedValue(P1->h,P2->h);
    alpha = 0.25;
    beta = 0.25;
}

inline void SPHInteracton::CalcForce(double dt)
{
    double di = P1->Density;
    double dj = P2->Density;
    double d0i = P1->Density0;
    double d0j = P2->Density0;
    Vec3_t vij = P2->v - P1->v;
    Vec3_t rij = P2->x - P1->x;
    double muij = h*dot(vij,rij)/(dot(rij,rij)+0.01*h*h);
    double cij = (SoundSpeed(di)+SoundSpeed(dj));
    double piij;
    //if (dot(vij,rij)<0) piij = (-alpha*cij*muij+beta*muij*muij)/(di+dj);
    //else                piij = 0.0;
    piij = -alpha*muij*cij/(di+dj);
    P1->a += d0j*(Pressure(di)/(di*di)+Pressure(dj)/(dj*dj)+piij)*rij*GradSPHKernel(norm(rij),h)/norm(rij);
    P2->a -= d0i*(Pressure(di)/(di*di)+Pressure(dj)/(dj*dj)+piij)*rij*GradSPHKernel(norm(rij),h)/norm(rij);
    P1->dDensity += d0j*dot(vij,rij)*GradSPHKernel(norm(rij),h)/norm(rij);
    P2->dDensity += d0i*dot(vij,rij)*GradSPHKernel(norm(rij),h)/norm(rij);
}

inline bool SPHInteracton::UpdateContacts (double alpha)
{
    if (Distance(P1->x,P2->x)<=P1->h+P2->h+2*alpha) return true;
    else return false;
}



#endif // MECHSYS_SPH_INTERACTON_H
