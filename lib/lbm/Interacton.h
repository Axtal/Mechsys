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


#ifndef MECHSYS_LBM_INTERACTON_H
#define MECHSYS_LBM_INTERACTON_H

// Mechsys
#include <mechsys/lbm/Dem.h>
#include <mechsys/dem/basic_functions.h>
namespace LBM
{
class Interacton
{
public:
    //Constructor
    Interacton () {};
    Interacton(Particle * D1, Particle * D2);

    //Methods
    void CalcForce      (double dt);
    bool UpdateContacts (double Alpha);

#ifdef USE_THREAD
    pthread_mutex_t lck;   ///< Lock to protect variables from race conditions.
#endif

    //Data
    Particle * D1;       ///< Pointer to first particle
    Particle * D2;       ///< Pointer to second particle
    double     Kn;       ///< Normal Spring constant 
    double     Kt;       ///< Tangential Spring constant 
    double     Gn;       ///< dissipation constant
    double     Mu;       ///< Friction coefficient
    Vec3_t    SFr;       ///< Vector of static friction
};

Interacton::Interacton(Particle * Dp1, Particle * Dp2)
{
    D1 = Dp1;
    D2 = Dp2;
    Kn = 2.0*ReducedValue(D1->Kn,D2->Kn);
    Kt = 2.0*ReducedValue(D1->Kt,D2->Kt);
    Gn = 2.0*ReducedValue(D1->Gn,D2->Gn)*ReducedValue(D1->M,D2->M);
    Mu = 2.0*ReducedValue(D1->Mu,D2->Mu);
    SFr= 0.0;
#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#endif
}

void Interacton::CalcForce(double dt)
{
    double dist  = norm(D2->X - D1->X);
    double delta = D1->R + D2->R - dist;
    if (delta>0)
    {
        Vec3_t n    = (D2->X - D1->X)/dist;
        Vec3_t Fn   = Kn*delta*n;
        Vec3_t Vrel = D1->V - D2->V;
        Vec3_t vt = Vrel - dot(n,Vrel)*n;
        SFr      += dt*vt;
        SFr      -= dot(SFr,n)*n;
        Vec3_t tan= SFr;
        if(norm(tan)>0.0) tan/=norm(tan);
        if(norm(SFr)>Mu*norm(Fn)/Kt)
        {
            SFr = Mu*norm(Fn)/Kt*tan;
        }
        Vec3_t F    = Fn + Gn*dot(n,Vrel)*n + Kt*SFr;
#ifdef USE_THREAD
        //pthread_mutex_lock(&D1->lck);
        //pthread_mutex_lock(&D2->lck);
        pthread_mutex_lock(&lck);
#endif
        D1->F      -= Kn*delta*n;
        D2->F      += Kn*delta*n;
#ifdef USE_THREAD
        //pthread_mutex_unlock(&D1->lck);
        //pthread_mutex_unlock(&D2->lck);
        pthread_mutex_unlock(&lck);
#endif
    }
}

bool Interacton::UpdateContacts(double Alpha)
{
    if (norm(D1->X-D2->X) <= D1->R + D2->R + 2*Alpha) return true;
    else                                              return false;
}
}
#endif
