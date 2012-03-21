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
    Particle * P1;       ///< Pointer to first particle
    Particle * P2;       ///< Pointer to second particle
    double     Kn;       ///< Normal Spring constant 
    double     Kt;       ///< Tangential Spring constant 
    double     Gn;       ///< dissipation constant
    double     Gt;       ///< dissipation constant
    double     Mu;       ///< Friction coefficient
    double     Beta;     ///< Rolling stiffness coeffcient
    double     Eta;      ///< Plastic moment coefficient
    Vec3_t    SFr;       ///< Vector of static friction
    Vec3_t    Fdr;       ///< Vector of rolling resistance
};

Interacton::Interacton(Particle * Dp1, Particle * Dp2)
{
    P1 = Dp1;
    P2 = Dp2;
    Kn = 2.0*ReducedValue(P1->Kn,P2->Kn);
    Kt = 2.0*ReducedValue(P1->Kt,P2->Kt);
    Gn = 2.0*ReducedValue(P1->Gn,P2->Gn)*ReducedValue(P1->M,P2->M);
    Gt = 2.0*ReducedValue(P1->Gt,P2->Gt)*ReducedValue(P1->M,P2->M);
    Mu = 2.0*ReducedValue(P1->Mu,P2->Mu);
    Eta= 2.0*ReducedValue(P1->Eta,P2->Eta);
    Beta= 2.0*ReducedValue(P1->Beta,P2->Beta);
    SFr= OrthoSys::O;
    Fdr= OrthoSys::O;
#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#endif
}

void Interacton::CalcForce(double dt)
{
    double dist  = norm(P2->X - P1->X);
    double delta = P1->R + P2->R - dist;
    if (delta>0)
    {
        //Force
        Vec3_t n    = (P2->X - P1->X)/dist;
        Vec3_t x    = P1->X+n*((P1->R*P1->R-P2->R*P2->R+dist*dist)/(2*dist));
        Vec3_t Fn   = Kn*delta*n;
        Vec3_t t1,t2,x1,x2;
        Rotation(P1->W,P1->Q,t1);
        Rotation(P2->W,P2->Q,t2);
        x1 = x - P1->X;
        x2 = x - P2->X;
        Vec3_t Vrel = -((P2->V-P1->V)+cross(t2,x2)-cross(t1,x1));
        Vec3_t vt = Vrel - dot(n,Vrel)*n;
        SFr      += dt*vt;
        SFr      -= dot(SFr,n)*n;
        Vec3_t tan= SFr;
        if(norm(tan)>0.0) tan/=norm(tan);
        if(norm(SFr)>Mu*norm(Fn)/Kt)
        {
            SFr = Mu*norm(Fn)/Kt*tan;
        }
        Vec3_t F    = Fn + Gn*dot(n,Vrel)*n + Kt*SFr + Gt*vt;
        //if (dot(F,n)<0) F-=dot(F,n)*n;
        
        //Rolling resistance
        Vec3_t Vr = P1->R*P2->R*cross(Vec3_t(t1 - t2),n)/(P1->R+P2->R);
        Fdr += Vr*dt;
        Fdr -= dot(Fdr,n)*n;
        tan = Fdr;
        if (norm(tan)>0.0) tan/=norm(tan);
        double Kr = Beta*Kt;
        if (norm(Fdr)>Eta*Mu*norm(Fn)/Kr)
        {
            Fdr = Eta*Mu*norm(Fn)/Kr*tan;
        }
        Vec3_t Ft = -Kr*Fdr;

        //Assigning torque values
        Vec3_t T1 = -cross(x1,F) + P1->R*cross(n,Ft);
        Vec3_t T2 =  cross(x2,F) - P2->R*cross(n,Ft);
        Quaternion_t q;
        Conjugate (P1->Q,q);
        Rotation  (T1,q,t1);
        Conjugate (P2->Q,q);
        Rotation  (T2,q,t2);

#ifdef USE_THREAD
        //pthread_mutex_lock(&D1->lck);
        //pthread_mutex_lock(&D2->lck);
        pthread_mutex_lock(&lck);
#endif
        P1->F      -= F;
        P2->F      += F;
        P1->T      += t1;
        P2->T      += t2;
#ifdef USE_THREAD
        //pthread_mutex_unlock(&D1->lck);
        //pthread_mutex_unlock(&D2->lck);
        pthread_mutex_unlock(&lck);
#endif
    }
}

bool Interacton::UpdateContacts(double Alpha)
{
    if (norm(P1->X-P2->X) <= P1->R + P2->R + 2*Alpha) return true;
    else                                              return false;
}
}
#endif
