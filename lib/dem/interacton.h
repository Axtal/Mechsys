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

#ifndef MECHSYS_DEM_INTERACTON_H
#define MECHSYS_DEM_INTERACTON_H

// Std Lib
#include <math.h>
#include <map>
#include <vector>
#include <utility>

// MechSys
#include "dem/particle.h"

class Interacton
{
public:
    // typedefs
    typedef std::map<std::pair<int,int>,Vec3_t> FrictionMap_t;

    // Constructor and destructor
    Interacton () {};                            ///< Default constructor
    Interacton (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles

    // Methods
    virtual void CalcForce (double dt = 0.0); ///< Calculates the contact force between particles

    // Data
    Particle    * P1;   ///< First particle
    Particle    * P2;   ///< Second particle
    double        Kn;   ///< Normal stiffness
    double        Kt;   ///< Tengential stiffness
    double        Gn;   ///< Normal viscous coefficient
    double        Gt;   ///< Tangential viscous coefficient
    double        Mu;   ///< Microscpic coefficient of friction
    double        Epot; ///< Potential elastic energy
    Vec3_t        Fn;   ///< Normal force between elements
    FrictionMap_t Fdee; ///< Static friction displacement for pair of edges
    FrictionMap_t Fdvf; ///< Static friction displacement for pair of vertex-face
    FrictionMap_t Fdfv; ///< Static friction displacement for pair of face-vertex
protected:
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, double dt);
};


class InteractonSphere: public Interacton
{
public:
    // Methods 
    InteractonSphere (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles
    void CalcForce (double dt = 0.0);                  ///< Calculates the contact force between particles
    
    // Data
    FrictionMap_t Fdvv;                                ///< Static Friction displacement for the vertex vertex pair
    Vec3_t        Fdr;                                 ///< Rolling displacement 
    double        beta;                                ///< Rolling stiffness coefficient
    double        eta;                                 ///< Plastic moment coefficient

protected:
    void _update_rolling_resistance(double dt);               ///< Calculates the rolling resistance torque
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Interacton::Interacton (Particle * Pt1, Particle * Pt2)
    : P1(Pt1), P2(Pt2), Kn(2*P1->Kn*P2->Kn/(P1->Kn + P2->Kn)), Kt(2*P1->Kt*P2->Kt/(P1->Kt + P2->Kt)), Gn(2*P1->Gn*P2->Gn/(P1->Gn + P2->Gn)), Gt(2*P1->Gt*P2->Gt/(P1->Gt + P2->Gt)), Mu(2*P1->Mu*P2->Mu/(P1->Mu + P2->Mu)), Epot(0.0)
{
    CalcForce(0.1);
}

inline void Interacton::CalcForce (double dt)
{
    Epot = 0.0;
    if (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax)
    {
        _update_disp_calc_force (P1->Edges,P2->Edges,Fdee,dt);
        _update_disp_calc_force (P1->Verts,P2->Faces,Fdvf,dt);
        _update_disp_calc_force (P1->Faces,P2->Verts,Fdfv,dt);
    }
}

template<typename FeatureA_T, typename FeatureB_T>
inline void Interacton::_update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, double dt)
{
    // update
    for (size_t i=0; i<A.Size(); ++i)
    for (size_t j=0; j<B.Size(); ++j)
    {
        Vec3_t xi, xf;
        Distance ((*A[i]), (*B[j]), xi, xf);
        double dist  = norm(xf-xi);
        double delta = P1->R + P2->R - dist;
        if (delta>=0)
        {
            // update force
            Vec3_t n = (xf-xi)/dist;
            Vec3_t x = xi+n*((P1->R*P1->R-P2->R*P2->R+dist*dist)/(2*dist));
            Vec3_t t1,t2,x1,x2;
            Rotation(P1->w,P1->Q,t1);
            Rotation(P2->w,P2->Q,t2);
            x1 = x - P1->x;
            x2 = x - P2->x;
            Vec3_t vrel = -((P2->v-P1->v)+cross(t2,x2)-cross(t1,x1));
            Vec3_t vt = vrel - dot(n,vrel)*n;
            pair<int,int> p;
            p = make_pair(i,j);
            FMap[p] += vt*dt;
            FMap[p] -= dot(FMap[p],n)*n;
            Vec3_t tan = FMap[p];
            if (norm(tan)>1.0e-22) tan/=norm(tan);
            Fn = Kn*delta*n;
            if (norm(FMap[p])>Mu*norm(Fn)/Kt)
            {
                FMap[p] = Mu*norm(Fn)/Kt*tan;
            }
            Vec3_t F = Fn + Kt*FMap[p] +Gn*dot(n,vrel)*n +Gt*vt;
            P1->F += -F;
            P2->F +=  F;

            // torque
            Vec3_t T, Tt, temp;
            temp = x - P1->x;
            Tt = cross (temp,F);
            Quaternion_t q;
            Conjugate (P1->Q,q);
            Rotation  (Tt,q,T);
            P1->T -= T;
            temp = x - P2->x;
            Tt = cross (temp,F);
            Conjugate (P2->Q,q);
            Rotation  (Tt,q,T);
            P2->T += T;

            // potential energy
            Epot += 0.5*Kn*delta*delta+0.5*Kt*dot(FMap[p],FMap[p]);
        }
    }
}

inline InteractonSphere::InteractonSphere (Particle * Pt1, Particle * Pt2)
{
    P1   = Pt1;
    P2   = Pt2;
    Kn   = 2*Pt1->Kn*Pt2->Kn/(Pt1->Kn+Pt2->Kn);
    Kt   = 2*Pt1->Kt*Pt2->Kt/(Pt1->Kn+Pt2->Kt);
    Gn   = 2*Pt1->Gn*Pt2->Gn/(Pt1->Gn+Pt2->Gn);
    Gt   = 2*Pt1->Gt*Pt2->Gt/(Pt1->Gn+Pt2->Gt);
    Mu   = 2*Pt1->Mu*Pt2->Mu/(Pt1->Mu+Pt2->Mu);
    beta = 0.12;
    eta  = 1.0;

    Epot = 0.0;
    Fdr  = 0.0, 0.0, 0.0;

    CalcForce(0.0);
}

inline void InteractonSphere::_update_rolling_resistance(double dt)
{
    //Vec3_t t1,t2;
    //Rotation(P1->w,P1->Q,t1);
    //Rotation(P2->w,P2->Q,t2);
    //Vec3_t Normal = Fn/norm(Fn);
    //Vec3_t Vr = P1->R*P2->R*cross(Vec3_t(t2 - t1),Normal)/(P1->R+P2->R);
    //Vec3_t Theta = cross(Vr,Normal);
    //if (norm(Theta)>1.0e-22) Theta /= norm(Theta);
    //double r = 0.5*(P1->R+P2->R);
    //Fdr += norm(Vr)*dt*Theta/r;
    //Vec3_t tan = Fdr;
    //if (norm(tan)>1.0e-22) tan/=norm(tan);
    //double Kr = beta*Kt*r*r;
    //if (norm(Fdr)>eta*r*norm(Fn)/Kr)
    //{
        //Fdr = eta*r*norm(Fn)/Kr*tan;
    //}
    //Vec3_t Tr = -Kr*Fdr;
    //Vec3_t T;
    //Quaternion_t q;
    //Conjugate (P1->Q,q);
    //Rotation  (Tr,q,T);
    //P2->T -= T;
    //Conjugate (P2->Q,q);
    //Rotation  (Tr,q,T);
    //P1->T += T;

    Vec3_t t1,t2;
    Rotation(P1->w,P1->Q,t1);
    Rotation(P2->w,P2->Q,t2);
    Vec3_t Normal = Fn/norm(Fn);
    Vec3_t Vr = P1->R*P2->R*cross(Vec3_t(t1 - t2),Normal)/(P1->R+P2->R);
    Fdr += Vr*dt;
    Fdr -= dot(Fdr,Normal)*Normal;
    Vec3_t tan = Fdr;
    if (norm(tan)>1.0e-22) tan/=norm(tan);
    double Kr = beta*Kt;
    if (norm(Fdr)>eta*norm(Fn)/Kr)
    {
        Fdr = eta*norm(Fn)/Kr*tan;
    }
    Vec3_t Ft = -Kr*Fdr;

    Vec3_t Tt = P1->R*cross(Normal,Ft);
    Vec3_t T;
    Quaternion_t q;
    Conjugate (P1->Q,q);
    Rotation  (Tt,q,T);
    P1->T += T;

    Tt = P2->R*cross(Normal,Ft);
    Conjugate (P2->Q,q);
    Rotation  (Tt,q,T);
    P2->T -= T;
}

inline void InteractonSphere::CalcForce(double dt)
{
    Epot = 0.0;
    _update_disp_calc_force (P1->Verts,P2->Verts,Fdvv,dt);
    if (Epot>0.0) _update_rolling_resistance(dt);
}


#endif //  MECHSYS_DEM_INTERACTON_H
