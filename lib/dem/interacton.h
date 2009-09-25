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
     Interacton (Particle * Pt1, Particle * Pt2); ///< Constructor requires pointers to both particles

    // Methods
    void CalcForce (double dt = 0.0); ///< Calculates the contact force between particles

    // Data
    double        Kn;   ///< Normal stiffness
    double        Kt;   ///< Tengential stiffness
    double        Gn;   ///< Normal viscous coefficient
    double        Gt;   ///< Tangential viscous coefficient
    double        Mu;   ///< Microscpic coefficient of friction
    double        Mur;  ///< Rolling resistance coefficient
    double        Epot; ///< Ptential elastic energy
    FrictionMap_t Fdee; ///< Static friction displacement for pair of edges
    FrictionMap_t Fdvf; ///< Static friction displacement for pair of vertex-face
    FrictionMap_t Fdfv; ///< Static friction displacement for pair of face-vertex
    Particle    * P1;   ///< First particle
    Particle    * P2;   ///< Second particle

private:
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap, double dt);
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Interacton::Interacton (Particle * Pt1, Particle * Pt2)
    : P1(Pt1), P2(Pt2), Kn(10000.0), Kt(5000.0), Gn(16.), Gt(8), Mu(0.4), Epot(0.0)
{
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
            Vec3_t tan = FMap[p]/norm(FMap[p]);
            Vec3_t F = Kn*delta*n;
            if (norm(FMap[p])>Mu*norm(F)/Kt)
            {
                FMap[p] = Mu*norm(F)/Kt*tan;
            }
            F += Kt*FMap[p]+Gn*dot(n,vrel)*n+Gt*vt;
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

#endif //  MECHSYS_DEM_INTERACTON_H
