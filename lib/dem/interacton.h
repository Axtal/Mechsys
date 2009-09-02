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

#include <math.h>
#include <map>
#include <vector>
#include <utility>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/particle.h"


class Interacton
{
public:
    // typedefs
    typedef std::map<std::pair<int,int>,Vec3_t> FrictionMap_t;

    // Constructor and destructor
     Interacton (Particle * Pt1, Particle * Pt2); ///< Constructor, it requires pointers to both particles
    ~Interacton ();                               ///< Destructor

    // Methods
    void CalcForce (double dt); ///< Calculates the contact force between particles

    // Data
    double        _Kn;   ///< Normal stiffness
    double        _Kt;   ///< Tengential stiffness
    double        _gn;   ///< Normal viscous coefficient
    double        _gt;   ///< Tangential viscous coefficient
    double        _mu;   ///< Microscpic coefficient of friction
    double        _mur;  ///< Rolling resistance coefficient
    double        _Epot; ///< Ptential elastic energy
    FrictionMap_t _fdee; ///< Static friction displacement for pair of edges
    FrictionMap_t _fdvf; ///< Static friction displacement for pair of vertex-face
    FrictionMap_t _fdfv; ///< Static friction displacement for pair of face-vertex
    Particle    * _p1;   ///< First particle
    Particle    * _p2;   ///< Second particle

private:
    template<typename FeatureA_T, typename FeatureB_T>
    void _update_disp_calc_force (FeatureA_T & A, FeatureB_T & B, FrictionMap_t & FMap)
    {
        for (size_t i=0; i<A.Size(); ++i)
        for (size_t j=0; j<B.Size(); ++j)
        {
            Vec3_t ri, rf;
            Distance ((*A[i]), (*B[j]), ri, rf);
            double dist = norm(rf-ri);
            double delta = _p1->Radius() + _p2->Radius() - dist;
            if(delta>=0)
            {
                Vec3_t n = (rf-ri)/dist;
                Vec3_t F = _Kn*delta*n;
                _p1->F()+=-F;
                _p2->F()+=F;
                Vec3_t r = ri+n*((_p1->Radius()*_p1->Radius()-_p2->Radius()*_p2->Radius()+dist*dist)/(2*dist));
                
                Vec3_t T,Tt,temp;
                temp = r - _p1->r();
                Tt = cross(temp,F);
                Quaternion_t q;
                Conjugate(_p1->Q(),q);
                Rotation(Tt,q,T);
                _p1->T()-=T;
                temp = r - _p2->r();
                Tt = cross(temp,F);
                Conjugate(_p2->Q(),q);
                Rotation(Tt,q,T);
                _p2->T()+=T;
                _Epot += 0.5*_Kn*delta*delta;
            }
            //std::pair<int,int> p = std::make_pair(i,j);
            //FMap[p] += 0.0;
        }
    }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Interacton::Interacton(Particle * Pt1, Particle * Pt2)
{
    _p1 = Pt1;
    _p2 = Pt2;
    _Kn = 1000;
}

inline void Interacton::CalcForce(double Dt)
{
    _Epot = 0;
    _update_disp_calc_force (_p1->_edges,  _p2->_edges,  _fdee);
    _update_disp_calc_force (_p1->_vertex, _p2->_faces,  _fdvf);
    _update_disp_calc_force (_p1->_faces,  _p2->_vertex, _fdfv);
}

#endif //  MECHSYS_DEM_INTERACTON_H
