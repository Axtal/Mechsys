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

#ifndef DEM_INTERACTON_H
#define DEM_INTERACTON_H

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
    void CalcForce (double Dt); ///< Calculates the contact force between particles

    // Data
    double        _Kn;   ///< Normal stiffness
    double        _Kt;   ///< Tengential stiffness
    double        _gn;   ///< Normal viscous coefficient
    double        _gt;   ///< Tangential viscous coefficient
    double        _mu;   ///< Microscpic coefficient of friction
    double        _mur;  ///< Rolling resistance coefficient
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
            Distance ((*A[i]), (*B[i]), ri, rf);
            std::pair<int,int> p = std::make_pair(i,j);
            FMap[p] += 0.0;
        }
    }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


//inline Interacton::Interacton(Particle * Pt1, Particle * Pt2)
    //: _p1(Pt2), _p2(Pt2)
//{
//}
//
//inline Interacton::~Interacton()
//{
//}
//
//inline void Interacton::CalcForce(double Dt)
//{
    //_update_disp_calc_force (_p1->_edges,  _p2->_edges,  _fdee);
    //_update_disp_calc_force (_p1->_vertex, _p2->_faces,  _fdvf);
    //_update_disp_calc_force (_p1->_faces,  _p2->_vertex, _fdfv);
//
    ///*
    //for(size_t i = 0; i < _p1->_edges.Size();i++)
    //for(size_t j = i; j < _p2->_edges.Size();j++)
    //{
        //Vec3_t rf,ri;
        //Distance(*_p1->_edges[i],*_p2->_edges[j],ri,rf);
        //double dist = norm(rf-ri);
        //Vec3_t n = (rf-ri)/dist;
        //double delta = dist-_p1->_R-_p2->_R;
        //if(delta<=0)
        //{
            //_edges[i].Others.Push(_edges[j]);
            //Vec3_t f = _Kn*delta*n;
            //Vec3_t vrel = _p2->_v-_p1->_v;
            //Vec3_t vt = vrel-dot(vrel,n)*n;
            //pair<int,int> p = make_pair(i,j);
            //_fdee[p] += Dt*vt;
            //if (norm(_fdee[p])>mu*f/Kt) _fdee=mu*f/Kt
            //if (_fdee.has_key(p)) asdfasdfa
            //else asdfasdfasd
            //f += -_Kt*_fdee[p];
        //}
        //else _fdee[p]=zero
    //}
//
    //for(size_t i = 0; i < _p1->_vertex.Size();i++)
    //for(size_t j = 0; j < _p2->_faces.Size();j++)
    //{
        //Vec3_t rf,ri;
        //Distance(*_p1->_vertex[i],*_p2->_faces[j],ri,rf);
        //double dist = norm(rf-ri);
        //Vec3_t n = (rf-ri)/dist;
        //double delta = dist-_p1->_R-_p2->_R;
        //if(delta<0)
        //{
            //Vec3_t f = _Kn*delta*n;
            //Vec3_t vrel = _p2->_v-_p1->_v;
            //Vec3_t vt = vrel-dot(vrel,n)*n;
            //pair<int,int> p = make_pair(i,j);
            //_fdvf[p] += Dt*vt;
        //}
    //}
//
    //for(size_t i = 0; i < _p1->_faces.Size();i++)
    //for(size_t j = 0; j < _p2->_vertex.Size();j++)
    //{
        //Vec3_t rf,ri;
        //Distance(*_p1->_faces[i],*_p2->_vertex[j],ri,rf);
        //double dist = norm(rf-ri);
        //Vec3_t n = (rf-ri)/dist;
        //double delta = dist-_p1->_R-_p2->_R;
        //if(delta<0)
        //{
            //Vec3_t f = _Kn*delta*n;
            //Vec3_t vrel = _p2->_v-_p1->_v;
            //Vec3_t vt = vrel-dot(vrel,n)*n;
            //pair<int,int> p = make_pair(i,j);
            //_fdfv[p] += Dt*vt;
        //}
    //}
    //*/
//
//}


#endif // DEM_INTERACTON_H
