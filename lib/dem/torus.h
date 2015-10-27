/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
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

#ifndef MECHSYS_DEM_TORUS_H
#define MECHSYS_DEM_TORUS_H

// Std lib
#include <cmath>

// MechSys
#include <mechsys/dem/face.h>

namespace DEM
{

class Torus
{
public:
    // Constructor
    Torus (Vec3_t const * X0, Vec3_t const * X1, Vec3_t const * X2); ///< X0: Center of the torus Xi: endpoints of edge
    Torus (Vec3_t const & X0, Vec3_t const & X1, Vec3_t const & X2); ///< X0: Center of the torus Xi: endpoints of edge

    // Methods
    void Draw      (std::ostream & os, double Radius=1.0, char const * Color="Blue", bool BPY=false);

    // Data
    Vec3_t const * X0; ///< Center of the Torus
    Vec3_t const * X1; ///< Right vector of the torus
    Vec3_t const * X2; ///< Left vector of the torus
    double          R; ///< Internal radius of the torus
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Torus::Torus (Vec3_t const * TheX0, Vec3_t const * TheX1, Vec3_t const * TheX2)
    : X0(TheX0), X1(TheX1), X2(TheX2)
{
    R = norm(*X0-*X1);
}

inline Torus::Torus (Vec3_t const & TheX0, Vec3_t const & TheX1, Vec3_t const & TheX2)
    : X0(&TheX0), X1(&TheX1), X2(&TheX2)
{
    R = norm(*X0-*X1);
}

inline void Torus::Draw (std::ostream & os, double Radius, char const * Color, bool BPY)
{
    if (BPY)
    {
        //throw new Fatal("Torus shape not implemented for blender");
    }
    else
    {
        //Vec3_t xp,yp,zp;
        //xp = (*X1)-(*X0);
        //zp = (*X2)-(*X0);
        //yp = cross(zp,xp);
        //xp/= norm(xp);
        //yp/= norm(yp);
        //zp/= norm(zp);
        //CheckDestroGiro(xp,yp,zp);
        //Quaternion_t Q;
        //Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
        //Q(1) = (yp(2)-zp(1))/(4*Q(0));
        //Q(2) = (zp(0)-xp(2))/(4*Q(0));
        //Q(3) = (xp(1)-yp(0))/(4*Q(0));
        //Q = Q/norm(Q);
        //double angle = 2*acos(Q(0));
        //double Rmax = norm((*X0)-(*X1));
        //double Rmin = Radius;
        //os << "torus  { \n" << Rmax << "," << Rmin << "\n Axis_Rotate_Trans (<" 
           //<< Q(1) << "," << Q(2) << "," << Q(3) << ">," << angle*180.0/M_PI << ") \n"
           //<< "translate <" << (*X0)(0) << "," << (*X0)(1) << "," << (*X0)(2) << ">"
           //<<"\n pigment { color "<<Color<<" } }\n";
    }
}
}
#endif // MECHSYS_DEM_TORUS_H
