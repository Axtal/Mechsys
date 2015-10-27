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

#ifndef MECHSYS_DEM_EDGE_H
#define MECHSYS_DEM_EDGE_H

// Std lib
#include <cmath>

// MechSys
#include <mechsys/linalg/quaternion.h>

namespace DEM
{

class Edge
{
public:
    // Constructor
    Edge (Vec3_t * X0, Vec3_t * X1); ///< Xi: endpoints of edge
    Edge (Vec3_t & X0, Vec3_t & X1); ///< Xi: endpoints of edge

    // Methods
    void UpdatedL  ();  
    void Draw      (std::ostream & os, double Radius=1.0, char const * Color="Blue", bool BPY=false);

    // Data
    Vec3_t * X0; ///< Left endpoint
    Vec3_t * X1; ///< Right endpoint
    Vec3_t   dL; ///< Delta(X) = X1 - X0. difference Vector
    double Dmax; ///< Maximun length from edge centre
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Edge::Edge (Vec3_t * TheX0, Vec3_t * TheX1)
    : X0(TheX0), X1(TheX1), dL(*X1-*X0)
{
    Dmax = 0.5*norm(dL);
}

inline Edge::Edge (Vec3_t & TheX0, Vec3_t & TheX1)
    : X0(&TheX0), X1(&TheX1), dL(*X1-*X0)
{
    Dmax = 0.5*norm(dL);
}

inline void Edge::UpdatedL()
{
    dL = (*X1) - (*X0);
}

inline void Edge::Draw (std::ostream & os, double Radius, char const * Color, bool BPY)
{
    if (BPY)
    {
        Vec3_t cen = ((*X0)+(*X1))/2.0;
        double len = norm(dL);
        double l   = sqrt(dL(0)*dL(0)+dL(1)*dL(1));
        double phi = atan2(dL(1), dL(0));
        double th  = atan2(l, dL(2));
        os << "bpy.ops.mesh.primitive_cylinder_add(vertices=32, radius="<<Radius<<", depth="<<len<<", location=("<<cen(0)<<","<<cen(1)<<","<<cen(2)<<"), rotation=(0,"<<th<<","<<phi<<"))\n";
    }
    else
    {
        os << "cylinder { <"<<(*X0)(0)<<","<<(*X0)(1)<<","<<(*X0)(2)<<">,<"<<(*X1)(0)<<","<<(*X1)(1)<<","<<(*X1)(2)<<">,"<<Radius<<"\n pigment { color "<<Color<<" } }\n";
    }
}

}
#endif // MECHSYS_DEM_EDGE_H
