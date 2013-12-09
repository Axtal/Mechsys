/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
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
        Vec3_t mid = ((*X0)+(*X1))/2.0;
        double L   = norm(dL);
        Vec3_t Bv(0.0,0.0,L/2.0);
        Vec3_t axis;
        axis = (*X1)-mid;
        axis = cross(Bv,axis);
        double angle;
        if (norm(axis)<1.0e-10)
        {
            axis =1.0,0.0,0.0;
            angle = 0.0;
        }
        else angle = acos(4.0*dot(Bv,(*X1)-mid)/(L*L));
        if (false) // does not work
        {
            Quaternion_t q(angle, axis(0), axis(1), axis(2));
            Vec3_t r(0,0,0);
            Rotation(Vec3_t(0,0,M_PI/2.0), q, r);
            os << "bpy.ops.mesh.primitive_cylinder_add(vertices=32, radius="<<Radius<<", depth="<<L<<", location=("<<mid(0)<<","<<mid(1)<<","<<mid(2)<<"), rotation=("<<r(0)<<","<<r(1)<<","<<r(2)<<"))\n";
        }
        if (false) // also does not work
        {
            os << "bpy.ops.mesh.primitive_cylinder_add(vertices=32, radius="<<Radius<<", depth="<<L<<", location=("<<mid(0)<<","<<mid(1)<<","<<mid(2)<<"))\n";
            os << "axis = ["<<axis(0)<<","<<axis(1)<<","<<axis(2)<<"]\n";
            os << "quat = mathutils.Quaternion(axis,"<<angle*180/M_PI<<")\n";
            os << "quat.normalize()\n";
            //os << "bpy.ops.object.editmode_togle()\n";
            os << "bpy.ops.transform.rotate(value=quat.angle, axis=quat.axis)\n";
        }
        if (false) // also does not work
        {
            os << "bpy.ops.mesh.primitive_cylinder_add(vertices=32, radius="<<Radius<<", depth="<<L<<")\n";
            os << "v  = mathutils.Vector(("<<dL(0)<<","<<dL(1)<<","<<dL(2)<<"))\n";
            os << "up = mathutils.Vector((0,0,1))\n";
            os << "if v != -up:\n";
            os << "    rot = up.rotation_difference(v)\n";
            os << "else:\n";
            os << "    rot = mathutils.Quaternion((1,0,0),pi)\n";
            os << "bpy.ops.transform.translate(value=(0,0,v.length/2.0))\n";
            //os << "bpy.ops.transform.translate(value=(0,0,v.length/2.0))\n";
            os << "bpy.ops.object.editmode_toggle()\n";
            os << "bpy.ops.transform.rotate(value=rot.angle, axis=rot.axis)\n";
            os << "bpy.ops.object.editmode_toggle()\n";
        }
    }
    else
    {
        os << "cylinder { <"<<(*X0)(0)<<","<<(*X0)(1)<<","<<(*X0)(2)<<">,<"<<(*X1)(0)<<","<<(*X1)(1)<<","<<(*X1)(2)<<">,"<<Radius<<"\n pigment { color "<<Color<<" } }\n";
    }
}

}
#endif // MECHSYS_DEM_EDGE_H
