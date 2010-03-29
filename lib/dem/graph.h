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
 * MERCHANTABILITY or finNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef MECHSYS_DEM_GRAPH_H
#define MECHSYS_DEM_GRAPH_H

// Std lib
#include <iostream>

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/linalg/matvec.h>


/////////////////////////////////////////////////////////////////////////////////////////// PovRay /////


inline void POVHeader (std::ostream & os)
{
    os << "#include \"colors.inc\" \n";
    os << "#include \"glass.inc\" \n";
    os << "background {color White} \n";
    os << "light_source{<10,0,0> color White shadowless}  \n";
    os << "light_source{<-10,0,0> color White shadowless}  \n";
    os << "light_source{<0,10,0> color White shadowless}  \n";
    os << "light_source{<0,-10,0> color White shadowless}  \n";
    os << "light_source{<0,0,10> color White shadowless}  \n";
    os << "light_source{<0,0,-10> color White shadowless}  \n";
}   

inline void POVSetCam (std::ostream & os, const Vec3_t & X, const Vec3_t & F)
{
    os << "camera { location <"<<X(0)<<","<<X(1)<<","<<X(2)<<"> sky <0,0,1> look_at <"<<F(0)<<","<<2<<","<<F(2)<<"> }\n";
}

inline void POVDrawVert (Vec3_t const & V, std::ostream & os, double Radius=1.0, char const * Color="Blue")
{
    os << "sphere  { <"<<V(0)<<","<<V(1)<<","<<V(2)<<">,"<<Radius<<"\n pigment { color "<<Color<<" } }\n";
}

inline void POVDrawPolygon (Array<Vec3_t> const & V, std::ostream & os, char const * Color="Blue")
{
    size_t N = V.Size();
    Vec3_t mid;
    mid = 0.0, 0.0, 0.0;
    for (size_t i=0; i<V.Size(); i++) mid += V[i];
    mid /= V.Size();
    for (size_t i=0; i<V.Size(); i++)
    {
        os << "polygon {"<<3<<", \n";
        os << "<"<<V[i](0)<<","<<V[i](1)<<","<<V[i](2)<<">";
        os << ",<"<<V[(i+1)%N](0)<<","<<V[(i+1)%N](1)<<","<<V[(i+1)%N](2)<<">";
        os << ",<"<<mid(0)<<","<<mid(1)<<","<<mid(2)<<">";
        os <<"\n pigment { color "<<Color<<" } \n } \n ";
    }
}

/////////////////////////////////////////////////////////////////////////////////////////// Blender /////


inline void BPYHeader (std::ostream & os)
{
    os << "from Blender import *\n";
    os << "from bpy import *\n";
    os << "from Blender import Mathutils\n";
    os << "from Blender.Mathutils import *\n";
    os << "s = data.scenes.active\n";
}

inline void BPYDrawVert (Vec3_t const & V, std::ostream & os, double Radius=1.0)
{
    os << "m = Mesh.Primitives.UVsphere(32,32,"<<Radius*2.0<<")\n";
    os << "o = s.objects.new(m,'Sphere')\n";
    os << "o.setLocation("<<V(0)<<","<<V(1)<<","<<V(2)<<")\n";
}

inline void BPYDrawPolygon (Array<Vec3_t> const & V, std::ostream & os)
{
    size_t N = V.Size();
    Vec3_t mid;
    mid = 0.0, 0.0, 0.0;
    for (size_t i=0; i<V.Size(); i++) mid += V[i];
    mid /= V.Size();
    for (size_t i=0; i<N; i++) 
    {   
        os << "m = data.meshes.new('Face') \no = s.objects.new(m,'Face') \nv = [["<<V[i](0)<<","<<V[i](1)<<","<<V[i](2)<<"],["<<V[(i+1)%N](0)<<","<<V[(i+1)%N](1)<<","<<V[(i+1)%N](2)<<"],["<<mid(0)<<","<<mid(1)<<","<<mid(2)<<"]]\n";
        os << "f = [[0,1,2]] \nm.verts.extend(v) \nm.faces.extend(f) \n";
    }
}

#endif // MECHSYS_DEM_GRAPH_H
