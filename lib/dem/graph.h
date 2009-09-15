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
#include "util/array.h"
#include "linalg/matvec.h"


/////////////////////////////////////////////////////////////////////////////////////////// PovRay /////


inline void PovHeader (std::ostream & os)
{
    os << "#include \"colors.inc\" \n";
    os << "background {color White} \n";
    os << "light_source{<10,0,0> color White shadowless}  \n";
    os << "light_source{<-10,0,0> color White shadowless}  \n";
    os << "light_source{<0,10,0> color White shadowless}  \n";
    os << "light_source{<0,-10,0> color White shadowless}  \n";
    os << "light_source{<0,0,10> color White shadowless}  \n";
    os << "light_source{<0,0,-10> color White shadowless}  \n";
}   

inline void PovSetCam (std::ostream & os, const Vec3_t & X, const Vec3_t & F)
{
    os << "camera { location <"<<X(0)<<","<<X(1)<<","<<X(2)<<"> sky <0,0,1> look_at <"<<F(0)<<","<<F(1)<<","<<F(2)<<"> }\n";
}

inline void PovDrawVert (Vec3_t const & V, std::ostream & os, double Radius=1.0, char const * Color="Blue")
{
    os << "sphere  { <"<<V(0)<<","<<V(1)<<","<<V(2)<<">,"<<Radius<<"\n pigment { color "<<Color<<" } }\n";
}

inline void PovDrawPolygon (Array<Vec3_t> const & V, std::ostream & os, char const * Color="Blue")
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
        os <<"\n pigment { color "<<Color<<" } }\n";
    }
}

/////////////////////////////////////////////////////////////////////////////////////////// Blender /////


inline void BlenderHeader (std::ostream & os)
{
    os << "from Blender import *\n";
    os << "from bpy import *\n";
    os << "from Blender import Mathutils\n";
    os << "from Blender.Mathutils import *\n";
    os << "s = data.scenes.active\n";
}

inline void BlenderDrawVert (Vec3_t const & V, std::ostream & os, double Radius=1.0)
{
    os << "m = Mesh.Primitives.UVsphere(32,32,"<<Radius*2.0<<")\n";
    os << "o = s.objects.new(m,'Sphere')\n";
    os << "o.setLocation("<<V(0)<<","<<V(1)<<","<<V(2)<<")\n";
}

inline void BlenderDrawPolygon (Array<Vec3_t> const & V, std::ostream & os)
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

/*
class Graph
{
public:
    // Constructor
    Graph(char const * FileName,           ///< Filename
          bool IsPovray=true);     ///< Flag for the user to choose between Povray or Blender for visualization

    // Methods
    void SetCamera (const Vec3_t & x,const Vec3_t & v);        ///< Put the Camera at position x looking at point v
    void DrawPoint (const Vec3_t & r,double R,char const *c);  ///< Draw a single point as a sphere at position r with radius R and color c
    void DrawEdge ( Edge & E,double R,char const *c);          ///< Draw a edge as a cylinder with radius R and color c
    void DrawFace ( Face & F,double R,char const *c);          ///< Draw a face as a polygonal face with radius R and color c
    void DrawPolygon (const Vec3_t *v,size_t N,char const *c); ///< Draw a polygon with N sides and color c whose vertices are stored in v
    void DrawParticle (Particle & P,char const *c);            ///< Draw the entire particle with all its geometric features
    void DrawEntireDomain (Domain & D,char const *c);          ///< Draw the entire Domain with all its particles
    void Close ();                                             ///< Flushes the working string into the final file

    void WritePov (char const * Filename);




protected:
    String _fn;              ///< The exporting File Name
    std::ostringstream _oss; ///< Working String
    bool _IsPovray;      ///< Povray-Blender Flag
};



/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void Graph::DrawFace (Face & F,double R,char const *c)
{
    Vec3_t * vs;
    Vec3_t * vi;
    Vec3_t n = cross(F.Edges(0)->dr(),F.Edges(1)->dr());
    n=n/norm(n);
    size_t ns = F.NumberofSides();
    vs = new Vec3_t [ns];
    vi = new Vec3_t [ns];
    for(size_t i=0;i<ns;i++) {
        vi[i]=F.Edges(i%ns)->ri()-R*n;
            vs[i]=F.Edges(i%ns)->ri()+R*n;
    }
    DrawPolygon(vi,ns,c);
        DrawPolygon(vs,ns,c);
    delete [] vs;
    delete [] vi;
}

inline void Graph::DrawParticle (Particle & P, const char *c)
{
    size_t nv = P.NumberVertices(),ne = P.NumberEdges(),nf = P.NumberFaces();
    for (size_t i = 0; i < nv; i++)
    {
        DrawPoint(*P.Vertex(i),P.Radius(),c);
    }

    for (size_t i = 0; i < ne; i++)
    {
        DrawEdge(*P.Edges(i),P.Radius(),c);
    }

    for (size_t i = 0; i < nf; i++)
    {
        DrawFace(*P.Faces(i),P.Radius(),c);
    }
}

inline void Graph::DrawEntireDomain(Domain & D,char const *c)
{
    for (size_t i = 0;i<D.NumberParticles();i++)
    {
        DrawParticle(*D.Particles(i),c);
    }
} 

inline void Graph::Close ()
{
    std::ofstream file(_fn.CStr());
    file << _oss.str();
    file.close();
}
*/

#endif // MECHSYS_DEM_GRAPH_H
