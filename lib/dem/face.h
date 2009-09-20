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

#ifndef MECHSYS_DEM_FACE_H
#define MECHSYS_DEM_FACE_H

// MechSys
#include "dem/edge.h"
#include "dem/graph.h"
#include "util/array.h"

class Face
{
public:
    // Constructor
    Face (Array<Vec3_t> const & V); ///< V: vertices of face
    Face () {};

    // Destructor
    ~Face ();

    // Methods
    void Rotate    (Quaternion_t const & Q, Vec3_t const & Xa); ///< Q: quaternion representing the rotation, Xa: position of the axis of rotation
    void Translate (Vec3_t const & dX);                         ///< Translate edge by dX
    void Draw      (std::ostream & os, double Radius=1.0, char const * Color="Blue", bool BPY=false);

    // Data
    Array<Edge*> Edges; ///< Edges
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Face::Face (Array<Vec3_t> const & V)
{
    if (V.Size()<3) throw new Fatal("Face::Face: Number of vertices must be greater than 2");
    Edges.Resize (V.Size());
    for (size_t i=0; i<Edges.Size(); i++) 
    {
        if (i==V.Size()-1) Edges[i] = new Edge (V[i], V[0]);
        else               Edges[i] = new Edge (V[i], V[i+1]);
    }
}

inline Face::~Face ()
{
    for (size_t i=0; i<Edges.Size(); i++) delete Edges[i];
}

inline void Face::Rotate (Quaternion_t const & Q, Vec3_t const & Xa)
{
    for (size_t i=0; i<Edges.Size(); i++) Edges[i]->Rotate (Q,Xa);
}

inline void Face::Translate (Vec3_t const & dX)
{
    for (size_t i=0; i<Edges.Size(); i++) Edges[i]->Translate (dX);
}

inline void Face::Draw (std::ostream & os, double Radius, char const * Color, bool BPY)
{
    Array<Vec3_t> vi, vs; // two "sandwich" faces due to the spheroradius (i:inferior, s:superior)
    Vec3_t n = cross(Edges[0]->dL, Edges[1]->dL);
    n = n/norm(n);
    for (size_t i=0; i<Edges.Size(); i++)
    {
        vi.Push(Edges[i]->X0 - Radius*n);
        vs.Push(Edges[i]->X0 + Radius*n);
    }
    if (BPY)
    {
        BPYDrawPolygon (vi,os);
        BPYDrawPolygon (vs,os);
    }
    else
    {
        POVDrawPolygon (vi,os,Color);
        POVDrawPolygon (vs,os,Color);
    }
}

#endif // MECHSYS_DEM_FACE_H
