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
#include <mechsys/dem/edge.h>
#include <mechsys/dem/graph.h>
#include <mechsys/util/array.h>

class Face
{
public:
    // Constructor
    Face (Array<Edge *> E); ///< E: Edges of the face
    Face (Array<Vec3_t> & V);
    Face (Array<Vec3_t*> & V);
    Face () {};

    // Destructor
    ~Face ();


    // Methods
    void UpdatedL (); ///< UdatedL for each edge
    void Draw      (std::ostream & os, double Radius=1.0, char const * Color="Blue", bool BPY=false);

    // Data
    Array<Edge*> Edges; ///< Edges
    bool         Allocate; ///< It allocates memory or not
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Face::Face (Array<Edge *> E)
{
    Edges = E;
    Allocate = false;
}

inline Face::Face (Array<Vec3_t> & V)
{
    for (size_t i = 0; i < V.Size() ; i++)
    {
        Edges.Push(new Edge(&V[i],&V[(i+1)%V.Size()]));
    }
    Allocate = true;
}

inline Face::Face(Array<Vec3_t*> & V)
{
    for (size_t i = 0; i < V.Size() ; i++)
    {
        Edges.Push(new Edge(V[i],V[(i+1)%V.Size()]));
    }
    Allocate = true;
}

inline Face::~Face ()
{
    if (Allocate) 
    {
        for (size_t i = 0; i<Edges.Size();i++)
        {
            delete Edges[i];
        }
    }
}

inline void Face::UpdatedL()
{
    for (size_t i = 0; i<Edges.Size();i++)
    {
        Edges[i]->UpdatedL();
    }
}

inline void Face::Draw (std::ostream & os, double Radius, char const * Color, bool BPY)
{
    Array<Vec3_t> vi, vs; // two "sandwich" faces due to the spheroradius (i:inferior, s:superior)
    Vec3_t n = cross(Edges[0]->dL, Edges[1]->dL);
    n = n/norm(n);
    for (size_t i=0; i<Edges.Size(); i++)
    {
        vi.Push(*Edges[i]->X0 - Radius*n);
        vs.Push(*Edges[i]->X0 + Radius*n);
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
