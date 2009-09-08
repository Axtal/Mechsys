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

#ifndef MECHSYS_DEM_DISTANCE_H
#define MECHSYS_DEM_DISTANCE_H

// MechSys
#include "dem/edge.h"
#include "dem/face.h"

inline void Distance (Vec3_t const & V, Edge const & E, Vec3_t & Xi, Vec3_t & Xf)
{
    double t = (dot(V,E.dL)-dot(E.X0,E.dL))/(dot(E.dL,E.dL));
    Xi = V;
    if      (t<0) Xf = E.X0;
    else if (t>1) Xf = E.X1;
    else          Xf = E.X0 + E.dL*t;
}

inline void Distance (Edge const & E, Vec3_t const & V, Vec3_t & Xi, Vec3_t & Xf)
{
    Distance (V,E,Xf,Xi);
}

inline void Distance (Edge const & E0, Edge const & E1, Vec3_t & Xi, Vec3_t & Xf)
{
    double a = dot(E0.dL, E0.X0-E1.X0);
    double b = dot(E1.dL, E0.X0-E1.X0);
    double c = dot(E0.dL, E0.dL);
    double d = dot(E1.dL, E1.dL);
    double e = dot(E0.dL, E1.dL);
    double t = (c*b-e*a)/(c*d-e*e);
    double s = (e*b-a*d)/(c*d-e*e);
    
    if ((s>0) && (s<1) && (t>0) && (t<1)) 
    {
        Xi = E0.X0+E0.dL*s;
        Xf = E1.X0+E1.dL*t;
    }
    else 
    {
        Vec3_t xi1,xf1;
        Vec3_t xi2,xf2;
        Vec3_t xi3,xf3;
        Vec3_t xi4,xf4;
        Distance (E0.X0,E1,xi1,xf1);
        Distance (E0.X1,E1,xi2,xf2);
        Distance (E1.X0,E0,xf3,xi3);
        Distance (E1.X1,E0,xf4,xi4);
        double l1 = norm(xf1-xi1);
        double l2 = norm(xf2-xi2);
        double l3 = norm(xf3-xi3);
        double l4 = norm(xf4-xi4);
        if ((l1<=l2) && (l1<=l3) && (l1<=l4))
        {   
            Xi = xi1;
            Xf = xf1;
        }
        if ((l2<=l1) && (l2<=l3) && (l2<=l4)) 
        {   
            Xi = xi2;
            Xf = xf2;
        }
        if ((l3<=l1) && (l3<=l2) && (l3<=l4)) 
        {   
            Xi = xi3;
            Xf = xf3;
        }
        if ((l4<=l1) && (l4<=l2) && (l4<=l3)) 
        {   
            Xi = xi4;
            Xf = xf4;
        }
    }
}

inline void Distance (Vec3_t const & V, Face const & F, Vec3_t & Xi, Vec3_t & Xf)
{
    // find the normal to face
    Vec3_t nor = cross(F.Edges[0]->dL, F.Edges[1]->dL);
    nor = nor/norm(nor);

    // find the projection
    double a   = dot(F.Edges[0]->X0-V, F.Edges[0]->dL);
    double b   = dot(F.Edges[0]->dL,   F.Edges[0]->dL);
    double c   = dot(F.Edges[0]->dL,   F.Edges[1]->dL);
    double d   = dot(F.Edges[0]->X0-V, F.Edges[1]->dL);
    double f   = dot(F.Edges[1]->dL,   F.Edges[1]->dL);
    double s   = (c*d-a*f)/(b*f-c*c);
    double t   = (a*c-b*d)/(b*f-c*c);
    Vec3_t pro = F.Edges[0]->X0 + s*F.Edges[0]->dL + t*F.Edges[1]->dL;

    // check if vertex is inside
    bool inside = true;
    for (size_t i=0; i<F.Edges.Size(); i++) 
    {
        Vec3_t tmp = pro-F.Edges[i]->X0;
        if (dot(cross(F.Edges[i]->dL,tmp),nor)<0) inside = false;
    }

    // get Xi, Xf
    if (inside)
    {
        Xi = V;
        Xf = pro;
    }
    else // compare V against each edge
    {
        Distance (V, (*F.Edges[0]), pro, nor); // pro=Xi, nor=Xf
        double lt = norm(nor-pro);
        double ld = lt;
        Xi = pro;
        Xf = nor;
        for (size_t i=1; i<F.Edges.Size(); i++) 
        {
            Distance (V, (*F.Edges[i]), pro, nor);
            lt = norm(nor-pro);
            if (lt<ld)
            {
                Xi = pro;
                Xf = nor;
                ld = lt;
            }
        }
    }
}

inline void Distance (Face const & F, Vec3_t const & V, Vec3_t & Xi, Vec3_t & Xf)
{
    Distance (V,F,Xf,Xi);
}

inline double Distance (Edge const & E, Vec3_t const & V)
{
    Vec3_t Xi,Xf;
    Distance (E,V,Xi,Xf);
    return norm(Xf-Xi);
}

inline double Distance (Vec3_t const & V, Edge const & E)
{
    return Distance (E,V);
}

inline double Distance (Edge const & E0, Edge const & E1)
{
    Vec3_t Xi,Xf;
    Distance (E0,E1,Xi,Xf);
    return norm(Xf-Xi);
}

inline double Distance (Face const & F, Vec3_t const & V)
{
    Vec3_t Xi,Xf;
    Distance (F,V,Xi,Xf);
    return norm(Xf-Xi);
}

inline double Distance (Vec3_t const & V, Face const & F)
{
    return Distance (F,V);
}

inline double Distance (Vec3_t const & V0, Vec3_t const & V1)
{
    return norm(V1-V0);
}

#endif // MECHSYS_DEM_DISTANCE_H
