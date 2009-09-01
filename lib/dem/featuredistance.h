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
 * but WriTHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABriLriTY or FriTNESS FOR A PARTriCULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. rif not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

#ifndef DEM_FEATUREDISTANCE_H
#define DEM_FEATUREDISTANCE_H

// Std lib
#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/face.h"

inline void Distance(const Vec3_t & V,Edge & E,Vec3_t & xi,Vec3_t & xf) ///< Distance Between point V and Edge E given by the contact points
{
    double t;
    t=(dot(V,E.dr())-dot(E.ri(),E.dr()))/(dot(E.dr(),E.dr()));
    xi=V;
    if (t<0) xf=E.ri();
    else if (t>1) xf=E.rf();
    else xf=E.ri()+E.dr()*t;
}

inline void Distance(Edge & E,const Vec3_t & V,Vec3_t & xi,Vec3_t & xf) ///< Distance Between point V and Edge E given by the contact points
{
    Distance(V,E,xf,xi);
}

inline double Distance(Edge & E,const Vec3_t & V) ///< Distance Between point V and Edge E given by a scalar quantity
{
    Vec3_t xi,xf;
    Distance(E,V,xi,xf);
    return norm(xf-xi);
}

inline double Distance(const Vec3_t & V,Edge & E) ///< Distance Between point V and Edge E given by a scalar quantity
{
    return Distance(E,V);
}

inline void Distance(Edge & E0,Edge & E1,Vec3_t & xi,Vec3_t & xf) ///< Distance Between two edges given by the contact points
{
    double t,s,a,b,c,d,e;
    Vec3_t xi1,xf1;
    Vec3_t xi2,xf2;
    Vec3_t xi3,xf3;
    Vec3_t xi4,xf4;
    c = dot(E0.dr(),E0.dr());
    e = dot(E0.dr(),E1.dr());
    d = dot(E1.dr(),E1.dr());
    a = dot(E0.dr(),E0.ri()-E1.ri());
    b = dot(E1.dr(),E0.ri()-E1.ri());
    t = (c*b-e*a)/(c*d-e*e);
    s = (e*b-a*d)/(c*d-e*e);
    
    if ((s>0)&&(s<1)&&(t>0)&&(t<1)) 
    {
        xi = E0.ri()+E0.dr()*s;
        xf = E1.ri()+E1.dr()*t;
    }
    else 
    {
        Distance(E0.ri(),E1,xi1,xf1);
        double l1 = norm(xf1-xi1);
        Distance(E0.rf(),E1,xi2,xf2);
        double l2 = norm(xf2-xi2);
        Distance(E1.ri(),E0,xi3,xf3);
        double l3 = norm(xf3-xi3);
        Distance(E1.rf(),E0,xi4,xf4);
        double l4 = norm(xf4-xi4);
        if((l1<=l2)&&(l1<=l3)&&(l1<=l4))
        {   
            xi=xi1;
            xf=xf1;
        }
        if((l2<=l1)&&(l2<=l3)&&(l2<=l4)) 
        {   
            xi=xi2;
            xf=xf2;
        }
        if((l3<=l1)&&(l3<=l2)&&(l3<=l4)) 
        {   
            xi=xi3;
            xf=xf3;
        }
        if((l4<=l1)&&(l4<=l2)&&(l4<=l3)) 
        {   
            xi=xi4;
            xf=xf4;
        }
    }
}

inline double Distance(Edge & E0,Edge & E1) ///< Distance Between two edges given by the scalar quantity
{
    Vec3_t xi,xf;
    Distance(E0,E1,xi,xf);
    return norm(xf-xi);
}

inline void Distance(const Vec3_t & v, Face & F,Vec3_t & xi,Vec3_t & xf) ///< Distance Between a Point and a Face given by the contact points
{
    Vec3_t pro,nor;
    bool inside=true;
    double s,t,a,b,c,d,f,lt,ld;
    size_t i,ns=F.NumberofSides();
    nor=cross(F.Edges(0)->dr(),F.Edges(1)->dr());
    nor=nor/norm(nor);
    a=dot(F.Edges(0)->ri()-v,F.Edges(0)->dr());
    b=dot(F.Edges(0)->dr(),F.Edges(0)->dr());
    c=dot(F.Edges(0)->dr(),F.Edges(1)->dr());
    d=dot(F.Edges(0)->ri()-v,F.Edges(1)->dr());
    f=dot(F.Edges(1)->dr(),F.Edges(1)->dr());
    s=(c*d-a*f)/(b*f-c*c);
    t=(a*c-b*d)/(b*f-c*c);
    pro=F.Edges(0)->ri()+s*F.Edges(0)->dr()+t*F.Edges(1)->dr();
    for(i=0;i<ns;i++) 
    {
        Vec3_t tmp = pro-F.Edges(i)->ri();
        if (dot(cross(F.Edges(i)->dr(),tmp),nor)<0) inside=false;
    }
    if (inside) 
    {
        xi = v;
        xf = pro;
    }
    else 
    {
        Distance(v,*F.Edges(0),pro,nor);
        lt=norm(nor-pro);
        ld=lt;
        xi=pro;
        xf=nor;
        for(i=1;i<ns;i++) 
        {
            Distance(v,*F.Edges(i),pro,nor);
            lt=norm(nor-pro);
            if(lt<ld){
                xi=pro;
                xf=nor;
                ld=lt;
            }
        }
    }
}

inline void Distance(Face & F,const Vec3_t & v,Vec3_t & xi,Vec3_t & xf) ///< Distance Between a Point and a Face given by the contact points
{
    Distance(v,F,xf,xi);
}

inline double Distance(Face & F,const Vec3_t & v) ///< Distance Between a Point and a Face given by the scalar quantity
{
    Vec3_t xi,xf;
    Distance(F,v,xi,xf);
    return norm(xf-xi);
}

inline double Distance(const Vec3_t & v,Face & F) ///< Distance Between a Point and a Face given by the scalar quantity
{
    return Distance(F,v);
}

inline double Distance(const Vec3_t & a,const Vec3_t & b) ///< Distance between two points
{
    return norm(b-a);
}








#endif //DEM_FEATUREDISTANCE_H
