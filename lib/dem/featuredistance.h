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
#include "dem/face3d.h"

inline void Distance(const Vec3_t & V,Edge3D & E,Vec3_t & xi,Vec3_t & xf) ///< Distance Between point V and Edge E
{
	double t;
	t=(dot(V,E.dr())-dot(E.ri(),E.dr()))/(dot(E.dr(),E.dr()));
	xi=V;
	if (t<0) xf=E.ri();
	else if (t>1) xf=E.rf();
	else xf=E.ri()+E.dr()*t;
}

inline void Distance(Edge3D & E0,Edge3D & E1,Vec3_t & xi,Vec3_t & xf) ///< Distance Between two edges
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

inline void Distance(Vec3_t & v,Face3D & F,Vec3_t & xi,Vec3_t & xf) ///< Distance Between a Point and a Face
{
	Vec3_t pro,nor;
	bool inside=true;
	double s,t,a,b,c,d,f,lt,ld;
	size_t i,ns=F.NumberofSides();
	nor=cross(F.Edge(0)->dr(),F.Edge(1)->dr());
	nor=nor/norm(nor);
	a=dot(F.Edge(0)->ri()-v,F.Edge(0)->dr());
	b=dot(F.Edge(0)->dr(),F.Edge(0)->dr());
	c=dot(F.Edge(0)->dr(),F.Edge(1)->dr());
	d=dot(F.Edge(0)->ri()-v,F.Edge(1)->dr());
	f=dot(F.Edge(1)->dr(),F.Edge(1)->dr());
	s=(c*d-a*f)/(b*f-c*c);
	t=(a*c-b*d)/(b*f-c*c);
	pro=F.Edge(0)->ri()+s*F.Edge(0)->dr()+t*F.Edge(1)->dr();
	for(i=0;i<ns;i++) 
	{
		Vec3_t tmp = pro-F.Edge(i)->ri();
		if (dot(cross(F.Edge(i)->dr(),tmp),nor)<0) inside=false;
	}
	if (inside) 
	{
		xi = v;
		xf = pro;
	}
	
	else 
	{
		Distance(v,*F.Edge(0),pro,nor);
		lt=norm(nor-pro);
		ld=lt;
		for(i=1;i<ns;i++) 
		{
			Distance(v,*F.Edge(i),pro,nor);
			lt=norm(nor-pro);
			if(lt<ld){
				xi=pro;
				xf=nor;
				ld=lt;
			}
		}
	}
	
}





#endif //DEM_FEATUREDISTANCE_H
