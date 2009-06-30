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

#ifndef MECHSYS_DEM3D_VEC3_H
#define MECHSYS_DEM3D_VEC3_H

#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

namespace DEM3D
{

typedef blitz::TinyVector<double,3> Vec3_t;

double dot(Vec3_t & a,Vec3_t & b)
{
	return a(0)*b(0)+a(1)*b(1)+a(2)*b(2);
}

double norm(Vec3_t & a) 
{
	return sqrt(dot(a,a));
}

void cross(Vec3_t & a,Vec3_t & b,Vec3_t & c) 
{
	c(0) = a(1)*b(2)-a(2)*b(1);
	c(1) = a(2)*b(0)-a(0)*b(2);
	c(2) = a(0)*b(1)-a(1)*b(0);
}

typedef blitz::TinyVector<double,4> Quaternion;

double norm(Quaternion & a)
{
	return sqrt(a(0)*a(0)+a(1)*a(1)+a(2)*a(2)+a(3)*a(3));
}

void normalize(Quaternion & a) 
{
	a = a/norm(a);
}

void normalize_rotation(double theta,Vec3_t & axis,Quaternion & c)
{
	Vec3_t a = axis/norm(axis);
	c(0)     = cos(theta/2);
	c(1)     = a(0)*sin(theta/2);
	c(2)     = a(1)*sin(theta/2);
	c(3)     = a(2)*sin(theta/2);
}

void conjugate(Quaternion & a,Quaternion & c)
{
	c(0) = a(0);
	c(1) = -a(1);
	c(2) = -a(2);
	c(3) = -a(3);
}

void getvector(Quaternion & a,Vec3_t & c)
{
	c(0) = a(1);
	c(1) = a(2);
	c(2) = a(3);
}

double getscalar(Quaternion & a) {
	return a(0);
}

void setQ(double scalar,Vec3_t & a,Quaternion & c)
{
	c(0) = scalar;
	c(1) = a(0);
	c(2) = a(1);
	c(3) = a(2);
}

void Qproduct(Quaternion & a,Quaternion & b,Quaternion & c)
{
	Vec3_t t1,t2;
	getvector(a,t1);
	getvector(b,t2);
	double scalar = a(0)*b(0)-dot(t1,t2);
	Vec3_t cp;
	cross(t1,t2,cp);
	Vec3_t vect = a(0)*t2+b(0)*t1+cp;
	setQ(scalar,vect,c);
}

void rotate(Vec3_t & a,Quaternion & b,Vec3_t & c)
{
	Quaternion t1,t2,t3;
	setQ(0,a,t1);
	Qproduct(b,t1,t2);
	conjugate(b,t3);
	Qproduct(t2,t3,t1);
	getvector(t1,c);
}



};

#endif


