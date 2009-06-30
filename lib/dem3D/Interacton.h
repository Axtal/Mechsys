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

#ifndef MECHSYS_DEM3D_INTERACTON_H
#define MECHSYS_DEM3D_INTERACTON_H

#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

#include "dem3D/sphere.h"

using namespace DEM3D;

class Interacton
{
public:
	Interacton(Sphere &p1,Sphere &p2); ///< Constructor, it requires pointers to both particles
	~Interacton ();                    ///< Destructor
	void calcForce(double dt);         ///< Calculates the contact force between particles
				
protected:
	double _Kn;     ///< Normal stiffness
	double _Kt;     ///< Tengential stiffness
	double _gn;     ///< Normal viscous coefficient
	double _gt;     ///< Tangential viscous coefficient
	double _mu;     ///< Microscpic coefficient of friction
	double _mur;    ///< Rolling resistance coefficient
	Sphere *p1,*p2; ///< Pointers to the pair of interacting particles
};

#endif

