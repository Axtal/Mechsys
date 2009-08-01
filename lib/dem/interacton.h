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

#ifndef DEM_INTERACTON_H
#define DEM_INTERACTON_H

#include <math.h>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "dem/particle.h"

class Interacton
{
public:
	// Constructor and destructor
	 Interacton (Particle * Pt1, Particle * Pt2); ///< Constructor, it requires pointers to both particles
	~Interacton ();                           ///< Destructor

	// Methods
	void CalcForce (double Dt); ///< Calculates the contact force between particles

protected:
	double     _Kn;  ///< Normal stiffness
	double     _Kt;  ///< Tengential stiffness
	double     _gn;  ///< Normal viscous coefficient
	double     _gt;  ///< Tangential viscous coefficient
	double     _mu;  ///< Microscpic coefficient of friction
	double     _mur; ///< Rolling resistance coefficient
	Particle * _p1;  ///< First particle
	Particle * _p2;  ///< Second particle
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Interacton::Interacton(Particle * Pt1, Particle * Pt2)
	: _p1(Pt2), _p2(Pt2)
{
}

inline Interacton::~Interacton()
{
}

inline void Interacton::CalcForce(double Dt)
{
}


#endif // DEM_INTERACTON_H