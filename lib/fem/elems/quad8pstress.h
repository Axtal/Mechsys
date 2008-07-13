/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

#ifndef MECHSYS_FEM_QUAD8PSTRESS_H
#define MECHSYS_FEM_QUAD8PSTRESS_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/quad8.h"

namespace FEM
{

class Quad8PStress : public Quad8, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int _geom() const { return 5;} ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

}; // class Quad8PStress

// Quad8PStress constants
char const * Quad8PStress::NAME = "Quad8PStress";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Quad8PStress element
Element * Quad8PStressMaker()
{
	return new Quad8PStress();
}

// Register a Quad8PStress element into ElementFactory array map
int Quad8PStressRegister()
{
	ElementFactory[Quad8PStress::NAME] = Quad8PStressMaker;
	return 0;
}

// Execute the autoregistration
int __Quad8PStress_dummy_int  = Quad8PStressRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD8PSTRESS_H
