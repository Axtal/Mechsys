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

#ifndef MECHSYS_FEM_QUAD8PSTRAIN_H
#define MECHSYS_FEM_QUAD8PSTRAIN_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/quad8.h"

namespace FEM
{

class Quad8PStrain : public Quad8, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int _geom() const { return 2;} ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

}; // class Quad8PStrain

// Quad8PStrain constants
char const * Quad8PStrain::NAME = "Quad8PStrain";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Quad8PStrain element
Element * Quad8PStrainMaker()
{
	return new Quad8PStrain();
}

// Register a Quad8PStrain element into ElementFactory array map
int Quad8PStrainRegister()
{
	ElementFactory[Quad8PStrain::NAME] = Quad8PStrainMaker;
	return 0;
}

// Execute the autoregistration
int __Quad8PStrain_dummy_int  = Quad8PStrainRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD8PSTRAIN_H
