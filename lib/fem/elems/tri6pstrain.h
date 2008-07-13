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

#ifndef MECHSYS_FEM_TRI6PSTRAIN_H
#define MECHSYS_FEM_TRI6PSTRAIN_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/tri6.h"

namespace FEM
{

class Tri6PStrain : public Tri6, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

private:
	// Private methods
	int _geom() const { return 2;} ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

}; // class Tri6PStrain

// Tri6PStrain constants
char const * Tri6PStrain::NAME = "Tri6PStrain";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri6PStrain element
Element * Tri6PStrainMaker()
{
	return new Tri6PStrain();
}

// Register a Tri6PStrain element into ElementFactory array map
int Tri6PStrainRegister()
{
	ElementFactory[Tri6PStrain::NAME] = Tri6PStrainMaker;
	return 0;
}

// Execute the autoregistration
int __Tri6PStrain_dummy_int  = Tri6PStrainRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI6PSTRAIN_H
