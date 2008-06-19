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

#ifndef MECHSYS_FEM_TRI6PSTRESS_H
#define MECHSYS_FEM_TRI6PSTRESS_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/tri6.h"

namespace FEM
{

class Tri6PStress : public Tri6, public EquilibElem
{
public:
	// Constants
	static String NAME;

	// Derived methods
	String Name() const { return NAME; };

private:
	// Private methods
	int _geom() const { return 5;} ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

}; // class Tri6PStress

// Tri6PStress constants
String Tri6PStress::NAME = "Tri6PStress";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri6PStress element
Element * Tri6PStressMaker()
{
	return new Tri6PStress();
}

// Register a Tri6PStress element into ElementFactory array map
int Tri6PStressRegister()
{
	ElementFactory[Tri6PStress::NAME] = Tri6PStressMaker;
	return 0;
}

// Execute the autoregistration
int __Tri6PStress_dummy_int  = Tri6PStressRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI6PSTRESS_H
