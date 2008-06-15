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

#ifndef MECHSYS_FEM_QUAD4EQUILIB_H
#define MECHSYS_FEM_QUAD4EQUILIB_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/quad4.h"

namespace FEM
{

class Quad4Equilib : public Quad4, public EquilibElem
{
public:
	// Constants
	static String NAME;

	// Derived methods
	String Name() const { return NAME; };

private:
}; // class Quad4Equilib

// Quad4Equilib constants
String Quad4Equilib::NAME = "Quad4Equilib";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Quad4Equilib element
Element * Quad4EquilibMaker()
{
	return new Quad4Equilib();
}

// Register a Quad4Equilib element into ElementFactory array map
int Quad4EquilibRegister()
{
	ElementFactory[Quad4Equilib::NAME] = Quad4EquilibMaker;
	return 0;
}

// Execute the autoregistration
int __Quad4Equilib_dummy_int  = Quad4EquilibRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4EQUILIB_H
