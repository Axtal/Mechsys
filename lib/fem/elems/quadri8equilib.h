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

#ifndef MECHSYS_FEM_QUADRI8EQUILIB_H
#define MECHSYS_FEM_QUADRI8EQUILIB_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/quadri8.h"

namespace FEM
{

class Quadri8Equilib : public Quadri8, public EquilibElem
{
public:
	// Constants
	static String NAME;

	// Derived methods
	String Name() const { return NAME; };

private:
}; // class Quadri8Equilib

// Quadri8Equilib constants
String Quadri8Equilib::NAME = "Quadri8Equilib";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Quadri8Equilib element
Element * Quadri8EquilibMaker()
{
	return new Quadri8Equilib();
}

// Register a Quadri8Equilib element into ElementFactory array map
int Quadri8EquilibRegister()
{
	ElementFactory[Quadri8Equilib::NAME] = Quadri8EquilibMaker;
	return 0;
}

// Execute the autoregistration
int __Quadri8Equilib_dummy_int  = Quadri8EquilibRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_QUADRI8EQUILIB_H
