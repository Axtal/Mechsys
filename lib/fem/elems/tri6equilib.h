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

#ifndef MECHSYS_FEM_TRI6EQUILIB_H
#define MECHSYS_FEM_TRI6EQUILIB_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/tri6.h"

namespace FEM
{

class Tri6Equilib : public Tri6, public EquilibElem
{
public:
	// Constants
	static String NAME;

	// Derived methods
	String Name() const { return NAME; };

private:
}; // class Tri6Equilib

// Tri6Equilib constants
String Tri6Equilib::NAME = "Tri6Equilib";

///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Tri6Equilib element
Element * Tri6EquilibMaker()
{
	return new Tri6Equilib();
}

// Register a Tri6Equilib element into ElementFactory array map
int Tri6EquilibRegister()
{
	ElementFactory[Tri6Equilib::NAME] = Tri6EquilibMaker;
	return 0;
}

// Execute the autoregistration
int __Tri6Equilib_dummy_int  = Tri6EquilibRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI6EQUILIB_H
