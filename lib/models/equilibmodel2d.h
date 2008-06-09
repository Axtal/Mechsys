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

#ifndef MECHSYS_EQUILIBMODEL2D_H
#define MECHSYS_EQUILIBMODEL2D_H

// MechSys
#include "util/string.h"
#include "util/util.h"
#include "linalg/laexpr.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class EquilibModel2D
{
public:
	// Destructor
	virtual ~EquilibModel2D () {}

	// Parameters and initial values
	virtual void SetPrms (String const & Prms) =0;
	virtual void SetInis (String const & Inis) =0;

	// Derived Methods
	virtual void TgStiffness  (Matrix<double> & D) const =0;
	virtual int  StressUpdate (Vector<double> const & DEps, Vector<double> & DSig) =0;
	virtual void BackupState  () =0;
	virtual void RestoreState () =0;

	// Access Methods
	virtual void Sig (Vector<double> & Stress ) const =0;
	virtual void Eps (Vector<double> & Strain ) const =0;
	virtual void Ivs (Array<double>  & IntVals) const =0;

}; // class EquilibModel2D


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new EquilibModel2D
typedef EquilibModel2D * (*EquilibModel2DMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes equilibmodels
typedef std::map<String, EquilibModel2DMakerPtr, std::less<String> > EquilibModel2DFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes equilibmodels
EquilibModel2DFactory_t EquilibModel2DFactory;

// Allocate a new equilibmodel according to a string giving the name of the equilibmodel
EquilibModel2D * AllocEquilibModel2D(String const & Name)
{
	// Check if there is Name model implemented
	EquilibModel2DMakerPtr ptr=NULL;
	ptr = EquilibModel2DFactory[Name];
	if (ptr==NULL)
		throw new Fatal(_("FEM::AllocEquilibModel2D: There is no < %s > implemented in this library"), Name.GetSTL().c_str());

	return (*ptr)();
}


#endif // MECHSYS_EQUILIBMODEL2D_H
