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

#ifndef MECHSYS_EQUILIBMODEL_H
#define MECHSYS_EQUILIBMODEL_H

// MechSys
#include "util/string.h"
#include "util/util.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class EquilibModel
{
public:
	// Constructor
	EquilibModel () : _geom(3) {}

	// Destructor
	virtual ~EquilibModel () {}

	// Set geometry type
	void SetGeom (int Type) { _geom=Type; } ///< Geometry type:  1:1D, 2:2D(plane-strain), 3:3D, 4:Axis-symmetric, 5:2D(plane-stress)

	// Derived Methods
	virtual void SetPrms      (String const & Prms) =0;
	virtual void SetInis      (String const & Inis) =0;
	virtual void TgStiffness  (Matrix<double> & D) const =0;
	virtual int  StressUpdate (Vector<double> const & DEps, Vector<double> & DSig) =0;
	virtual void BackupState  () =0;
	virtual void RestoreState () =0;

	// Access Methods
	virtual void Sig (Vector<double> & Stress ) const =0;
	virtual void Eps (Vector<double> & Strain ) const =0;
	virtual void Ivs (Array<double>  & IntVals) const =0;

protected:
	int _geom; ///< Geometry type:  1:1D, 2:2D(plane-strain), 3:3D, 4:Axis-symmetric, 5:2D(plane-stress)

}; // class EquilibModel


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new EquilibModel
typedef EquilibModel * (*EquilibModelMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes equilibmodels
typedef std::map<String, EquilibModelMakerPtr, std::less<String> > EquilibModelFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes equilibmodels
EquilibModelFactory_t EquilibModelFactory;

// Allocate a new equilibmodel according to a string giving the name of the equilibmodel
EquilibModel * AllocEquilibModel(String const & Name)
{
	// Check if there is Name model implemented
	EquilibModelMakerPtr ptr=NULL;
	ptr = EquilibModelFactory[Name];
	if (ptr==NULL)
		throw new Fatal(_("FEM::AllocEquilibModel: There is no < %s > implemented in this library"), Name.GetSTL().c_str());

	return (*ptr)();
}


#endif // MECHSYS_EQUILIBMODEL_H
