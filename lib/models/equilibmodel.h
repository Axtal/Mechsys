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
#include "models/model.h"
#include "util/string.h"
#include "util/util.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class EquilibModel : public Model
{
public:
	// Constructor
	EquilibModel () : _geom(-1) {}

	// Destructor
	virtual ~EquilibModel () {}

	// Set geometry type
	void SetGeom (int Type) { _geom=Type; } ///< Geometry type:  1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

	// Derived Methods
	virtual void SetPrms      (char const * Prms) =0;
	virtual void SetInis      (char const * Inis) =0;
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

#endif // MECHSYS_EQUILIBMODEL_H
