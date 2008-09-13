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

#ifndef MECHSYS_DIFFUSIONMODEL_H
#define MECHSYS_DIFFUSIONMODEL_H

// MechSys
#include "models/model.h"
#include "util/string.h"
#include "util/util.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class DiffusionModel : public Model
{
public:
	// Destructor
	virtual ~DiffusionModel () {}

	// Methods
	void SetGeom      (int Type) { _geom = Type; } ///< Geometry type:  1:1D, 2:2D, 3:3D
	void BackupState  () {}
	void RestoreState () {}

	// Methods that MUST be derived
	virtual void TgConductivity (Matrix<double> & D) const =0;

	// Access methods
	double Val (char const * Name  ) const;

protected:
	// Private methods that MUST be derived
	virtual double _val (char const * Name) const =0; ///< Return internal values

}; // class DiffusionModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline double DiffusionModel::Val(char const * Name) const
{
	return _val(Name);
}

#endif // MECHSYS_DIFFUSIONMODEL_H
