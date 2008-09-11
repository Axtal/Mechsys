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

#ifndef MECHSYS_HEATMODEL_H
#define MECHSYS_HEATMODEL_H

// MechSys
#include "models/model.h"
#include "util/string.h"
#include "util/util.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class HeatModel : public Model
{
public:
	// Destructor
	virtual ~HeatModel () {}

	// Methods
	void SetGeom      (int Type) { _geom = Type; } ///< Geometry type:  1:1D, 2:2D, 3:3D
	void BackupState  () { _u_bkp=_u; _f_bkp=_f; }
	void RestoreState () { _u=_u_bkp; _f=_f_bkp; }

	// Methods that MUST be derived
	virtual void TgConductivity (Matrix<double> & C) const =0;
	//virtual int  UpdateState    (Vector<double> const & DEps, Vector<double> & DSig) =0;

	// Access methods
	double Val (char const * Name  ) const;

protected:
	// Data
	double _u;     ///< primary variable == Temperature
	double _f;     ///< source variable == Heat source
	double _u_bkp; ///< backup: primary variable
	double _f_bkp; ///< backup: source variable

	// Private methods that MUST be derived
	virtual double _val  (char const * Name) const =0; ///< Return internal values
	// Private methods that MUST be derived
	//virtual void   _cond (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const =0; ///< Tangent or secant stiffness

}; // class HeatModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline double HeatModel::Val(char const * Name) const
{
	     if (strcmp(Name,"u" )==0) return _u;
	else if (strcmp(Name,"f" )==0) return _f;
	else                           return _val(Name);
}

#endif // MECHSYS_HEATMODEL_H
