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

using Tensors::Tensor2;
using Tensors::Tensor2ToVector;
using Util::SQ2;
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

	// Methods to be derived
	virtual void SetPrms      (char const * Prms) =0;
	virtual void SetInis      (char const * Inis) =0;
	virtual void TgStiffness  (Matrix<double> & D) const =0;
	virtual int  StressUpdate (Vector<double> const & DEps, Vector<double> & DSig) =0;
	virtual void BackupState  () =0;
	virtual void RestoreState () =0;

	// Access Methods to be derived
	virtual void Ivs (Array<double> & IntVals) const =0;

	// Access methods
	void   Sig (Vector<double> & Stress) const { Tensor2ToVector(_geom, _sig, Stress); }
	void   Eps (Vector<double> & Strain) const { Tensor2ToVector(_geom, _eps, Strain); }
	double Val (char const     * Name  ) const;

protected:
	// Data
	int     _geom;    ///< Geometry type:  1:1D, 2:2D(plane-strain), 3:3D, 4:Axis-symmetric, 5:2D(plane-stress)
	Tensor2 _sig;     ///< Stress
	Tensor2 _eps;     ///< Strain
	Tensor2 _sig_bkp; ///< Backup stress
	Tensor2 _eps_bkp; ///< Backup strain

	// Private methods to be derived
	virtual double _val (char const * Name) const =0;

}; // class EquilibModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////
//

inline double EquilibModel::Val(char const * Name) const
{
	     if (strncmp(Name,"Sx" ,2)==0)                             return _sig(0);
	else if (strncmp(Name,"Sy" ,2)==0)                             return _sig(1);
	else if (strncmp(Name,"Sz" ,2)==0)                             return _sig(2);
	else if (strncmp(Name,"Sxy",2)==0 || strncmp(Name,"Syx",2)==0) return _sig(3)/SQ2;
	else if (strncmp(Name,"Syz",2)==0 || strncmp(Name,"Szy",2)==0) return _sig(4)/SQ2;
	else if (strncmp(Name,"Szx",2)==0 || strncmp(Name,"Sxz",2)==0) return _sig(5)/SQ2;
	else if (strncmp(Name,"Ex" ,2)==0)                             return _eps(0);
	else if (strncmp(Name,"Ey" ,2)==0)                             return _eps(1);
	else if (strncmp(Name,"Ez" ,2)==0)                             return _eps(2);
	else if (strncmp(Name,"Exy",2)==0 || strncmp(Name,"Eyx",2)==0) return _eps(3)/SQ2;
	else if (strncmp(Name,"Eyz",2)==0 || strncmp(Name,"Ezy",2)==0) return _eps(4)/SQ2;
	else if (strncmp(Name,"Ezx",2)==0 || strncmp(Name,"Exz",2)==0) return _eps(5)/SQ2;
	else                                                           return _val(Name);
}

#endif // MECHSYS_EQUILIBMODEL_H
