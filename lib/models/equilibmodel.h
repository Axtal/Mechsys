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
	// Destructor
	virtual ~EquilibModel () {}

	// Methods
	void SetGeom (int Type) { _geom = Type; } ///< Geometry type:  1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

	// Methods that MUST be derived
	virtual void TgStiffness  (Matrix<double> & D) const =0;
	virtual int  StressUpdate (Vector<double> const & DEps, Vector<double> & DSig) =0;
	virtual void BackupState  () =0;
	virtual void RestoreState () =0;

	// Access methods
	void   Sig (Vector<double> & Stress) const { Tensor2ToVector(_geom, _sig, Stress); }
	void   Eps (Vector<double> & Strain) const { Tensor2ToVector(_geom, _eps, Strain); }
	double Val (char const     * Name  ) const;

protected:
	// Data
	int     _geom;    ///< Geometry type (in FEM must be the ELEMENT geometry):  1:1D, 2:2D(plane-strain), 3:3D, 4:Axis-symmetric, 5:2D(plane-stress)
	Tensor2 _sig;     ///< Stress
	Tensor2 _eps;     ///< Strain
	Tensor2 _sig_bkp; ///< Backup stress
	Tensor2 _eps_bkp; ///< Backup strain

	// Private methods that MUST be derived
	virtual double _val (char const * Name) const =0;

}; // class EquilibModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline double EquilibModel::Val(char const * Name) const
{
	     if (strcmp(Name,"Sx" )==0)                          return _sig(0);
	else if (strcmp(Name,"Sy" )==0)                          return _sig(1);
	else if (strcmp(Name,"Sz" )==0)                          return _sig(2);
	else if (strcmp(Name,"Sxy")==0 || strcmp(Name,"Syx")==0) return _sig(3)/SQ2;
	else if (strcmp(Name,"Syz")==0 || strcmp(Name,"Szy")==0) return _sig(4)/SQ2;
	else if (strcmp(Name,"Szx")==0 || strcmp(Name,"Sxz")==0) return _sig(5)/SQ2;
	else if (strcmp(Name,"p"  )==0)                          return (_sig(0)+_sig(1)+_sig(2))/3.0;
	else if (strcmp(Name,"q"  )==0)                          return sqrt(((_sig(0)-_sig(1))*(_sig(0)-_sig(1)) + (_sig(1)-_sig(2))*(_sig(1)-_sig(2)) + (_sig(2)-_sig(0))*(_sig(2)-_sig(0)) + 3.0*(_sig(3)*_sig(3) + _sig(4)*_sig(4) + _sig(5)*_sig(5)))/2.0);
	else if (strcmp(Name,"Ex" )==0)                          return _eps(0);
	else if (strcmp(Name,"Ey" )==0)                          return _eps(1);
	else if (strcmp(Name,"Ez" )==0)                          return _eps(2);
	else if (strcmp(Name,"Exy")==0 || strcmp(Name,"Eyx")==0) return _eps(3)/SQ2;
	else if (strcmp(Name,"Eyz")==0 || strcmp(Name,"Ezy")==0) return _eps(4)/SQ2;
	else if (strcmp(Name,"Ezx")==0 || strcmp(Name,"Exz")==0) return _eps(5)/SQ2;
	else if (strcmp(Name,"Ev" )==0)                          return _eps(0)+_eps(1)+_eps(2); 
	else if (strcmp(Name,"Ed" )==0)                          return sqrt(2.0*((_eps(0)-_eps(1))*(_eps(0)-_eps(1)) + (_eps(1)-_eps(2))*(_eps(1)-_eps(2)) + (_eps(2)-_eps(0))*(_eps(2)-_eps(0)) + 3.0*(_eps(3)*_eps(3) + _eps(4)*_eps(4) + _eps(5)*_eps(5))))/3.0;
	else                                                     return _val(Name);
}

#endif // MECHSYS_EQUILIBMODEL_H
