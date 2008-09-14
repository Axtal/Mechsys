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

#ifndef MECHSYS_LINDIFFUSION_H
#define MECHSYS_LINDIFFUSION_H

// MechSys
#include "models/diffusionmodel.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class LinDiffusion : public DiffusionModel
{
public:
	// Destructor
	virtual ~LinDiffusion () {}

	// Derived Methods
	void         SetPrms (char const * Prms);
	void         SetInis (char const * Inis);
	char const * Name    () const { return "LinDiffusion"; }

private:
	// Data
	TinyMat _K; ///< Conductivity

	// Private methods
	void   _cond (TinyVec const & DuDx, TinyVec const & Vel, IntVals const & Ivs,  TinyMat & D, Array<TinyVec> & B) const;
	double _val  (char const * Name) const { throw new Fatal("LinDiffusion::_val: The Name==%s is invalid",Name); }

}; // class LinDiffusion


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void LinDiffusion::SetPrms(char const * Prms)
{
	if (_geom<0) throw new Fatal("LinDiffusion::SetPrms: Geometry type:\n\t[1:1D, 2:2D, 3:3D] must be set via SetGeom before calling this method");

	/* "kxx=1.0 kxy=0.0 kxz=0.0
	            kyy=1.0 kxz=0.0
	                    kzz=1.0"
	    or "k=1.0" => isotropic
	*/
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Conductivity matrix
	if (_geom<1 || _geom>3) new Fatal("LinDiffusion::SetPrms: Geometry==%d is invalid. The valid one must be one of [1:1D, 2:2D, 3:3D]",_geom);
	_K = 0.0;

	// Set
	if (names.Size()==1)
	{
		if (names[0]=="k")
		{
			_K = values[0],       0.0,       0.0,
			           0.0, values[0],       0.0,
			           0.0,       0.0, values[0];
		}
		else throw new Fatal("LinDiffusion::SetPrms: Parameter key==%s for isotropic models is invalid. It must be equal to 'k'. Ex.: k=1.0",names[0].CStr());
	}
	else
	{
		if (_geom==1) throw new Fatal("LinDiffusion::SetPrms: For unidimensional problems, only one parameter key (equal to 'k') must be used. Ex.: k=1.0 (%s is invalid)");
		for (size_t i=0; i<names.Size(); ++i)
		{
			      if (names[i]=="kxx")                                { _K(0,0) = values[i];                       }
			 else if (names[i]=="kxy" || names[i]=="kyx")             { _K(0,1) = values[i];  _K(1,0) = values[i]; }
			 else if (names[i]=="kxz" || names[i]=="kzx" && _geom==3) { _K(0,2) = values[i];  _K(2,0) = values[i]; }
			 else if (names[i]=="kyy")                                { _K(1,1) = values[i];                       }
			 else if (names[i]=="kyz" || names[i]=="kzy" && _geom==3) { _K(1,2) = values[i];  _K(2,1) = values[i]; }
			 else if (names[i]=="kzz"                    && _geom==3) { _K(2,2) = values[i];                       }
			 else throw new Fatal("LinDiffusion::SetPrms: Parameter key==%s is invalid. It must be: kxx, kxy, kxz,  kyy, kyz,  kzz  (or kyx, kzx, kzy), where the 'z-coefficients' are valid only for 3D problems.",names[i].CStr());
		}
	}
}

inline void LinDiffusion::SetInis(char const * Inis)
{
	/* "ZERO" */
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	for (size_t i=0; i<names.Size(); i++)
	{
		if (names[i]!="ZERO") throw new Fatal("LinDiffusion::SetInis: Initial value key==%s is invalid.",names[i].CStr());
	}
}


/* private */

inline void LinDiffusion::_cond(TinyVec const & DuDx, TinyVec const & Vel, IntVals const & Ivs,  TinyMat & D, Array<TinyVec> & B) const
{
	D = _K;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new LinDiffusion model
Model * LinDiffusionMaker()
{
	return new LinDiffusion();
}

// Register an LinDiffusion model into ModelFactory array map
int LinearDiffusionRegister()
{
	ModelFactory["LinDiffusion"] = LinDiffusionMaker;
	return 0;
}

// Execute the autoregistration
int __LinearDiffusion_dummy_int = LinearDiffusionRegister();


#endif // MECHSYS_LINDIFFUSION_H
