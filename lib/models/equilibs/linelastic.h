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

#ifndef MECHSYS_LINELASTIC_H
#define MECHSYS_LINELASTIC_H

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

using Tensors::Tensor2;
using Tensors::Tensor4;
using Tensors::VectorToTensor2;
using Tensors::Tensor2ToVector;
using Tensors::Tensor4ToMatrix;
using Util::SQ2;
using LinAlg::Vector;
using LinAlg::Matrix;

class LinElastic : public EquilibModel
{
public:
	// Destructor
	virtual ~LinElastic () {}

	// Derived Methods
	void SetPrms      (char const * Prms);
	void SetInis      (char const * Inis);
	void TgStiffness  (Matrix<double> & D) const { Tensor4ToMatrix (_geom, _De, D); }
	int  StressUpdate (Vector<double> const & DEps, Vector<double> & DSig);
	void BackupState  ();
	void RestoreState ();

private:
	// Data
	Tensor4 _De;

	// Private methods
	double _val (char const * Name) const { throw new Fatal("LinElastic::_val The Name==%s is invalid",Name); }

}; // class LinElastic


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void LinElastic::SetPrms(char const * Prms)
{
	if (_geom<0) throw new Fatal("LinElastic::SetPrms: Geometry type:\n\t[1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)] must be set via SetGeom before calling this method");

	/* "E=20000.0 nu=0.2" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	double E  = 0.0;
	double nu = 0.0;
	for (size_t i=0; i<names.Size(); ++i)
	{
			 if (names[i]=="E" ) E  = values[i];
		else if (names[i]=="nu") nu = values[i];
	}
	if (_geom==1) _De(0,0) = E;
	else
	{
		double c  = (_geom==5 ? E/(1.0-nu*nu)  : E/((1.0+nu)*(1.0-2.0*nu)) ); // plane-stress != (plane-strain=3D)
		double c1 = (_geom==5 ? c*1.0          : c*(1.0-nu)                ); // plane-stress != (plane-strain=3D)
		double c2 = (_geom==5 ? c*0.5*(1.0-nu) : c*(1.0-2.0*nu)/2.0        ); // plane-stress != (plane-strain=3D)
		double c3 = c*nu;
		_De = c1     , c3     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
		      c3     , c1     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
		      c3     , c3     , c1     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
		      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, c2 *2.0, 0.0*2.0, 0.0*2.0,
		      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, c2 *2.0, 0.0*2.0,
		      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, 0.0*2.0, c2 *2.0; // In Mandel's basis
	}
}

inline void LinElastic::SetInis(char const * Inis)
{
	/* "Sx=0.0 Sy=0.0 Sxy=0.0" */
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	_sig = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	_eps = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	for (size_t i=0; i<names.Size(); i++)
	{
			 if (names[i]=="Sx")                     _sig(0) = values[i];
		else if (names[i]=="Sy")                     _sig(1) = values[i];
		else if (names[i]=="Sz")                     _sig(2) = values[i];
		else if (names[i]=="Sxy" || names[i]=="Syx") _sig(3) = values[i]*SQ2;
		else if (names[i]=="Syz" || names[i]=="Szy") _sig(4) = values[i]*SQ2;
		else if (names[i]=="Szx" || names[i]=="Sxz") _sig(5) = values[i]*SQ2;
	}
}

inline int LinElastic::StressUpdate(Vector<double> const & DEps, Vector<double> & DSig)
{
	if (_geom==1)
	{
		DSig.Resize(1);
		DSig(0) = _De(0,0) * DEps(0);
		_sig(0) += DSig(0);
		_eps(0) += DEps(0);
	}
	else
	{
		Tensor2 deps;  deps = 0.0;
		Tensor2 dsig;  dsig = 0.0;
		VectorToTensor2 (_geom, DEps, deps);
		dsig = blitz::product(_De, deps);
		Tensor2ToVector (_geom, dsig, DSig);
		_sig += dsig;
		_eps += deps;
	}
	return 1;
}

inline void LinElastic::BackupState()
{
	_sig_bkp = _sig;
	_eps_bkp = _eps;
}

inline void LinElastic::RestoreState()
{
	_sig = _sig_bkp;
	_eps = _eps_bkp;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new LinElastic model
Model * LinElasticMaker()
{
	return new LinElastic();
}

// Register an LinElastic model into ModelFactory array map
int LinearElasticRegister()
{
	ModelFactory["LinElastic"] = LinElasticMaker;
	return 0;
}

// Execute the autoregistration
int __LinearElastic_dummy_int = LinearElasticRegister();


#endif // MECHSYS_LINELASTIC_H
