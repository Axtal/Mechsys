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
	// Constructor
	LinElastic ();
	
	// Destructor
	virtual ~LinElastic () {}

	// Derived Methods
	void SetPrms      (String const & Prms);
	void SetInis      (String const & Inis);
	void TgStiffness  (Matrix<double> & D) const { Tensor4ToMatrix (_geom, _De, D); }
	int  StressUpdate (Vector<double> const & DEps, Vector<double> & DSig);
	void BackupState  ();
	void RestoreState ();

	// Access Methods
	void Sig (Vector<double> & Stress)  const { Tensor2ToVector(_geom, _sig, Stress); }
	void Eps (Vector<double> & Strain)  const { Tensor2ToVector(_geom, _eps, Strain); }
	void Ivs (Array<double>  & IntVals) const { IntVals.Resize(0); } ///< No internal values

private:
	// Data
	Tensor2 _sig;
	Tensor2 _eps;
	Tensor2 _sig_bkp;
	Tensor2 _eps_bkp;
	Tensor4 _De;

}; // class LinElastic


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline LinElastic::LinElastic()
{
	// Parameters
	SetPrms ("E=10000.0 nu=0.25");

	// Initial values
	SetInis ("Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0 Syz=0.0 Szx=0.0");
}

inline void LinElastic::SetPrms(String const & Prms)
{
	/* "E=20000.0 nu=0.2" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	if (names.Size()==2 && values.Size()==2)
	{
		int    count = 0;
		double E     = 0;
		double nu    = 0;
		for (size_t i=0; i<names.Size(); ++i)
		{
			     if (names[i]=="E" ) { E  = values[0];  count++; }
			else if (names[i]=="nu") { nu = values[1];  count++; }
		}
		if (count==2)
		{
			double c  = E/((1.0+nu)*(1.0-2.0*nu));
			double c1 = c*(1.0-nu);
			double c2 = c*(1.0-2.0*nu)/2.0;
			double c3 = c*nu;
			_De = c1     , c3     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
			      c3     , c1     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
			      c3     , c3     , c1     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
			      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, c2 *2.0, 0.0*2.0, 0.0*2.0,
			      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, c2 *2.0, 0.0*2.0,
			      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, 0.0*2.0, c2 *2.0; // In Mandel's basis
			return;
		}
	}

	// Wrong parameters
	throw Fatal("LinElastic::SetPrms: Parameters definition is incorrect. The syntax must be as in:\n\t E=10000.0 nu=0.25\n");
}

inline void LinElastic::SetInis(String const & Inis)
{
	/* "Sx=0.0 Sy=0.0 Sxy=0.0" */
	LineParser lp(Inis);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Check
	if (names.Size()==values.Size() && names.Size()>0)
	{
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
		return;
	}

	// Wrong parameters
	throw Fatal("LinElastic::SetInis: Initial values definition is incorrect. The syntax must be as in:\n\t Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0 Syz=0.0 Szx=0.0\n");
}

inline int LinElastic::StressUpdate(Vector<double> const & DEps, Vector<double> & DSig)
{
	Tensor2 deps;  deps = 0.0;
	Tensor2 dsig;  dsig = 0.0;
	VectorToTensor2 (_geom, DEps, deps);
	dsig = blitz::product(_De, deps);
	Tensor2ToVector (_geom, dsig, DSig);
	_sig += dsig;
	_eps += deps;
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
EquilibModel * LinElasticMaker()
{
	return new LinElastic();
}

// Register an LinElastic model into ModelFactory array map
int LinearElasticRegister()
{
	EquilibModelFactory["LinElastic"] = LinElasticMaker;
	return 0;
}

// Execute the autoregistration
int __LinearElastic_dummy_int = LinearElasticRegister();


#endif // MECHSYS_LINELASTIC_H
