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

// STL
#include <cmath> // for fabs, sqrt, etc.

// MechSys
#include "models/model.h"
#include "util/string.h"
#include "util/util.h"
#include "tensors/functions.h"

using Tensors::Tensor2;
using Tensors::Tensor4;
using Tensors::Tensor2ToVector;
using Tensors::VectorToTensor2;
using Tensors::Tensor4ToMatrix;
using Util::SQ2;
using LinAlg::Vector;
using LinAlg::Matrix;

class EquilibModel : public Model
{
public:
	// Typedefs
	typedef Array<double> IntVals; ///< Internal values (specific volume, yield surface size, etc.)

	// Constructor
	EquilibModel () { STOL().dTini().mMin().mMax().maxSS(); }

	// Destructor
	virtual ~EquilibModel () {}

	// Methods
	void SetGeom      (int Type) { _geom = Type; }                          ///< Geometry type:  1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void TgStiffness  (Matrix<double> & Dmat) const;                        ///< Tangent stiffness tensor
	int  StressUpdate (Vector<double> const & DEps, Vector<double> & DSig); ///< Update stress/strain state for given strain increment
	void BackupState  () { _sig_bkp=_sig; _eps_bkp=_eps; _ivs_bkp=_ivs; };  ///< Backup internal state (sig, eps, ivs; stress, strain, internal values)
	void RestoreState () { _sig=_sig_bkp; _eps=_eps_bkp; _ivs=_ivs_bkp; };  ///< Restore internal state (sig, eps, ivs; stress, strain, internal values)

	// Set methods
	EquilibModel & STOL  (double Val=1.0e-5) { _STOL =Val; return (*this); }
	EquilibModel & dTini (double Val=1.0   ) { _dTini=Val; return (*this); }
	EquilibModel & mMin  (double Val=0.1   ) { _mMin =Val; return (*this); }
	EquilibModel & mMax  (double Val=10.0  ) { _mMax =Val; return (*this); }
	EquilibModel & maxSS (size_t Val=2000  ) { _maxSS=Val; return (*this); }

	// Access methods
	void   Sig (Vector<double> & Stress) const { Tensor2ToVector(_geom, _sig, Stress); }
	void   Eps (Vector<double> & Strain) const { Tensor2ToVector(_geom, _eps, Strain); }
	double Val (char const     * Name  ) const;

protected:
	// Data
	Tensor2 _sig;     ///< Stress
	Tensor2 _eps;     ///< Strain
	IntVals _ivs;     ///< Internal values
	Tensor2 _sig_bkp; ///< Backup stress
	Tensor2 _eps_bkp; ///< Backup strain
	IntVals _ivs_bkp; ///< Backup internal values

	// Constants for the stress update algorithm
	double _STOL;
	double _dTini;
	double _mMin;
	double _mMax;
	size_t _maxSS;

	// Private methods that MUST be derived
	virtual void   _stiff (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const =0; ///< Tangent or secant stiffness
	virtual double _val   (char const * Name) const =0; ///< Return internal values

private:
	// Private methods
	void _tg_incs (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor2 & DSig, IntVals & DIvs);

}; // class EquilibModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void EquilibModel::TgStiffness(Matrix<double> & Dmat) const
{
	Tensor2 deps; deps = 0.0;
	Tensor4        D;
	Array<Tensor2> B;
	_stiff (deps,_sig,_eps,_ivs, D,B);
	Tensor4ToMatrix (_geom,D, Dmat);
}

inline int EquilibModel::StressUpdate(Vector<double> const & DEps, Vector<double> & DSig)
{
	// Auxiliar variables
	Tensor2  sig_i; sig_i = _sig;  // Initial stress state
	Tensor2  deps;                 // Strain increment tensor
	Tensor2  deps_T;               // Driver increment of strain
	Tensor2  dsig_1;               // FE increment of stress
	IntVals  divs_1(_ivs.Size());  // FE increment of internal values
	Tensor2  eps_1;                // FE strain state
	Tensor2  sig_1;                // FE stress state
	IntVals  ivs_1(_ivs.Size());   // FE internal values
	Tensor2  dsig_2;               // Intermediary increment of stress evaluated with a tangent computed at the FE state
	IntVals  divs_2(_ivs.Size());  // Intermediary increment of internal velues evaluated with a tangent computed at the FE state
	Tensor2  sig_ME;               // ME stress state
	IntVals  ivs_ME(_ivs.Size());  // ME internal values
	double   error;                // Estimated error

	// Read input vector
	VectorToTensor2 (_geom,DEps, deps);

	// Solve
	double T  = 0.0;
	double dT = _dTini;
	for (size_t k=0; k<=_maxSS; ++k)
	{
		// Exit point
		if (T>=1.0)
		{
			// Write output vector
			Tensor2 dsig; dsig = _sig - sig_i;
			Tensor2ToVector (_geom,dsig, DSig);
			return k;
		}

		// Driver increment
		deps_T = deps*dT;

		// FE increments (dsig_1 & divs_1)
		_tg_incs (deps_T,_sig,_eps,_ivs, dsig_1,divs_1);

		// FE state
		eps_1    = _eps    + deps_T;
		sig_1    = _sig    + dsig_1;   for (size_t i=0; i<_ivs.Size(); ++i)
		ivs_1[i] = _ivs[i] + divs_1[i];

		// Intermediary increments (dsig_2 & divs_2)
		_tg_incs (deps_T,sig_1,eps_1,ivs_1, dsig_2,divs_2);

		// ME state
		sig_ME    = _sig    + 0.5*(dsig_1   +dsig_2   );  for (size_t i=0; i<_ivs.Size(); ++i)
		ivs_ME[i] = _ivs[i] + 0.5*(divs_1[i]+divs_2[i]);

		// Local error estimate
		error = Tensors::Norm(sig_ME-sig_1)/(1.0+Tensors::Norm(sig_ME));
		for (size_t i=0; i<_ivs.Size(); ++i)
		{
			double err = fabs(ivs_ME[i]-ivs_1[i])/(1.0+ivs_ME[i]);
			if (err>error) error = err;
		}

		// Step multiplier
		double m = 0.9*sqrt(_STOL/error);

		// Update
		if (error<_STOL) // step accepted
		{
			T   += dT;
			_eps = eps_1;
			_sig = sig_ME;
			_ivs = ivs_ME;
			if (m>_mMax) m = _mMax;
		}
		else // step rejected
			if (m<_mMin) m = _mMin;

		// Change next step size
		dT = m * dT;
		if (dT>1.0-T) dT = 1.0-T; // last step
	}
	throw new Fatal("EquilibModel::StressUpdate did not converge for %d steps",_maxSS);
}

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


/* private */

inline void EquilibModel::_tg_incs(Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor2 & DSig, IntVals & DIvs)
{
	// Stiffness
	Tensor4        D;
	Array<Tensor2> B;
	_stiff (DEps,Sig,Eps,Ivs, D,B);

	// Increments
	Tensors::Dot (D,DEps, DSig);      // DSig    = D:DEps
	for (size_t i=0; i<_ivs.Size(); ++i) DIvs[i] = blitz::dot(B[i],DEps);
}


#endif // MECHSYS_EQUILIBMODEL_H
