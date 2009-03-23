/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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
#include <cstring> // for strcmp

// MechSys
#include "models/model.h"
#include "util/string.h"
#include "util/util.h"
#include "tensors/tensors.h"
#include "tensors/functions.h"

using Tensors::Tensor2;
using Tensors::Tensor4;
using Tensors::Tensor2ToVector;
using Tensors::VectorToTensor2;
using Tensors::Tensor4ToMatrix;
using Util::SQ2;

typedef LinAlg::Vector<double> Vec_t;
typedef LinAlg::Matrix<double> Mat_t;
typedef char const           * Str_t;

class EquilibModel : public Model
{
public:
	// Typedefs
	typedef Array<double> IntVals; ///< Internal values (specific volume, yield surface size, etc.)

	// Constructor
	EquilibModel () { _deps = 0.0,0.0,0.0, 0.0,0.0,0.0; }

	// Destructor
	virtual ~EquilibModel () {}

	/* Tangent stiffness. */
	void TgStiffness (Tensor2 const & Sig,
	                  Tensor2 const & Eps,
	                  IntVals const & Ivs,
	                  Mat_t         & Dmat,
	                  bool            First) const;

	/* State update. */
	int StateUpdate (Vec_t   const & DEps,
	                 Tensor2       & Sig,
	                 Tensor2       & Eps,
	                 IntVals       & Ivs,
	                 Vec_t         & DSig);

	/* State update. */
	int StateUpdate (Vec_t   const & DEps,
	                 Vec_t   const & DM,
	                 Tensor2       & Sig,
	                 Tensor2       & Eps,
	                 IntVals       & Ivs,
	                 Vec_t         & DSig);

protected:
	/* Tangent or secant stiffness. */
	virtual void _stiff (Tensor2 const  & DEps,
	                     Tensor2 const  & Sig,
	                     Tensor2 const  & Eps,
	                     IntVals const  & Ivs,
	                     Tensor4        & D,
	                     Array<Tensor2> & B,
	                     bool             First) const =0;

private:
	// Data
	Tensor2 _deps;

	/* Tangent increments. */
	void _tg_incs (Tensor2 const & DEps,
	               Tensor2 const & Sig,
	               Tensor2 const & Eps,
	               IntVals const & Ivs,
	               Tensor2       & DSig,
	               IntVals       & DIvs) const;

}; // class EquilibModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void EquilibModel::TgStiffness(Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs, Mat_t & Dmat, bool First) const
{
	Tensor4        D;
	Array<Tensor2> B;
	_stiff          (_deps,Sig,Eps,Ivs, D,B, First);
	Tensor4ToMatrix (_gi,D, Dmat);
}

inline int EquilibModel::StateUpdate(Vec_t const & DEps, Tensor2 & Sig, Tensor2 & Eps, IntVals & Ivs, Vec_t & DSig)
{
	// Auxiliar variables
	Tensor2  sig_i; sig_i = Sig;  // Initial stress state
	Tensor2  deps_T;              // Driver increment of strain
	Tensor2  dsig_1;              // FE increment of stress
	IntVals  divs_1(Ivs.Size());  // FE increment of internal values
	Tensor2  eps_1;               // FE strain state
	Tensor2  sig_1;               // FE stress state
	IntVals  ivs_1(Ivs.Size());   // FE internal values
	Tensor2  dsig_2;              // Intermediary increment of stress evaluated with a tangent computed at the FE state
	IntVals  divs_2(Ivs.Size());  // Intermediary increment of internal velues evaluated with a tangent computed at the FE state
	Tensor2  sig_ME;              // ME stress state
	IntVals  ivs_ME(Ivs.Size());  // ME internal values
	double   error;               // Estimated error

	// Read input vector
	VectorToTensor2 (_gi,DEps, _deps);

	// Solve
	double T  = 0.0;
	double dT = _dTini;
	for (size_t k=0; k<=_maxSS; ++k)
	{
		// Exit point
		if (T>=1.0)
		{
			// Write output vector
			Tensor2 dsig; dsig = Sig - sig_i;
			Tensor2ToVector (_gi,dsig, DSig);
			return k;
		}

		// Driver increment
		deps_T = _deps*dT;

		// FE increments (dsig_1 & divs_1)
		_tg_incs (deps_T,Sig,Eps,Ivs, dsig_1,divs_1);

		// FE state
		eps_1    = Eps    + deps_T;
		sig_1    = Sig    + dsig_1;   for (size_t i=0; i<Ivs.Size(); ++i)
		ivs_1[i] = Ivs[i] + divs_1[i];

		// Intermediary increments (dsig_2 & divs_2)
		_tg_incs (deps_T,sig_1,eps_1,ivs_1, dsig_2,divs_2);

		// ME state
		sig_ME    = Sig    + 0.5*(dsig_1   +dsig_2   );  for (size_t i=0; i<Ivs.Size(); ++i)
		ivs_ME[i] = Ivs[i] + 0.5*(divs_1[i]+divs_2[i]);

		// Local error estimate
		error = Tensors::Norm(sig_ME-sig_1)/(1.0+Tensors::Norm(sig_ME));
		for (size_t i=0; i<Ivs.Size(); ++i)
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
			Eps = eps_1;
			Sig = sig_ME;
			Ivs = ivs_ME;
			if (m>_mMax) m = _mMax;
		}
		else // step rejected
			if (m<_mMin) m = _mMin;

		// Change next step size
		dT = m * dT;
		if (dT>1.0-T) dT = 1.0-T; // last step
	}
	throw new Fatal("EquilibModel::StateUpdate: %s:Tag=%d: Update did not converge for %d steps",Name(),_tag,_maxSS);
}

inline int EquilibModel::StateUpdate(Vec_t const & DEps, Vec_t const & DM, Tensor2 & Sig, Tensor2 & Eps, IntVals & Ivs, Vec_t & DSig)
{
	// Auxiliar variables
	Tensor2  sig_i; sig_i = Sig;  // Initial stress state
	Tensor2  deps_T;              // Driver increment of strain
	Tensor2  dsig_1;              // FE increment of stress
	IntVals  divs_1(Ivs.Size());  // FE increment of internal values
	Tensor2  eps_1;               // FE strain state
	Tensor2  sig_1;               // FE stress state
	IntVals  ivs_1(Ivs.Size());   // FE internal values
	Tensor2  dsig_2;              // Intermediary increment of stress evaluated with a tangent computed at the FE state
	IntVals  divs_2(Ivs.Size());  // Intermediary increment of internal velues evaluated with a tangent computed at the FE state
	Tensor2  sig_ME;              // ME stress state
	IntVals  ivs_ME(Ivs.Size());  // ME internal values
	double   error;               // Estimated error

	Tensor2  dM;                  // Term due to viscosity
	Tensor2  dM_T;                // Term due to viscosity

	// Read input vector
	VectorToTensor2 (_gi,DEps, _deps);
	VectorToTensor2 (_gi,DM,    dM);

	// Solve
	double T  = 0.0;
	double dT = _dTini;
	for (size_t k=0; k<=_maxSS; ++k)
	{
		// Exit point
		if (T>=1.0)
		{
			// Write output vector
			Tensor2 dsig; dsig = Sig - sig_i;
			Tensor2ToVector (_gi,dsig, DSig);
			return k;
		}

		// Driver increment
		deps_T = _deps*dT;
		dM_T   =    dM*dT;

		// FE increments (dsig_1 & divs_1)
		_tg_incs (deps_T,Sig,Eps,Ivs, dsig_1,divs_1);
		dsig_1 -= dM_T;

		// FE state
		eps_1    = Eps    + deps_T;
		sig_1    = Sig    + dsig_1;   for (size_t i=0; i<Ivs.Size(); ++i)
		ivs_1[i] = Ivs[i] + divs_1[i];

		// Intermediary increments (dsig_2 & divs_2)
		_tg_incs (deps_T,sig_1,eps_1,ivs_1, dsig_2,divs_2);
		dsig_2 -= dM_T;

		// ME state
		sig_ME    = Sig    + 0.5*(dsig_1   +dsig_2   );  for (size_t i=0; i<Ivs.Size(); ++i)
		ivs_ME[i] = Ivs[i] + 0.5*(divs_1[i]+divs_2[i]);

		// Local error estimate
		error = Tensors::Norm(sig_ME-sig_1)/(1.0+Tensors::Norm(sig_ME));
		for (size_t i=0; i<Ivs.Size(); ++i)
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
			Eps = eps_1;
			Sig = sig_ME;
			Ivs = ivs_ME;
			if (m>_mMax) m = _mMax;
		}
		else // step rejected
			if (m<_mMin) m = _mMin;

		// Change next step size
		dT = m * dT;
		if (dT>1.0-T) dT = 1.0-T; // last step
	}
	throw new Fatal("EquilibModel::StateUpdate: %s:Tag=%d: Update did not converge for %d steps",Name(),_tag,_maxSS);
}


/* private */

inline void EquilibModel::_tg_incs(Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor2 & DSig, IntVals & DIvs) const
{
	// Stiffness
	Tensor4        D;
	Array<Tensor2> B;
	_stiff (DEps,Sig,Eps,Ivs, D,B, false);

	// Increments
	Tensors::Dot (D,DEps, DSig);     // DSig    = D:DEps
	for (size_t i=0; i<Ivs.Size(); ++i) DIvs[i] = blitz::dot(B[i],DEps);
}


#endif // MECHSYS_EQUILIBMODEL_H
