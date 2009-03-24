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
	// Constructor
	EquilibModel () { _deps = 0.0,0.0,0.0, 0.0,0.0,0.0; }

	// Destructor
	virtual ~EquilibModel () {}

	// Derived methods
	void TgStiffness (MechState const & State, Mat_t & Dmat) const;
	int  StateUpdate (double Time, double Dt, Vec_t const & DEps, MechState & State, Vec_t & DSig);

protected:
	/** Tangent or secant stiffness. */
	virtual void _stiff (Tensor2 const  & Sig,
	                     Tensor2 const  & Eps,
	                     IntVals const  & Ivs,
	                     Tensor2 const  & DEps,
	                     Tensor4        & D,
	                     Array<Tensor2> & B) const =0;

	/** Calculate delta sigma star (viscosity). */
	virtual void _dsig_star (double           Time,
	                         double           Dt,
	                         Tensor2 const  & Sig,
	                         Tensor2 const  & Eps,
	                         IntVals const  & Ivs,
	                         Tensor2        & DSigStar) const {}

private:
	// Data
	Tensor2 _deps; ///< Current strain increment. Used to decide between De or Dep (unloading/loading matrices)

	/** Tangent increments. */
	void _tg_incs (double          Time,
	               Tensor2 const & Sig,
	               Tensor2 const & Eps,
	               IntVals const & Ivs,
	               double          Dt,
	               Tensor2 const & DEps,
	               Tensor2       & DSig,
	               IntVals       & DIvs) const;

}; // class EquilibModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void EquilibModel::TgStiffness(MechState const & State, Mat_t & Dmat) const
{
	Tensor4        D;
	Array<Tensor2> B;
	_stiff          (State.Sig,State.Eps,State.Ivs, _deps, D,B);
	Tensor4ToMatrix (_gi,D, Dmat);
}

inline int EquilibModel::StateUpdate(double Time, double Dt, Vec_t const & DEps, MechState & State, Vec_t & DSig)
{
	// Number of internal values
	size_t nivs = State.Ivs.Size();

	// Auxiliar variables
	Tensor2  sig_i;         // Initial stress state
	double   dt_T;          // Driver increment of time
	Tensor2  deps_T;        // Driver increment of strain
	Tensor2  dsig_1;        // FE increment of stress
	IntVals  divs_1(nivs);  // FE increment of internal values
	double   t_1;           // FE time
	Tensor2  eps_1;         // FE strain state
	Tensor2  sig_1;         // FE stress state
	IntVals  ivs_1(nivs);   // FE internal values
	Tensor2  dsig_2;        // Intermediary increment of stress evaluated with a tangent computed at the FE state
	IntVals  divs_2(nivs);  // Intermediary increment of internal velues evaluated with a tangent computed at the FE state
	Tensor2  sig_ME;        // ME stress state
	IntVals  ivs_ME(nivs);  // ME internal values
	double   error;         // Estimated error

	// Initial stress state
   	sig_i = State.Sig;

	// Read input vector
	VectorToTensor2 (_gi,DEps, _deps);

	// Time
	double t = Time;

	// Solve
	double T  = 0.0;
	double dT = _dTini;
	for (size_t k=0; k<=_maxSS; ++k)
	{
		// Exit point
		if (T>=1.0)
		{
			// Write output vector
			Tensor2 dsig;  dsig = State.Sig - sig_i;
			Tensor2ToVector (_gi,dsig, DSig);
			return k;
		}

		// Driver increment
		dt_T   =    Dt*dT;
		deps_T = _deps*dT;

		// FE increments (dsig_1 & divs_1)
		_tg_incs (t, State.Sig,State.Eps,State.Ivs, dt_T,deps_T, dsig_1,divs_1);

		// FE state
		t_1      = t            + dt_T;
		eps_1    = State.Eps    + deps_T;
		sig_1    = State.Sig    + dsig_1;   for (size_t i=0; i<nivs; ++i)
		ivs_1[i] = State.Ivs[i] + divs_1[i];

		// Intermediary increments (dsig_2 & divs_2)
		_tg_incs (t_1, sig_1,eps_1,ivs_1, dt_T,deps_T, dsig_2,divs_2);

		// ME state
		sig_ME    = State.Sig    + 0.5*(dsig_1   +dsig_2   );  for (size_t i=0; i<nivs; ++i)
		ivs_ME[i] = State.Ivs[i] + 0.5*(divs_1[i]+divs_2[i]);

		// Local error estimate
		error = Tensors::Norm(sig_ME-sig_1)/(1.0+Tensors::Norm(sig_ME));
		for (size_t i=0; i<nivs; ++i)
		{
			double err = fabs(ivs_ME[i]-ivs_1[i])/(1.0+ivs_ME[i]);
			if (err>error) error = err;
		}

		// Step multiplier
		double m = 0.9*sqrt(_STOL/error);

		// Update
		if (error<_STOL) // step accepted
		{
			T        += dT;
			t         = t_1;
			State.Eps = eps_1;
			State.Sig = sig_ME;
			State.Ivs = ivs_ME;
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

inline void EquilibModel::_tg_incs(double Time, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs, double Dt, Tensor2 const & DEps, Tensor2 & DSig, IntVals & DIvs) const
{
	// Stiffness
	Tensor4        D;
	Array<Tensor2> B;
	_stiff (Sig,Eps,Ivs,DEps,D,B);

	// Increments
	Tensors::Dot (D,DEps, DSig);     // DSig    = D:DEps
	for (size_t i=0; i<Ivs.Size(); ++i) DIvs[i] = blitz::dot(B[i],DEps);

	// Delta sigma star
	if (HasDSigStar())
	{
		Tensor2 dss;
		_dsig_star (Time,Dt,Sig,Eps,Ivs,dss);
		DSig -= dss;
	}
}


#endif // MECHSYS_EQUILIBMODEL_H
