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

// STL
#include <cmath> // for fabs, sqrt, etc.
#include <cstring> // for strcmp

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "models/model.h"
#include "util/string.h"
#include "util/util.h"
#include "tensors/functions.h"

typedef LinAlg::Vector<double> Vec_t;
typedef LinAlg::Matrix<double> Mat_t;
typedef char const           * Str_t;

class DiffusionModel : public Model
{
public:
	// Typedefs
	typedef Array<double>                 IntVals; ///< Internal values
	typedef blitz::TinyVector<double,3>   Vec3_t;
	typedef blitz::TinyMatrix<double,3,3> Mat3_t;

	// Constructor
	DiffusionModel () { _dgra = 0.0,0.0,0.0; }

	// Destructor
	virtual ~DiffusionModel () {}

	/* Tangent conductivity. */
	void TgConductivity (Vec3_t  const & Vel,
	                     Vec3_t  const & Gra,
	                     IntVals const & Ivs,
	                     Mat_t         & Dmat) const;

	/* State update. */
	int StateUpdate (Vec_t   const & DGra,
	                 Vec3_t        & Vel,
	                 Vec3_t        & Gra,
	                 IntVals       & Ivs,
	                 Vec_t         & DVel);

	static void Vec3ToVec (int GeomIdx, Vec3_t const & D, Vec_t  & Dvec);
	static void VecToVec3 (int GeomIdx, Vec_t  const & Dvec, Vec3_t & D);
	static void Mat3ToMat (int GeomIdx, Mat3_t const & D, Mat_t  & Dmat);

protected:
	/* Tangent or secant conductivity. */
	virtual void _cond (Vec3_t  const & DGra,
	                    Vec3_t  const & Vel,
	                    Vec3_t  const & Gra,
	                    IntVals const & Ivs,
	                    Mat3_t        & D,
	                    Array<Vec3_t> & B) const =0;

private:
	// Data
	Vec3_t _dgra;

	/* Tangent increments. */
	void _tg_incs (Vec3_t  const & DGra,
	               Vec3_t  const & Vel,
	               Vec3_t  const & Gra,
	               IntVals const & Ivs,
	               Vec3_t        & DVel,
	               IntVals       & DIvs) const;

	double _norm        (Vec3_t const & V) const { return sqrt(blitz::dot(V,V)); }

}; // class DiffusionModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void DiffusionModel::TgConductivity(Vec3_t const & Vel, Vec3_t const & Gra, IntVals const & Ivs, Mat_t & Dmat) const
{
	Mat3_t        D;
	Array<Vec3_t> B;
	_cond        (_dgra,Vel,Gra,Ivs, D,B);
	Mat3ToMat    (_gi, D, Dmat);
}

inline int DiffusionModel::StateUpdate(Vec_t const & DGra, Vec3_t & Vel, Vec3_t & Gra, IntVals & Ivs, Vec_t & DVel)
{
	// Auxiliar variables
	Vec3_t  vel_i; vel_i = Vel;  // Initial stress state
	Vec3_t  dgra_T;              // Driver increment of strain
	Vec3_t  dvel_1;              // FE increment of stress
	IntVals  divs_1(Ivs.Size());  // FE increment of internal values
	Vec3_t  gra_1;               // FE strain state
	Vec3_t  vel_1;               // FE stress state
	IntVals  ivs_1(Ivs.Size());   // FE internal values
	Vec3_t  dvel_2;              // Intermediary increment of stress evaluated with a tangent computed at the FE state
	IntVals  divs_2(Ivs.Size());  // Intermediary increment of internal velues evaluated with a tangent computed at the FE state
	Vec3_t  vel_ME;              // ME stress state
	IntVals  ivs_ME(Ivs.Size());  // ME internal values
	double   error;               // Estimated error

	// Read input vector
	VecToVec3 (_gi, DGra, _dgra);

	// Solve
	double T  = 0.0;
	double dT = _dTini;
	for (size_t k=0; k<=_maxSS; ++k)
	{
		// Exit point
		if (T>=1.0)
		{
			// Write output vector
			Vec3_t dvel; dvel = Vel - vel_i;
			Vec3ToVec (_gi, dvel, DVel);
			return k;
		}

		// Driver increment
		dgra_T = _dgra*dT;

		// FE increments (dvel_1 & divs_1)
		_tg_incs (dgra_T,Vel,Gra,Ivs, dvel_1,divs_1);

		// FE state
		gra_1    = Gra    + dgra_T;
		vel_1    = Vel    + dvel_1;   for (size_t i=0; i<Ivs.Size(); ++i)
		ivs_1[i] = Ivs[i] + divs_1[i];

		// Intermediary increments (dvel_2 & divs_2)
		_tg_incs (dgra_T,vel_1,gra_1,ivs_1, dvel_2,divs_2);

		// ME state
		vel_ME    = Vel    + 0.5*(dvel_1   +dvel_2   );  for (size_t i=0; i<Ivs.Size(); ++i)
		ivs_ME[i] = Ivs[i] + 0.5*(divs_1[i]+divs_2[i]);

		// Local error estimate
		error = _norm(vel_ME-vel_1)/(1.0+_norm(vel_ME));
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
			Gra = gra_1;
			Vel = vel_ME;
			Ivs = ivs_ME;
			if (m>_mMax) m = _mMax;
		}
		else // step rejected
			if (m<_mMin) m = _mMin;

		// Change next step size
		dT = m * dT;
		if (dT>1.0-T) dT = 1.0-T; // last step
	}
	throw new Fatal("DiffusionModel::StateUpdate: %s:Tag=%s: Update did not converge for %d stgra",Name(),_tag,_maxSS);
}

inline void DiffusionModel::Vec3ToVec(int GeomIdx, Vec3_t const & D, Vec_t & Dvec)
{
	if (GeomIdx==0) // 3D
	{
		Dvec.Resize(3);
		Dvec = D(0), D(1), D(2);
	}
	else if (GeomIdx==1) // 2D
	{
		Dvec.Resize(2);
		Dvec = D(0), D(1);
	}
	else throw new Fatal("DiffusionModel::Vec3ToVec: GeomIdx==%d is invalid",GeomIdx);
}

inline void DiffusionModel::VecToVec3(int GeomIdx, Vec_t const & Dvec, Vec3_t & D)
{
	if (GeomIdx==0) // 3D
	{
		D = Dvec(0), Dvec(1), Dvec(2);
	}
	else if (GeomIdx==1) // 2D
	{
		D = Dvec(0), Dvec(1), 0.0;
	}
	else throw new Fatal("DiffusionModel::VecToVec3: GeomIdx==%d is invalid",GeomIdx);
}

inline void DiffusionModel::Mat3ToMat(int GeomIdx, Mat3_t const & D, Mat_t & Dmat)
{
	if (GeomIdx==0) // 3D
	{
		Dmat.Resize(3,3);
		Dmat = D(0,0), D(0,1), D(0,2),
		       D(1,0), D(1,1), D(1,2),
		       D(2,0), D(2,1), D(2,2);
	}
	else if (GeomIdx==1) // 2D
	{
		Dmat.Resize(2,2);
		Dmat = D(0,0), D(0,1),
		       D(1,0), D(1,1);
	}
	else throw new Fatal("DiffusionModel::Mat3ToMat: GeomIdx==%d is invalid",GeomIdx);
}


/* private */

inline void DiffusionModel::_tg_incs(Vec3_t const & DGra, Vec3_t const & Vel, Vec3_t const & Gra, IntVals const & Ivs,  Vec3_t & DVel, IntVals & DIvs) const
{
	// Conductivity
	Mat3_t        D;
	Array<Vec3_t> B;
	_cond (DGra,Vel,Gra,Ivs, D,B);

	// Increments
	DVel = blitz::product(D,DGra);   // DVel    = D * DGra
	for (size_t i=0; i<Ivs.Size(); ++i) DIvs[i] = blitz::dot(B[i],DGra);
}

#endif // MECHSYS_DIFFUSIONMODEL_H
