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

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "models/model.h"
#include "util/string.h"
#include "util/util.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"

using LinAlg::Vector;
using LinAlg::Matrix;

class DiffusionModel : public Model
{
public:
	// Typedefs
	typedef blitz::TinyVector<double,3>   TinyVec;
	typedef blitz::TinyMatrix<double,3,3> TinyMat;

	// Destructor
	virtual ~DiffusionModel () {}

	// Methods
	void SetGeom        (int Type) { _geom = Type; if (_geom<1 || _geom>3) throw new Fatal("DiffusionModel::SetGeom: Geometry type must be: 1:1D, 2:2D, 3:3D."); } ///< Geometry type:  1:1D, 2:2D, 3:3D
	void TgConductivity (Matrix<double> & Dmat) const;                        ///< Tangent conductivity
	int  StateUpdate    (Vector<double> const & DuDx, Vector<double> & DVel); ///< Update internal state for given du_dx

	// Access methods
	double Val (char const * Name  ) const;

protected:
	// Data
	TinyVec _vel;     ///< Diffusion velocity: {vel} = [K]{du_dx}
	TinyVec _vel_bkp; ///< Backup velocity

	// Private methods that MUST be derived
	virtual void   _cond (TinyVec const & DuDx, TinyVec const & Vel, IntVals const & Ivs,  TinyMat & D, Array<TinyVec> & B) const =0; ///< Tangent or secant conductivity
	virtual double _val  (char const * Name) const =0; ///< Return internal values

private:
	// Derived methods
	void _backup_state  () { _vel_bkp=_vel; } ///< Backup internal state
	void _restore_state () { _vel=_vel_bkp; } ///< Restore internal state

	// Private methods
	void _tg_incs (TinyVec const & DuDx, TinyVec const & Vel, IntVals const & Ivs,  TinyVec & DVel, IntVals & DIvs) const; ///< Tangent increments

}; // class DiffusionModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void DiffusionModel::TgConductivity(Matrix<double> & Dmat) const
{
	// Conductivity
	TinyVec dudx; dudx = 0.0;
	TinyMat        D;
	Array<TinyVec> B;
	_cond (dudx,_vel,_ivs, D,B);
	switch (_geom)
	{
		case 1: // 1D
		{
			Dmat.Resize(1,1);
			Dmat(0,0) = D(0,0);
			return;
		}
		case 2: // 2D
		{
			Dmat.Resize(2,2);
			Dmat = D(0,0), D(0,1),
			       D(1,0), D(1,1);
			return;
		}
		case 3: // 3D
		{
			Dmat.Resize(3,3);
			Dmat = D(0,0), D(0,1), D(0,2),
			       D(1,0), D(1,1), D(1,2),
			       D(2,0), D(2,1), D(2,2);
			return;
		}
	}
}

inline int DiffusionModel::StateUpdate(Vector<double> const & DuDx, Vector<double> & DVel)
{
	// Gradient
	TinyVec dudx;
	switch (_geom)
	{
		case 1: { dudx = DuDx(0);                   break; }
		case 2: { dudx = DuDx(0), DuDx(1);          break; }
		case 3: { dudx = DuDx(0), DuDx(1), DuDx(2); break; }
	}

	// Forward-Euler update
	TinyVec dvel;
	IntVals divs;
	_tg_incs (dudx,_vel,_ivs, dvel,divs);
	_vel    = _vel    + dvel;   for (size_t i=0; i<_ivs.Size(); ++i)
	_ivs[i] = _ivs[i] + divs[i];

	return 0;
}

inline double DiffusionModel::Val(char const * Name) const
{
	     if (strcmp(Name,"Vx" )==0) return _vel(0);
	else if (strcmp(Name,"Vy" )==0) return _vel(1);
	else if (strcmp(Name,"Vz" )==0) return _vel(2);
	return _val(Name);
}


/* private */

inline void DiffusionModel::_tg_incs(TinyVec const & DuDx, TinyVec const & Vel, IntVals const & Ivs,  TinyVec & DVel, IntVals & DIvs) const
{
	// Conductivity
	TinyMat        D;
	Array<TinyVec> B;
	_cond (DuDx,Vel,Ivs, D,B);

	// Increments
	DVel = blitz::product (D,DuDx);   // DVel    = D*DuDx
	for (size_t i=0; i<_ivs.Size(); ++i) DIvs[i] = blitz::dot(B[i],DuDx);
}

#endif // MECHSYS_DIFFUSIONMODEL_H
