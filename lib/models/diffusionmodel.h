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

// Std lib
#include <cstring> // for strcmp

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

	// Constructor
	DiffusionModel () { _vel=0.0; _gra=0.0; _dgra=0.0; }

	// Destructor
	virtual ~DiffusionModel () {}

	// Methods
	void SetGeom        (int Type) { _geom = Type; if (_geom<1 || _geom>3) throw new Fatal("DiffusionModel::SetGeom: Geometry type must be: 1:1D, 2:2D, 3:3D."); } ///< Geometry type:  1:1D, 2:2D, 3:3D
	void TgConductivity (Matrix<double> & Dmat) const;                        ///< Tangent conductivity
	int  StateUpdate    (Vector<double> const & DGra, Vector<double> & DVel); ///< Update internal state for given DGra = Delta(du_dx)

	// Access methods
	void   CalcDepVars () const {}                     ///< Calculate dependent variables (to be called before Val() or OutNodes() for example).
	double Val         (char const * Name) const;      ///< Return stress/strain components, internal values, or principal components of stress/strain
	void   Vel         (Vector<double> & Veloc) const; ///< Return stress tensor

protected:
	// Data
	TinyVec _vel;      ///< Diffusion velocity: {vel} = [K]{du_dx}
	TinyVec _gra;      ///< Flow gradient:      {gra} = du_dx
	TinyVec _dgra;     ///< Delta gradient
	TinyVec _vel_bkp;  ///< Backup velocity
	TinyVec _gra_bkp;  ///< Backup gradient
	TinyVec _dgra_bkp; ///< Backup delta gradient

	// Private methods that MUST be derived
	virtual void   _cond (TinyVec const & DGra, TinyVec const & Vel, TinyVec const & Gra, IntVals const & Ivs,  TinyMat & D, Array<TinyVec> & B) const =0; ///< Tangent or secant conductivity
	virtual double _val  (char const * Name) const =0; ///< Return internal values

private:
	// Derived methods
	void _backup_state  () { _vel_bkp=_vel; _gra_bkp=_gra; _dgra_bkp=_dgra; } ///< Backup internal state
	void _restore_state () { _vel=_vel_bkp; _gra=_gra_bkp; _dgra=_dgra_bkp; } ///< Restore internal state

	// Private methods
	void _tg_incs  (TinyVec const & DGra, TinyVec const & Vel, TinyVec const & Gra, IntVals const & Ivs,  TinyVec & DVel, IntVals & DIvs) const; ///< Tangent increments
	void _tvec2vec (TinyVec const & TVec, Vector<double> & Vect) const;
	void _tmat2mat (TinyMat const & TMat, Matrix<double> & Matr) const;

}; // class DiffusionModel


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void DiffusionModel::TgConductivity(Matrix<double> & Dmat) const
{
	// Conductivity
	TinyMat        D;
	Array<TinyVec> B;
	_cond (_dgra,_vel,_gra,_ivs, D,B);

	// Result
	_tmat2mat (D, Dmat);
}

inline int DiffusionModel::StateUpdate(Vector<double> const & DGra, Vector<double> & DVel)
{
	// Gradient
	_dgra = 0.0;
	switch (_geom)
	{
		case 1: { _dgra = DGra(0);                   break; }
		case 2: { _dgra = DGra(0), DGra(1);          break; }
		case 3: { _dgra = DGra(0), DGra(1), DGra(2); break; }
	}

	// Forward-Euler update
	TinyVec dvel;
	IntVals divs;
	_tg_incs (_dgra,_vel,_gra,_ivs, dvel,divs);
	_vel    = _vel    +  dvel;
	_gra    = _gra    + _dgra;   for (size_t i=0; i<_ivs.Size(); ++i)
	_ivs[i] = _ivs[i] +  divs[i];

	// Result
	_tvec2vec (dvel, DVel);

	// 1 substep
	return 1;
}

inline double DiffusionModel::Val(char const * Name) const
{
	     if (strcmp(Name,"Vx" )==0) return _vel(0); // velocity
	else if (strcmp(Name,"Vy" )==0) return _vel(1);
	else if (strcmp(Name,"Vz" )==0) return _vel(2);
	else if (strcmp(Name,"Ix" )==0) return _gra(0); // gradient
	else if (strcmp(Name,"Iy" )==0) return _gra(1);
	else if (strcmp(Name,"Iz" )==0) return _gra(2);
	return _val(Name);
}

inline void DiffusionModel::Vel(Vector<double> & Veloc) const
{
	_tvec2vec(_vel, Veloc);
}


/* private */

inline void DiffusionModel::_tg_incs(TinyVec const & DGra, TinyVec const & Vel, TinyVec const & Gra, IntVals const & Ivs,  TinyVec & DVel, IntVals & DIvs) const
{
	// Conductivity
	TinyMat        D;
	Array<TinyVec> B;
	_cond (DGra,Vel,Gra,Ivs, D,B);

	// Increments
	DVel = -blitz::product (D,DGra);  // DVel    = -D*Dgra
	for (size_t i=0; i<_ivs.Size(); ++i) DIvs[i] = blitz::dot(B[i],DGra);
}

inline void DiffusionModel::_tvec2vec(TinyVec const & TVec, Vector<double> & Vect) const
{
	switch (_geom)
	{
		case 1: // 1D
		{
			Vect.Resize(1);
			Vect(0) = TVec(0);
			return;
		}
		case 2: // 2D
		{
			Vect.Resize(2);
			Vect = TVec(0), TVec(1);
			return;
		}
		case 3: // 3D
		{
			Vect.Resize(3);
			Vect = TVec(0), TVec(1), TVec(2);
			return;
		}
	}
}

inline void DiffusionModel::_tmat2mat(TinyMat const & TMat, Matrix<double> & Matr) const
{
	switch (_geom)
	{
		case 1: // 1D
		{
			Matr.Resize(1,1);
			Matr(0,0) = TMat(0,0);
			return;
		}
		case 2: // 2TMat
		{
			Matr.Resize(2,2);
			Matr = TMat(0,0), TMat(0,1),
			       TMat(1,0), TMat(1,1);
			return;
		}
		case 3: // 3TMat
		{
			Matr.Resize(3,3);
			Matr = TMat(0,0), TMat(0,1), TMat(0,2),
			       TMat(1,0), TMat(1,1), TMat(1,2),
			       TMat(2,0), TMat(2,1), TMat(2,2);
			return;
		}
	}
}

#endif // MECHSYS_DIFFUSIONMODEL_H
