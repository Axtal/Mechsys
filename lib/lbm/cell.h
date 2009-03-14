/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo, Fernando Alonso                   *
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

#ifndef MECHSYS_LBM_CELL_H
#define MECHSYS_LBM_CELL_H

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "util/array.h"
#include "util/exception.h"

typedef blitz::TinyVector<double,3> Vec3_t;

namespace LBM
{

class Cell
{
public:
	// Enums
	enum BCType_T { NONE_T, VELOCITY_T, DENSITY_T };
	enum BCSide_T { LEFT_T, RIGHT_T, BOTTOM_T, TOP_T, FRONT_T, BACK_T };

	// Typedefs
	typedef double LVeloc_T[3]; ///< Local velocities type

	// Constants
	static const double   WEIGHTS2D     [ 9]; ///< Weights for the equilibrium distribution functions (2D)
	static const double   WEIGHTS3D     [27]; ///< Weights for the equilibrium distribution functions (3D)
	static const LVeloc_T LOCAL_VELOC2D [ 9]; ///< Local velocities (2D)
	static const LVeloc_T LOCAL_VELOC3D [27]; ///< Local velocities (3D)

	// Constructor
	Cell (bool Is3D);

	// Initialization
	void Initialize (double Tau, double Rho0, double Vx, double Vy, double Vz=0.0); ///< Tau: Relaxation time, V0: Initial velocity, Rho0: Initial density

	// Access methods
	bool     IsSolid  () const       { return _is_solid; }     ///< Is solid or fluid cell?
	double & F        (size_t Index) { return _f[Index]; }     ///< Return the current value of the distribution function
	double & TmpF     (size_t Index) { return _f_tmp[Index]; } ///< Return the current value of the distribution function
	double   Density  ()           const;                      ///< Calculate the current density of the fluid in this cell
	void     Velocity (Vec3_t & V) const;                      ///< Calculate the current velocity of the fluid in this cell

	// Set methods
	void SetSolid      () { _is_solid = true; }                             ///< TODO
	void SetVelocityBC (BCSide_T Side,double Vx, double Vy, double Vz=0.0); ///< Set velocity boundary conditions (Neumann)
	void SetDensityBC  (BCSide_T Side, double Rho);                         ///< Set density boundary conditions (Dirichlet)
	void ApplyBC       ();                                                  ///< Apply bonudary conditions

	// Methods
	double EqFun    (size_t Index, Vec3_t const & V, double Rho); ///< Calculate the equilibrium distribution function in this cell for the Index direction. Note V/Rho may not be the velocity in this cell.
	void Collide    ();                                           ///< Calculate the collision of the fluid in this cell
	void BounceBack ();                                           ///< 
	void ApplyForce (double Fx, double Fy, double Fz=0);          ///< TODO
	void ApplyGravity(double Gx, double Gy, double Gz=0);         ///< TODO
	void ResetForce () { _bforce = 0.0, 0.0, 0.0; }               ///< TODO

protected:
	// Data
	bool             _is_3D;    ///< Is 3D?
	bool             _is_solid; ///< Solid cell or fluid cell?
	double   const * _w;        ///< Weights for the equilibrium distribution function. 2D: size==9, 3D: size==27
	LVeloc_T const * _c;        ///< Local velocities
	size_t           _nv;       ///< Number of velocities
	double           _tau;      ///< Relaxation time
	Array<double>    _f;        ///< Distribution functions. 2D: size==9, 3D: size==27
	Array<double>    _f_tmp;    ///< Distribution functions. 2D: size==9, 3D: size==27
	Vec3_t           _V0;       ///< Initial velocity from boundary condition
	double           _rho0;     ///< Initial density from boundary condition
	BCSide_T         _bc_side;  ///< Side by where the boundary condiction is applied
	BCType_T         _bc_type;  ///< Type of boundary condiction
	Vec3_t           _bforce;   ///< Body force (gravity)

}; // class Cell


// Set constants
const double Cell::WEIGHTS2D [ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const double Cell::WEIGHTS3D [27] = { 1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9 }; // TODO: Correct this
const Cell::LVeloc_T Cell::LOCAL_VELOC2D[ 9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
const Cell::LVeloc_T Cell::LOCAL_VELOC3D[27] =
{
	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Cell::Cell(bool Is3D)
	: _is_3D    (Is3D),
	  _is_solid (false),
	  _bc_type  (NONE_T)
{
	if (_is_3D)
	{
		_w  = WEIGHTS3D;
		_c  = LOCAL_VELOC3D;
		_nv = 27;
	}
	else
	{
		_w  = WEIGHTS2D;
		_c  = LOCAL_VELOC2D;
		_nv = 9;
	}
	_f.Resize     (_nv);
	_f_tmp.Resize (_nv);
}

inline void Cell::Initialize(double Tau, double Rho0, double Vx, double Vy, double Vz)
{
	// Initial velocity
	Vec3_t v0;  v0 = Vx, Vy, Vz;

	// Initilize distribution function with initial velocity V0 and initial density Rho0
	_tau = Tau;
	for (size_t i=0; i<_nv; i++)
		_f[i] = EqFun(i, v0, Rho0);
}

inline double Cell::Density() const
{
	// Skip if it is solid
	if (_is_solid) return 0.0;

	// Calulate current density
	double density = 0.0;
	for (size_t i=0; i<_nv; i++)
		density += _f[i];
	return density;
}

inline void Cell::Velocity(Vec3_t & V) const
{
	// Skip if it is solid
	if (_is_solid)
	{
		V = 0.0, 0.0, 0.0;
		return;
	}

	double vx      = 0.0;
	double vy      = 0.0;
	double vz      = 0.0;
	double density = Density();

	for (size_t i=0; i<_nv; i++)
		vx += _f[i]*_c[i][0];
	vx /= density;

	for (size_t i=0; i<_nv; i++)
		vy += _f[i]*_c[i][1];
	vy /= density;

	if (_is_3D)
	{
		for (size_t i=0; i<_nv; i++)
			vz += _f[i]*_c[i][2];
		vz /= density;
	}

	// Return curren tvelocity
	V = vx, vy, vz;
}

inline double Cell::EqFun(size_t Index, Vec3_t const & V, double Rho) 
{
	double f_eq;

	if (_is_3D)
	{
		f_eq = 0.0; // TODO: correct this
	}
	else
	{
		Vec3_t v = V + _bforce*_tau/Rho;
		double vxy  = v(0)*_c[Index][0] + v(1)*_c[Index][1];
		double vsqr = v(0)*v(0) + v(1)*v(1);
		f_eq = _w[Index] * Rho * (1.0 + 3.0*vxy + 4.5*vxy*vxy - 1.5*vsqr);
	}

	return f_eq;
}

inline void Cell::SetVelocityBC(BCSide_T Side, double Vx, double Vy, double Vz)
{
	_V0      = Vx, Vy, Vz;
	_bc_side = Side;
	_bc_type = VELOCITY_T;
}

inline void Cell::SetDensityBC(BCSide_T Side, double Rho)
{
	_rho0    = Rho;
	_bc_side = Side;
	_bc_type = DENSITY_T;
}

inline void Cell::Collide()
{
	Vec3_t V;      Velocity (V);
	double   rho = Density  ();
	double omega = 1.0/_tau;

	// Updating _f
	for (size_t k=0; k<_nv; k++)
	{
		double f_eq = EqFun(k, V, rho);
		_f[k] = (1.0-omega)*_f[k] + omega*f_eq;
		if (_f[k]<0) throw new Fatal("F[%d]<0 detected",k);
	}
}

inline void Cell::ApplyBC()
{
	if (_bc_type==VELOCITY_T)
	{
		// Set velocity BC at left side
		if (_bc_side==LEFT_T)
		{
			if (_is_3D)
			{
				// TODO: Implement this
			}
			else
			{
				double vx  = _V0(0);
				double vy  = _V0(1);
				double rho = (_f[0]+_f[2]+_f[4] + 2.0*(_f[3]+_f[6]+_f[7]))/(1.0 - vx);

				_f[1] = _f[3] + 2.0/3.0*rho*vx;
				_f[5] = _f[7] + 1.0/6.0*rho*vx + 0.5*rho*vy - 0.5*(_f[2]-_f[4]);
				_f[8] = _f[6] + 1.0/6.0*rho*vx - 0.5*rho*vy + 0.5*(_f[2]-_f[4]);
			}
		}
		else throw new Fatal("Cell::ApplyBC: This feature is not available yet for BCType_T==%d, BCSide==%d",_bc_type,_bc_side);
		// TODO: Add other directions
	}
	
	if (_bc_type==DENSITY_T)
	{
		// Set density boundary condition at right side
		if (_bc_side==RIGHT_T)
		{
			if (_is_3D)
			{
				// TODO: Implement this
			}
			else
			{
				double rho = _rho0;
				double vx  = -1.0 + (_f[0]+_f[2]+_f[4] + 2.0*(_f[1]+_f[5]+_f[8]))/rho;
				double vy  = 0.0;

				_f[3] = _f[1] - 2.0/3.0*rho*vx; 
				_f[7] = _f[5] - 1.0/6.0*rho*vx - 0.5*rho*vy + 0.5*(_f[2]-_f[4]);
				_f[6] = _f[8] - 1.0/6.0*rho*vx + 0.5*rho*vy - 0.5*(_f[2]-_f[4]);
			}
		}
		else throw new Fatal("Cell::ApplyBC: This feature is not available yet for BCType==%d, BCSide==%d",_bc_type,_bc_side);
		// TODO: Add other directions
	}
}

inline void Cell::BounceBack()
{
	if (_is_3D)
	{
		//TODO: Implement this
	}
	else
	{
		// Copy values to the temporary f_tmp
		for (size_t i=1; i<9; i++) _f_tmp[i] = _f[i];

		// Bonce back
		_f[1] = _f_tmp[3];  _f[2] = _f_tmp[4];
		_f[3] = _f_tmp[1];  _f[4] = _f_tmp[2];
		_f[5] = _f_tmp[7];  _f[6] = _f_tmp[8];
		_f[7] = _f_tmp[5];  _f[8] = _f_tmp[6];
	}
}

inline void Cell::ApplyForce (double Fx, double Fy, double Fz)
{
	_bforce(0) = Fx;
	_bforce(1) = Fy;
	_bforce(2) = Fz;
}

inline void Cell::ApplyGravity (double Gx, double Gy, double Gz)
{
	double Rho = Density();
	_bforce(0) += Gx*Rho;
	_bforce(1) += Gy*Rho;
	_bforce(2) += Gz*Rho;
}

/** Outputs a Cell. */
std::ostream & operator<< (std::ostream & os, Cell * C)
{
	for(size_t i=0; i<9; i++)
		os << "F[" << i << "]=" << C->F(i) << "  ";

	os << "\n";

	return os;
}

}; // namespace LBM

#endif // MECHSYS_LBM_CELL_H