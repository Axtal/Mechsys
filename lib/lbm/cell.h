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
	// Typedefs
	typedef double LVeloc_T[3]; ///< Local velocities type

	// Constants
	static const double   WEIGHTS2D     [ 9]; ///< Weights for the equilibrium distribution functions (2D)
	static const double   WEIGHTS3D     [27]; ///< Weights for the equilibrium distribution functions (3D)
	static const LVeloc_T LOCAL_VELOC2D [ 9]; ///< Local velocities (2D)
	static const LVeloc_T LOCAL_VELOC3D [27]; ///< Local velocities (3D)
	static const size_t   OPPOSITE2D    [ 9]; ///< Opposite directions

	// Constructor
	Cell (bool Is3D, double Tau, long i, long j, long k, size_t Nx, size_t Ny, size_t Nz=1);

	// Set methods
	void Initialize (double Rho0, Vec3_t const & V0);                   ///< V0: Initial velocity, Rho0: Initial density
	void SetSolid   (bool IsSolid=true)     { _is_solid = IsSolid; }    ///< Set solid cell
	void SetSolid   (double Vx, double Vy, double Vz=0.0) { _is_solid = true; _vel_bc=Vx,Vy,Vz; }  ///< Set solid cell
	void SetRhoBC   (double Rho)            { _rho_bc = Rho; } ///< Set density boundary condition
	void SetVelBC   (Vec3_t const &  V)     { _vel_bc = V;   } ///< Set velocity boundary condition

	// Access methods
	bool     IsSolid ()                   const { return _is_solid;  } ///< Is solid or fluid cell?
	double   RhoBC   ()                   const { return _rho_bc;    } ///< Initial density
	double   VelBC   (size_t i)           const { return _vel_bc(i); } ///< Component of initial velocity
	double   W       (size_t k)           const { return _w[k];      } ///< Component of weight
	double   C       (size_t k, size_t i) const { return _c[k][i];   } ///< Component of local velocity
	size_t   Opp     (size_t k)           const { return OPPOSITE2D[k];      } ///< Calculate the opposite direction
	size_t   Neigh   (size_t k)           const { return _neigh[k];  } ///< Returns the index of neighbour k
	bool     Left    () const { return static_cast<bool>(_side[0]);  } ///< Is this cell on the left side of lattice ?
	bool     Right   () const { return static_cast<bool>(_side[1]);  } ///< Is this cell on the right side of lattice ?
	bool     Bottom  () const { return static_cast<bool>(_side[2]);  } ///< Is this cell on the bottom side of lattice ?
	bool     Top     () const { return static_cast<bool>(_side[3]);  } ///< Is this cell on the top side of lattice ?
	double & F       (size_t Index)         { return _f[Index];      } ///< Return the current value of the distribution function
	double & TmpF    (size_t Index)         { return _f_tmp[Index];  } ///< Return the current value of the distribution function
	Vec3_t & BForce  ()                     { return _bforce;        } ///< Body force

	// Methods
	double   Density  () const;                                       ///< Calculate the current density of the fluid in this cell
	void     Velocity (Vec3_t & V) const;                             ///< Calculate the current velocity of the fluid in this cell
	double   EqFun    (size_t k, Vec3_t const & V, double Rho) const; ///< Calculate the equilibrium distribution function in this cell for the Index direction. Note V/Rho may not be the velocity in this cell.
	Vec3_t & MixVelocity() { return _mix_vel; }


protected:
	// Data
	bool             _is_3d;    ///< Is 3D?
	bool             _is_solid; ///< Solid cell or fluid cell?
	double           _tau;      ///< TODO
	size_t           _nneigh;   ///< Number of neighbours
	double           _rho_bc;   ///< Initial density from boundary condition
	Vec3_t           _vel_bc;   ///< Initial velocity from boundary condition
	double   const * _w;        ///< Weights for the equilibrium distribution function. 2D: size==9, 3D: size==27
	LVeloc_T const * _c;        ///< Local velocities
	Array<double>    _f;        ///< Distribution functions. 2D: size==9, 3D: size==27
	Array<double>    _f_tmp;    ///< Distribution functions. 2D: size==9, 3D: size==27
	Vec3_t           _bforce;   ///< Body force (gravity)
	Array<size_t>    _side;     ///< Which side on the boundary this cell is located
	Array<size_t>    _neigh;    ///< Indices of neighbour cells
	Vec3_t           _mix_vel;  ///< Mixed velocity for multicomponent analysis

}; // class Cell


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Cell::Cell(bool Is3D, double Tau, long i, long j, long k, size_t Nx, size_t Ny, size_t Nz)
	: _is_3d    (Is3D),
	  _is_solid (false),
	  _tau      (Tau)
{
	_rho_bc = 0.0;
	_vel_bc = 0.0, 0.0, 0.0;
	_bforce = 0.0, 0.0, 0.0;
	if (_is_3d)
	{
		_nneigh = 27;
		_w      = WEIGHTS3D;
		_c      = LOCAL_VELOC3D;
		_f    .Resize    (_nneigh);
		_f_tmp.Resize    (_nneigh);
		_neigh.Resize    (_nneigh);
		_side .Resize    (6);
		_side .SetValues (false);
		throw new Fatal("Cell::Cell: 3D simulation is not implemented yet");
	}
	else
	{
		// Set constants
		_nneigh = 9;
		_w      = WEIGHTS2D;
		_c      = LOCAL_VELOC2D;
		_f    .Resize    (_nneigh);
		_f_tmp.Resize    (_nneigh);
		_neigh.Resize    (_nneigh);
		_side .Resize    (4);
		_side .SetValues (false);

		// Find side in lattice
		long nx = static_cast<long>(Nx);
		long ny = static_cast<long>(Ny);
		if (i==0)    _side[0] = true; // left
		if (i==nx-1) _side[1] = true; // right
		if (j==0)    _side[2] = true; // bottom
		if (j==ny-1) _side[3] = true; // top

		// Set neighbours
		for (size_t k=0; k<_nneigh; ++k)
		{
			long p = i + _c[k][0]; // neighbour i
			long q = j + _c[k][1]; // neighbour j
			if (p==-1) p = nx-1;
			if (p==nx) p = 0;
			if (q==-1) q = ny-1;
			if (q==ny) q = 0;
			_neigh[k] = p + q*nx;
		}
	}
}

inline void Cell::Initialize(double Rho0, Vec3_t const & V0)
{
	for (size_t k=0; k<_nneigh; k++) _f[k] = EqFun(k, V0, Rho0);
}

inline double Cell::Density() const
{
	// Skip if it is solid
	if (_is_solid) return 0.0;
	else
	{
		// Calulate current density
		double rho = 0.0;
		for (size_t k=0; k<_nneigh; k++) rho += _f[k];
		return rho;
	}
}

inline void Cell::Velocity(Vec3_t & V) const
{
	// Skip if it is solid
	if (_is_solid)
	{
		V = 0.0, 0.0, 0.0;
		return;
	}
	else
	{
		double rho = Density();
		V = 0.0, 0.0, 0.0;
		for (size_t k=0; k<_nneigh; ++k)
		{
			            V(0) += _f[k]*_c[k][0]/rho;
			            V(1) += _f[k]*_c[k][1]/rho;
			if (_is_3d) V(2) += _f[k]*_c[k][2]/rho;
		}
	}
}

inline double Cell::EqFun(size_t k, Vec3_t const & V, double Rho) const
{
	if (_is_3d) throw new Fatal("Cell::EqFun: 3D simulation is not implemented yet");
	else
	{
		Vec3_t v;  v = V + _bforce*_tau/Rho;
		double vdotc = v(0)*_c[k][0] + v(1)*_c[k][1];
		double vdotv = v(0)*v(0) + v(1)*v(1);
		return _w[k]*Rho*(1.0 + 3.0*vdotc + 4.5*vdotc*vdotc - 1.5*vdotv);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////// Constants /////


const double Cell::WEIGHTS2D [ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const double Cell::WEIGHTS3D [27] = { 1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9 }; // TODO: Correct this
const Cell::LVeloc_T Cell::LOCAL_VELOC2D[ 9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
const Cell::LVeloc_T Cell::LOCAL_VELOC3D[27] =
{
	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}
};
const size_t Cell::OPPOSITE2D[ 9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

}; // namespace LBM

#endif // MECHSYS_LBM_CELL_H
