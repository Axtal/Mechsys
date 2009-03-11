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

#ifndef MECHSYS_LBM_CELL_H
#define MECHSYS_LBM_CELL_H

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// MechSys
#include "fem/node.h"

typedef blitz::TinyVector<double,3> Vec3_t;

namespace LBM
{

class Cell
{
public:
	// Enums
	enum BCSide_T { NORTH_T=0, EAST_T=1, SOUTH_T=2, WEST_T=3 };
	enum BCType_T { NEUMANN_T=0, DIRICHLET_T=1 };

	// Typedefs
	typedef double LVeloc_T[3]; ///< Local velocities type

	// Constants
	static const double   WEIGHTS2D     [ 9]; ///< Weights for the equilibrium distribution functions (2D)
	static const double   WEIGHTS3D     [27]; ///< Weights for the equilibrium distribution functions (3D)
	static const LVeloc_T LOCAL_VELOC2D [ 9]; ///< Local velocities (2D)
	static const LVeloc_T LOCAL_VELOC3D [27]; ///< Local velocities (3D)

	// Initialization
	void Initialize (double Tau, Vec3_t const & V0, double Rho0); ///< Tau: Relaxation time, V0: Initial velocity, Rho0: Initial density

	// Access methods
	bool   IsSolid () const       { return _is_solid; }     ///< Is solid or fluid cell?
	size_t NNeigh  () const       { return _neigh.Size(); } ///< Number of neighbours
	Cell * Neigh   (size_t Index) { return _neigh[Index]; } ///< Return the access to a neighbour cell
	double & F     (size_t Index) { return _f[Index]; }     ///< Return the current value of the distribution function
	double & TmpF  (size_t Index) { return _f_tmp[Index]; } ///< Return the current value of the distribution function
	double Density ()           const;                      ///< Calculate the density of the fluid in this cell
	void   Veloc   (Vec3_t & V) const;                      ///< Calculate the velocity of the fluid in this cell

	// Set methods
	void SetVelocBC    (BCSide_T Side,double Vx, double Vy, double Vz); ///< Set velocity boundary conditions (Neumann)
	void SetDensityBC  (BCSide_T Side, double Rho);                     ///< Set density boundary conditions (Dirichlet)
	void ApplyBC       ();                                              ///< Apply bonudary conditions
	void NeumannWest   ();                                              ///< Set velocity BC by west
	void DirichletEast ();                                              ///< Set density boundary condition by east

	// Methods
	double EqFun    (size_t Index, Vec3_t const & V, double Rho); ///< Calculate the equilibrium distribution function in this cell for the Index direction. Note V/Rho may not be the velocity in this cell.
	//void   Stream  ();                                          ///< Calculate the movement of the fluid in this cell TODO:implement in Lattice class
	void Collide    ();                                           ///< Calculate the collision of the fluid in this cell
	void BounceBack ();                                           ///< 

protected:
	// Data
	double           _tau;      ///< Relaxation time
	bool             _is_solid; ///< Solid cell or fluid cell?
	Array<Cell*>     _neigh;    ///< Neighbours. Note: First index [0] referes to itself
	Array<double>    _f;        ///< Distribution functions. 2D: size==9, 3D: size==27
	Array<double>    _f_tmp;    ///< Distribution functions. 2D: size==9, 3D: size==27
	double   const * _w;        ///< Weights for the equilibrium distribution function. 2D: size==9, 3D: size==27
	LVeloc_T const * _c;        ///< Local velocities
	int              _nv;       ///< Number of velocities
	bool             _is_3D;    ///< Is 3D?
	Vec3_t           _V0;       ///< Initial velocity from boundary condition
	double           _rho0;     ///< Initial density from boundary condition

	BCSide_T         _bc_side;  ///< Side by where the boundary condiction is applied
	BCType_T         _bc_type;  ///< Type of boundary condiction

}; // class Cell

const double Cell::WEIGHTS2D[ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const double Cell::WEIGHTS3D[27] = { 1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9 }; // TODO: Correct this
const Cell::LVeloc_T Cell::LOCAL_VELOC2D[9] = {
	{0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0}
};
const Cell::LVeloc_T Cell::LOCAL_VELOC3D[27] = {
   	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
   	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
   	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}
};

class Cell2D : public Cell
{
public:
	Cell2D () { _w = WEIGHTS2D; _c = LOCAL_VELOC2D; _nv=9; _f.Resize(_nv); _f_tmp.Resize(_nv); _is_3D=false; }
private:
}; // class Cell2D

class Cell3D : public Cell
{
public:
	Cell3D () { _w = WEIGHTS3D; _c = LOCAL_VELOC3D; _nv=27; _f.Resize(_nv); _f_tmp.Resize(_nv); _is_3D=true;  }
private:
}; // class Cell3D


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////



void Cell::Initialize(double Tau, Vec3_t const & V0, double Rho0)
{
	_tau = Tau;
	for (int i=0; i<_nv; i++)
		_f[i] = EqFun(i, V0, Rho0);

}

double Cell::Density() const
{
	double density=0.0;
	for (int i=0; i<_nv; i++)
		density += _f[i];
	return density;
}

void Cell::Veloc(Vec3_t & V) const
{
	double vx=0;
	double vy=0;
	double vz=0;

	double density=Density();

	for (int i=0; i<_nv; i++)
		vx += _f[i]*_c[i][0];
	vx /= density;

	for (int i=0; i<_nv; i++)
		vy += _f[i]*_c[i][1];
	vy /= density;

	if (_is_3D)
	{
		for (int i=0; i<_nv; i++)
			vz += _f[i]*_c[i][2];
		vz /= density;
	}

	V = vx, vy, vz;
}

double Cell::EqFun(size_t Index, Vec3_t const & V, double Rho) 
{
	double f_eq;

	if (_is_3D) {} //TODO 
	else
	{
		double vxy  = V(0)*_c[Index][0] + V(1)*_c[Index][1];
		double vsqr = V(0)*V(0) + V(1)*V(1);
		f_eq = _w[Index] * Rho * ( 1.0 + 3.0*vxy + 4.5*vxy*vxy - 1.5*vsqr);
	}

	return f_eq;
}

void Cell::SetVelocBC(BCSide_T Side,double Vx, double Vy, double Vz)
{
	_V0 = Vx, Vy, Vz;
	_bc_side = Side;
	_bc_type = NEUMANN_T;
}

void Cell::SetDensityBC(BCSide_T Side, double Rho)
{
	_rho0 = Rho;
	_bc_side = Side;
	_bc_type = DIRICHLET_T;
}

void Cell::Collide()
{
	Vec3_t V;      Veloc(V);
	double   rho = Density();
	double omega = 1.0/_tau;

	for (int k=0; k<_nv; k++) 
	{
		double f_eq = EqFun(k, V, rho);
		_f[k] = (1.0-omega)*_f[k] + omega*f_eq;
	}
}

void Cell::ApplyBC()
{
	if (_bc_type==NEUMANN_T)
	{
		if (_bc_side==WEST_T) NeumannWest();
	}
	
	if (_bc_type==DIRICHLET_T)
	{
		if (_bc_side==EAST_T) DirichletEast();
	}

}

void Cell::NeumannWest()
{
	if (_is_3D); //TODO

	double vx = _V0(0);
	double vy = _V0(1);

	double rho = (_f[0]+_f[2]+_f[4] + 2.0*(_f[3]+_f[6]+_f[7]))/(1.0 - vx);

	_f[1] = _f[3] + 2.0/3.0*rho*vx;
	_f[5] = _f[7] + 1.0/6.0*rho*vx + 0.5*rho*vy - 0.5*(_f[2]-_f[4]);
	_f[8] = _f[6] + 1.0/6.0*rho*vx - 0.5*rho*vy + 0.5*(_f[2]-_f[4]);
}

void Cell::DirichletEast()
{
	if (_is_3D); //TODO

	double rho = _rho0;

	double vx = -1.0 + (_f[0]+_f[2]+_f[4] + 2.0*(_f[1]+_f[5]+_f[8]))/rho;
	double vy = 0.0;

	_f[3] = _f[1] - 2.0/3.0*rho*vx; 
	_f[7] = _f[5] - 1.0/6.0*rho*vx - 0.5*rho*vy + 0.5*(_f[2]-_f[4]);
	_f[6] = _f[8] - 1.0/6.0*rho*vx + 0.5*rho*vy - 0.5*(_f[2]-_f[4]);
}

void Cell::BounceBack()
{
	if (_is_3D); //TODO

	for (int i=1; i<9; i++) _f_tmp[i] = _f[i]; // copy values to the temporary f_tmp

	// Bonce back
	_f[1] = _f_tmp[3]; _f[2] = _f_tmp[4];
	_f[3] = _f_tmp[1]; _f[4] = _f_tmp[2];
	_f[5] = _f_tmp[7]; _f[6] = _f_tmp[8];
	_f[7] = _f_tmp[5]; _f[8] = _f_tmp[6];

}


}; // namespace LBM

#endif // MECHSYS_LBM_CELL_H
