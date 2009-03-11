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

	// Constructor
	Cell (double Tau, SETNEIGHS) : _tau(Tau), _is_solid(false) {}

	// Initialization
	void Initialize (Vec3_t const & V0, double Rho0); ///< Tau: Relaxation time, V0: Initial velocity, Rho0: Initial density

	// Access methods
	bool   IsSolid () const       { return _is_solid; }     ///< Is solid or fluid cell?
	size_t NNeigh  () const       { return _neigh.Size(); } ///< Number of neighbours
	Cell * Neigh   (size_t Index) { return _neigh[Index]; } ///< Return the access to a neighbour cell
	double F       (size_t Index) { return _f[Index]; }     ///< Return the current value of the distribution function
	void   Veloc   (Vec3_t & V) const;                      ///< Calculate the velocity of the fluid in this cell
	double Density ()           const;                      ///< Calculate the density of the fluid in this cell

	// Set methods
	void SetVelocBC   (); ///< Set velocity boundary conditions (Neumann)
	void SetDensityBC (); ///< Set density boundary conditions (Dirichlet)

	// Methods
	double EqFun   (size_t Index, Vec3_t const & V, double Rho); ///< Calculate the equilibrium distribution function in this cell for the Index direction. Note V/Rho may not be the velocity in this cell.
	void   Stream  ();                                           ///< Calculate the movement of the fluid in this cell
	void   Collide ();                                           ///< Calculate the collision of the fluid in this cell

protected:
	// Data
	double           _tau;      ///< Relaxation time
	bool             _is_solid; ///< Solid cell or fluid cell?
	Array<Cell*>     _neigh;    ///< Neighbours. Note: First index [0] referes to itself
	Array<double>    _f;        ///< Distribution functions. 2D: size==9, 3D: size==27
	double   const * _w;        ///< Weights for the equilibrium distribution function. 2D: size==9, 3D: size==27
	LVeloc_T const * _c;        ///< Local velocities

}; // class Cell

const double Cell::WEIGHTS2D[ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const double Cell::WEIGHTS3D[27] = { 1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9, 1,2,3,4,5,6,7,8,9 }; // TODO: Correct this
const Cell::LVeloc_T Cell::LOCAL_VELOC2D[9] = {
   	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
};
const Cell::LVeloc_T Cell::LOCAL_VELOC3D[27] = {
   	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
   	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0},
   	{0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}
};

class Cell2D : public Cell
{
public:
	Cell2D () { _w = WEIGHTS2D; _c = LOCAL_VELOC2D; }
private:
}; // class Cell2D


class Cell3D : public Cell
{
public:
	Cell3D () { _w = WEIGHTS3D; _c = LOCAL_VELOC3D; }
private:
}; // class Cell3D



/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////



}; // namespace LBM

#endif // MECHSYS_LBM_CELL_H
