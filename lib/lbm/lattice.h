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

#ifndef MECHSYS_LBM_LATTICE_H
#define MECHSYS_LBM_LATTICE_H

// MechSys
#include "lbm/cell.h"

namespace LBM
{

class Lattice
{
public:
	// Constructor
	Lattice (bool Is3D, double Tau, size_t Nx, size_t Ny, double dL, size_t Nz=0); ///< Nx,Ny,Nz number of cells along each direction. dL = Length of each size

	// Destructor
	~Lattice ();

	// Access methods
	Array<Cell*> const & Bottom () { return _bottom; }
	Array<Cell*> const & Left   () { return _left;   }
	Array<Cell*> const & Top    () { return _top;    }
	Array<Cell*> const & Right  () { return _right;  }
	Array<Cell*> const & Front  () { return _front;  }
	Array<Cell*> const & Back   () { return _back;   }

	// Access method
	Cell & GetCells (size_t i, size_t j, size_t k=0);

	// Solve
	void Solve (double tIni, double tFin, double dt, double dtOut); ///< Solve

protected:
	Array<Cell*> _cells;
	Array<Cell*> _bottom;
	Array<Cell*> _top;
	Array<Cell*> _left;
	Array<Cell*> _right;
	Array<Cell*> _front;
	Array<Cell*> _back;

}; // class Lattice


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Lattice::Lattice(bool Is3D, double Tau, size_t Nx, size_t Ny, double dL, size_t Nz)
{
}

inline Lattice::~Lattice()
{
}

inline Cell & Lattice::GetCells(size_t i, size_t j, size_t k)
{
}

inline void Lattice::Solve(double tIni, double tFin, double dt, double dtOut)
{
}


}; // namespace LBM

#endif // MECHSYS_LBM_LATTICE_H
