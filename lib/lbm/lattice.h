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
	Lattice (bool Is3D, double Tau, double dL, size_t Nx, size_t Ny, size_t Nz=0); ///< Nx,Ny,Nz number of cells along each direction. dL = Length of each size

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
	void Stream  ();
	void ApplyBC ();
	void Collide ();
	void BounceBack ();

protected:
	bool   _is_3d;
	double _tau;
	int    _nx;
	int    _ny;
	int    _nz;
	double _dL;
	int    _size;

	Array<Cell*> _cells;
	Array<Cell*> _bottom;
	Array<Cell*> _top;
	Array<Cell*> _left;
	Array<Cell*> _right;
	Array<Cell*> _front;
	Array<Cell*> _back;

}; // class Lattice


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Lattice::Lattice(bool Is3D, double Tau, double dL, size_t Nx, size_t Ny, size_t Nz)
{
	_is_3d = Is3D;
	_tau = Tau;
	_nx  = Nx;
	_ny  = Ny;
	_nz  = Nz;
	_size = Nx*Ny; //TODO: update for 3D
	_cells.Resize(Nx*Ny);

	for (int i=0; i<_size; i++)
	{
		if (Is3D)
			_cells[i] = new Cell3D;
		else
			_cells[i] = new Cell2D;
	}

}

inline Lattice::~Lattice()
{
	for (int i=0; i<_size; i++) delete _cells[i];
}

inline Cell & Lattice::GetCells(size_t i, size_t j, size_t k)
{
	return *_cells[i+_nx*j]; //TODO: update for 3D
}

inline void Lattice::Solve(double tIni, double tFin, double dt, double dtOut)
{
}

inline void Lattice::Stream()
{
	if (_is_3d); //TODO

	// Streaming
	Cell::LVeloc_T const * c = Cell::LOCAL_VELOC2D; // Local velocities
	for (int i=0; i<_nx; i++)
		for (int j=0; j<_ny; j++)
		{
			for (int k=0; k<9; k++)
			{
				int next_i = i + c[k][0];
				int next_j = j + c[k][1];
				if (next_i==-1)  next_i = _nx-1;
				if (next_i==_nx) next_i = 0;
				if (next_j==-1)  next_j = _ny-1;
				if (next_j==_ny) next_j = 0;
				_cells[next_i+_nx*next_j]->TmpF(k) = _cells[i+_nx*j]->F(k);
			}
		}
	
	// Swap the distribution function values
	for (int i=0; i<_size; i++) 
		for (int k=0; k<9; k++)
			_cells[i]->F(k) = _cells[i]->TmpF(k);
}

inline void Lattice::ApplyBC()
{
	for (int i=0; i<_size; i++) _cells[i]->ApplyBC();
}

inline void Lattice::Collide()
{
	for (int i=0; i<_size; i++)
	{
		if (_cells[i]->IsSolid()==false)
			_cells[i]->Collide();
	}

}

inline void Lattice::BounceBack()
{
	for (int i=0; i<_size; i++)
		if (_cells[i]->IsSolid())
			_cells[i]->BounceBack();
}


}; // namespace LBM

#endif // MECHSYS_LBM_LATTICE_H
