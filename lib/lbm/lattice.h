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
	Array<Cell*> const & Left   () { return _left;   }
	Array<Cell*> const & Right  () { return _right;  }
	Array<Cell*> const & Bottom () { return _bottom; }
	Array<Cell*> const & Top    () { return _top;    }
	Array<Cell*> const & Front  () { return _front;  }
	Array<Cell*> const & Back   () { return _back;   }

	// Access method
	size_t Nx() const { return _nx; }
	size_t Ny() const { return _ny; }
	size_t Nz() const { return _nz; }
	Cell * GetCell (size_t i, size_t j, size_t k=0);

	// Solve
	void Solve      (double tIni, double tFin, double dt, double dtOut); ///< Solve
	void Stream     ();
	void ApplyBC    ();
	void Collide    ();
	void BounceBack ();

protected:
	bool   _is_3d;        ///< TODO:
	double _tau;          ///< TODO:
	size_t _nx;           ///< TODO:
	size_t _ny;           ///< TODO:
	size_t _nz;           ///< TODO:
	double _dL;           ///< TODO:
	size_t _size;         ///< TODO:

	Array<Cell*> _cells;  ///< TODO:
	Array<Cell*> _bottom; ///< TODO:
	Array<Cell*> _top;    ///< TODO:
	Array<Cell*> _left;   ///< TODO:
	Array<Cell*> _right;  ///< TODO:
	Array<Cell*> _front;  ///< TODO:
	Array<Cell*> _back;   ///< TODO:

}; // class Lattice


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Lattice::Lattice(bool Is3D, double Tau, double dL, size_t Nx, size_t Ny, size_t Nz)
{
	// Set data
	_is_3d = Is3D;
	_tau   = Tau;
	_nx    = Nx;
	_ny    = Ny;
	_nz    = Nz;
	_size  = (Is3D ? Nx*Ny*Nz : Nx*Ny);
	_cells.Resize(_size);

	// Allocate memory
	for (size_t i=0; i<_size; i++)
		_cells[i] = new Cell(_is_3d);
}

inline Lattice::~Lattice()
{
	for (size_t i=0; i<_size; i++) delete _cells[i];
}

inline Cell * Lattice::GetCell(size_t i, size_t j, size_t k)
{
	if (_is_3d)
	{
		// TODO: Implement this
	}
	return _cells[i+_nx*j];
}

inline void Lattice::Solve(double tIni, double tFin, double dt, double dtOut)
{
	double t    = tIni;
	double tout = t + dtOut;
	//WriteState (t);
	while (t<tFin)
	{
		Collide    ();
		BounceBack ();
		Stream     ();
		ApplyBC    ();
		t += dt;
		if (t>=tout)
		{
			std::cout << "t = " << t << std::endl;
			//WriteState (t);
			tout += dtOut;
		}
	}
}

inline void Lattice::Stream()
{
	if (_is_3d)
	{
	}
	else
	{
		// Streaming
		Cell::LVeloc_T const * c = Cell::LOCAL_VELOC2D; // Local velocities
		for (size_t i=0; i<_nx; i++)
		for (size_t j=0; j<_ny; j++)
		{
			for (size_t k=0; k<9; k++)
			{
				int next_i = i + c[k][0];
				int next_j = j + c[k][1];
				if (next_i==-1)                    next_i = _nx-1;
				if (next_i==static_cast<int>(_nx)) next_i = 0;
				if (next_j==-1)                    next_j = _ny-1;
				if (next_j==static_cast<int>(_ny)) next_j = 0;
				_cells[next_i+_nx*next_j]->TmpF(k) = _cells[i+_nx*j]->F(k);
			}
		}
		
		// Swap the distribution function values
		for (size_t i=0; i<_size; i++) 
		for (size_t k=0; k<9;     k++)
			_cells[i]->F(k) = _cells[i]->TmpF(k);
	}
}

inline void Lattice::ApplyBC()
{
	for (size_t i=0; i<_size; i++) _cells[i]->ApplyBC();
}

inline void Lattice::Collide()
{
	for (size_t i=0; i<_size; i++)
	{
		if (_cells[i]->IsSolid()==false)
			_cells[i]->Collide();
	}
}

inline void Lattice::BounceBack()
{
	for (size_t i=0; i<_size; i++)
	{
		if (_cells[i]->IsSolid())
			_cells[i]->BounceBack();
	}
}


}; // namespace LBM

#endif // MECHSYS_LBM_LATTICE_H
