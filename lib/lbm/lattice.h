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

// Std Lib
#include <fstream>
#include <sstream> // for std::ostringstream

// MechSys
#include "lbm/cell.h"

namespace LBM
{

typedef char const * Str_t;

class Lattice
{
public:
	// Constructor
	Lattice (Str_t FileKey, bool Is3D, double Tau, double dL, size_t Nx, size_t Ny, size_t Nz=1); ///< Nx,Ny,Nz number of cells along each direction. dL = Length of each size

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

	// Methods
	void Solve      (double tIni, double tFin, double dt, double dtOut); ///< Solve
	void Stream     ();
	void ApplyBC    ();
	void Collide    ();
	void BounceBack ();
	void ApplyForce (double Fx, double Fy, double Fz=0);
	void WriteState (size_t TimeStep); ///< TODO

protected:
	String _file_key;     ///< TODO:
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


inline Lattice::Lattice(Str_t FileKey, bool Is3D, double Tau, double dL, size_t Nx, size_t Ny, size_t Nz)
	: _file_key (FileKey),
	  _is_3d    (Is3D),
	  _tau      (Tau),
	  _nx       (Nx),
	  _ny       (Ny),
	  _nz       (Nz),
	  _dL       (dL),
	  _size     (Nx*Ny*Nz)
{
	// Allocate memory for cells
	_cells.Resize (_size);
	for (size_t i=0; i<_size; i++)
		_cells[i] = new Cell(_is_3d);

	// Set auxiliary arrays
	if (_is_3d) throw new Fatal("Lattice::Lattice: 3D Lattice is not available yet :-(");
	else
	{
		_bottom.Resize (_nx);
		_top   .Resize (_nx);
		_left  .Resize (_ny);
		_right .Resize (_ny);
		for (size_t i=0; i<_nx; ++i)
		{
			_bottom[i] = GetCell (i,0);
			_top   [i] = GetCell (i,_ny-1);
		}
		for (size_t i=0; i<_ny; ++i)
		{
			_left [i] = GetCell (0,i);
			_right[i] = GetCell (_nx-1,i);
		}
	}
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
	size_t T    = 0;         // Normalized time (integers) used for output only
	double t    = tIni;      //
	double tout = t + dtOut; //
	WriteState (T);
	while (t<tFin)
	{
		ApplyForce (0.0001,0);
		Collide    ();
		BounceBack ();
		Stream     ();
		ApplyBC    ();
		t += dt;
		if (t>=tout)
		{
			std::cout << "[1;34mMechSys[0m::LBM::Lattice::Solve: [1;31mt = " << t << "[0m\n";
			T++;
			WriteState (T);
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

inline void Lattice::ApplyForce(double Fx, double Fy, double Fz)
{
	for (size_t i=0; i<_size; i++)
		_cells[i]->ApplyForce(Fx, Fy, Fz);
}

inline void Lattice::WriteState(size_t TimeStep)
{
	// Open/create file
	String         fn;  fn.Printf("%s_%d.vtk",_file_key.CStr(),TimeStep);
	std::ofstream  of;
	of.open(fn.CStr(), std::ios::out);

	// Header
	std::ostringstream oss;
	oss << "# vtk DataFile Version 2.0\n";
	oss << "TimeStep = " << TimeStep << "\n";
	oss << "ASCII\n";
	oss << "DATASET STRUCTURED_POINTS\n";
	oss << "DIMENSIONS "   << _nx << " " << _ny << " " << _nz << "\n";
	oss << "ORIGIN "       << 0   << " " << 0   << " " << 0   << "\n";
	oss << "SPACING "      << _dL << " " << _dL << " " << _dL << "\n";
	oss << "POINT_DATA "   << _size << "\n";

	// Solid cells
	oss << "SCALARS Geom float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<_size; ++i)
	{
		if (_cells[i]->IsSolid()) oss << "1.0\n";
		else                      oss << "0.0\n";
	}

	// Density field
	oss << "SCALARS Density float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<_size; ++i)
		oss << _cells[i]->Density() << "\n";

	// Velocity field
	oss << "VECTORS Velocity float\n";
	for (size_t i=0; i<_size; ++i)
	{
		Vec3_t v;  _cells[i]->Velocity(v);
		oss << v(0) << " " << v(1) << " " << v(2) << "\n";
	}

	// Write to file and close file
	of << oss.str();
	of.close();
}

}; // namespace LBM

#endif // MECHSYS_LBM_LATTICE_H
