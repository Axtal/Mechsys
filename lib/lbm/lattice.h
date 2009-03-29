/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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
	Lattice (Str_t  FileKey, ///< Key such as "mytest" to be used when generating output files: Ex.: mytest_11.vtk
	         bool   Is3D,    ///< 
	         size_t Nx,      ///< Number of cells along x direction
	         size_t Ny,      ///< Number of cells along y direction
	         size_t Nz=1);   ///< Number of cells along z direction

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
	size_t   Nx        () const { return _nx; }
	size_t   Ny        () const { return _ny; }
	size_t   Nz        () const { return _nz; }
	size_t   NNeigh    () const { return _nneigh; }
	double   Tau       () const { return _tau; }
	double   TotalMass () const;
	Cell   * GetCell   (size_t i, size_t j, size_t k=0);
	Cell   * GetCell   (size_t i) { return _cells[i]; }
	double   Curl      (size_t i, size_t j);

	// Set constants
	Lattice * SetTau       (double Val) { _tau     = Val;  return this; } ///< Set the dynamic viscocity
	Lattice * SetG         (double Val) { _G       = Val;  return this; } ///< Set the interaction strenght
	Lattice * SetGSolid    (double Val) { _G_solid = Val;  return this; } ///< Set the interaction strenght with solids
	Lattice * SetRhoRef    (double Val) { _rho_ref = Val;  return this; } ///< Set the density reference value
	Lattice * SetPsiRef    (double Val) { _psi_ref = Val;  return this; } ///< Set 
	Lattice * SetMultiComp (double Val) { _is_mc   = true; return this; } ///< Set the flag that defines if the analysis is multicomponent

	// Set methods
	void SetGravity    (double Gx, double Gy, double Gz=0.0);
	void SetVelocityBC (size_t i, size_t j, Vec3_t const & V);
	void SetDensityBC  (size_t i, size_t j, double Rho);

	// Methods
	void Solve        (double tIni, double tFin, double dt, double dtOut); ///< Solve
	void Homogenize   ();
	void Stream       ();
	void ApplyBC      ();
	void Collide      ();
	void BounceBack   ();
	void ApplyForce   ();
	void ApplyGravity ();
	void WriteState   (size_t TimeStep); ///< Write the current state to a vtk output file
	double Psi(double Rho) const;        ///< Interaction potential function

protected:
	String       _file_key; ///< File key used as part of the output filenames
	bool         _is_3d;    ///< Flag that defines if the analysis is 3D
	bool         _is_mc;    ///< Flag to define if the analysis is multi-component
	double       _tau;      ///< 
	double       _G;        ///< Interaction strenght
	double       _G_solid;  ///< Interaction strenght with solids
	double       _rho_ref;  ///< Density reference value
	double       _psi_ref;  ///< Interaction potential reference value
	size_t       _nx;       ///< Number of cells in the x direction
	size_t       _ny;       ///< Number of cells in the y direction
	size_t       _nz;       ///< Number of cells in the z direction
	size_t       _size;     ///< Total number of cells
	size_t       _T;        ///< Normalized time
	Vec3_t       _gravity;  ///< Value of the applied gravity
	size_t       _nneigh;   ///< Number of neighbors
	Array<Cell*> _cells;    ///< Array of cell pointers
	Array<Cell*> _bottom;   ///< Array of cell pointers of the bottom line
	Array<Cell*> _top;      ///< Array of cell pointers of the left side
	Array<Cell*> _left;     ///< 
	Array<Cell*> _right;    ///< 
	Array<Cell*> _front;    ///< 
	Array<Cell*> _back;     ///< 
	Array<Cell*> _cpveloc;  ///< Cells with prescribed velocity
	Array<Cell*> _cpdens;   ///< Cells with prescribed density


}; // class Lattice


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Lattice::Lattice(Str_t FileKey, bool Is3D, size_t Nx, size_t Ny, size_t Nz)
	: _file_key (FileKey),
	  _is_3d    (Is3D),
	  _is_mc    (false),
	  _tau      (1.0),
	  _G        (0.0),
	  _G_solid  (0.0),
	  _rho_ref  (1.0),
	  _psi_ref  (1.0),
	  _nx       (Nx),
	  _ny       (Ny),
	  _nz       (Nz),
	  _size     (Nx*Ny*Nz),
	  _T        (0)
{
	// Gravity
	_gravity = 0.0, 0.0, 0.0;

	// Number of neighbours
	_nneigh = (_is_3d ? 27 : 9);

	// Allocate memory for cells
	_cells.Resize (_size);
	for (size_t i=0; i<Nx; ++i)
	for (size_t j=0; j<Ny; ++j)
	{
		size_t n = i+j*Nx;
		_cells[n] = new Cell(n,_is_3d,_tau,i,j,1,Nx,Ny,Nz);
	}

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

inline double Lattice::TotalMass() const
{
	double mass = 0.0;
	for (size_t i=0; i<_size; i++) mass += _cells[i]->Density();
	return mass;
}

inline Cell * Lattice::GetCell(size_t i, size_t j, size_t k)
{
	if (_is_3d)
	{
		// TODO: Implement this
	}
	return _cells[i+_nx*j];
}

inline double Lattice::Curl(size_t i, size_t j)
{
	if (GetCell(i,j)->IsSolid()) return 0.0;
	size_t i1 = (i+1+_nx) % _nx;
	size_t i3 = (i-1+_nx) % _nx;
	size_t j2 = (j+1+_ny) % _ny;
	size_t j4 = (j-1+_ny) % _ny;

	Vec3_t v1; GetCell(i1,j)->Velocity(v1);
	Vec3_t v3; GetCell(i3,j)->Velocity(v3);
	Vec3_t v2; GetCell(i,j2)->Velocity(v2);
	Vec3_t v4; GetCell(i,j4)->Velocity(v4);

	double dvydx = (v1(1)-v3(1))/2.0; // dVy/dx
	double dvxdy = (v2(0)-v4(0))/2.0; // dVx/dy
	
	return dvydx - dvxdy;
}

inline void Lattice::SetGravity(double Gx, double Gy, double Gz)
{
	_gravity = Gx, Gy, Gz;
}

inline void Lattice::SetVelocityBC(size_t i, size_t j, Vec3_t const & V)
{
	LBM::Cell * c = GetCell(i,j);
	_cpveloc.Push (c);
	c->SetVelBC   (V);
}

inline void Lattice::SetDensityBC(size_t i, size_t j, double Rho)
{
	LBM::Cell * c = GetCell(i,j);
	_cpdens.Push (c);
	c->SetRhoBC  (Rho);
}

inline void Lattice::Solve(double tIni, double tFin, double dt, double dtOut)
{
	double t    = tIni;
	double tout = t + dtOut;
	WriteState (_T);
	while (t<tFin)
	{
		ApplyForce   ();
		ApplyGravity ();
		Collide      ();
		BounceBack   ();
		Stream       ();
		ApplyBC      ();
		t += dt;
		if (t>=tout)
		{
			String buf;
			buf.Printf("[1;34mMechSys[0m::LBM::Lattice::Solve: [1;31mt = %g   TotalMass = %g[0m\n", t, TotalMass());
			std::cout << buf;
			_T++;
			tout += dtOut;
			WriteState (_T);
		}
	}
}

inline void Lattice::Stream()
{
	// Stream to temp
	for (size_t i=0; i<_size; ++i)
	{
		LBM::Cell * c = _cells[i];
		for (size_t k=0; k<_nneigh; ++k)
		{
			LBM::Cell * nb = _cells[c->Neigh(k)];
			nb->TmpF(k) = c->F(k);
		}
	}

	// Swap the distribution function values
	for (size_t i=0; i<_size;   i++) 
	for (size_t k=0; k<_nneigh; k++)
		_cells[i]->F(k) = _cells[i]->TmpF(k);
}

inline void Lattice::Homogenize()
{
	// Homogenize to temp
	for (size_t i=0; i<_size; ++i)
	{
		LBM::Cell * c = _cells[i];
		if (c->IsSolid()) continue;
		for (size_t k=0; k<_nneigh; ++k)
		{
			double       f = 0.0;
			for (size_t kk=0; kk<_nneigh; ++kk)
			{
				LBM::Cell * nb = _cells[c->Neigh(kk)];
				f += nb->F(k);
			}
			c->TmpF(k) = f/_nneigh;
		}
	}

	// Swap the distribution function values
	for (size_t i=0; i<_size;   i++)
	for (size_t k=0; k<_nneigh; k++)
		_cells[i]->F(k) = _cells[i]->TmpF(k);
}

inline void Lattice::ApplyBC()
{
	// Cells with prescribed velocity
	for (size_t i=0; i<_cpveloc.Size(); ++i)
	{
		LBM::Cell * c = _cpveloc[i];
		if (c->IsSolid()) continue;
		if (_is_3d) throw new Fatal("Lattice::ApplyBC: 3D simulation is not implemented yet");
		if (c->Left()) // Cell is on the left side
		{
			double rho = (c->F(0)+c->F(2)+c->F(4) + 2.0*(c->F(3)+c->F(6)+c->F(7)))/(1.0-c->VelBC(0));
			c->F(1) = c->F(3) + (2.0/3.0)*rho*c->VelBC(0);
			c->F(5) = c->F(7) + (1.0/6.0)*rho*c->VelBC(0) + 0.5*rho*c->VelBC(1) - 0.5*(c->F(2)-c->F(4));
			c->F(8) = c->F(6) + (1.0/6.0)*rho*c->VelBC(0) - 0.5*rho*c->VelBC(1) + 0.5*(c->F(2)-c->F(4));
		}
	}

	// Cells with prescribed density
	for (size_t i=0; i<_cpdens.Size(); ++i)
	{
		LBM::Cell * c = _cpdens[i];
		if (c->IsSolid()) continue;
		if (_is_3d) throw new Fatal("Lattice::ApplyBC: 3D simulation is not implemented yet");
		if (c->Right()) // Cell is on the right side
		{
			double vx = -1.0 + (c->F(0)+c->F(2)+c->F(4) + 2.0*(c->F(1)+c->F(5)+c->F(8)))/c->RhoBC();
			c->F(3) = c->F(1) - (2.0/3.0)*c->RhoBC()*vx; 
			c->F(7) = c->F(5) - (1.0/6.0)*c->RhoBC()*vx + 0.5*(c->F(2)-c->F(4));
			c->F(6) = c->F(8) - (1.0/6.0)*c->RhoBC()*vx - 0.5*(c->F(2)-c->F(4));
		}
	}
}

inline void Lattice::Collide()
{
	double om = 1.0/_tau;
	for (size_t i=0; i<_size; i++)
	{
		LBM::Cell * c = _cells[i];
		if (c->IsSolid()==false) // not solid
		{
			Vec3_t v;
			if (_is_mc)	v = c->MixVelocity(); // For multi-components analysis
			else            c->Velocity(v);   // Fore one component analysis

			double rho = c->Density  ();
			for (size_t k=0; k<_nneigh; ++k)
			{
				double feq = c->EqFun (k,v,rho);
				c->F(k) = (1.0-om)*c->F(k) + om*feq;
				if (c->F(k)<0.0) throw new Fatal("Lattice::Collide: Cell(%d)->F(%d)<0.0 detected",i,k);
			}
		}
	}
}

inline void Lattice::BounceBack()
{
	if (_is_3d) throw new Fatal("Lattice::BounceBack: 3D simulation is not implemented yet");
	for (size_t i=0; i<_size; i++)
	{
		LBM::Cell * c = _cells[i];
		if (c->IsSolid())
		{
			// Copy values to the temporary f_tmp
			for (size_t i=1; i<_nneigh; i++) c->TmpF(i) = c->F(i);

			// Performing bounce back
			c->F(1) = c->TmpF(3);  c->F(2) = c->TmpF(4);
			c->F(3) = c->TmpF(1);  c->F(4) = c->TmpF(2);
			c->F(5) = c->TmpF(7);  c->F(6) = c->TmpF(8);
			c->F(7) = c->TmpF(5);  c->F(8) = c->TmpF(6);
                                                       
			// Adding the component related with the solids surface velocity ( alpha = 6*w*rho/Cs^2 )
			double   rho = 0.0;
			for (size_t i=0; i<_nneigh; i++) rho += c->F(i);

			c->F(1) += -(6.0*c->W(3)*rho) * (c->C(3,0)*c->VelBC(0) + c->C(3,1)*c->VelBC(1));
			c->F(2) += -(6.0*c->W(4)*rho) * (c->C(4,0)*c->VelBC(0) + c->C(4,1)*c->VelBC(1));
			c->F(3) += -(6.0*c->W(1)*rho) * (c->C(1,0)*c->VelBC(0) + c->C(1,1)*c->VelBC(1));
			c->F(4) += -(6.0*c->W(2)*rho) * (c->C(2,0)*c->VelBC(0) + c->C(2,1)*c->VelBC(1));
			c->F(5) += -(6.0*c->W(7)*rho) * (c->C(7,0)*c->VelBC(0) + c->C(7,1)*c->VelBC(1));
			c->F(6) += -(6.0*c->W(8)*rho) * (c->C(8,0)*c->VelBC(0) + c->C(8,1)*c->VelBC(1));
			c->F(7) += -(6.0*c->W(5)*rho) * (c->C(5,0)*c->VelBC(0) + c->C(5,1)*c->VelBC(1));
			c->F(8) += -(6.0*c->W(6)*rho) * (c->C(6,0)*c->VelBC(0) + c->C(6,1)*c->VelBC(1));

		}
	}
}

inline void Lattice::ApplyForce()
{
	if (_is_3d) throw new Fatal("Lattice::BounceBack: 3D simulation is not implemented yet");
	for (size_t i=0; i<_size; i++)
	{
		LBM::Cell * c = _cells[i];
		Vec3_t F; F   = 0.0, 0.0, 0.0;
		double    psi = Psi(c->Density());
		for (size_t k=1; k<_nneigh; ++k)
		{
			LBM::Cell * nb = _cells[c->Neigh(k)];
			double  nb_psi = (nb->IsSolid() ? 1.0      : Psi(nb->Density()));
			double       G = (nb->IsSolid() ? _G_solid : _G);
			//double  nb_psi =  Psi(nb->Density());
			//double       G =  _G;
			F(0) += -G*psi*c->W(k)*nb_psi*c->C(k,0);
			F(1) += -G*psi*c->W(k)*nb_psi*c->C(k,1);
		}
		c->BForce() = F;
	}
}

inline void Lattice::ApplyGravity()
{
	for (size_t i=0; i<_size; ++i)
	{
		double rho = _cells[i]->Density();
		_cells[i]->BForce() += _gravity*rho;
	}
}

inline void Lattice::WriteState(size_t TimeStep)
{
	// Header
	std::ostringstream oss;
	oss << "# vtk DataFile Version 2.0\n";
	oss << "TimeStep = " << TimeStep << "\n";
	oss << "ASCII\n";
	oss << "DATASET STRUCTURED_POINTS\n";
	oss << "DIMENSIONS "   << _nx << " " << _ny << " " << _nz << "\n";
	oss << "ORIGIN "       << 0   << " " << 0   << " " << 0   << "\n";
	oss << "SPACING "      << 1   << " " << 1   << " " << 1   << "\n";
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
	for (size_t i=0; i<_size; i++)
		oss << _cells[i]->Density() << "\n";

	// Curl field
	oss << "SCALARS Curl float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<_size; ++i)
		oss << Curl(i % _nx, i / _nx) << "\n";

	// Velocity field
	oss << "VECTORS Velocity float\n";
	for (size_t i=0; i<_size; ++i)
	{
		Vec3_t v;  _cells[i]->Velocity(v);
		oss << v(0) << " " << v(1) << " " << v(2) << "\n";
	}

	// Open/create file, write to file, and close file
	String         fn;  fn.Printf("%s_%d.vtk",_file_key.CStr(),TimeStep);
	std::ofstream  of;
	of.open (fn.CStr(), std::ios::out);
	of << oss.str();
	of.close();
}

/* private */
inline double Lattice::Psi(double Rho) const
{
	return _psi_ref*(1.0-exp(-Rho/_rho_ref));
}

}; // namespace LBM

#endif // MECHSYS_LBM_LATTICE_H
