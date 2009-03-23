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

#ifndef MECHSYS_LBM_MIXTURE_H
#define MECHSYS_LBM_MIXTURE_H

// Std Lib
#include <fstream>
#include <sstream> // for std::ostringstream

// MechSys
#include "lbm/cell.h"

namespace LBM
{

typedef char const * Str_t;

class Mixture
{
public:
	// Constructor
	Mixture (Str_t FileKey, bool Is3D, size_t NComp, size_t Nx, size_t Ny, size_t Nz=1); ///< Nx,Ny,Nz number of cells along each direction. dL = Length of each size

	// Destructor
	~Mixture ();

	// Access method
	size_t Nx() const { return _nx; }
	size_t Ny() const { return _ny; }
	size_t Nz() const { return _nz; }


	// Methods
	Lattice * GetLattice(size_t Index) { return _latts[Index]; }
	void SetMixG        (double G) { _G_mix=G; }
	void MixVelocity    (long i, long j, Vec3_t & V);
	void SetMixVelocity ();
	void ApplyMixForce  ();
	void SetGravity     (double Gx, double Gy, double Gz);
	void Solve          (double tIni, double tFin, double dt, double dtOut); ///< Solve
	void WriteState     (size_t TimeStep); ///< TODO
	
protected:
	String _file_key;     ///< TODO:
	bool   _is_3d;        ///< TODO:
	size_t _n_comp;       ///< TODO: MultiPhase ?
	size_t _nx;           ///< TODO:
	size_t _ny;           ///< TODO:
	size_t _nz;           ///< TODO:
	size_t _size;         ///< TODO:
	size_t _T;            ///< Normalized time
	Array<Lattice *> _latts;  ///< Component Lattices

	double _G_mix;       ///< TODO:

private:
	double _psi(double Density) const;

}; // class Mixture


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Mixture::Mixture(Str_t FileKey, bool Is3D, size_t NComp, size_t Nx, size_t Ny, size_t Nz)
	: _file_key (FileKey),
	  _is_3d    (Is3D),
	  _n_comp   (NComp),
	  _nx       (Nx),
	  _ny       (Ny),
	  _nz       (Nz),
	  _size     (Nx*Ny*Nz),
	  _T        (0),
	  _G_mix    (0.0)
{
	_latts.Resize(_n_comp);
	for (size_t n=0; n<_n_comp; n++)
	{
		_latts[n] = new Lattice(FileKey,/*Is3D*/false, Nx, Ny);
		_latts[n]->SetMultiComp(true);
	}
}

inline Mixture::~Mixture()
{
	for (size_t n=0; n<_n_comp; n++) delete _latts[n];
}

inline void Mixture::MixVelocity(long i, long j, Vec3_t & V)
{
	Vec3_t num; num = 0.0, 0.0, 0.0;
	double      den = 0.0; // denominator
	//size_t nneigh=_latts[0]->NNeigh();

	for (size_t n=0; n<_n_comp; n++)
	{
		Cell   * c = _latts[n]->GetCell(i,j);
		double tau = _latts[n]->Tau();
		double rho = c->Density();
		Vec3_t   v;  c->Velocity(v);
		num += v*rho/tau;
		//for (size_t k=0; k<nneigh; k++)
		//{
		//				V(0) += c->F(k)*c->C(k,0)/tau;
		//				V(1) += c->F(k)*c->C(k,1)/tau;
		//	if (_is_3d) V(2) += c->F(k)*c->C(k,2)/tau;
		//}

		den += rho/tau;
	}

	V = num/den;
}

inline void Mixture::ApplyMixForce()
{
	for (size_t n=0; n<_n_comp; n++)
	for (size_t i=0; i<_size; i++)
	{
		// n: means the index for the current lattice 
		// m: means the index for the another lattice 

		size_t m = (n==0)?1:0;

		LBM::Cell * c = _latts[n]->GetCell(i);
		Vec3_t F; F   = 0.0, 0.0, 0.0;
		double    psi = c->Density();
		for (size_t k=1; k<_latts[m]->NNeigh(); ++k)
		{
			LBM::Cell * nb = _latts[m]->GetCell(c->Neigh(k));
			double  nb_psi = (nb->IsSolid() ? 0.0 : nb->Density()); // TODO Check the value of nb_psi	

			F(0) += -_G_mix*psi*c->W(k)*nb_psi*c->C(k,0);
			F(1) += -_G_mix*psi*c->W(k)*nb_psi*c->C(k,1);
		}
		c->BForce() += F;
	}
}

inline void Mixture::SetMixVelocity()
{
	for (size_t n=0; n<_n_comp; n++)
	{
		for (size_t i=0; i<Nx(); i++)
		for (size_t j=0; j<Ny(); j++)
		{
			Vec3_t V; MixVelocity(i,j,V);
			_latts[n]->GetCell(i,j)->MixVelocity() = V;
		}
	}
}

inline void Mixture::SetGravity(double Gx, double Gy, double Gz)
{
	for (size_t n=0; n<_n_comp; n++)
		_latts[n]->SetGravity(Gx, Gy, Gz);
}

inline void Mixture::Solve(double tIni, double tFin, double dt, double dtOut)
{
	double t    = tIni;      //
	double tout = t + dtOut; //
	_latts[0]->Homogenize();
	_latts[1]->Homogenize();
	WriteState (_T);
	while (t<tFin)
	{

		_latts[0]->ApplyForce ();
		_latts[1]->ApplyForce ();

		_latts[0]->ApplyGravity();
		_latts[1]->ApplyGravity();

		ApplyMixForce();
		SetMixVelocity();

		_latts[0]->Collide    ();
		_latts[1]->Collide    ();

		_latts[0]->BounceBack ();
		_latts[1]->BounceBack ();

		_latts[0]->Stream     ();
		_latts[1]->Stream     ();

		_latts[0]->ApplyBC    ();
		_latts[1]->ApplyBC    ();

		t += dt;
		if (t>=tout)
		{
			std::cout << "[1;34mMechSys[0m::LBM::Mixture::Solve: [1;31mT = " << _T << "[0m";
			std::cout << " Total mass:  l0 = " << _latts[0]->TotalMass();
			std::cout << "  l1 = " << _latts[1]->TotalMass() << "\n";
			_T++;
			WriteState (_T);
			tout += dtOut;
			//std::cout << GetCell(25,0) << std::endl;
		}
	}
}

inline void Mixture::WriteState(size_t TimeStep)
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
	double _dL = 1.0;
	oss << "SPACING "      << _dL << " " << _dL << " " << _dL << "\n";
	oss << "POINT_DATA "   << _size << "\n";

	// Solid cells
	oss << "SCALARS Geom float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<_size; ++i)
	{
		if (_latts[0]->GetCell(i)->IsSolid()) oss << "1.0\n";
		else                      oss << "0.0\n";
	}

	// Density field
	for (size_t n=0; n<_n_comp; n++)
	{
		oss << "SCALARS Density_" << n << " float 1\n";
		oss << "LOOKUP_TABLE default\n";
		for (size_t i=0; i<_size; ++i)
			oss <<_latts[n]->GetCell(i)->Density() << "\n";
	}

	// Velocity field
	for (size_t n=0; n<_n_comp; n++)
	{
		oss << "VECTORS Velocity_" << n << " float\n";
		for (size_t i=0; i<_size; ++i)
		{
			Vec3_t v; _latts[n]->GetCell(i)->Velocity(v);
			oss << v(0) << " " << v(1) << " " << v(2) << "\n";
		}
	}

	// Mass transfer field
	for (size_t n=0; n<_n_comp; n++)
	{
		oss << "VECTORS Mass_flux_" << n << " float\n";
		for (size_t i=0; i<_size; ++i)
		{
			Vec3_t v; _latts[n]->GetCell(i)->Velocity(v);
			v *= _latts[n]->GetCell(i)->Density();
			oss << v(0) << " " << v(1) << " " << v(2) << "\n";
		}
	}

	// Total Mass field
	oss << "VECTORS Total_mass_flux float\n";
	for (size_t i=0; i<_size; ++i)
	{
		Vec3_t v0; _latts[0]->GetCell(i)->Velocity(v0);
		Vec3_t v1; _latts[1]->GetCell(i)->Velocity(v1);
		v0 *= _latts[0]->GetCell(i)->Density();
		v1 *= _latts[1]->GetCell(i)->Density();
		Vec3_t v; v = v0+v1;
		oss << v(0) << " " << v(1) << " " << v(2) << "\n";
	}

	// Write to file and close file
	of << oss.str();
	of.close();
}

}; // namespace LBM

#endif // MECHSYS_LBM_MIXTURE_H
