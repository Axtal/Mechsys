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

#ifndef MECHSYS_FEM_BEAM_H
#define MECHSYS_FEM_BEAM_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/lin2.h"
#include "util/exception.h"

namespace FEM
{

class Beam : public Lin2, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name() const { return NAME; };

	// Methods
	void Order1Matrix__ (size_t Index, LinAlg::Matrix<double> & Ke) const; ///< Stiffness
	void B_Matrix__     (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

private:
	// Private methods
	int  _geom () const { return 1; }    ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	bool _beam () const { return true; } ///< This is a beam element

	// Private methods
	void _set_ndim                    (int nDim); ///< Set space dimension
	void _calc_initial_internal_state__ ();         ///< Calculate initial internal state

}; // class Beam

// Beam constants
char const * Beam::NAME = "Beam";


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void Beam::Order1Matrix__(size_t Index, LinAlg::Matrix<double> & Ke) const
{
	if (_ndim==2)
	{
		double E  = 100;
		double A  = 125;
		double I  = 250;
		double dx = _connects[1]->X()-_connects[0]->X();
		double dy = _connects[1]->Y()-_connects[0]->Y();
		double LL = dx*dx+dy*dy;
		double L  = sqrt(LL);
		double c  = dx/L;
		double s  = dy/L;
		double c1 = E*(A*c*c+12.0*I*s*s/LL)/L;
		double c2 = E*((A-12.0*I/LL)*c*s)/L;
		double c3 = E*(6.0*I*s/L)/L;
		double c4 = E*(A*s*s+12.0*I*c*c/LL)/L;
		double c5 = E*(6.0*I*c/L)/L;
		double c6 = E*(4.0*I)/L;
		double c7 = E*(2.0*I)/L;
		Ke.Resize(_nd*_n_nodes, _nd*_n_nodes);
		Ke =  c1,  c2, -c3, -c1, -c2, -c3,
		      c2,  c4,  c5, -c2, -c4,  c5,
		     -c3,  c5,  c6,  c3, -c5,  c7,
		     -c1, -c2,  c3,  c1,  c2,  c3,
		     -c2, -c4, -c5,  c2,  c4, -c5,
		     -c3,  c5,  c7,  c3, -c5,  c6;
	}
	else throw new Fatal("Beam::Order1Matrix: Feature no available for nDim==%d",_ndim);
}

inline void Beam::B_Matrix__(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	throw new Fatal("Beam::B_Matrix: Feature not available");
}


/* private */

inline void Beam::_set_ndim(int nDim)
{
	if (nDim<1) throw new Fatal("Beam::_set_ndim: For this element, nDim must be greater than or equal to 1 (%d is invalid)",nDim);
	_ndim = nDim;
	_d    = _ndim-1;
	_nd   = EquilibElem::NDB[_d];
}

inline void Beam::_calc_initial_internal_state__()
{
	//throw new Fatal("Beam::_calc_initial_internal_state: Feature not available");
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Beam element
Element * BeamMaker()
{
	return new Beam();
}

// Register a Beam element into ElementFactory array map
int BeamRegister()
{
	ElementFactory[Beam::NAME] = BeamMaker;
	return 0;
}

// Execute the autoregistration
int __Beam_dummy_int  = BeamRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_BEAM_H
