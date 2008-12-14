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

#ifndef MECHSYS_BIOTELASTIC_H
#define MECHSYS_BIOTELASTIC_H

// MechSys
#include "models/equilibmodel.h"
#include "tensors/tensors.h"
#include "util/string.h"
#include "util/util.h"
#include "util/lineparser.h"

typedef char const * Str_t;

class BiotElastic : public EquilibModel
{
public:
	// Constants
	static const char BIOTELASTIC_PN[3][8];

	// Destructor
	virtual ~BiotElastic () {}

	// Methods
	void TgPermeability (Mat_t & Kmat) const { Kmat = _K; }

	// Derived methods
	Str_t Name () const { return "BiotElastic"; }

private:
	// Data
	Tensor4 _De; ///< Constant tangent stiffness
	Mat_t   _K;  ///< Constant tangent permeability

	// Private methods
	void _set_ctes   () { _np=3;  PRMS=BIOTELASTIC_PN; }
	void _initialize ();
	void _stiff      (Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const;

}; // class BiotElastic

const char BiotElastic::BIOTELASTIC_PN[3][8] = {"E", "nu", "k"};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void BiotElastic::_initialize()
{
	// 3D:gi=0, 2D:gi=1
	if (_gi<0 || _gi>1) throw new Fatal("BiotElastic::_initialize: Tag=%d. This model is only available for 3D(gi==0) and 2D(plane-strain, gi==1) problems. (gi==%d is invalid)",_tag,_gi);

	// Parameters
	double E  = _prms[0];
	double nu = _prms[1];
	double k  = _prms[2];

	// Check
	if (E<=0.0)             throw new Fatal("BiotElastic::_initialize: Tag=%d. Young modulus (E) must positive. E==%f is invalid",_tag,E);
	if (nu<0.0 || nu>0.499) throw new Fatal("BiotElastic::_initialize: Tag=%d. Poisson ratio (nu) must be in the range: 0 < nu < 0.5). nu==%f is invalid",_tag,nu);
	if (k<0.0)              throw new Fatal("BiotElastic::_initialize: Tag=%d. Isotropic permeability must positive. k=%f is invalid",_tag,k);

	// Set stiffness
	double c  = E/((1.0+nu)*(1.0-2.0*nu));
	double c1 = c*(1.0-nu);
	double c2 = c*(1.0-2.0*nu)/2.0;
	double c3 = c*nu;
	_De = c1     , c3     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	      c3     , c1     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	      c3     , c3     , c1     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, c2 *2.0, 0.0*2.0, 0.0*2.0,
	      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, c2 *2.0, 0.0*2.0,
	      0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, 0.0*2.0, c2 *2.0; // In Mandel's basis

	// Set permeability
	double kx = k;
	double ky = k;
	double kz = k;
	if (_gi==1) // 2D
	{
		_K.Resize(2,2);
		_K =  kx,  0.0,
		     0.0,   ky;
	}
	else
	{
		_K.Resize(3,3);
		_K =  kx,  0.0,  0.0,
		     0.0,   ky,  0.0,
		     0.0,  0.0,   kz;
	}
}

inline void BiotElastic::_stiff(Tensor2 const & DEps, Tensor2 const & Sig, Tensor2 const & Eps, IntVals const & Ivs,  Tensor4 & D, Array<Tensor2> & B) const
{
	D = _De;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new model
Model * BiotElasticMaker() { return new BiotElastic(); }

// Register model
int BiotElasticRegister() { ModelFactory["BiotElastic"]=BiotElasticMaker;  return 0; }

// Call register
int __BiotElastic_dummy_int = BiotElasticRegister();


#endif // MECHSYS_BIOTELASTIC_H
