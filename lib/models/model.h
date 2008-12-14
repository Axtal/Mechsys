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

#ifndef MECHSYS_MODEL_H
#define MECHSYS_MODEL_H

// MechSys
#include "util/string.h"
#include "util/array.h"
#include "util/util.h"
#include "util/lineparser.h"
#include "tensors/tensors.h"

using Tensors::Tensor2;
using Tensors::Tensor4;

typedef LinAlg::Vector<double>  Vec_t;
typedef LinAlg::Matrix<double>  Mat_t;
typedef char const            * Str_t;
typedef std::map<String,double> Prm_t;        ///< Parameters type
typedef const char              PrmName_t[8]; ///< Parameters names. Ex: "E", "nu", "lam", "kap", ...
typedef std::map<String,double> Ini_t;        ///< Initial values. Ex.: Sx=0.0

class Model
{
public:
	// Typedefs
	typedef Array<double> IntVals;      ///< Internal values (specific volume, yield surface size, etc.)

	// Constructor
	Model () : _gi(-1), _tag(0) { STOL().dTini().mMin().mMax().maxSS(); }

	// Destructor
	virtual ~Model () {}

	// Methods
	void   SetGeomIdx (int GeomIdx);            ///< Set geometry index. MUST be called before Initialize
	void   Initialize (int Tag, Str_t StrPrms); ///< Initialize this model
	double Prm        (Str_t Key) const;        ///< Value of a parameter

	// Methods to be derived
	virtual int         NPrms () const =0; ///< Number of parameters
	virtual PrmName_t * Prms  () const =0; ///< Parameters names. Ex: "E", "nu", "lam", "kap", ...
	virtual Str_t       Name  () const =0; ///< Model name

	/* Initialize internal values. */
	virtual void InitIVS (Ini_t const & Ini, Tensor2 const & Sig, Tensor2 const & Eps, IntVals & Ivs) const {}

	/* Tangent stiffness. */
	virtual void TgStiffness (Tensor2 const & Sig,
	                          Tensor2 const & Eps,
	                          IntVals const & Ivs,
	                          Mat_t         & Dmat) const {}

	/* Tangent stiffness. */
	virtual void TgPermeability (Mat_t & Kmat) const {}

	/* State update. */
	virtual int StateUpdate (Vec_t   const & DEps,
	                         Tensor2       & Sig,
	                         Tensor2       & Eps,
	                         IntVals       & Ivs,
	                         Vec_t         & DSig) { return -1; }

	// Integration constants
	Model & STOL  (double Val=1.0e-5) { _STOL =Val; return (*this); }
	Model & dTini (double Val=1.0   ) { _dTini=Val; return (*this); }
	Model & mMin  (double Val=0.1   ) { _mMin =Val; return (*this); }
	Model & mMax  (double Val=10.0  ) { _mMax =Val; return (*this); }
	Model & maxSS (size_t Val=2000  ) { _maxSS=Val; return (*this); }

protected:
	// Data
	int   _gi;   ///< Geometry index: Equilib: 3D=0, PStrain=1, PStress=2, Axis=3. Others: 3D=0, 2D=1
	int   _tag;  ///< The tag of the model
	Prm_t _prms; ///< Parameters

	// Constants for the stress update algorithm
	double _STOL;
	double _dTini;
	double _mMin;
	double _mMax;
	size_t _maxSS;

	// Methods
	virtual void _initialize () =0; ///< Initialize the model

}; // class Model


////////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void Model::SetGeomIdx(int GeomIdx)
{
	if (_gi<0) _gi = GeomIdx; // first time set
	else if (GeomIdx!=_gi) throw new Fatal("Model::SetGeomIdx: GeomIdx==%d was already set. All other calls to this method must be with the same GeomIdx (%d is invalid).\nTo use different models with different GeomIdx, different tags should be provided in eatts array.",_gi,GeomIdx);
	if (_gi<0) throw new Fatal("Model::SetGeomIdx: Geometry index==%d is invalid\n It must be in: 0==3D, 1==2D(plane-strain), 2==2D(plane-stress), 3==2D(axis-symmetric)",_gi);
}

inline void Model::Initialize(int Tag, Str_t StrPrms)
{
	// Set PRMS and resize _prms
	_tag = Tag;

	// Read parameters
	LineParser lp(StrPrms);
	lp.ReadVariables (NPrms(), Prms(), _prms, "parameters", "Model", _tag);

	// Initialize model
	_initialize ();
}

inline double Model::Prm(Str_t Key) const
{
	Prm_t::const_iterator it = _prms.find(Key);
	if (it==_prms.end()) throw new Fatal("Model::Prm: Could not find parameter < %d > in _prms array",Key);
	return it->second;
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new Model
typedef Model * (*ModelMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes Models
typedef std::map<String, ModelMakerPtr, std::less<String> > ModelFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes Models
ModelFactory_t ModelFactory;

// Allocate a new Models according to a string giving the name of the Model
Model * AllocModel(Str_t Name)
{
	// Check if there is Name model implemented
	ModelMakerPtr ptr=NULL;
	ptr = ModelFactory[Name];
	if (ptr==NULL) throw new Fatal(_("FEM::AllocModel: There is no < %s > implemented in this library"), Name);
	return (*ptr)();
}


#endif // MECHSYS_MODEL_H
