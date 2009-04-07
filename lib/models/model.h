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

typedef LinAlg::Vector<double>        Vec_t;
typedef LinAlg::Matrix<double>        Mat_t;
typedef char const                  * Str_t;
typedef std::map<String,double>       Prm_t;        ///< Parameters type
typedef const char                    PrmName_t[8]; ///< Parameters names. Ex: "E", "nu", "lam", "kap", ...
typedef std::map<String,double>       Ini_t;        ///< Initial values. Ex.: Sx=0.0
typedef blitz::TinyVector<double,3>   Vec3_t;
typedef blitz::TinyMatrix<double,3,3> Mat3_t;
typedef Array<double>                 IntVals;      ///< Internal values (specific volume, yield surface size, etc.)

/** Mechanical state. */
struct MechState
{
	Tensor2 Sig; ///< Stress
	Tensor2 Eps; ///< Strain
	IntVals Ivs; ///< Internal values
};

/** Diffusion state. */
struct DiffState
{
	Vec3_t  Vel; ///< Velocity
	Vec3_t  Gra; ///< Gradient
	IntVals Ivs; ///< Internal values
};

class Model
{
public:
	// Typedefs

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
	virtual double    * DefPrms() const =0; ///< Default parameters values
	virtual Str_t       Name  () const =0; ///< Model name

	// Viscosity
	virtual bool HasDSigStar  () const { return false; } ///< Has delta sigma star (to be subtract to dsig)
	virtual void CalcDSigStar (double t, double Dt, MechState const & State, Vec_t & DSigStar) const {}

	// Methods that may be derived
	virtual void SetPyName (Str_t ScriptFileName) {} ///< Set script file name for models that can be extended using Python

	/* Initialize internal values. */
	virtual void InitIVS (Ini_t const & Ini, MechState & State) const {} ///< In/Out: MechState
	virtual void InitIVS (Ini_t const & Ini, DiffState & State) const {} ///< In/Out: DiffState

	/* Tangent stiffness, conductivity and permeability. */
	virtual void TgStiffness    (MechState const & State, Mat_t & Dmat) const {}
	virtual void TgConductivity (DiffState const & State, Mat_t & Dmat) const {}
	virtual void TgPermeability (                         Mat_t & Kmat) const {}

	/* State update. */
	virtual int StateUpdate (double t, double Dt, Vec_t const & DEps, MechState & State, Vec_t & DSig) { return -1; }
	virtual int StateUpdate (double t, double Dt, Vec_t const & DGra, DiffState & State, Vec_t & DVel) { return -1; }

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
	lp.ReadSomeVariables (NPrms(), Prms(), DefPrms(), _prms, "Parameter", "Model", _tag);

	// Initialize model
	_initialize ();
}

inline double Model::Prm(Str_t Key) const
{
	Prm_t::const_iterator it = _prms.find(Key);
	if (it==_prms.end()) throw new Fatal("Model::Prm: Could not find parameter < %s > in _prms array",Key);
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
	// Check if ModelName has a comma (separating the name of a script)
	// Ex.: "PyEquilib,mymodel.py"
	Array<String> res;
	LineParser lp(Name);
	lp.SplitLine (",", res);
	if (res.Size()==0) throw new Fatal("AllocModel: Model name cannot be empty. Use, for instance: LinElastic");
	
	// Check if there is Name model implemented
	ModelMakerPtr ptr=NULL;
	ptr = ModelFactory[res[0].CStr()];
	if (ptr==NULL) throw new Fatal(_("AllocModel: There is no < %s > implemented in this library"), Name);
	Model * mdl = (*ptr)();

	// Set script name
	if (res.Size()==2) mdl->SetPyName (res[1].CStr());
	
	return mdl;
}


#endif // MECHSYS_MODEL_H
