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

#ifndef MECHSYS_FEM_PROBELEM_H
#define MECHSYS_FEM_PROBELEM_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "util/string.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "linalg/laexpr.h"

namespace FEM
{

typedef LinAlg::Matrix<double> Mat_t;
typedef LinAlg::Vector<double> Vec_t;
typedef char const *           Str_t;

class ProbElem
{
public:

	// Constructor
	ProbElem() : IsActive(true) {}

	// Methods related to PROBLEM
	virtual void      AddVolForces ()                                           {}
	virtual void      ClearDisp    ()                                           {}
	virtual void      SetActive    (bool Activate)                             =0;
	virtual bool      CheckMdl     () const                                    =0;
	virtual void      EdgeBry      (Str_t Key, double Val, int iEdge)           {}
	virtual void      EdgeBry      (Str_t Key, double V0, double V1, int iEdge) {}
	virtual void      FaceBry      (Str_t Key, double Val, int iFace)           {}
	virtual void      CalcDeps     () const                                     {}
	virtual Str_t     ModelName    () const                                    =0;
	virtual double    Val          (int iNod, Str_t Key) const                 =0;
	virtual double    Val          (          Str_t Key) const                 =0;
	virtual bool      IsEssen      (Str_t Key) const                           =0;
	virtual void      SetProps     (Str_t Properties)                          =0;
	virtual void      SetModel     (Str_t ModelName, Str_t Prms, Str_t Inis)   =0;
	virtual void      SetConn      (int iNod, FEM::Node * ptNode)              =0;
	virtual void      Update       (double h, Vec_t const & dU, Vec_t & dFint) =0;
	virtual void      Backup       ()                                          =0;
	virtual void      Restore      ()                                          =0;
	virtual void      GetLbls      (Array<String> & Lbls) const                =0;
	virtual void      OutNodes     (Mat_t & Vals, Array<String> & Lbls) const  =0;
	virtual void      OutInfo      (std::ostream & os) const                   =0;
	virtual bool      HasExtra     () const                      { return false; }
	virtual void      OutExtra     (Mat_t & Coords, Vec_t & Norm,
	virtual                         Mat_t & Vals, Array<String> & Lbls) const   {}
	virtual size_t    NCMats       () const                          { return 0; }
	virtual size_t    NHMats       () const                          { return 0; }
	virtual size_t    NUVecs       () const                          { return 0; }
	virtual void      CMatrix      (size_t Idx, Mat_t & M) const                {}
	virtual void      HMatrix      (size_t Idx, Mat_t & M) const                {}
	virtual void      UVector      (size_t Idx, Vec_t & V) const                {}
	virtual void      CMatMap      (size_t Idx,
	                                Array<size_t> & RMap,
	                                Array<size_t> & CMap,
	                                Array<bool> & RUPresc,
	                                Array<bool> & CUPresc) const                {}
	virtual void      HMatMap      (size_t Idx,
	                                Array<size_t> & RMap,
	                                Array<size_t> & CMap,
	                                Array<bool> & RUPresc,
	                                Array<bool> & CUPresc) const                {}
	virtual void      UVecMap      (size_t Idx, Array<size_t> & RMap) const     {}

	// Methods
	void Initialize (GeomElem * GE, bool IsAct); ///< Initialize the element

	// Public data (read only)
	bool IsActive;
	
protected:
	// Data
	GeomElem * _ge; ///< Geometry element

}; // class ProbElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void ProbElem::Initialize(GeomElem * GE, bool IsAct)
{
	_ge      = GE;
	IsActive = IsAct;
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new ProbElement
typedef ProbElem * (*ProbElemMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes elements
typedef std::map<String, ProbElemMakerPtr, std::less<String> > ProbElemFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes elements
ProbElemFactory_t ProbElemFactory;

// Allocate a new ProbElem according to a string giving the name of the element
ProbElem * AllocProbElem(char const * ProbElemName)
{
	ProbElemMakerPtr ptr=NULL;
	ptr = ProbElemFactory[ProbElemName];
	if (ptr==NULL) throw new Fatal("FEM::AllocProbElem: There is no < %s > implemented in this library", ProbElemName);
	return (*ptr)();
}

}; // namespace FEM

#endif // MECHSYS_FEM_PROBELEM
