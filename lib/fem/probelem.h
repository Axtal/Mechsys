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
//#include "fem/probelem.h"
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
	ProbElem() :  IsActive(true) {}

	// Methods
	void  Initialize (GeomElem * GE, int NDim, bool IsActive);                       ///< Initialize the element
	virtual void      AddVolForces ()                                           { }  ///< Method to apply volumetric (body) forces as boundary condition
	virtual void      ClearDisp    ()                                           { }  ///< Clear displacements and strains (for equilibrium/coupled problems)
	virtual void      SetActive    (bool Activate)=0;                                ///< Activate element (construction/excavation)
	virtual bool      CheckMdl     () const=0;                                       ///< Check if constitutive models are OK
	virtual void      EdgeBry      (Str_t Key, double Val, int iEdge)           { }  ///< Set edge boundary conditions (Initialize MUST be called first)
	virtual void      EdgeBry      (Str_t Key, double V0, double V1, int iEdge) { }  ///< Set edge boundary conditions (Initialize MUST be called first)
	virtual void      FaceBry      (Str_t Key, double Val, int iFace)           { }  ///< Set face boundary conditions (Initialize MUST be called first)
	virtual void      CalcDeps     () const                                     { }  ///< Calculate dependent variables (to be called before Val() or OutNodes() for example). Necessary for output of principal stresses, for example.
	virtual Str_t     ModelName    () const =0;                                      ///< Return the name of the model of the first IP of this element
	virtual double    Val          (int iNod, Str_t Key) const =0;                   ///< Return computed values at the Nodes of the element. Ex.: Key="ux", "fx", "Sx", "Sxy", "Ex", etc.
	virtual double    Val          (          Str_t Key) const =0;                   ///< Return computed values at the CG of the element. Ex.: Key="Sx", "Sxy", "Ex", etc.
	virtual bool      IsEssen      (Str_t Key) const =0;                             ///< Is the correspondent DOFKey (Degree of Freedom, such as "Dux") essential (such displacements)?
	virtual void      SetProps     (Str_t Properties) =0;                            ///< Set element properties such as body forces, internal heat source, water pumping, etc.
	virtual void      SetModel     (Str_t ModelName, Str_t Prms, Str_t Inis) =0;     ///< (Re)allocate model with parameters and initial values
	virtual void      SetConn      (int iNod, FEM::Node * ptNode) =0;                ///< Set connectivity, by linking the local node ID with the pointer to the connection node
	virtual void      Update       (double h, Vec_t const & dU, Vec_t & dFint) =0;   ///< Update the internal state of this element for given dU and update the DOFs related to this element inside dFint (internal forces increment vector)
	virtual void      Backup       () =0;                                            ///< Backup internal state
	virtual void      Restore      () =0;                                            ///< Restore internal state from a previously backup state
	virtual void      GetLbls      (Array<String> & Lbls) const =0;                  ///< Get the labels of all values to be output
	virtual void      OutNodes     (Mat_t & Vals, Array<String> & Lbls) const =0;    ///< Output values at nodes
	virtual void      OutInfo      (std::ostream & os) const =0;                     ///< Output extra info of the derived element
	virtual bool      HasExtra     () const =0                                       ///< Has extra output ?
	virtual void      OutExtra     (Mat_t & Coords, Vec_t & Norm,                
	virtual                         Mat_t & Vals, Array<String> & Lbls) const =0;    ///< Output extra information
	virtual size_t    NCMats       () const                                     {return 0; } ///< Number of C matrices such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix
	virtual size_t    NHMats       () const                                     {return 0; } ///< Number of H matrices such as H:Permeability
	virtual size_t    NUVecs       () const                                     {return 0; } ///< Number of U vectors such as U:Displacements, P:Pore-pressures
	virtual void      CMatrix      (size_t Idx, Mat_t & M) const                { }  ///< C matrix such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix
	virtual void      HMatrix      (size_t Idx, Mat_t & M) const                { }  ///< H matrix such as H:Permeability
	virtual void      UVector      (size_t Idx, Vec_t & V) const                { }  ///< U vector such as U:Displacement, P:Pore-pressure
	virtual void      CMatMap      (size_t Idx,
	                                Array<size_t> & RMap,
	                                Array<size_t> & CMap,
	                                Array<bool> & RUPresc,
	                                Array<bool> & CUPresc) const                { } ///< CMatrix map to convert local DOFs into global equation positions
	virtual void     HMatMap       (size_t Idx,
	                                Array<size_t> & RMap,
	                                Array<size_t> & CMap,
	                                Array<bool> & RUPresc,
	                                Array<bool> & CUPresc) const                { } ///< HMatrix map to convert local DOFs into global equation positions
	virtual void     UVecMap       (size_t Idx, Array<size_t> & RMap) const     { } ///< UVector map to convert local DOFs into global equation positions


	bool             IsActive;
	
private:
	// Data
	GeomElem * _ge;   ///< Geometry element

}; // class ProbElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void  Initialize (GeomElem * GE, int NDimVal, bool IsActiveVal)
{
	_ge      = GE;
	NDim     = NDimVal;
	IsActive = IsActiveVal;
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new ProbElement
typedef ProbElem * (*ProbElemMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes elements
typedef std::map<String, ProbElemMakerPtr, std::less<String> > ProbElemFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes elements
ProbElemFactory_t ProbElemFactory;

// Allocate a new element according to a string giving the name of the element
ProbElem * AllocProbElem(char const * ProbElemName)
{
	ProbElemMakerPtr ptr=NULL;
	ptr = ProbElemFactory[ProbElemName];
	if (ptr==NULL)
		throw new Fatal(_("FEM::AllocProbElem: There is no < %s > implemented in this library"), ProbElemName);

	return (*ptr)();
}

}; // namespace FEM

#endif // MECHSYS_FEM_PROBELEM
