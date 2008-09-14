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


/* __ Element for diffusion transport simulations __

   Primary variale:    u == temperature/head                    ~~ displacements
   Source variable:    f == heat source/wate recharge (pumping) ~~ body forces
   Secondary variable: q == heat flow/seepage                   ~~ tractions

   Heat transfer:
      u = temperature
      f = heat source
	  q = heat flux due to conduction

   Groundwater flow:
      u = piezometric head
      f = recharge (pumping => -f)
	  q = water flux (seepage)

*/


#ifndef MECHSYS_FEM_DIFFUSION_H
#define MECHSYS_FEM_DIFFUSION_H

// MechSys
#include "fem/element.h"
#include "models/diffusionmodel.h"
#include "util/string.h"
#include "linalg/laexpr.h"

namespace FEM
{

class DiffusionElem : public virtual Element
{
public:
	// Constructor
	DiffusionElem () : _source(0.0), _has_source(false) {}

	// Destructor
	virtual ~DiffusionElem() {}

	// Derived methods
	bool         IsReady         () const;
	bool         IsEssential     (char const * DOFName) const;
	void         SetModel        (char const * ModelName, char const * Prms, char const * Inis);
	void         SetProps        (Array<double> const & P);
	Element    * Connect         (int iNodeLocal, FEM::Node * ptNode);
	void         UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	bool         AddVolForces    (LinAlg::Vector<double> & dFext) const;
	void         BackupState     ();
	void         RestoreState    ();
	void         GetLabels       (Array<String> & Labels) const;
	void         Deactivate      ();
	char const * ModelName       () const { return (_a_model.Size()>0 ? _a_model[0]->Name() : "__no_model__"); }

	// Derived methods to assemble DAS matrices
	size_t nOrder1Matrices () const { return 1; }
	void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; ///< Permeability/Conductivity

	// Methods
	void B_Matrix (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

	// Access methods
	double Val (int iNodeLocal, char const * Name) const;
	double Val (                char const * Name) const;

private:
	// Data
	Array<DiffusionModel*> _a_model;    ///< Array of pointers to diffusion models
	double                 _source;     ///< Source (heat) or recharge/pumping (water)
	bool                   _has_source; ///< Has source?

	// Private methods
	void _calc_initial_internal_state (); ///< Calculate initial internal state

	// Private methods that MUST be derived
	virtual int _geom() const =0; ///< Geometry of the element: 1:1D, 2:2D, 3:3D

}; // class DiffusionElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
	
// Derived methods

inline bool DiffusionElem::IsReady() const
{
	if (_a_model.Size()==static_cast<size_t>(_n_int_pts) && _connects.Size()==static_cast<size_t>(_n_nodes)) return true;
	else return false;
}

inline bool DiffusionElem::IsEssential(char const * DOFName) const
{
	if (strcmp(DOFName,"u")==0) return true;
	return false;
}

inline void DiffusionElem::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("DiffusionElem::SetModel: The space dimension (SetDim) must be set before calling this method");

	// If pointers to model was not already defined => No model was allocated
	if (_a_model.Size()==0)
	{
		// Resize the array of model pointers
		_a_model.Resize(_n_int_pts);

		// Loop along integration points
		for (size_t i=0; i<_n_int_pts; ++i)
		{
			// Allocate a new model and set parameters
			_a_model[i] = static_cast<DiffusionModel*>(AllocModel(ModelName));
			_a_model[i]->SetGeom (_geom());
			_a_model[i]->SetPrms (Prms);
			_a_model[i]->SetInis (Inis);
		}

		// Calculate initial internal forces
		_calc_initial_internal_state();
	}
	else throw new Fatal("DiffusionElem::SetModel: Feature not implemented.");
}

inline void DiffusionElem::SetProps(Array<double> const & P)
{
	if (P.Size()==1)
	{
		_source     = P[0];
		_has_source = true;
	}
	else throw new Fatal("DiffusionElem::SetProps: ElemProp(P)==Source must have size equal to 1 (%d is invalid)",P.Size());
}

inline Element * DiffusionElem::Connect(int iNodeLocal, FEM::Node * ptNode)
{
	// Connects
	_connects[iNodeLocal] = ptNode;

	// Add Degree of Freedom to a node (Essential, Natural)
	_connects[iNodeLocal]->AddDOF("u", "q");

	// Shared
	_connects[iNodeLocal]->SetSharedBy(_my_id);

	return this;
}

inline void DiffusionElem::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) temperature/head vector
	LinAlg::Vector<double> du(_n_nodes); // Delta temperature/head of this element
	
	// Assemble (local/element) temperature/head vector
	for (size_t i=0; i<_n_nodes; ++i)
		du(i) = dUglobal(_connects[i]->DOFVar("u").EqID);
	
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> dF(_n_nodes); // Delta internal force of this element
	dF.SetValues(0.0);
	
	// Allocate entities used for every integration point
	LinAlg::Vector<double> shape;   // shape functions
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> du_dx;   // gradient
	LinAlg::Vector<double> vel;     // velocity

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		Shape    (r,s,t, shape);   // Calculate shape functions
		Derivs   (r,s,t, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs, J, B);   // Calculate B matrix for i Integration Point

		// Calculate gradient
		du_dx = B*du;
		
		// Update model
		_a_model[i]->StateUpdate(du_dx, vel);

		// Calculate internal force vector
		//dF += trn(B)*vel*det(J)*w;
	}

	// Return internal forces
	for (size_t i=0; i<_n_nodes; ++i)
		dFint(_connects[i]->DOFVar("q").EqID) += dF(i);
}

inline bool DiffusionElem::AddVolForces(LinAlg::Vector<double> & FVol) const
{
	if (_has_source)
	{
		// Allocate (local/element) external volume force vector
		LinAlg::Vector<double> fvol(_n_nodes);
		fvol.SetValues(0.0);

		// Allocate entities used for every integration point
		LinAlg::Vector<double> shape;
		LinAlg::Matrix<double> derivs;
		LinAlg::Matrix<double> J;

		// Loop along integration points
		for (size_t i=0; i<_n_int_pts; ++i)
		{
			// Temporary Integration Points
			double r = _a_int_pts[i].r;
			double s = _a_int_pts[i].s;
			double t = _a_int_pts[i].t;
			double w = _a_int_pts[i].w;

			Shape    (r,s,t, shape);   // Calculate shape functions for i IP
			Derivs   (r,s,t, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
			Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i Integration Point

			// Calculate external volume force
			for (size_t j=0; j<_n_nodes; j++)
				fvol(j) += _source*shape(j)*det(J)*w;
		}

		// Add to external force vector
		for (size_t i=0; i<_n_nodes; ++i)
			FVol(_connects[i]->DOFVar("q").EqID) += fvol(i);

		// Flag that there are volumetric forces
		return true;
	}
	else return false; // there aren't volumetric forces
}

inline void DiffusionElem::BackupState()
{
	for (size_t i=0; i<_n_int_pts; ++i)
		_a_model[i]->BackupState();
}

inline void DiffusionElem::RestoreState()
{
	for (size_t i=0; i<_n_int_pts; ++i)
		_a_model[i]->RestoreState();
}

inline void DiffusionElem::GetLabels(Array<String> & Labels) const
{
	// Get labels of all values to output
	Labels.Resize(2);
	Labels[0]="u";
	Labels[1]="q";
}

inline double DiffusionElem::Val(int iNodeLocal, char const * Name) const
{
	// Essential
	if (strcmp(Name,"u")==0)
		return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Natural
	else if (strcmp(Name,"q")==0)
		return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	// Velocities, internal values, etc.
	else
	{
		// Vectors for extrapolation
		LinAlg::Vector<double>    ip_values (_n_int_pts);
		LinAlg::Vector<double> nodal_values (_n_nodes);

		// Get integration point values
		if (_a_model.Size()==_n_int_pts)
			for (size_t i=0; i<_n_int_pts; i++)
				ip_values(i) = _a_model[i]->Val(Name);
		else throw new Fatal("DiffusionElem::Val: Constitutive models for this element (ID==%d) were not set yet", _my_id);

		// Extrapolation
		Extrapolate (ip_values, nodal_values);

		// Output single value
		return nodal_values (iNodeLocal);
	}
}

inline double DiffusionElem::Val(char const * Name) const
{
	// Get integration point values
	double sum = 0.0;
	for (size_t i=0; i<_n_int_pts; i++)
		sum += _a_model[i]->Val(Name);

	// Output single value at CG
	return sum/_n_int_pts;
}

inline void DiffusionElem::Deactivate()
{
	throw new Fatal("DiffusionElem::Deactivate: Feature not implemented in this element");
}

// Derived methods to assemble DAS matrices

inline void DiffusionElem::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	// Size of Ke
	int n_rows = _n_nodes; // == n_cols

	// Mounting a map of positions from Ke to Global
	int idx_Ke = 0;                // position (idx) inside Ke matrix
	RowsMap       .Resize(n_rows); // size=Ke.Rows()=Ke.Cols()
	RowsEssenPresc.Resize(n_rows); // size=Ke.Rows()=Ke.Cols()

	// Fill map of Ke position to K position of DOFs components
	for (size_t i_node=0; i_node<_n_nodes; ++i_node)
	{
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("u").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("u").IsEssenPresc; 
		idx_Ke++;
	}
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
}

inline void DiffusionElem::Order1Matrix(size_t index, LinAlg::Matrix<double> & Ke) const
{
	/* Conductivity:
	   ============
	
	                 /    T
	        [Ke]  =  | [B]  * [D] * [B]  * dV
	                 /
	*/

	// Resize Ke
	Ke.Resize(_n_nodes, _n_nodes); // sum(Bt*D*B*det(J)*w)
	Ke.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs; // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;      // Jacobian matrix
	LinAlg::Matrix<double> B;      // B matrix
	LinAlg::Matrix<double> D;      // Conductivity matrix

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs,J, B);   // Calculate B matrix for i Integration Point

		// Conductivity
		_a_model[i]->TgConductivity(D); 

		// Calculate Tangent Conductivity
		Ke += trn(B)*D*B*det(J)*w;
	}
}
	
// Methods

inline void DiffusionElem::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	// Calculate Bp matrix (NDIM x _n_nodes)
	B = inv(J)*derivs; // equal to the cartesian derivatives matrix for diffusion elements
}


/* private */

inline void DiffusionElem::_calc_initial_internal_state()
{
}


}; // namespace FEM

#endif // MECHSYS_FEM_DIFFUSIO_H
