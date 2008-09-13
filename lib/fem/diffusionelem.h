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

   Primary variale:    u
   Source variable:    f
   Secondary variable: q

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
	// Destructor
	virtual ~DiffusionElem() {}

	// Derived methods
	bool         IsReady         () const;
	bool         IsEssential     (char const * DOFName) const;
	void         SetModel        (char const * ModelName, char const * Prms, char const * Inis);
	Element    * Connect         (int iNodeLocal, FEM::Node * ptNode);
	void         UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void         BackupState     ();
	void         RestoreState    ();
	void         SetProperties   (Array<double> const & EleProps) { _unit_weight=EleProps[0]; }
	void         GetLabels       (Array<String> & Labels) const;
	void         Deactivate      ();
	char const * ModelName       () const { return (_a_model.Size()>0 ? _a_model[0]->Name() : "__no_model__"); }

	// Derived methods to assemble DAS matrices
	size_t nOrder0Matrices () const { return 1; }
	void   Order0MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void   Order0Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; ///< Permeability/Conductivity

	// Methods
	void B_Matrix (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

	// Access methods
	double Val (int iNodeLocal, char const * Name) const;
	double Val (                char const * Name) const;

private:
	// Data
	int                    _n_stress;
	Array<DiffusionModel*> _a_model;
	double                 _unit_weight;

	// Private methods
	void _calc_initial_internal_state ();

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
	else
		throw new Fatal("DiffusionElem::Val: This key==%s is invalid.",Name);
}

inline double DiffusionElem::Val(char const * Name) const
{
	// Get integration point values
	double sum = 0.0;
	for (size_t i_ip=0; i_ip<_n_int_pts; i_ip++)
		sum += _a_model[i_ip]->Val(Name);

	// Output single value at CG
	return sum/_n_int_pts;
}

inline void DiffusionElem::Deactivate()
{
	throw new Fatal("DiffusionElem::Deactivate: Feature not implemented in this element");
}

// Derived methods to assemble DAS matrices

inline void DiffusionElem::Order0MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
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

inline void DiffusionElem::Order0Matrix(size_t index, LinAlg::Matrix<double> & Ke) const
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
	for (size_t i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double t = _a_int_pts[i_ip].t;
		double w = _a_int_pts[i_ip].w;

		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix (derivs,J, B);   // Calculate B matrix for i_ip Integration Point

		// Conductivity
		_a_model[i_ip]->TgConductivity(D); 

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
