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


/* __ Equilibrium capable element __

T = Nodal temperature [ temp ]     ex. C°
F = Nodal heat flow   [ work/time] ex. Joules/sec

*/

#ifndef MECHSYS_FEM_HEAT_H
#define MECHSYS_FEM_HEAT_H

// MechSys
#include "fem/element.h"
#include "models/heatmodel.h"
#include "util/string.h"
#include "linalg/laexpr.h"

namespace FEM
{

class HeatElem : public virtual Element
{
public:
	// Destructor
	virtual ~HeatElem() {}

	// Derived methods
	bool      IsReady         () const;
	bool      IsEssential     (char const * DOFName) const;
	void      SetModel        (char const * ModelName, char const * Prms, char const * Inis);
	Element * SetNode         (int iNodeLocal, int iNodeGlobal);
	void      UpdateState     (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void      BackupState     ();
	void      RestoreState    ();
	void      SetGeometryType (int Geom) {};  
	void      SetProperties   (Array<double> const & EleProps) { _unit_weight=EleProps[0]; }
	void      GetLabels       (Array<String> & Labels) const;
	void      Deactivate      ();

	// Derived methods to assemble DAS matrices
	size_t nOrder1Matrices () const { return 1; }
	void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & Ke) const; // Stiffness

	// Methods
	void B_Matrix (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;

	// Access methods
	double Val (int iNodeLocal, char const * Name) const;
	double Val (                char const * Name) const;

private:
	// Data
	int                  _n_stress;
	Array<HeatModel*> _a_model;
	double               _unit_weight;

	// Private methods
	void _calc_initial_internal_forces ();

}; // class HeatElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */
	
// Derived methods

inline bool HeatElem::IsReady() const
{
	if (_a_model.Size()==static_cast<size_t>(_n_int_pts) && _connects.Size()==static_cast<size_t>(_n_nodes)) return true;
	else return false;
}

inline bool HeatElem::IsEssential(char const * DOFName) const
{
	if (strcmp(DOFName,"T")==0) return true; // Temperature is the essential variable
	return false;
}

inline void HeatElem::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// If pointers to model was not already defined => No model was allocated
	if (_a_model.Size()==0)
	{
		// Resize the array of model pointers
		_a_model.Resize(_n_int_pts);

		// Loop along integration points
		for (size_t i=0; i<_n_int_pts; ++i)
		{
			// Allocate a new model and set parameters
			_a_model[i] = static_cast<HeatModel*>(AllocModel(ModelName));
			_a_model[i]->SetPrms(Prms);
			_a_model[i]->SetInis(Inis);
		}

		// Calculate initial internal forces
		//_calc_initial_internal_forces(); TODO
	}
	else throw new Fatal("HeatElem::SetModel: Feature not implemented.");
}

inline Element * HeatElem::SetNode(int iNodeLocal, int iNodeGlobal)
{
	// Connects
	_connects[iNodeLocal] = Nodes[iNodeGlobal];

	// Add Degree of Freedom to a node (Essential, Natural)
	Nodes[iNodeGlobal]->AddDOF("T", "F"); // T stands for temperature and F for nodal heat source

	// Shared
	Nodes[iNodeGlobal]->SetSharedBy(_my_id);

	return this;
}

inline void HeatElem::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> dT(_n_nodes); // Delta temp. of this element
	
	// Assemble (local/element) temperatures vector
	for (size_t i=0; i<_n_nodes; ++i)
	   dT(i) = dUglobal(_connects[i]->DOFVar("T").EqID);
	
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> dF(_n_nodes); // Delta internal heat source of this element
	dF.SetValues(0.0);
	
	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> DGrad;   // Gradient vector 
	LinAlg::Vector<double> DFlow;   // Stress vector 

	// Loop along integration points
	for (size_t i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double t = _a_int_pts[i_ip].t; // only for 3D cases
		double w = _a_int_pts[i_ip].w;

		Derivs   (r,s,t, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix (derivs, J, B);   // Calculate B matrix for i_ip Integration Point

		// Calculate the gradient
		DGrad = B*dT;
		
		// Update model
		_a_model[i_ip]->UpdateState(DGrad, DFlow);

		// Calculate internal force vector;
		dF += trn(B)*DFlow*det(J)*w;
	}

	// Return internal F
	for (size_t i=0; i<_n_nodes; ++i)
	{
		// Sum up contribution to internal forces vector
		dFint(_connects[i]->DOFVar("F").EqID) += dF(i);
	}
}

inline void HeatElem::BackupState()
{
	for (size_t i=0; i<_n_int_pts; ++i)
		_a_model[i]->BackupState();
}

inline void HeatElem::RestoreState()
{
	for (size_t i=0; i<_n_int_pts; ++i)
		_a_model[i]->RestoreState();
}

inline void HeatElem::GetLabels(Array<String> & Labels) const
{
	// Get labels of all values to output
	Labels.Resize(2);
	Labels[0]="T"; // Nodal temperature
	Labels[1]="F"; // Nodal heat source
}

inline double HeatElem::Val(int iNodeLocal, char const * Name) const
{
	// Essential
	if (strcmp(Name,"T")==0)
		return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Natural
	else if (strcmp(Name,"F")==0)
		return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;
	else
		return 0.0;
}

inline double HeatElem::Val(char const * Name) const
{
	// Get integration point values
	double sum = 0.0;
	for (size_t i_ip=0; i_ip<_n_int_pts; i_ip++)
		sum += _a_model[i_ip]->Val(Name);

	// Output single value at CG
	return sum/_n_int_pts;
}

inline void HeatElem::Deactivate()
{
	throw new Fatal("HeatElem::Deactivate: Feature not implemented in this element");
}

// Derived methods to assemble DAS matrices

inline void HeatElem::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
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
		RowsMap        [idx_Ke] = _connects[i_node]->DOFVar("T").EqID; 
		RowsEssenPresc [idx_Ke] = _connects[i_node]->DOFVar("T").IsEssenPresc; 
		idx_Ke++;
	}
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
}

inline void HeatElem::Order1Matrix(size_t index, LinAlg::Matrix<double> & Ke) const
{
	/* Stiffness:
	   ==========
	
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
	LinAlg::Matrix<double> B;      // strain-displacement matrix
	LinAlg::Matrix<double> D;      // Constitutive matrix

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

		// Constitutive tensor 
		_a_model[i_ip]->TgStiffness(D); 

		// Calculate Tangent Stiffness
		Ke += trn(B)*D*B*det(J)*w;
	}
}
	
// Methods

inline void HeatElem::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	// Calculate Bp matrix (NDIM x _n_nodes)
	B = inv(J)*derivs;     // equal to the cartesian derivatives matrix for flow elements
}


/* private */

inline void HeatElem::_calc_initial_internal_forces() // TODO
{
	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> F(_n_nodes);
	F.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix
	LinAlg::Matrix<double> B;       // strain-displacement matrix
	LinAlg::Vector<double> sig;     // Stress vector in Mandel's notation 

	// Loop along integration points
	for (size_t i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double t = _a_int_pts[i_ip].t; // only for 3D
		//double w = _a_int_pts[i_ip].w;

		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i_ip Integration Point
		B_Matrix (derivs, J, B);  // Calculate B matrix for i_ip Integration Point

		//_a_model[i_ip]->Sig(sig); 

		// Calculate internal force vector;
		//F += trn(B)*sig*det(J)*w;
	}

	// Update nodal NaturVals
//	for (int i_node=0; i_node<_n_nodes; ++i_node)
//	{
//		// Assemble (local/element) displacements vector.
//		if (_ndim_prob==3) _connects[i_node]->DOFVar("F").NaturalVal += F(i_node*_ndim_prob+2);
//	}
}


/////////////////////////////////////////////////////////////////////////////////////////// Map /////


// Register the DOF information of HeatElement into DOFInfoMap
int EquilibDOFInfoRegister()
{
	// Temporary 
	DOFInfo D; 

	// Nodal

	// Insert into DOFInfoMap
	DOFInfoMap["Equilibrium"] = D;

	return 0;
}

// Execute the autoregistration
int __HeatElemDOFInfo_dummy_int  = EquilibDOFInfoRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIB_H
