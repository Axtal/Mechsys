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

#ifndef MECHSYS_FEM_ELASTICBEAM_H
#define MECHSYS_FEM_ELASTICBEAM_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/elems/lin2.h"

namespace FEM
{

class ElasticBeam : public Lin2, public EquilibElem
{
public:
	// Constants
	static char const * NAME;

	// Derived methods
	char const * Name        () const { return NAME; };
	bool         IsEssential (char const * DOFName) const;
	void         SetModel    (char const * ModelName, char const * Prms, char const * Inis);
	Element    * Connect     (int iNodeLocal, FEM::Node * ptNode);
	void         UpdateState (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void         AddVolForces(LinAlg::Vector<double> & FVol) const;
	void         GetLabels   (Array<String> & Labels) const;

	// Derived methods to assemble DAS matrices
	void Order1MatMap (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const;
	void Order1Matrix (size_t Index, LinAlg::Matrix<double> & Ke) const;

	// Access methods
	double Val (int iNodeLocal, char const * Name) const;

private:
	// Data
	double _E;
	double _A;
	double _I3;
	double _sig;
	double _eps;

	// Private methods
	int  _geom () const { return 1; }     ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _calc_initial_internal_state (); ///< Calculate initial internal state

}; // class ElasticBeam

// ElasticBeam constants
char const * ElasticBeam::NAME = "ElasticBeam";


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


/* public */

inline bool ElasticBeam::IsEssential(char const * DOFName) const
{
	if (strcmp(DOFName,"ux")==0 || strcmp(DOFName,"uy")==0 || strcmp(DOFName,"rz"))                         return true;
	if (_ndim==3                && strcmp(DOFName,"uz")==0 || strcmp(DOFName,"rx") || strcmp(DOFName,"ry")) return true;
	return false;
}

inline void ElasticBeam::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("ElasticBeam::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("ElasticBeam::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	// If pointers to model was not already defined => No model was allocated
	if (_a_model.Size()==0)
	{
		/* "E=20000.0 nu=0.2" */
		LineParser lp(Prms);
		Array<String> names;
		Array<double> values;
		lp.BreakExpressions(names,values);

		// Set parameters
		for (size_t i=0; i<names.Size(); ++i)
		{
				 if (names[i]=="E" ) _E  = values[i];
			else if (names[i]=="A" ) _A  = values[i];
			else if (names[i]=="I3") _I3 = values[i];
		}

		/* "Sx=0.0 Sy=0.0 Sxy=0.0 ..." or "ZERO" */
		lp.Reset(Inis);
		lp.BreakExpressions(names,values);

		// Parse input
		_sig = 0.0;
		_eps = 0.0;
		for (size_t i=0; i<names.Size(); i++)
		{
				 if (names[i]=="ZERO") break;
			else if (names[i]=="Sa")   _sig = values[i];
			else throw new Fatal("ElasticBeam::SetModel: '%s' component of stress is invalid",names[i].CStr());
		}

		// Calculate initial internal state
		_calc_initial_internal_state ();
	}
	else throw new Fatal("EquilibElem::SetModel: Feature not implemented.");
}


inline Element * ElasticBeam::Connect(int iNodeLocal, FEM::Node * ptNode)
{
	// Connects
	_connects[iNodeLocal] = ptNode;

	// Add Degree of Freedom to a node (Essential, Natural)
	if (_ndim==3)
	{
		_connects[iNodeLocal]->AddDOF("ux", "fx");
		_connects[iNodeLocal]->AddDOF("uy", "fy");
		_connects[iNodeLocal]->AddDOF("uz", "fz");
		_connects[iNodeLocal]->AddDOF("rx", "mx");
		_connects[iNodeLocal]->AddDOF("ry", "my");
		_connects[iNodeLocal]->AddDOF("rz", "mz");
	}
	else
	{
		_connects[iNodeLocal]->AddDOF("ux", "fx");
		_connects[iNodeLocal]->AddDOF("uy", "fy");
		_connects[iNodeLocal]->AddDOF("rz", "mz");
	}

	// Shared
	_connects[iNodeLocal]->SetSharedBy(_my_id);

	return this;
}

inline void ElasticBeam::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	std::cout << "ElasticBeam::UpdateState: HI <<<<<<<<<<" << std::endl;
}

inline void ElasticBeam::AddVolForces(LinAlg::Vector<double> & FVol) const
{
	std::cout << "ElasticBeam::AddVolForces: HI <<<<<<<<<<" << std::endl;
}

inline void ElasticBeam::GetLabels(Array<String> & Labels) const
{
	switch (_ndim)
	{
		case 2:
		{
			Labels.Resize(8);
			Labels[0] = "ux"; Labels[1] = "uy"; Labels[2] = "rz";
			Labels[3] = "fx"; Labels[4] = "fy"; Labels[5] = "mz";
			Labels[6] = "Ea"; // axial strain
			Labels[7] = "Sa"; // axial stress
			return;
		}
		default: throw new Fatal("ElasticBeam::GetLabels: GeometryType==%d & NDim==%d is not implemented yet",_geom(),_ndim);
	}
}

inline void ElasticBeam::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	// Size of Ke
	int n_rows = _ndim*_n_nodes; // == n_cols

	// Mounting a map of positions from Ke to Global
	int idx_Ke = 0;                // position (idx) inside Ke matrix
	RowsMap       .Resize(n_rows); // size=Ke.Rows()=Ke.Cols()
	RowsEssenPresc.Resize(n_rows); // size=Ke.Rows()=Ke.Cols()

	// Fill map of Ke position to K position of DOFs components
	for (size_t i=0; i<_n_nodes; ++i)
	{
		if (_ndim==3)
		{
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("ux").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("ux").IsEssenPresc; 
			idx_Ke++;
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("uy").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("uy").IsEssenPresc; 
			idx_Ke++;
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("uz").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("uz").IsEssenPresc; 
			idx_Ke++;
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("rx").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("rx").IsEssenPresc; 
			idx_Ke++;
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("ry").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("ry").IsEssenPresc; 
			idx_Ke++;
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("rz").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("rz").IsEssenPresc; 
			idx_Ke++;
		}
		else
		{
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("ux").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("ux").IsEssenPresc; 
			idx_Ke++;
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("uy").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("uy").IsEssenPresc; 
			idx_Ke++;
			RowsMap        [idx_Ke] = _connects[i]->DOFVar("rz").EqID; 
			RowsEssenPresc [idx_Ke] = _connects[i]->DOFVar("rz").IsEssenPresc; 
			idx_Ke++;
		}
	}
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
}

inline void ElasticBeam::Order1Matrix(size_t index, LinAlg::Matrix<double> & Ke) const
{
	if (_ndim==3)
	{
		throw new Fatal("ElasticBeam:Order1Matrix: Feature not available yet (ndim=3D)");
	}
	else
	{
		double l   = sqrt(pow(_connects[0]->X()-_connects[1]->X(),2.0)+pow(_connects[0]->Y()-_connects[1]->Y(),2.0));
		double c   = (_connects[1]->X()-_connects[0]->X())/l;
		double s   = (_connects[1]->Y()-_connects[0]->Y())/l;
		double cc  = c*c;
		double ss  = s*s;
		double ll  = l*l;
		double lll = l*l*l;
		double a1  = (12.0*ss*_E*_I3)/lll+(cc*_A*_E)/l;
		double a2  = (c*s*_A*_E)/l-(12.0*c*s*_E*_I3)/lll;
		double a3  = -(6.0*s*_E*_I3)/ll;
		double a4  = (12.0*cc*_E*_I3)/lll+(ss*_A*_E)/l;
		double a5  = (6.0*c*_E*_I3)/ll;
		double a6  = (4.0*_E*_I3)/l;
		double a7  = (2.0*_E*_I3)/l;
		Ke.Resize(6,6);
		Ke =  a1,  a2,  a3, -a1, -a2,  a3,
		      a2,  a4,  a5, -a2, -a4,  a5,
		      a3,  a5,  a6, -a3, -a5,  a7,
		     -a1, -a2, -a3,  a1,  a2, -a3,
		     -a2, -a4, -a5,  a2,  a4, -a5,
		      a3,  a5,  a7, -a3, -a5,  a6;
	}
}

inline double ElasticBeam::Val(int iNodeLocal, char const * Name) const
{
	// Displacements and rotations
	if (strcmp(Name,"ux")==0 || strcmp(Name,"uy")==0 || strcmp(Name,"uz")==0 ||
	    strcmp(Name,"rx")==0 || strcmp(Name,"ry")==0 || strcmp(Name,"rz")==0)
		return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces and momenta
	else if (strcmp(Name,"fx")==0 || strcmp(Name,"fy")==0 || strcmp(Name,"fz")==0 ||
	         strcmp(Name,"mx")==0 || strcmp(Name,"my")==0 || strcmp(Name,"mz")==0)
		return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	else
	{
		// Vectors for extrapolation
		LinAlg::Vector<double>    ip_values (_a_int_pts.Size());
		LinAlg::Vector<double> nodal_values (_n_nodes);

		// Get integration point values
		if (_a_model.Size()==_a_int_pts.Size())
			for (size_t i=0; i<_a_int_pts.Size(); i++)
				ip_values(i) = _a_model[i]->Val(Name);
		else throw new Fatal("ElasticBeam::Val: Constitutive models for this element (ID==%d) were not set yet", _my_id);

		// Extrapolation
		Extrapolate (ip_values, nodal_values);

		// Output single value
		return nodal_values (iNodeLocal);
	}
}


/* private */

inline void ElasticBeam::_calc_initial_internal_state()
{
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new ElasticBeam element
Element * ElasticBeamMaker()
{
	return new ElasticBeam();
}

// Register a ElasticBeam element into ElementFactory array map
int ElasticBeamRegister()
{
	ElementFactory[ElasticBeam::NAME] = ElasticBeamMaker;
	return 0;
}

// Execute the autoregistration
int __ElasticBeam_dummy_int  = ElasticBeamRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_ELASTICBEAM_H
