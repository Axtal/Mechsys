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
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Beam : public EquilibElem
{
public:
	// Constants
	static const size_t NDIV;        ///< Number of points for extra output
	
	// Constructor
	Beam () : _E(-1), _A(-1), _Izz(-1), _q0(0.0), _q1(0.0), _has_q(false) { _n_nodes=2; _connects.Resize(_n_nodes); _connects.SetValues(NULL); }

	// Derived methods
	char const * Name() const { return "Beam"; }

	// Derived methods
	bool   CheckModel   () const;
	void   SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void   UpdateState  (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void   CalcDepVars  () const;
	double Val          (int iNodeLocal, char const * Name) const;
	double Val          (char const * Name) const;
	void   Order1Matrix (size_t Index, LinAlg::Matrix<double> & Ke) const; ///< Stiffness
	void   B_Matrix     (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;
	int    VTKCellType  () const { return VTK_LINE; }
	void   VTKConnect   (String & Nodes) const { Nodes.Printf("%d %d",_connects[0]->GetID(),_connects[1]->GetID()); }
	void   OutExtra     (LinAlg::Matrix<double> & Coords, LinAlg::Vector<double> & Norm, LinAlg::Matrix<double> & Values, Array<String> & Labels) const;

	// Methods
	Beam * EdgeBry (char const * Key, double q, int EdgeLocalID) { return EdgeBry(Key,q,q,EdgeLocalID); } ///< Set distributed load with key = q0 or q1
	Beam * EdgeBry (char const * Key, double q0, double q1, int EdgeLocalID);                             ///< Set distributed load with key = q0 or q1

	// Methods
	double N(double l) const; ///< Axial force      (0 < l < 1) (Must be used after CalcDepVars())
	double M(double l) const; ///< Bending momentum (0 < l < 1) (Must be used after CalcDepVars())
	double V(double l) const; ///< Shear force      (0 < l < 1) (Must be used after CalcDepVars())

private:
	// Data
	double _E;     ///< Young modulus
	double _A;     ///< Cross-sectional area
	double _Izz;   ///< Cross-sectional inertia
	double _q0;    ///< Normal distributed load (value at node # 0. Or constant q)
	double _q1;    ///< Normal distributed load (value at node # 1. Or constant q)
	bool   _has_q; ///< Has distributed load (q0 and/or q1)

	// Depedent variables (calculated by CalcDepVars)
	mutable double         _L;  ///< Beam length
	mutable Vector<double> _uL; ///< Beam-Local displacements/rotations

	// Private methods
	int  _geom                        () const { return 1; }              ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _set_ndim                    (int nDim);                         ///< Set space dimension
	void _calc_initial_internal_state ();                                 ///< Calculate initial internal state
	void _transf_mat                  (LinAlg::Matrix<double> & T) const; ///< Calculate transformation matrix

}; // class Beam

const size_t Beam::NDIV = 10;  ///< Number of points for extra output


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline Beam * Beam::EdgeBry(char const * Key, double q0, double q1, int EdgeLocalID)
{
	// Transformation matrix and calculate _L
	LinAlg::Matrix<double> T;
	_transf_mat(T);

	_has_q = true;
	_q0    = q0;
	_q1    = q1;
	double LL = _L*_L;

	// Beam-Local increment of force
	LinAlg::Vector<double> f(_nd*_n_nodes);
	f = 0.0, _L*(7.0*_q0+3.0*_q1)/20.0,  LL*(3.0*_q0+2.0*_q1)/60.0,
	    0.0, _L*(3.0*_q0+7.0*_q1)/20.0, -LL*(2.0*_q0+3.0*_q1)/60.0;

	// Increment of force in global coordinates
	f = inv(T)*f;

	// Add to nodes Brys
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_connects[i]->Bry (FD[_d][j], f(i*_nd+j));

	return this;
}

inline bool Beam::CheckModel() const
{
	if (_E<0.0 || _A<0.0 || _Izz<0.0) return false;
	return true;
}

inline void Beam::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1) throw new Fatal("Beam::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("Beam::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	/* "E=1 A=1 Izz=1" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		     if (names[i]=="E")    _E   = values[i];
		else if (names[i]=="A")    _A   = values[i];
		else if (names[i]=="Izz" ) _Izz = values[i];
		else throw new Fatal("Beam::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}
}

inline void Beam::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> du(_nd*_n_nodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dUglobal(_connects[i]->DOFVar(UD[_d][j]).EqID);

	// Allocate (local/element) internal force vector
	LinAlg::Vector<double> df(_nd*_n_nodes); // Delta internal force of this element
	df.SetValues(0.0);

	LinAlg::Matrix<double> Ke;
	Order1Matrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_connects[i]->DOFVar(UD[_d][j]).EqID) += df(i*_nd+j);
}

inline void Beam::CalcDepVars() const
{
	// Element displacements vector
	_uL.Resize(_nd*_n_nodes);
	for (size_t i=0; i<_n_nodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_uL(i*_nd+j) = _connects[i]->DOFVar(UD[_d][j]).EssentialVal;

	// Transform to beam-local coordinates
	LinAlg::Matrix<double> T;
	_transf_mat(T);
	_uL = T * _uL;
}

inline double Beam::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	if (_uL.Size()<1) throw new Fatal("Beam::Val: Please, call CalcDepVars() before calling this method");
	double l = (iNodeLocal==0 ? 0 : 1.0);
	     if (strcmp(Name,"N" )==0) return N(l);
	else if (strcmp(Name,"M" )==0) return M(l);
	else if (strcmp(Name,"V" )==0) return V(l);
	else if (strcmp(Name,"Ea")==0) return    (_uL(_nd)-_uL(0))/_L;
	else if (strcmp(Name,"Sa")==0) return _E*(_uL(_nd)-_uL(0))/_L;
	else throw new Fatal("Beam::Val: This element does not have a Val named %s",Name);
}

inline double Beam::Val(char const * Name) const
{
	throw new Fatal("Beam::Val: Feature not available");
}

inline void Beam::Order1Matrix(size_t Index, LinAlg::Matrix<double> & Ke) const
{
	if (_ndim==2)
	{
		double dx = _connects[1]->X()-_connects[0]->X();
		double dy = _connects[1]->Y()-_connects[0]->Y();
		double LL = dx*dx+dy*dy;
		      _L  = sqrt(LL);
		double c  = dx/_L;
		double s  = dy/_L;
		double c1 = _E*(_A*c*c+12.0*_Izz*s*s/LL)/_L;
		double c2 = _E*((_A-12.0*_Izz/LL)*c*s)/_L;
		double c3 = _E*(6.0*_Izz*s/_L)/_L;
		double c4 = _E*(_A*s*s+12.0*_Izz*c*c/LL)/_L;
		double c5 = _E*(6.0*_Izz*c/_L)/_L;
		double c6 = _E*(4.0*_Izz)/_L;
		double c7 = _E*(2.0*_Izz)/_L;
		Ke.Resize(_nd*_n_nodes, _nd*_n_nodes);
		Ke =  c1,  c2, -c3, -c1, -c2, -c3,
		      c2,  c4,  c5, -c2, -c4,  c5,
		     -c3,  c5,  c6,  c3, -c5,  c7,
		     -c1, -c2,  c3,  c1,  c2,  c3,
		     -c2, -c4, -c5,  c2,  c4, -c5,
		     -c3,  c5,  c7,  c3, -c5,  c6;
	}
	else throw new Fatal("Beam::Order1Matrix: Feature not available for nDim==%d",_ndim);
}

inline void Beam::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	throw new Fatal("Beam::B_Matrix: Feature not available");
}

inline double Beam::N(double l) const
{
	return _E*_A*(_uL(_nd)-_uL(0))/_L;
}

inline double Beam::M(double l) const
{
	double M = 0.0;
	if (_ndim==2)
	{
		double s   = l*_L;
		double LL  = _L*_L;
		double LLL = LL*_L;
		M = _E*_Izz*(_uL(5)*((6*s)/LL-2/_L)+_uL(2)*((6*s)/LL-4/_L)+_uL(4)*(6/LL-(12*s)/LLL)+_uL(1)*((12*s)/LLL-6/LL));
		if (_has_q)
		{
			double ss  = s*s;
			double sss = ss*s;
			M += (2.0*_q1*LLL+3.0*_q0*LLL-9.0*_q1*s*LL-21.0*_q0*s*LL+30.0*_q0*ss*_L+10.0*_q1*sss-10.0*_q0*sss)/(60.0*_L);
		}
	}
	else throw new Fatal("Beam::M: Feature not available for nDim==%d",_ndim);
	return M;
}

inline double Beam::V(double l) const
{
	double V = 0.0;
	if (_ndim==2)
	{
		double LL  = _L*_L;
		double LLL = LL*_L;
		V = _E*_Izz*((6*_uL(5))/LL+(6*_uL(2))/LL-(12*_uL(4))/LLL+(12*_uL(1))/LLL);
		if (_has_q)
		{
			double s  = l*_L;
			double ss = s*s;
			V += -(3.0*_q1*LL+7.0*_q0*LL-20.0*_q0*s*_L-10.0*_q1*ss+10.0*_q0*ss)/(20.0*_L);
		}
	}
	else throw new Fatal("Beam::V: Feature not available for nDim==%d",_ndim);
	return V;
}

inline void Beam::OutExtra(LinAlg::Matrix<double> & Coords, LinAlg::Vector<double> & Norm, LinAlg::Matrix<double> & Values, Array<String> & Labels) const
{
	if (_ndim==2)
	{
		// Generate coordinates for the extra points
		double x0  = _connects[0]->X();
		double y0  = _connects[0]->Y();
		double x1  = _connects[1]->X();
		double y1  = _connects[1]->Y();
		double len = sqrt(pow(x1-x0,2) + pow(y1-y0,2));
		Coords.Resize(NDIV+1, 2);
		for (size_t i=0; i<NDIV+1; i++)
		{
			Coords(i,0) = x0 + i*(x1-x0)/NDIV;
			Coords(i,1) = y0 + i*(y1-y0)/NDIV;
		}

		// Normal vector
		Norm.Resize(2);
		double v = (x1-x0)/len;
		double w = (y1-y0)/len;
		Norm(0) = -w;
		Norm(1) =  v;

		// Labels
		Labels.Resize(3);
		Labels[0] = "M"; Labels[1] = "N"; Labels[2] = "V";

		// Values
		Values.Resize(NDIV+1, Labels.Size());
		for (size_t i=0; i<NDIV+1; i++)
		{
			double l = static_cast<double>(i)/NDIV;
			Values(i,0) = M(l);
			Values(i,1) = N(l);
			Values(i,2) = V(l);
		}
	}
	else throw new Fatal("Beam::OutExtra: Feature not available for nDim==%d",_ndim);
}

/* private */

inline void Beam::_set_ndim(int nDim)
{
	if (nDim<1) throw new Fatal("Beam::_set_ndim: For this element, nDim must be greater than or equal to 1 (%d is invalid)",nDim);
	_ndim = nDim;
	_d    = _ndim-1;
	_nd   = EquilibElem::NDB[_d];
}

inline void Beam::_calc_initial_internal_state()
{
	throw new Fatal("Beam::_calc_initial_internal_state: Feature not available");
}

inline void Beam::_transf_mat(LinAlg::Matrix<double> & T) const
{
	// Transformation matrix
	if (_ndim==2)
	{
		double dx = _connects[1]->X()-_connects[0]->X();
		double dy = _connects[1]->Y()-_connects[0]->Y();
		double LL = dx*dx+dy*dy;
		      _L  = sqrt(LL);
		double c  = dx/_L;
		double s  = dy/_L;
		T.Resize(6,6);
		T =    c,   s, 0.0, 0.0, 0.0, 0.0,
		      -s,   c, 0.0, 0.0, 0.0, 0.0,
		     0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
		     0.0, 0.0, 0.0,   c,   s, 0.0,
		     0.0, 0.0, 0.0,  -s,   c, 0.0,
		     0.0, 0.0, 0.0, 0.0, 0.0, 1.0;
	}
	else throw new Fatal("Beam::_transf_mat: Feature not available for nDim==%d",_ndim);
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
	ElementFactory["Beam"] = BeamMaker;
	return 0;
}

// Execute the autoregistration
int __Beam_dummy_int  = BeamRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_BEAM_H
