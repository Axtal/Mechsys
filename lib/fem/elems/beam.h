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
#include "fem/elems/lin2.h"
#include "fem/equilibelem.h"
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Beam : public EquilibElem
{
public:
	// Constants
	static const size_t NDIV_EXTRA;        ///< Number of points for extra output
	
	// Constructor
	Beam () : _E(-1), _A(-1), _Izz(-1), _q0(0.0), _q1(0.0), _has_q(false), _use_cor(true) {}

	// Derived methods
	Str_t Name() const { return "Beam"; }

	// Derived methods
	void   AddVolForces ();
	void   ClearDisp    ();
	void   EdgeBry      (Str_t Key, double q, int iEdge) { return EdgeBry(Key,q,q,iEdge); } ///< Set distributed load with key = q0 or q1
	void   EdgeBry      (Str_t Key, double q0, double q1, int iEdge);                       ///< Set distributed load with key = q0 or q1
	void   CalcDeps     () const;
	Str_t  ModelName    () const { return "__no_model__"; }
	double Val          (int iNod, Str_t Name) const;
	void   SetProps     (Str_t Properties);
	void   SetModel     (Str_t ModelName, Str_t Prms, Str_t Inis);
	void   Update       (double h, Vec_t const & dU, Vec_t & dFint);
	bool   HasExtra     () const { return true; }
	void   OutExtra     (Mat_t & Coords, Vec_t & Norm, Mat_t & Values, Array<String> & Labels) const;
	void   CMatrix      (size_t Idx, Mat_t & Ke) const;
	bool   CheckModel   () const;

	// Methods
	double N(double l) const; ///< Axial force      (0 < l < 1) (Must be used after CalcDeps())
	double M(double l) const; ///< Bending momentum (0 < l < 1) (Must be used after CalcDeps())
	double V(double l) const; ///< Shear force      (0 < l < 1) (Must be used after CalcDeps())

private:
	// Data
	double _E;       ///< Young modulus
	double _A;       ///< Cross-sectional area
	double _Izz;     ///< Cross-sectional inertia
	double _q0;      ///< Normal distributed load (value at node # 0. Or constant q)
	double _q1;      ///< Normal distributed load (value at node # 1. Or constant q)
	bool   _has_q;   ///< Has distributed load (q0 and/or q1)
	bool   _use_cor; ///< Use correction for distributed (q) load?

	// Depedent variables (calculated by CalcDeps)
	mutable double _L;  ///< Beam length
	mutable Vec_t  _uL; ///< Beam-Local displacements/rotations

	// Private methods
	void _initialize ();                ///< Initialize the element
	void _transf_mat (Mat_t & T) const; ///< Calculate transformation matrix

}; // class Beam

const size_t Beam::NDIV_EXTRA = 10;  ///< Number of points for extra output


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void Beam::EdgeBry(Str_t Key, double q0, double q1, int iEdge)
{
	// Skip if key is not "Qb"
	if (strcmp(Key,"Qb")!=0) return;

	// Check which node is the left-most and adjust q0 and q1
	bool is_node0_leftmost = true;
	double x0 = _ge->Conn[0]->X();  double y0 = _ge->Conn[0]->Y();
	double x1 = _ge->Conn[1]->X();  double y1 = _ge->Conn[1]->Y();
	if (fabs(x1-x0)<1.0e-5) { if (y1<y0) is_node0_leftmost = false; } // vertical segment
	else                    { if (x1<x0) is_node0_leftmost = false; } 
	if (is_node0_leftmost==false)
	{
		q0 = -q0;
		q1 = -q1;
	}

	// Transformation matrix and calculate _L
	Mat_t T;
	_transf_mat(T);

	_has_q = true;
	_q0    = q0;
	_q1    = q1;
	double LL = _L*_L;

	// Beam-Local increment of force
	Vec_t f(_nd*_ge->NNodes);
	f = 0.0, _L*(7.0*_q0+3.0*_q1)/20.0,  LL*(3.0*_q0+2.0*_q1)/60.0,
	    0.0, _L*(3.0*_q0+7.0*_q1)/20.0, -LL*(2.0*_q0+3.0*_q1)/60.0;

	// Increment of force in global coordinates
	f = inv(T)*f;

	// Add to nodes Brys
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->Bry (FD[j], f(i*_nd+j));
}

inline bool Beam::CheckModel() const
{
	if (_E<0.0 || _A<0.0 || _Izz<0.0) return false;
	return true;
}

inline void Beam::SetModel(Str_t ModelName, Str_t Prms, Str_t Inis)
{
	// Check _ge->NDim
	if (_ge->NDim<1)             throw new Fatal("Beam::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (_ge->CheckConn()==false) throw new Fatal("Beam::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");
	if (_ge->NNodes!=2)          throw new Fatal("Beam::SetModel: The number of nodes (%d) of Geometry element is incorret. It must be equal to 2",_ge->NNodes);

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

inline void Beam::SetProps(Str_t Properties)
{
	/* "cq=1 gam=20 */
	LineParser lp(Properties);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		      if (names[i]=="cq")  _use_cor = static_cast<bool>(values[i]); // cq == correct q
		 else if (names[i]=="gam") _gam     = values[i];
	}
}

inline void Beam::Update(double h, Vec_t const & dU, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*_ge->NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		du(i*_nd+j) = dU(_ge->Conn[i]->DOFVar(UD[j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df(_nd*_ge->NNodes); // Delta internal force of this element
	df.SetValues(0.0);

	Mat_t Ke;
	CMatrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[j]).EqID) += df(i*_nd+j);
}

inline void Beam::AddVolForces() 
{
	// Verify if element is active
	if (IsActive==false) return;

	// Weight
	double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
	double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
	double L  = sqrt(dx*dx+dy*dy);
	double W  = _A*L*_gam;

	// Set boundary conditions
	if (_ge->NDim==2)
	{
		_ge->Conn[0]->Bry("fy", -W/2.0);
		_ge->Conn[1]->Bry("fy", -W/2.0);
	}
	else if (_ge->NDim==3)
	{
		_ge->Conn[0]->Bry("fz", -W/2.0);
		_ge->Conn[1]->Bry("fz", -W/2.0);
	}
}

inline void Beam::CalcDeps() const
{
	if (IsActive==false) throw new Fatal("Beam::CalcDeps: This element is inactive");

	// Element displacements vector
	_uL.Resize(_nd*_ge->NNodes);
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_uL(i*_nd+j) = _ge->Conn[i]->DOFVar(UD[j]).EssentialVal;

	// Transform to beam-local coordinates
	Mat_t T;
	_transf_mat(T);
	_uL = T * _uL;
}

inline double Beam::Val(int iNod, Str_t Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[j])==0) return _ge->Conn[iNod]->DOFVar(Name).NaturalVal;

	if (_uL.Size()<1) throw new Fatal("Beam::Val: Please, call CalcDeps() before calling this method");
	double l = (iNod==0 ? 0 : 1.0);
	     if (strcmp(Name,"N" )==0) return      N(l);
	else if (strcmp(Name,"M" )==0) return fabs(M(l));
	else if (strcmp(Name,"V" )==0) return fabs(V(l));
	else if (strcmp(Name,"Ea")==0) return    (_uL(_nd)-_uL(0))/_L;
	else if (strcmp(Name,"Sa")==0) return _E*(_uL(_nd)-_uL(0))/_L;
	else throw new Fatal("Beam::Val: This element does not have a Val named %s",Name);
}

inline void Beam::ClearDisp()
{
	if (IsActive==false) return;

	// Clear displacements
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->DOFVar(UD[j]).EssentialVal = 0.0;
}

inline void Beam::CMatrix(size_t Idx, Mat_t & Ke) const
{
	if (_ge->NDim==2)
	{
		double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
		double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
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
		Ke.Resize(_nd*_ge->NNodes, _nd*_ge->NNodes);
		Ke =  c1,  c2, -c3, -c1, -c2, -c3,
		      c2,  c4,  c5, -c2, -c4,  c5,
		     -c3,  c5,  c6,  c3, -c5,  c7,
		     -c1, -c2,  c3,  c1,  c2,  c3,
		     -c2, -c4, -c5,  c2,  c4, -c5,
		     -c3,  c5,  c7,  c3, -c5,  c6;
	}
	else throw new Fatal("Beam::CMatrix: Feature not available for nDim==%d",_ge->NDim);
}

inline double Beam::N(double l) const
{
	return _E*_A*(_uL(_nd)-_uL(0))/_L;
}

inline double Beam::M(double l) const
{
	double M = 0.0;
	if (_ge->NDim==2)
	{
		double s   = l*_L;
		double LL  = _L*_L;
		double LLL = LL*_L;
		M = _E*_Izz*(_uL(5)*((6*s)/LL-2/_L)+_uL(2)*((6*s)/LL-4/_L)+_uL(4)*(6/LL-(12*s)/LLL)+_uL(1)*((12*s)/LLL-6/LL));
		if (_has_q && _use_cor)
		{
			double ss  = s*s;
			double sss = ss*s;
			M += (2.0*_q1*LLL+3.0*_q0*LLL-9.0*_q1*s*LL-21.0*_q0*s*LL+30.0*_q0*ss*_L+10.0*_q1*sss-10.0*_q0*sss)/(60.0*_L);
		}
	}
	else throw new Fatal("Beam::M: Feature not available for nDim==%d",_ge->NDim);
	return M;
}

inline double Beam::V(double l) const
{
	double V = 0.0;
	if (_ge->NDim==2)
	{
		double LL  = _L*_L;
		double LLL = LL*_L;
		V = _E*_Izz*((6*_uL(5))/LL+(6*_uL(2))/LL-(12*_uL(4))/LLL+(12*_uL(1))/LLL);
		if (_has_q && _use_cor)
		{
			double s  = l*_L;
			double ss = s*s;
			V += -(3.0*_q1*LL+7.0*_q0*LL-20.0*_q0*s*_L-10.0*_q1*ss+10.0*_q0*ss)/(20.0*_L);
		}
	}
	else throw new Fatal("Beam::V: Feature not available for nDim==%d",_ge->NDim);
	return V;
}

inline void Beam::OutExtra(Mat_t & Coords, Vec_t & Norm, Mat_t & Values, Array<String> & Labels) const
{
	if (_uL.Size()<1) throw new Fatal("Beam::OutExtra: Please, call CalcDeps() before calling this method");
	if (_ge->NDim==2)
	{
		// Generate coordinates for the extra points
		double x0  = _ge->Conn[0]->X();
		double y0  = _ge->Conn[0]->Y();
		double x1  = _ge->Conn[1]->X();
		double y1  = _ge->Conn[1]->Y();
		double len = sqrt(pow(x1-x0,2) + pow(y1-y0,2));
		Coords.Resize(NDIV_EXTRA+1, 2);
		for (size_t i=0; i<NDIV_EXTRA+1; i++)
		{
			Coords(i,0) = x0 + i*(x1-x0)/NDIV_EXTRA;
			Coords(i,1) = y0 + i*(y1-y0)/NDIV_EXTRA;
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
		Values.Resize(NDIV_EXTRA+1, Labels.Size());
		for (size_t i=0; i<NDIV_EXTRA+1; i++)
		{
			double l = static_cast<double>(i)/NDIV_EXTRA;
			Values(i,0) = M(l);
			Values(i,1) = N(l);
			Values(i,2) = V(l);
		}
	}
	else throw new Fatal("Beam::OutExtra: Feature not available for nDim==%d",_ge->NDim);
}


/* private */

inline void Beam::_initialize()
{
	if (_ge->NDim==2)
	{
		_gi = 1;
		_nd = ND_BEAM_2D;
		UD  = UD_BEAM_2D;
		FD  = FD_BEAM_2D;
		_nl = NL_BEAM_2D;
		LB  = LB_BEAM_2D;
	}
	else if (_ge->NDim==3)
	{
		_gi = 0;
		_nd = ND_BEAM_3D;
		UD  = UD_BEAM_3D;
		FD  = FD_BEAM_3D;
		_nl = NL_BEAM_3D;
		LB  = LB_BEAM_3D;
	}
	else throw new Fatal("Beam::_initialize: nDim==%d is invalid",_ge->NDim);
}

inline void Beam::_transf_mat(Mat_t & T) const
{
	// Transformation matrix
	if (_ge->NDim==2)
	{
		double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
		double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
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
	else throw new Fatal("Beam::_transf_mat: Feature not available for nDim==%d",_ge->NDim);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new Beam element
ProbElem * BeamMaker() { return new Beam(); }

// Register element
int BeamRegister() { ProbElemFactory["Beam"]=BeamMaker;  return 0; }

// Call register
int __Beam_dummy_int = BeamRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_BEAM_H
