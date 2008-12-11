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

#ifndef MECHSYS_FEM_ROD_H
#define MECHSYS_FEM_ROD_H

// MechSys
#include "fem/equilibelem.h"
#include "fem/probelem.h"
#include "util/exception.h"
#include "fem/elems/vtkCellType.h"

namespace FEM
{

class Rod: public ProbElem
{
public:
	// Constants. Note: _di==dimension index, _gi==geometry index
	static const char   UD [2][3][4];  ///< Essential DOF names. Access: UD[_di][iDOF]
	static const char   FD [2][3][4];  ///< Natural DOF names.   Access: FD[_di][iDOF]
	static const char   LB [3][4];     ///< Additional lbls (exceed. those from UD/FD). Access: LB[_gi][iLbl]

	// Constructor
	Rod () : _gam(0.0), _E(-1), _A(-1) { };

	// Derived methods
	char const * Name() const { return "Rod"; }

	// Derived methods
	bool   CheckModel   () const;
	void   SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void   SetProps     (char const * Properties);
	void    Update       (double h, Vec_t const & dU, Vec_t & dFint) { }
	void    Backup       () { }
	void    Restore      () { }
	void   AddVolForces();
	void   CalcDepVars  () const;
	double Val          (int iNodeLocal, char const * Name) const;
	double Val          (char const * Name) const;
	bool   IsEssen      (Str_t Name) const;
	void   CMatrix      (size_t Index, Mat_t & Ke) const; ///< Stiffness
	void   B_Matrix     (Mat_t const & derivs, Mat_t const & J, Mat_t & B) const;
	void   SetConn(int iNod, FEM::Node * ptNode, int ID);
	int    VTKCellType  () const { return VTK_LINE; }
	void   VTKConn   (String & Nodes) const { Nodes.Printf("%d %d", _ge->Conn[0]->GetID(), _ge->Conn[1]->GetID()); }
	void   ClearDisp();
	Str_t   ModelName    () const { return "__no_model__"; }
	
	void   SetActive(bool Activate, int ID);
	

	// Methods
	double N(double l) const; ///< Axial force (0 < l < 1) (Must be used after CalcDepVars())

private:
	// Data
	double _gam; ///< Specific weigth
	double _E; ///< Young modulus
	double _A; ///< Cross-sectional area

	int    _nd;
	int    _di;
	int    _nl;

	// Depedent variables (calculated by CalcDepVars)
	mutable double         _L;  ///< Rod length
	mutable Vector<double> _uL; ///< Rod-Local displacements/rotations

	// Private methods
	int  _geom                        () const { return 1; }              ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	void _initialize                  ();                                 ///< Initialize the element
	void _calc_initial_internal_state ();                                 ///< Calculate initial internal state
	void _transf_mat                  (Mat_t & T) const; ///< Calculate transformation matrix

}; // class Rod

// UD[_di][iDOF]                         2D               3D
const char Rod:: UD [2][3][4] = {{"ux","uy",""},  {"ux","uy","uz"}};

const char Rod:: FD [2][3][4] = {{"fx","fy",""},  {"fx","fy","fz"}};
const char Rod:: LB [3][4]    = {"Ea", "Sa", "N" }; // 2D and 3D


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void Rod::ClearDisp()
{
	if (IsActive==false) return;

	// Clear displacements
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;         ++j)
		_ge->Conn[i]->DOFVar(UD[_di][j]).EssentialVal = 0.0;

}

inline bool Rod::IsEssen(Str_t Name) const
{
	for (int i=0; i<_nd; ++i) if (strcmp(Name,UD[_di][i])==0) return true;
	return false;
}

inline void Rod::SetActive(bool Activate, int ID)
{
	if (IsActive==false && Activate)
	{
		// Set active
		IsActive = true;

		for (size_t i=0; i<_ge->Conn.Size(); ++i)
		{
			// Add Degree of Freedom to a node (Essential, Natural)
			for (int j=0; j<_nd; ++j) _ge->Conn[i]->AddDOF (UD[_di][j], FD[_di][j]);

			// Set SharedBy
			_ge->Conn[i]->SetSharedBy (ID);
		}

	}
	if (IsActive && Activate==false)
	{
		// Set active
		IsActive = false;

		for (size_t i=0; i<_ge->Conn.Size(); ++i)
		{
			// Remove SharedBy
			_ge->Conn[i]->RemoveSharedBy (ID);

			// Remove Degree of Freedom to a node (Essential)
			if (_ge->Conn[i]->nSharedBy()==0) 
				for (int j=0; j<_nd; ++j) _ge->Conn[i]->RemoveDOF (UD[_di][j]);
		}

	}
}

inline bool Rod::CheckModel() const
{
	if (_E<0.0 || _A<0.0) return false;
	return true;
}

inline void Rod::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check NDim
	if (_ge->NDim<1) throw new Fatal("Rod::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (_ge->CheckConn()==false) throw new Fatal("Rod::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	/* "E=1 A=1" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		     if (names[i]=="E") _E = values[i];
		else if (names[i]=="A") _A = values[i];
		else throw new Fatal("Rod::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}
}

inline void Rod::SetProps(char const * Properties)
{
	/* "gam=20 */
	LineParser lp(Properties);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	for (size_t i=0; i<names.Size(); ++i)
	{
		 if (names[i]=="gam") _gam = values[i];
	}
}

inline void Rod::SetConn(int iNod, FEM::Node * ptNode, int ID)
{
	// Check
	if (_ge->NNodes<1)              throw new Fatal("Rod::Connect: __Internal Error__: There is a problem with the number of nodes: maybe derived elemet did not set _ge->NNodes");
	if (_ge->Conn.Size()<1)         throw new Fatal("Rod::Connect: __Internal Error__: There is a problem with connectivity array: maybe derived elemet did not allocate _connect");
	if (_di<0||_nd<0||_nl<0) throw new Fatal("Rod::Connect: __Internal Error__: There is a problem with _di=%d, _nd=%d, or _nd=%d\n (_di=dimension index, _gi=geometry index, _nd=number of degrees of freedom, _nl=number of additional labels)",_di,_nd,_nl);

	// Connectivity
	_ge->Conn[iNod] = ptNode;

	if (IsActive)
	{
		// Add Degree of Freedom to a node (Essential, Natural)
		for (int i=0; i<_nd; ++i) _ge->Conn[iNod]->AddDOF (UD[_di][i], FD[_di][i]);

		// Set shared
		_ge->Conn[iNod]->SetSharedBy (ID);
	}
}


inline void Rod::UpdateState(double TimeInc, Vec_t const & dUglobal, Vec_t & dFint)
{
	// Allocate (local/element) displacements vector
	Vec_t du(_nd*_ge->NNodes); // Delta disp. of this element

	// Assemble (local/element) displacements vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		du(i*_nd+j) = dUglobal(_ge->Conn[i]->DOFVar(UD[_di][j]).EqID);

	// Allocate (local/element) internal force vector
	Vec_t df(_nd*_ge->NNodes); // Delta internal force of this element
	df.SetValues(0.0);

	Mat_t Ke;
	CMatrix(0,Ke);
	df = Ke * du;

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		dFint(_ge->Conn[i]->DOFVar(UD[_di][j]).EqID) += df(i*_nd+j);
}

inline void Rod::AddVolForces() 
{
	// Verify if element is active
	if (IsActive==false) return;

	// Weight
	double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
	double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
	double L  = sqrt(dx*dx+dy*dy);
	double W  = _A*L*_gam;

	// Set boundary conditions
	if (_ge->NDim==1) throw new Fatal("Rod::ApplyBodyForces: feature not available for NDim==1");
	else if (_ge->NDim==2)
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

inline void Rod::CalcDepVars() const
{
	// Element displacements vector
	_uL.Resize(_nd*_ge->NNodes);
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (int    j=0; j<_nd;      ++j)
		_uL(i*_nd+j) = _ge->Conn[i]->DOFVar(UD[_di][j]).EssentialVal;

	// Transform to rod-local coordinates
	Mat_t T;
	_transf_mat(T);
	_uL = T * _uL;
}

inline double Rod::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_di][j])==0) return _ge->Conn[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_di][j])==0) return _ge->Conn[iNodeLocal]->DOFVar(Name).NaturalVal;

	if (_uL.Size()<1) throw new Fatal("Rod::Val: Please, call CalcDepVars() before calling this method");
	double l = (iNodeLocal==0 ? 0 : 1.0);
	     if (strcmp(Name,"N" )==0) return N(l);
	else if (strcmp(Name,"Ea")==0) return    (_uL(_nd)-_uL(0))/_L;
	else if (strcmp(Name,"Sa")==0) return _E*(_uL(_nd)-_uL(0))/_L;
	else throw new Fatal("Rod::Val: This element does not have a Val named %s",Name);
}

inline double Rod::Val(char const * Name) const
{
	throw new Fatal("Rod::Val: Feature not available");
}

inline void Rod::CMatrix(size_t Index, Mat_t & Ke) const
{
	if (_ge->NDim==2)
	{
		double dx = _ge->Conn[1]->X()-_ge->Conn[0]->X();
		double dy = _ge->Conn[1]->Y()-_ge->Conn[0]->Y();
		double LL = dx*dx+dy*dy;
		      _L  = sqrt(LL);
		double c  = dx/_L;
		double s  = dy/_L;
		double c1 = _E*_A*c*c/_L;
		double c2 = _E*_A*s*c/_L;
		double c3 = _E*_A*s*s/_L;
		Ke.Resize(_nd*_ge->NNodes, _nd*_ge->NNodes);
		Ke =   c1,  c2, -c1, -c2,
		       c2,  c3, -c2, -c3,
		      -c1, -c2,  c1,  c2,
		      -c2, -c3,  c2,  c3;
	}
	else throw new Fatal("Rod::CMatrix: Feature not available for nDim==%d",_ge->NDim);
}

inline void Rod::B_Matrix(Mat_t const & derivs, Mat_t const & J, Mat_t & B) const
{
	throw new Fatal("Rod::B_Matrix: Feature not available");
}

inline double Rod::N(double l) const
{
	return _E*_A*(_uL(_nd)-_uL(0))/_L;
}


/* private */

inline void Rod::_initialize()
{
	if (_ge->NDim<1) throw new Fatal("Rod::_initialize: For this element, _ge->NDim must be greater than or equal to 1 (%d is invalid)",_ge->NDim);
	_di  = _ge->NDim-2;
	_nd = _ge->NDim;
	_nl = 3;
}

inline void Rod::_calc_initial_internal_state()
{
	throw new Fatal("Rod::_calc_initial_internal_state: Feature not available");
}

inline void Rod::_transf_mat(Mat_t & T) const
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
		T.Resize(4,4);
		T =    c,   s, 0.0, 0.0,
		      -s,   c, 0.0, 0.0,
		     0.0, 0.0,   c,   s,
		     0.0, 0.0,  -s,   c;
	}
	else throw new Fatal("Rod::_transf_mat: Feature not available for nDim==%d",_ge->NDim);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
ProbElem * RodMaker()
{
	return new Rod();
}

// Register element
int RodRegister()
{
	ProbElemFactory["Rod"] = RodMaker;
	return 0;
}

// Call register
int __Rod_dummy_int  = RodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_ROD_H
