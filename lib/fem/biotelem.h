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

#ifndef MECHSYS_FEM_BIOTELEM_H
#define MECHSYS_FEM_BIOTELEM_H

// MechSys
#include "fem/element.h"
#include "util/string.h"
#include "util/util.h"
#include "linalg/laexpr.h"
#include "tensors/tensors.h"
#include "tensors/functions.h"

using Util::SQ2;
using Tensors::Tensor2;
using Tensors::Tensor4;

namespace FEM
{

class BiotElem : public virtual Element
{
public:
	// Constants
	static const size_t NDE[3];        ///< Number of DOFs of Equilibrium only. NDF(Flow) = ND - NDE
	static const size_t ND [3];        ///< Number of DOFs to add to a node == ND[_ndim-1]
	static const char   UD [3][4][4];  ///< Essential DOF names == UD[_ndim-1][iDOF]
	static const char   FD [3][4][4];  ///< Natural DOF names
	static const size_t NL [3];        ///< Number of additional labels (exceeding ND)
	static const char   LB [3][22][4]; ///< Additional labels

	// Constructor
	BiotElem () : _gam(0.0), _d(-1), _nd(-1), _gw(-1) {}

	// Destructor
	virtual ~BiotElem() {}

	// Derived methods
	virtual bool CheckModel   () const;
	bool         IsEssential  (char const * Name) const;
	virtual void SetModel     (char const * ModelName, char const * Prms, char const * Inis);
	void         SetProps     (char const * Properties);
	Element    * Connect      (int iNodeLocal, FEM::Node * ptNode);
	virtual void UpdateState  (double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint);
	void         ApplyBodyForces ();
	void         GetLabels    (Array<String> & Labels) const;
	void         Deactivate   ();
	char const * ModelName    () const { return "LinElastic/LinFlow"; }

	// Derived methods to assemble DAS matrices
	size_t nOrder0Matrices () const { return 1; }                                                                                                              ///< Number of zero order matrices: H:Permeability.
	size_t nOrder1Matrices () const { return 3; }                                                                                                              ///< Number of first order matrices: K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix.
	size_t nOrder0Vectors  () const { return 2; }                                                                                                              ///< Number of zero order vectors: H*P and Qb.
	void   Order0MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const; ///< Order0Matrix' map to convert local DOFs into global equation positions.
	void   Order0VecMap    (size_t Index, Array<size_t> & RowsMap)                                                                                      const; ///< Order0Vector' map to convert local DOFs into global equation positions.
	void   Order0Matrix    (size_t Index, LinAlg::Matrix<double> & M)                                                                                   const; ///< Zero order matrix: H:Permeability.
	void   Order0Vector    (size_t Index, LinAlg::Vector<double> & V)                                                                                   const; ///< Zero order vector: [U P]^T, where U is displacement and P, porepressure.
	void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const; ///< Order0Matrix' map to convert local DOFs into global equation positions.
	void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & M)                                                                                   const; ///< First order matrix: K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix.

	// Methods
	virtual void B_Matrix  (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const;
	virtual void Bp_Matrix (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & Bp) const;

	// Access methods
	virtual void   CalcDepVars () const;                                  ///< Calculate dependent variables (to be called before Val() or OutNodes() for example). Necessary for output of principal stresses, for example.
	virtual double Val         (int iNodeLocal, char const * Name) const; ///< Return values at nodes
	virtual double Val         (                char const * Name) const; ///< Return values at the CG of the element

protected:
	// Data
	double           _gam;    ///< Specific weight
	int              _d;      ///< Dimension index == _ndim-1
	int              _nd;     ///< Number of DOFs == ND[_d]
	double           _gw;     ///< Water specific weight (gamma W)
	Matrix<double>   _De;     ///< Constant tangent stiffness
	Matrix<double>   _Ke;     ///< Constant tangent permeability
	Array<Vector<double> >   _stress; ///< Total stress at each integration point. Size==_n_int_pts
	Array<Vector<double> >   _strain; ///< Strain at each integration point. Size==_n_int_pts

	// Private methods that MUST be derived
	virtual int  _geom() const =0; ///< Geometry of the element: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)

private:
	void _dist_to_face_nodes  (char const * Key, double const FaceValue, Array<Node*> const & FaceConnects) const;
	void _compute_K           (LinAlg::Matrix<double> & Ke) const;
	void _compute_C           (LinAlg::Matrix<double> & Ce) const;
	void _compute_L           (LinAlg::Matrix<double> & Le) const;
	void _compute_H           (LinAlg::Matrix<double> & He) const;
	void _compute_Qb          (LinAlg::Vector<double> & Qb) const;
	void _compute_equilib_map (Array<size_t> & RowsMap, Array<bool> & RowsEssenPresc) const;
	void _compute_flow_map    (Array<size_t> & RowsMap, Array<bool> & RowsEssenPresc) const;

    double _val_ip (size_t iIP, char const * Name) const; ///< Output values at a specific Integration Point (IP)

}; // class BiotElem

// NDE[_ndim-1]                  1D  2D  3D
const size_t BiotElem::NDE[3] = { 1,  2,  3};

// UD[_ndim-1][iDOF]                  1D                   2D                    3D
const size_t BiotElem::ND[3]       = { 2,                   3,                    4};
const char   BiotElem::UD[3][4][4] = {{"ux","pwp","",""},  {"ux","uy","pwp",""},  {"ux","uy","uz","pwp"}};
const char   BiotElem::FD[3][4][4] = {{"fx","vol","",""},  {"fx","fy","vol",""},  {"fx","fy","fz","vol"}};

// LB[_geom-1][iLabel]
const size_t BiotElem::NL[3]        = { 3, 15, 22 };
const char   BiotElem::LB[3][22][4] = {
	{"Ea", "Sa", "Vx",  ""   , ""   , ""   , ""  , ""   , ""  , ""   , ""   , ""   , ""  , ""  , ""  , ""  , ""  , ""  , ""  , ""  , ""  , "" }, // 1D
	{"Ex", "Ey", "Ez",  "Exy", "Sx" , "Sy" , "Sz", "Sxy", "E1", "E2" , "S1" , "S2" , "Vx", "Vy", "H" , ""  , ""  , ""  , ""  , ""  , ""  , "" }, // 2D (plane-strain)
	{"Ex", "Ey", "Ez",  "Exy", "Eyz", "Ezx", "Sx", "Sy" , "Sz", "Sxy", "Syz", "Szx", "E1", "E2", "E3", "S1", "S2", "S3", "Vx", "Vy", "Vz", "H"}, // 3D
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

// Derived methods

inline bool BiotElem::CheckModel() const
{
	if (_gw<0.0 || _De.Rows()==0 || _Ke.Rows()==0) return false;
	return true;
}

inline bool BiotElem::IsEssential(char const * Name) const
{
	for (int i=0; i<_nd; ++i) if (strcmp(Name,UD[_d][i])==0) return true;
	return false;
}

inline void BiotElem::SetModel(char const * ModelName, char const * Prms, char const * Inis)
{
	// Check _ndim
	if (_ndim<1)               throw new Fatal("BiotElem::SetModel: The space dimension (SetDim) must be set before calling this method");
	if (CheckConnect()==false) throw new Fatal("BiotElem::SetModel: Connectivity is not correct. Connectivity MUST be set before calling this method");

	/* "gw=10 E=200 nu=0.2 k=1.0e-5" */
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions(names,values);

	// Set
	double E  = -1.0;
	double nu = -1.0;
	double k  = -1.0;
	for (size_t i=0; i<names.Size(); ++i)
	{
		     if (names[i]=="gw") _gw = values[i];
		else if (names[i]=="E" )  E  = values[i];
		else if (names[i]=="nu")  nu = values[i];
		else if (names[i]=="k" )  k  = values[i];
		else throw new Fatal("BiotElem::SetModel: Parameter name (%s) is invalid",names[i].CStr());
	}

	// Check
	if (_gw<0.0)            throw new Fatal("BiotElem::SetModel: GammaW must be provided (and positive). gw==%f is invalid",_gw);
	if (E<=0.0)             throw new Fatal("BiotElem::SetModel: Young modulus (E) must be provided (and positive). E==%f is invalid",E);
	if (nu<0.0 || nu>0.499) throw new Fatal("BiotElem::SetModel: Poisson ratio (nu) must be provided (and in the range: 0 <= nu < 0.5). nu==%f is invalid",nu);
	if (k<0.0)              throw new Fatal("BiotElem::SetModel: Isotropic permeability must be provided (and positive). k=%f is invalid",k);

	// Set stiffness
	double c  = (_geom()==5 ? E/(1.0-nu*nu)  : E/((1.0+nu)*(1.0-2.0*nu)) ); // plane-stress != (plane-strain=3D)
	double c1 = (_geom()==5 ? c*1.0          : c*(1.0-nu)                ); // plane-stress != (plane-strain=3D)
	double c2 = (_geom()==5 ? c*0.5*(1.0-nu) : c*(1.0-2.0*nu)/2.0        ); // plane-stress != (plane-strain=3D)
	double c3 = c*nu;
	Tensor4 De;
	De = c1     , c3     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	     c3     , c1     , c3     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	     c3     , c3     , c1     , 0.0*SQ2, 0.0*SQ2, 0.0*SQ2,
	     0.0*SQ2, 0.0*SQ2, 0.0*SQ2, c2 *2.0, 0.0*2.0, 0.0*2.0,
	     0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, c2 *2.0, 0.0*2.0,
	     0.0*SQ2, 0.0*SQ2, 0.0*SQ2, 0.0*2.0, 0.0*2.0, c2 *2.0; // In Mandel's basis
	Tensors::Tensor4ToMatrix (_geom(),De, _De);

	// Set permeability
	double kx = k;
	double ky = k;
	double kz = k;
	if (_ndim==1)
	{
		_Ke.Resize(1,1);
		_Ke = kx;
	}
	else if (_ndim==2)
	{
		_Ke.Resize(2,2);
		_Ke =  kx,  0.0,
		      0.0,   ky;
	}
	else if (_ndim==3)
	{
		_Ke.Resize(3,3);
		_Ke =  kx,  0.0,  0.0,
		      0.0,   ky,  0.0,
		      0.0,  0.0,   kz;
	}

	// Set arrays of stress/strain
	_stress.Resize(_n_int_pts);
	_strain.Resize(_n_int_pts);
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		_stress[i].Resize(6);
		_stress[i] = 0.0,0.0,0.0, 0.0,0.0,0.0;
		_strain[i].Resize(6);
		_strain[i] = 0.0,0.0,0.0, 0.0,0.0,0.0;
	}
}

inline void BiotElem::SetProps(char const * Properties)
{
	/* "gam=20.0 */
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

inline Element * BiotElem::Connect(int iNodeLocal, FEM::Node * ptNode)
{
	// Check
	if (_n_nodes<1)               throw new Fatal("BiotElem::Connect: __Internal Error__: There is a problem with the number of nodes: maybe derived elemet did not set _n_nodes");
	if (_connects.Size()<1)       throw new Fatal("BiotElem::Connect: __Internal Error__: There is a problem with connectivity array: maybe derived elemet did not allocate _connect");
	if (_ndim<0 || _d<0 || _nd<0) throw new Fatal("BiotElem::Connect: __Internal Error__: There is a problem with _ndim=%d, _d=%d, or _nd=%d\n (_ndim=space dimension, _d=dimension index==_ndim-1, and _nd=number of degrees of freedom)",_ndim,_d,_nd);

	// Connectivity
	_connects[iNodeLocal] = ptNode;

	// Add Degree of Freedom to a node (Essential, Natural)
	for (int i=0; i<_nd; ++i) _connects[iNodeLocal]->AddDOF (UD[_d][i], FD[_d][i]);

	// Set shared
	_connects[iNodeLocal]->SetSharedBy(_my_id);

	return this;
}

inline void BiotElem::UpdateState(double TimeInc, LinAlg::Vector<double> const & dUglobal, LinAlg::Vector<double> & dFint)
{
	// Allocate (local/element) displacements vector
	LinAlg::Vector<double> du(_ndim*_n_nodes); // Delta disp. of this element
	LinAlg::Vector<double> dp(      _n_nodes); // Delta pore-water pressure of this element
	LinAlg::Vector<double>  p(      _n_nodes); // Total pore-water pressure of this element

	// Assemble (local/element) displacements vector
	size_t nde =       NDE[_d]; // nDOFs Equilibrium
	size_t ndf = _nd - NDE[_d]; // nDOFs Flow
	for (size_t i=0; i<_n_nodes; ++i)
	{
		for (size_t j=0; j<nde; ++j) du(i*nde+j) = dUglobal(_connects[i]->DOFVar(UD[_d][j    ]).EqID);
		for (size_t j=0; j<ndf; ++j) dp(i*ndf+j) = dUglobal(_connects[i]->DOFVar(UD[_d][j+nde]).EqID);
		for (size_t j=0; j<ndf; ++j)  p(i*ndf+j) =          _connects[i]->DOFVar(UD[_d][j+nde]).EssentialVal;
	}

	// Calculate increments of internal forces/volumes
	LinAlg::Vector<double> df  (_ndim*_n_nodes); // Delta internal force of this element
	LinAlg::Vector<double> dvol(_n_nodes);       // Delta internal volumes of this element
	LinAlg::Matrix<double> Ke;                   // stiffness matrix
	LinAlg::Matrix<double> Ce;                   // coupling matrix 1
	LinAlg::Matrix<double> Le;                   // coupling matrix 2
	LinAlg::Matrix<double> He;                   // permeability matrix
	Order1Matrix (0,Ke);
	Order1Matrix (1,Ce);
	Order1Matrix (2,Le);
	Order0Matrix (0,He);
	df   = Ke*du + Ce*dp;
	dvol = Le*du + TimeInc*He*p;

	// Stress update
	LinAlg::Vector<double> shape;
	LinAlg::Matrix<double> derivs;
	LinAlg::Matrix<double> J;
	LinAlg::Matrix<double> B;
	LinAlg::Vector<double> deps(6);
	LinAlg::Vector<double> dsig(6);
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;

		// Derivatives and another variables
		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs,J, B);   // Calculate B matrix for i Integration Point
		
		LinAlg::Vector<double> deps4;
		LinAlg::Vector<double> dsig4;
		deps4 = B*du;
		dsig4 = _De*deps4;
		deps  = deps4(0), deps4(1), deps4(2), deps4(3), 0.0, 0.0;
		dsig  = dsig4(0), dsig4(1), dsig4(2), dsig4(3), 0.0, 0.0;
		_stress[i] += dsig;
		_strain[i] += deps;
	}

	// Sum up contribution to internal forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	{
		for (size_t j=0; j<nde; ++j) dFint(_connects[i]->DOFVar(UD[_d][j    ]).EqID) += df  (i*nde+j);
		for (size_t j=0; j<ndf; ++j) dFint(_connects[i]->DOFVar(UD[_d][j+nde]).EqID) += dvol(i*ndf+j);
	}
}

inline void BiotElem::ApplyBodyForces()
{
	// Allocate (local/element) external volume force vector
	LinAlg::Matrix<double> fvol(_n_nodes, _ndim);
	fvol.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Vector<double> shape;
	LinAlg::Matrix<double> derivs;
	LinAlg::Matrix<double> J;
	LinAlg::Vector<double> b;

	// Mounting the body force vector
	     if (_ndim==3) { b.Resize(3); b = 0.0, 0.0, -_gam; }
	else if (_ndim==2) { b.Resize(2); b = 0.0, -_gam; }
	else if (_ndim==1) { b.Resize(1); b = -_gam; }

	// Calculate local external volume force
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		Shape    (r,s,t, shape);   // Calculate shape functions for i IP
		Derivs   (r,s,t, derivs);  // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);      // Calculate J (Jacobian) matrix for i Integration Point

		fvol += shape*trn(b)*det(J)*w;
	}

	// Sum up contribution to external forces vector
	for (size_t i=0; i<_n_nodes; ++i)
	{
		              _connects[i]->Bry("fx",fvol(i,0));
		if (_ndim==2) _connects[i]->Bry("fy",fvol(i,1));
		if (_ndim==3) _connects[i]->Bry("fz",fvol(i,2));
	}
}

inline void BiotElem::GetLabels(Array<String> & Labels) const
{
	const int p  = _geom()-1;
	const int q  = NL[p];   // number of additional labels
	const int nl = 2*_nd+q; // total number of labels
	Labels.Resize(nl);
	size_t k = 0;
	for (int i=0; i<_nd; ++i)
	{
		Labels[k] = UD[_d][i];  k++;
		Labels[k] = FD[_d][i];  k++;
	}
	for (int i=0; i<q; ++i)
	{
		Labels[k] = LB[p][i];  k++;
	}
}

inline void BiotElem::CalcDepVars() const
{
}

inline double BiotElem::Val(int iNodeLocal, char const * Name) const
{
	// Displacements
	for (int j=0; j<_nd; ++j) if (strcmp(Name,UD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).EssentialVal;

	// Forces
	for (int j=0; j<_nd; ++j) if (strcmp(Name,FD[_d][j])==0) return _connects[iNodeLocal]->DOFVar(Name).NaturalVal;

	// Stress, strains, internal values, etc.
	LinAlg::Vector<double>    ip_values (_n_int_pts); // Vectors for extrapolation
	LinAlg::Vector<double> nodal_values (_n_nodes);

	// Get integration point values
	for (size_t i=0; i<_n_int_pts; i++) ip_values(i) = _val_ip(i,Name);

	Extrapolate (ip_values, nodal_values);
	return nodal_values (iNodeLocal);
}

inline double BiotElem::Val(char const * Name) const
{
	throw new Fatal("BiotElem::Val: Feature not implemented yet");
}

inline void BiotElem::Deactivate()
{
	throw new Fatal("BiotElem::Deactivate: Feature not implemented yet");
}


// Derived methods to assemble DAS matrices

inline void BiotElem::Order0MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	_compute_flow_map (RowsMap, RowsEssenPresc);
	ColsMap        = RowsMap;
	ColsEssenPresc = RowsEssenPresc;
}

inline void BiotElem::Order0VecMap(size_t Index, Array<size_t> & RowsMap) const
{
	Array<bool> rows_essen_presc;
	_compute_flow_map (RowsMap, rows_essen_presc);
}

inline void BiotElem::Order0Matrix(size_t index, LinAlg::Matrix<double> & M) const
{
	_compute_H(M);
}

inline void BiotElem::Order0Vector(size_t Index, LinAlg::Vector<double> & V) const
{
	if (Index==0)
	{
		V.Resize(_n_nodes);
		size_t nde =       NDE[_d]; // nDOFs Equilibrium
		size_t ndf = _nd - NDE[_d]; // nDOFs Flow
		for (size_t i=0; i<_n_nodes; ++i)
		for (size_t j=0; j<ndf;      ++j)
			V(i*ndf+j) = _connects[i]->DOFVar(UD[_d][j+nde]).EssentialVal;

		LinAlg::Matrix<double> H;
		_compute_H(H);
		V = H*V;
	}
	else if (Index==1)
	{
		_compute_Qb(V);
	}
}

inline void BiotElem::Order1MatMap(size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const
{
	if (Index==0)
	{
		_compute_equilib_map (RowsMap, RowsEssenPresc);
		ColsMap        = RowsMap;
		ColsEssenPresc = RowsEssenPresc;
	}
	else if (Index==1)
	{
		_compute_equilib_map (RowsMap, RowsEssenPresc);
		_compute_flow_map    (ColsMap, ColsEssenPresc);
	}
	else if (Index==2)
	{
		_compute_flow_map    (RowsMap, RowsEssenPresc);
		_compute_equilib_map (ColsMap, ColsEssenPresc);
	}
}

inline void BiotElem::Order1Matrix(size_t Index, LinAlg::Matrix<double> & M) const
{
	     if (Index==0) _compute_K(M);
	else if (Index==1) _compute_C(M);
	else if (Index==2) _compute_L(M);
}


// Methods

inline void BiotElem::B_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & B) const
{
	/* OBS.:
	 *          This B matrix considers Solid Mechanics sign convention of stress and strains
	 *          Ex.: Compressive stresses/strains are negative
	 *          The B Matrix returns strains in Mandel notation
	 *
	 *          Traction    => Positive
	 *          Compression => Negative
	 */

	// Cartesian derivatives
	LinAlg::Matrix<double> dN;
	dN = inv(J)*derivs;

	// geometry type: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	switch (_geom())
	{
		case 2: // 2D(plane-strain)
		{
			const int n_scomps = 4; // number of stress compoments
			B.Resize (n_scomps,_ndim*_n_nodes);
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				B(0,0+i*_ndim) =     dN(0,i);  B(0,1+i*_ndim) =         0.0;
				B(1,0+i*_ndim) =         0.0;  B(1,1+i*_ndim) =     dN(1,i);
				B(2,0+i*_ndim) =         0.0;  B(2,1+i*_ndim) =         0.0;
				B(3,0+i*_ndim) = dN(1,i)/SQ2;  B(3,1+i*_ndim) = dN(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 3: // 3D
		{
			const int n_scomps = 6; // number of stress compoments
			B.Resize (n_scomps,_ndim*_n_nodes);
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				B(0,0+i*_ndim) =     dN(0,i);  B(0,1+i*_ndim) =         0.0;  B(0,2+i*_ndim) =         0.0;
				B(1,0+i*_ndim) =         0.0;  B(1,1+i*_ndim) =     dN(1,i);  B(1,2+i*_ndim) =         0.0;
				B(2,0+i*_ndim) =         0.0;  B(2,1+i*_ndim) =         0.0;  B(2,2+i*_ndim) =     dN(2,i);
				B(3,0+i*_ndim) = dN(1,i)/SQ2;  B(3,1+i*_ndim) = dN(0,i)/SQ2;  B(3,2+i*_ndim) =         0.0; // SQ2 => Mandel representation
				B(4,0+i*_ndim) =         0.0;  B(4,1+i*_ndim) = dN(2,i)/SQ2;  B(4,2+i*_ndim) = dN(1,i)/SQ2; // SQ2 => Mandel representation
				B(5,0+i*_ndim) = dN(2,i)/SQ2;  B(5,1+i*_ndim) =         0.0;  B(5,2+i*_ndim) = dN(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 5: // 2D(plane-stress)
		{
			const int n_scomps = 3; // number of stress compoments
			B.Resize(n_scomps,_ndim*_n_nodes);
			for (size_t i=0; i<_n_nodes; ++i) // i row of B
			{
				B(0,0+i*_ndim) =      dN(0,i);   B(0,1+i*_ndim) =         0.0;
				B(1,0+i*_ndim) =          0.0;   B(1,1+i*_ndim) =     dN(1,i);
				B(2,0+i*_ndim) =  dN(1,i)/SQ2;   B(2,1+i*_ndim) = dN(0,i)/SQ2; // SQ2 => Mandel representation
			}
			return;
		}
		case 1: // 1D
		case 4: // 2D(axis-symmetric)
		default: throw new Fatal("BiotElem::B_Matrix: B_Matrix() method is not available for GeometryType==%d",_geom());
	}
}

inline void BiotElem::Bp_Matrix(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> const & J, LinAlg::Matrix<double> & Bp) const
{
	Bp = inv(J)*derivs;
}


/* private */

inline void BiotElem::_dist_to_face_nodes(char const * Key, double const FaceValue, Array<Node*> const & FaceConnects) const
{
	// Call parent Element method for face boundary condition
	if (strcmp(Key,"Q")!=0)
	{
		Element::_dist_to_face_nodes(Key, FaceValue, FaceConnects);
		return;
	}

	// Normal traction boundary condition
	LinAlg::Matrix<double> values;  values.Resize(_n_face_nodes, _ndim);  values.SetValues(0.0);
	LinAlg::Matrix<double> J;                         // Jacobian matrix. size = [1,2] x 3
	LinAlg::Vector<double> face_shape(_n_face_nodes); // Shape functions of a face/edge. size = _n_face_nodes
	LinAlg::Matrix<double> F;                         // Shape function matrix
	LinAlg::Vector<double> P;                         // Vector perpendicular to the face 
	for (size_t i=0; i<_n_face_int_pts; i++)
	{
		double r = _a_face_int_pts[i].r;
		double s = _a_face_int_pts[i].s;
		double w = _a_face_int_pts[i].w;
		FaceShape    (r, s, face_shape);
		F = trn(trn(face_shape)); // trick just to convert Vector face_shape to a col Matrix

		// Calculate perpendicular vector
		if (_ndim==3)
		{
			FaceJacobian (FaceConnects, r, s, J);
			LinAlg::Vector<double> V(3); V = J(0,0), J(0,1), J(0,2);
			LinAlg::Vector<double> W(3); W = J(1,0), J(1,1), J(1,2);
			P.Resize(3);
			P = V(1)*W(2) - V(2)*W(1),      // vectorial product
			    V(2)*W(0) - V(0)*W(2),
			    V(0)*W(1) - V(1)*W(0);
		}
		else
		{
			FaceJacobian (FaceConnects, r, J);
			P.Resize(2);
			P = J(0,1), -J(0,0);  
		}
		values += FaceValue*F*trn(P)*w;
	}

	// Set nodes Brys
	for (size_t i=0; i<_n_face_nodes; ++i)
	{
		              FaceConnects[i]->Bry("fx",values(i,0));
		              FaceConnects[i]->Bry("fy",values(i,1));
		if (_ndim==3) FaceConnects[i]->Bry("fz",values(i,2));
	}
}

inline void BiotElem::_compute_K(LinAlg::Matrix<double> & Ke) const
{
	/* Stiffness K:
	   ===========
	
	                 /    T
	        [Ke]  =  | [B]  * [D] * [B]  * dV
	                 /
	*/

	// Resize Ke
	Ke.Resize    (_ndim*_n_nodes, _ndim*_n_nodes);
	Ke.SetValues (0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs; // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;      // Jacobian matrix
	LinAlg::Matrix<double> B;      // strain-displacement matrix

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		// Derivatives and another variables
		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs,J, B);   // Calculate B matrix for i Integration Point

		// Calculate tangent stiffness
		Ke += trn(B)*_De*B*det(J)*w;
	}
}

inline void BiotElem::_compute_C(LinAlg::Matrix<double> & Ce) const
{
	/*  Coupling Matrix L1:  ( n_dim*_n_nodes x _n_nodes )
	    ==================================================
	
	                 /    T      T
	        [C]   =  | [B]  * {m} * {N}  * dV    // coupled saturated
	                 /    u            p
	
	    Note:   p => pore-pressure
	    =====   u => displacement
	
	               T
	            {m} = [ 1 1 1 0 0 0 ]
	*/
	
	// Resize Ce
	Ce.Resize    (_ndim*_n_nodes, _n_nodes); 
	Ce.SetValues (0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs; // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Vector<double> shape;  // size = _n_nodes
	LinAlg::Matrix<double> J;      // Jacobian matrix
	LinAlg::Matrix<double> B;      // strain-displacement matrix
	LinAlg::Vector<double> m;      // vector with ones in the first three positions ~ 6D representation of the 2nd order identity tensor I
	     if (_ndim==2) { m.Resize(4); m = 1.0,1.0,1.0, 0.0; }
	else if (_ndim==3) { m.Resize(6); m = 1.0,1.0,1.0, 0.0,0.0,0.0; }

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		// Derivatives and another variables
		Shape    (r,s,t, shape);  // Calculate Np vector (equal to shape functions vector)
		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		B_Matrix (derivs,J, B);   // Calculate B matrix for i Integration Point

		// Calculate tangent coupling matrix
		Ce += trn(B)*m*trn(shape)*det(J)*w;
	}
}

inline void BiotElem::_compute_L(LinAlg::Matrix<double> & Le) const
{
	Matrix<double> Ce;
	_compute_C (Ce);
	Le = trn(Ce);
}

inline void BiotElem::_compute_H(LinAlg::Matrix<double> & He) const
{
	/*  Permeability Matrix H:
	    ======================
	
	                   /    T
	         [H]   =  -| [Bp]  * [K] * [Bp]  * dV
	                   /
	
	    Note:   p => pore-pressure
	    =====   u => displacement
	
	
	               T
	            {m} = [ 1 1 1 0 0 0 ]
	*/
	
	// Resize He
	He.Resize    (_n_nodes, _n_nodes);
	He.SetValues (0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs; // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Vector<double> shape;  // size = _n_nodes
	LinAlg::Matrix<double> J;      // Jacobian matrix
	LinAlg::Matrix<double> Bp;     // pore-pressure - gradient matrix

	// Loop along integration points
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		double w = _a_int_pts[i].w;

		// Derivatives and another variables
		Shape     (r,s,t, shape);  // Calculate Np vector (equal to shape functions vector)
		Derivs    (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian  (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point
		Bp_Matrix (derivs,J, Bp);  // Calculate Bp matrix for i Integration Point

		// Calculate tangent stiffness
		He += -trn(Bp)*_Ke*Bp*det(J)*w/_gw;
	}
}

inline void BiotElem::_compute_Qb(LinAlg::Vector<double> & Qb) const // {{{
{
	//	
	//	 Permeability Matrix Qh:
	//	 ============================
	//       
	//                    1   /    T                   
	//         [Qb]   =  ---  | [B]  * [K] * {b}  * dV
	//                    gw  /    p            w  
	//           
	//   OBS.: [K] = [k]/gammaW
	
	// Resize Qb
	Qb.Resize(_n_nodes); 
	Qb.SetValues(0.0);

	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs; // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Vector<double> shape;  // size = _n_nodes
	LinAlg::Matrix<double> J;      // Jacobian matrix
	LinAlg::Matrix<double> Bp;     // pore-pressure - gradient matrix
	LinAlg::Vector<double> b;   // fluid mass vector
	     if (_ndim==2) { b.Resize(2); b = 0.0, _gw; }
	else if (_ndim==3) { b.Resize(3); b = 0.0, 0.0, _gw; }

	// Loop along integration points
	for (size_t i_ip=0; i_ip<_n_int_pts; ++i_ip)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double t = _a_int_pts[i_ip].t;
		double w = _a_int_pts[i_ip].w;
		Shape(r,s,t, shape);              // Calculate Np vector (equal to shape functions vector)
		Derivs(r,s,t, derivs);            // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian(derivs, J);              // Calculate J (Jacobian) matrix for i_ip Integration Point
		Bp_Matrix(derivs,J, Bp);          // Calculate Bp matrix for i_ip Integration Point

		// Calculate the body mass vector for flow
		Qb += trn(Bp)*_Ke*b*det(J)*w/_gw;
	}
} // }}}

inline void BiotElem::_compute_equilib_map(Array<size_t> & RowsMap, Array<bool> & RowsEssenPresc) const
{
	RowsMap       .Resize(_ndim*_n_nodes);
	RowsEssenPresc.Resize(_ndim*_n_nodes);
	int    p   = 0;             // position (idx) inside Ke matrix
	size_t nde = NDE[_d];       // nDOFs Equilibrium
	for (size_t i=0; i<_n_nodes; ++i)
	{
		for (size_t j=0; j<nde; ++j)
		{
			RowsMap        [p] = _connects[i]->DOFVar(UD[_d][j]).EqID;
			RowsEssenPresc [p] = _connects[i]->DOFVar(UD[_d][j]).IsEssenPresc;
			p++;
		}
	}
}

inline void BiotElem::_compute_flow_map(Array<size_t> & RowsMap, Array<bool> & RowsEssenPresc) const
{
	RowsMap       .Resize(_n_nodes);
	RowsEssenPresc.Resize(_n_nodes);
	int    p   = 0;             // position (idx) inside Ce matrix
	size_t nde = NDE[_d];       // nDOFs Equilibrium
	size_t ndf = _nd - NDE[_d]; // nDOFs Flow
	for (size_t i=0; i<_n_nodes; ++i)
	{
		for (size_t j=0; j<ndf; ++j)
		{
			RowsMap        [p] = _connects[i]->DOFVar(UD[_d][j+nde]).EqID;
			RowsEssenPresc [p] = _connects[i]->DOFVar(UD[_d][j+nde]).IsEssenPresc;
			p++;
		}
	}
}

inline double BiotElem::_val_ip(size_t iIP, char const * Name) const
{
	     if (strcmp(Name,"Sx" )==0)                          return _stress[iIP](0);
	else if (strcmp(Name,"Sy" )==0)                          return _stress[iIP](1);
	else if (strcmp(Name,"Sz" )==0)                          return _stress[iIP](2);
	else if (strcmp(Name,"Sxy")==0 || strcmp(Name,"Syx")==0) return _stress[iIP](3)/SQ2;
	else if (strcmp(Name,"Syz")==0 || strcmp(Name,"Szy")==0) return _stress[iIP](4)/SQ2;
	else if (strcmp(Name,"Szx")==0 || strcmp(Name,"Sxz")==0) return _stress[iIP](5)/SQ2;
	else if (strcmp(Name,"p"  )==0)                          return (_stress[iIP](0)+_stress[iIP](1)+_stress[iIP](2))/3.0;
	else if (strcmp(Name,"q"  )==0)                          return sqrt(((_stress[iIP](0)-_stress[iIP](1))*(_stress[iIP](0)-_stress[iIP](1)) + (_stress[iIP](1)-_stress[iIP](2))*(_stress[iIP](1)-_stress[iIP](2)) + (_stress[iIP](2)-_stress[iIP](0))*(_stress[iIP](2)-_stress[iIP](0)) + 3.0*(_stress[iIP](3)*_stress[iIP](3) + _stress[iIP](4)*_stress[iIP](4) + _stress[iIP](5)*_stress[iIP](5)))/2.0);
	else if (strcmp(Name,"Ex" )==0)                          return _strain[iIP](0);
	else if (strcmp(Name,"Ey" )==0)                          return _strain[iIP](1);
	else if (strcmp(Name,"Ez" )==0)                          return _strain[iIP](2);
	else if (strcmp(Name,"Exy")==0 || strcmp(Name,"Eyx")==0) return _strain[iIP](3)/SQ2;
	else if (strcmp(Name,"Eyz")==0 || strcmp(Name,"Ezy")==0) return _strain[iIP](4)/SQ2;
	else if (strcmp(Name,"Ezx")==0 || strcmp(Name,"Exz")==0) return _strain[iIP](5)/SQ2;
	else if (strcmp(Name,"Ev" )==0)                          return _strain[iIP](0)+_strain[iIP](1)+_strain[iIP](2); 
	else if (strcmp(Name,"Ed" )==0)                          return sqrt(2.0*((_strain[iIP](0)-_strain[iIP](1))*(_strain[iIP](0)-_strain[iIP](1)) + (_strain[iIP](1)-_strain[iIP](2))*(_strain[iIP](1)-_strain[iIP](2)) + (_strain[iIP](2)-_strain[iIP](0))*(_strain[iIP](2)-_strain[iIP](0)) + 3.0*(_strain[iIP](3)*_strain[iIP](3) + _strain[iIP](4)*_strain[iIP](4) + _strain[iIP](5)*_strain[iIP](5))))/3.0;

	// Principal components of stress
	else if (strcmp(Name,"S1" )==0) { double sigp[3];/* Tensors::Eigenvals(_stress[iIP], sigp);*/ return sigp[2]; }
	else if (strcmp(Name,"S2" )==0) { double sigp[3];/* Tensors::Eigenvals(_stress[iIP], sigp);*/ return sigp[1]; }
	else if (strcmp(Name,"S3" )==0) { double sigp[3];/* Tensors::Eigenvals(_stress[iIP], sigp);*/ return sigp[0]; }

	// Principal components of strain
	else if (strcmp(Name,"E1" )==0) { double epsp[3];/* Tensors::Eigenvals(_strain[iIP], epsp);*/ return epsp[2]; }
	else if (strcmp(Name,"E2" )==0) { double epsp[3];/* Tensors::Eigenvals(_strain[iIP], epsp);*/ return epsp[1]; }
	else if (strcmp(Name,"E3" )==0) { double epsp[3];/* Tensors::Eigenvals(_strain[iIP], epsp);*/ return epsp[0]; }

	// Flow velocities
	else if (strcmp(Name,"Vx" )==0) return 0.0;
	else if (strcmp(Name,"Vy" )==0) return 0.0;
	else if (strcmp(Name,"Vz" )==0) return 0.0;

	// Total head
	else if (strcmp(Name,"H" )==0) return 0.0;

	else throw new Fatal("BiotElem::_val_ip: Name==%s if not available for this element",Name);
}

}; // namespace FEM

#endif // MECHSYS_FEM_BIOTELEM_H
