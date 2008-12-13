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
#include "util/lineparser.h"
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
	// Typedefs
	typedef const char VarName_t[4]; ///< Variables names. Ex: "ux", "uy", "Sx", ...
	typedef const char PrmName_t[8]; ///< Parameters names. Ex: "E", "nu", "lam", "kap", ...
	typedef const char ProName_t[8]; ///< Properties names. Ex: "gam", "s", "cq", ...

	// Constructor
	ProbElem() : IsActive(true), _nd(-1), _nl(-1), _gi(-1) {}

	// Destructor
	virtual ~ProbElem() {}

	// Methods related to PROBLEM
	virtual void    AddVolForces ()                                           {}
	virtual void    ClearDisp    ()                                           {}
	virtual void    SetActive    (bool Activate, int ID)                     =0;
	virtual void    CalcDeps     () const                                     {}
	virtual Str_t   ModelName    () const                                    =0;
	virtual double  Val          (int iNod, Str_t Key) const                 =0;
	virtual double  Val          (          Str_t Key) const                 =0;
	virtual void    EdgeBry      (Str_t Key, double V0, double V1, int iEdge) {}
	virtual void    Update       (double h, Vec_t const & dU, Vec_t & dFint) =0;
	virtual void    Backup       ()                                          =0;
	virtual void    Restore      ()                                          =0;
	virtual void    OutInfo      (std::ostream & os) const                   =0;
	virtual bool    HasExtra     () const                      { return false; }
	virtual void    OutExtra     (Mat_t & Coords, Vec_t & Norm,
	                              Mat_t & Vals, Array<String> & Lbls) const   {}
	virtual size_t  NCMats       () const                          { return 0; }
	virtual size_t  NHMats       () const                          { return 0; }
	virtual size_t  NUVecs       () const                          { return 0; }
	virtual void    CMatrix      (size_t Idx, Mat_t & M) const                {}
	virtual void    HMatrix      (size_t Idx, Mat_t & M) const                {}
	virtual void    UVector      (size_t Idx, Vec_t & V) const                {}
	virtual void    CMatMap      (size_t Idx,
	                              Array<size_t> & RMap,
	                              Array<size_t> & CMap,
	                              Array<bool> & RUPresc,
	                              Array<bool> & CUPresc) const                {}
	virtual void    HMatMap      (size_t Idx,
	                              Array<size_t> & RMap,
	                              Array<size_t> & CMap,
	                              Array<bool> & RUPresc,
	                              Array<bool> & CUPresc) const                {}
	virtual void    UVecMap      (size_t Idx, Array<size_t> & RMap) const     {}

	/* Initialize the element. */
	void Initialize (GeomElem * GE, int ID, Array<Node*> const & CONN, Str_t Model, Str_t Prms, Str_t Inis, Str_t Props, bool IsAct); ///< Initialize the element

	// Methods
	virtual bool CheckModel () const =0; ///< Check constitutive model

	// Methods related to PROBLEM implemented here
	bool  IsEssen  (Str_t Key) const;
	void  GetLbls  (Array<String> & Lbls) const;
	void  EdgeBry  (Str_t Key, double Val, int iEdge);
	void  FaceBry  (Str_t Key, double Val, int iFace);

	// Public data (read only)
	bool IsActive;
	
protected:
	// Data (set in _initialize)
	int         _nd;  ///< (set in _initialize) Number of DOFs
	int         _nl;  ///< (set in _initialize) Number of labels == NL[_gi]
	int         _gi;  ///< (set in _initialize) Geometry index: 3D=0, PStrain=1, PStress=2, Axis=3. Others: 3D=0, 2D=1
	VarName_t * UD;   ///< (set in _initialize) Essential DOF names. Access: UD[iDOF]
	VarName_t * FD;   ///< (set in _initialize) Natural   DOF names. Access: FD[iDOF]
	VarName_t * LB;   ///< (set in _initialize) Additional lbls (exceed. those from UD/FD).  Access: LB[iLbl]
	PrmName_t * PRMS; ///< (set in _initialize) Parameters names. Ex: "E", "nu", "lam", "kap", ...
	ProName_t * PROP; ///< (set in _initialize) Properties names. Ex: "gam", "s", "cq", ...

	// Data
	GeomElem      * _ge;    ///< Geometry element
	double          _gam;   ///< Specific weight Ex.: [kN/m3]
	Array<double>   _prms;  ///< Parameters
	Array<double>   _props; ///< Properties

	// Methods
	virtual void _initialize         (Str_t Model)                                          =0; ///< Initialize derived element. Set a pointer to Prms and return _nprms and _nprops
	virtual void _init_model         (Str_t Inis)                                            {} ///< Initialize model for already set _prms and given, i.e., Sx=1.0, Sy=2.0, etc.
	virtual void _dist_to_face_nodes (Str_t Key, double FaceValue, Array<Node*> const & FConn); ///< Distribute values from face/edges to nodes

}; // class ProbElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void ProbElem::Initialize(GeomElem * GE, int ID, Array<Node*> const & CONN, Str_t Model, Str_t Prms, Str_t Inis, Str_t Props, bool IsAct)
{
	// Data
	_ge      = GE;
	IsActive = IsAct;

	// Check GeomElement
	if (_ge->NNodes<1)      throw new Fatal("ProbElem::Initialize: Element # %d: There is a problem with the number of nodes: maybe derived elemet did not set _ge->NNodes",ID);
	if (_ge->Conn.Size()<1) throw new Fatal("ProbElem::Initialize: Element # %d: There is a problem with connectivity array: maybe derived elemet did not allocate _ge->Conn (connecitvity) array",ID);

	// Initialize: _nd, _nl, _gi, UD, FD, LB, PRMS and PROP. Also resize: _prms, _props
	_initialize (Model);

	// Check variables
	if (_nd<0 || _nl<0 || _prms.Size()<0 || _props.Size()<0) throw new Fatal("ProbElem::Initialize: Element # %d: There is a problem with nd=%d or nl=%d or nprms=%d or nprops=%d\n (nd=number of degrees of freedom, nl=number of additional labels, nprms=number of parameters, nprops=number of properties)",ID,_nd,_nl,_prms.Size(),_props.Size());

	// Set connectivity, by linking the local node ID with the pointer to the connected node
	for (size_t i=0; i<CONN.Size(); ++i)
	{
		_ge->Conn[i] = CONN[i];
		if (IsActive)
		{
			// Add Degree of Freedom to a node (Essential, Natural)
			for (int j=0; j<_nd; ++j) _ge->Conn[i]->AddDOF (UD[j], FD[j]);

			// Set shared
			_ge->Conn[i]->SetSharedBy (ID);
		}
	}

	// Check connectivity
	if (_ge->CheckConn()==false) throw new Fatal("EquilibElem::Initialize: Element # %d: Connectivity was not set properly",ID);

	// Read parameters
	LineParser lp(Prms);
	Array<String> names;
	Array<double> values;
	lp.BreakExpressions (names,values);
	_prms.SetValues(-1);
	for (size_t i=0; i<names.Size(); ++i)
	for (size_t j=0; j<_prms.Size(); ++j)
	{
		if (names[i]==PRMS[j]) // found parameter
		{
			_prms[j] = values[i];
			break;
		}
	}

	// Initialize model
	_init_model (Inis);

	// Read properties
	lp.Reset            (Props);
	lp.BreakExpressions (names,values);
	 _props.SetValues(0.0);
	for (size_t i=0; i<names.Size();  ++i)
	for (size_t j=0; j<_props.Size(); ++j)
	{
		if (names[i]==PROP[j]) _props[j] = values[i];
		else throw new Fatal("ProbElem::Initialize: Element # %d: Property name < %s > is invalid",ID,names[i].CStr());
	}
}

inline bool ProbElem::IsEssen(Str_t Name) const
{
	for (int i=0; i<_nd; ++i) if (strcmp(Name,UD[i])==0) return true;
	return false;
}

inline void ProbElem::GetLbls(Array<String> & Lbls) const
{
	const int nl = 2*_nd+_nl; // total number of labels
	Lbls.Resize (nl);
	size_t k = 0;
	for (int i=0; i<_nd; ++i)
	{
		Lbls[k] = UD[i];  k++;
		Lbls[k] = FD[i];  k++;
	}
	for (int i=0; i<_nl; ++i)
	{
		Lbls[k] = LB[i];  k++;
	}
}

inline void ProbElem::EdgeBry(Str_t Key, double Value, int iEdge)
{
	// Check
	if (_ge->NDim==3) throw new Fatal("ProbElem::EdgeBry: This method is not available for 3D meshes");

	// Skip if key is "Qb", Beam Normal Loading
	if (strcmp(Key,"Qb")==0) return;

	// For 2D meshes, edges correspond to faces
	Array<Node*> fnodes;
	_ge->GetFNodes      (iEdge, fnodes);
	_dist_to_face_nodes (Key, Value, fnodes);
}

inline void ProbElem::FaceBry(Str_t Key, double Value, int iFace)
{
	// Check
	if (_ge->NDim==2) throw new Fatal("ProbElem::FaceBry: This method is not available for 2D meshes and must be called only with 3D meshes");

	// Skip if key is "Qb", Beam Normal Loading
	if (strcmp(Key,"Qb")==0) return;

	// For 3D meshes, faces are faces ;-)
	Array<Node*> fnodes;
	_ge->GetFNodes      (iFace, fnodes);
	_dist_to_face_nodes (Key, Value, fnodes);
}


/* private */

inline void ProbElem::_dist_to_face_nodes(Str_t Key, double FaceValue, Array<Node*> const & FConn)
{
	// Check if the element is active
	if (IsActive==false) return;

	// Assing essential values directly and integrate for natural values
	if (IsEssen(Key)) for (size_t i=0; i<_ge->NFNodes; ++i) FConn[i]->Bry (Key,FaceValue);
	else
	{
		// Compute face nodal values (integration along the face)
		Vec_t values;  values.Resize(_ge->NFNodes);  values.SetValues(0.0);
		Mat_t J;
		Vec_t FN(_ge->NFNodes);
		for (size_t i=0; i<_ge->NFIPs; i++)
		{
			_ge->FaceShape (_ge->FIPs[i].r, _ge->FIPs[i].s, FN);
			_ge->FaceJacob (FConn, _ge->FIPs[i].r, _ge->FIPs[i].s, J);
			values += FaceValue*FN*det(J)*_ge->FIPs[i].w;
		}

		// Set nodes Brys
		for (size_t i=0; i<_ge->NFNodes; ++i) FConn[i]->Bry (Key,values(i));
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new ProbElem
typedef ProbElem * (*ProbElemMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes elements
typedef std::map<String, ProbElemMakerPtr, std::less<String> > ProbElemFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes elements
ProbElemFactory_t ProbElemFactory;

// Allocate a new ProbElem according to a string giving the name of the element
ProbElem * AllocProbElem(Str_t ProbElemName)
{
	ProbElemMakerPtr ptr=NULL;
	ptr = ProbElemFactory[ProbElemName];
	if (ptr==NULL) throw new Fatal("FEM::AllocProbElem: There is no < %s > implemented in this library", ProbElemName);
	return (*ptr)();
}

}; // namespace FEM

#endif // MECHSYS_FEM_PROBELEM
