/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo, Fernando Alonso                   *
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
#include "models/model.h"
#include "util/string.h"
#include "util/lineparser.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "linalg/laexpr.h"

namespace FEM
{

typedef LinAlg::Matrix<double>  Mat_t;
typedef LinAlg::Vector<double>  Vec_t;
typedef char const            * Str_t;
typedef std::map<String,double> Prop_t;       ///< Properties type
typedef const char              ProName_t[8]; ///< Properties names. Ex: "gam", "s", "cq", ...

typedef double(*Fun_t)(double);

class ProbElem
{
public:
	// Typedefs
	typedef const char VarName_t[4]; ///< Variables names. Ex: "ux", "uy", "Sx", ...
	typedef const char PrmName_t[8]; ///< Parameters names. Ex: "E", "nu", "lam", "kap", ...

	// Constructor
	ProbElem() : IsActive(true), _nd(-1), _nl(-1), _gi(-1), UD(NULL), FD(NULL), LB(NULL), _mdl(NULL), _srcfun(NULL) {}

	// Destructor
	virtual ~ProbElem() {}

	// Methods related to PROBLEM
	virtual int         InitCtes     (int nDim)        { return (nDim==3 ? 0 : 1); } ///< Return geometry index: 3D==>0, 2D==>1
	virtual int         NProps       () const                                    =0;
	virtual ProName_t * Props        () const                                    =0;
	virtual void        AddVolForces ()                                           {}
	virtual void        ClearDisp    ()                                           {}
	virtual void        SetActive    (bool Activate, int ID)                     =0;
	virtual void        CalcDeps     () const                                     {}
	virtual double      Val          (int iNod, Str_t Key) const                 =0;
	virtual double      Val          (          Str_t Key) const                 =0;
	virtual void        EdgeBry      (Str_t Key, double Val, int iEdge);
	virtual void        EdgeBry      (Str_t Key, double V0, double V1, int iEdge) {}
	virtual void        Update       (double h, Vec_t const & dU, Vec_t & dFint) =0;
	virtual void        Backup       ()                                          =0;
	virtual void        Restore      ()                                          =0;
	virtual void        OutInfo      (std::ostream & os) const                   =0;
	virtual bool        HasExtra     () const                      { return false; }
	virtual void        OutState     (std::ostream & os, bool OnlyCaption) const  {}
	virtual void        OutExtra     (Mat_t & Coords, Vec_t & Norm,
	                                  Mat_t & Vals, Array<String> & Lbls) const   {}
	virtual size_t      NCMats       () const                          { return 0; }
	virtual size_t      NHMats       () const                          { return 0; }
	virtual size_t      NUVecs       () const                          { return 0; }
	virtual void        CMatrix      (size_t Idx, Mat_t & M) const                {}
	virtual void        HMatrix      (size_t Idx, Mat_t & M) const                {}
	virtual void        UVector      (size_t Idx, Vec_t & V) const                {}
	virtual void        CMatMap      (size_t Idx,
	                                  Array<size_t> & RMap,
	                                  Array<size_t> & CMap,
	                                  Array<bool> & RUPresc,
	                                  Array<bool> & CUPresc) const                {}
	virtual void        HMatMap      (size_t Idx,
	                                  Array<size_t> & RMap,
	                                  Array<size_t> & CMap,
	                                  Array<bool> & RUPresc,
	                                  Array<bool> & CUPresc) const                {}
	virtual void        UVecMap      (size_t Idx, Array<size_t> & RMap) const     {}

	/* Initialize the element. */
	void Initialize (GeomElem * GE, int ID, Array<Node*> const & CONN, Model * Mdl, Str_t Inis, Prop_t * Prp, Fun_t SrcFun, bool IsAct); ///< Initialize the element

	// Methods related to PROBLEM implemented here
	Str_t MdlName  () const { return (_mdl==NULL ? "__no_model__" : _mdl->Name()); }
	bool  IsEssen  (Str_t Key) const;
	void  GetLbls  (Array<String> & Lbls) const;
	void  FaceBry  (Str_t Key, double Val, int iFace);

	// Methods
	double Prop (Str_t Key) const; ///< Value of a property

	// Public data (read only)
	bool IsActive;
	
protected:
	// Data (set by InitCtes)
	int         _nd;  ///< (set by InitCtes) Number of DOFs
	int         _nl;  ///< (set by InitCtes) Number of labels == NL[_gi]
	int         _gi;  ///< (set by InitCtes) Geometry index: 3D=0, PStrain=1, PStress=2, Axis=3. Others: 3D=0, 2D=1
	VarName_t * UD;   ///< (set by InitCtes) Essential DOF names. Access: UD[iDOF]
	VarName_t * FD;   ///< (set by InitCtes) Natural   DOF names. Access: FD[iDOF]
	VarName_t * LB;   ///< (set by InitCtes) Additional lbls (exceed. those from UD/FD).  Access: LB[iLbl]

	// Data
	GeomElem * _ge;     ///< Geometry element
	Model    * _mdl;    ///< Constitutive model
	Prop_t   * _prp;    ///< Properties
	Fun_t      _srcfun; ///< Source function

	// Methods
	virtual void _initialize         (Str_t Inis)                                           =0; ///< Initialize the element
	virtual void _dist_to_face_nodes (Str_t Key, double FaceValue, Array<Node*> const & FConn); ///< Distribute values from face/edges to nodes

}; // class ProbElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void ProbElem::Initialize(GeomElem * GE, int ID, Array<Node*> const & CONN, Model * Mdl, Str_t Inis, Prop_t * Prp, Fun_t SrcFun, bool IsAct)
{
	// Data
	_ge      = GE;
	_prp     = Prp;
	_srcfun  = SrcFun;
	IsActive = IsAct;

	// Check GeomElement
	if (_ge->NNodes<1)      throw new Fatal("ProbElem::Initialize: Element # %d: There is a problem with the number of nodes: maybe derived elemet did not set _ge->NNodes",ID);
	if (_ge->Conn.Size()<1) throw new Fatal("ProbElem::Initialize: Element # %d: There is a problem with connectivity array: maybe derived elemet did not allocate _ge->Conn (connecitvity) array",ID);

	// Check constants
	if (_nd<0 || _nl<0 || _gi<0)          throw new Fatal("ProbElem::Initialize: Element # %d: There is a problem with nd=%d or nl=%d\n (nd=number of degrees of freedom, nl=number of additional labels)",ID,_nd,_nl);
	if (UD==NULL || FD==NULL || LB==NULL) throw new Fatal("ProbElem::Initialize: Element # %d: There is a problem with UD, FD or LB\n (UD=Essential DOFs, FD=Corresponding Natural variables, LB=Labels)",ID);
	if (_prp==NULL)                       throw new Fatal("ProbElem::Initialize: Element # %d: There is a problem with _prp\n (_prp=properties)",ID);

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
	if (_ge->CheckConn()==false) throw new Fatal("ProbElem::Initialize: Element # %d: Connectivity was not properly set",ID);

	// Initialize
	_mdl = Mdl;
	_initialize (Inis);
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

inline double ProbElem::Prop(Str_t Key) const
{
	Prop_t::const_iterator it = _prp->find(Key);
	if (it==_prp->end()) throw new Fatal("ProbElem::Prop: Could not find property < %d > in _prp array",Key);
	return it->second;
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
