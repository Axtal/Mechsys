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
	typedef const char VarName_t[4];

	// Constructor
	ProbElem() : IsActive(true), _nd(-1), _nl(-1), _gam(0.0) {}

	// Methods related to PROBLEM
	virtual void    AddVolForces ()                                           {}
	virtual void    ClearDisp    ()                                           {}
	virtual void    SetActive    (bool Activate, int ID)                     =0;
	virtual void    CalcDeps     () const                                     {}
	virtual Str_t   ModelName    () const                                    =0;
	virtual double  Val          (int iNod, Str_t Key) const                 =0;
	virtual double  Val          (          Str_t Key) const                 =0;
	virtual void    SetModel     (Str_t ModelName, Str_t Prms, Str_t Inis)   =0;
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

	// Methods
	        void Initialize (GeomElem * GE, bool IsAct); ///< Initialize the element
	virtual bool CheckModel () const =0;                 ///< Check constitutive model

	// Methods related to PROBLEM implemented here
	bool  IsEssen  (Str_t Key) const;
	void  SetProps (Str_t Properties);
	void  SetConn  (int iNod, FEM::Node * ptNode, int ID);
	void  GetLbls  (Array<String> & Lbls) const;
	void  EdgeBry  (Str_t Key, double Val, int iEdge);
	void  FaceBry  (Str_t Key, double Val, int iFace);

	// Public data (read only)
	bool IsActive;
	
protected:
	// Data
	int         _nd;  ///< Number of DOFs
	int         _nl;  ///< Number of lbls == NL[_gi]
	double      _gam; ///< Specific weight Ex.: [kN/m3]
	GeomElem  * _ge;  ///< Geometry element
	VarName_t * UD;   ///< Essential DOF names. Access: UD[iDOF]
	VarName_t * FD;   ///< Natural   DOF names. Access: FD[iDOF]
	VarName_t * LB;   ///< Additional lbls (exceed. those from UD/FD).  Access: LB[iLbl]

	// Methods
	virtual void _initialize         ()                                                     =0; ///< Initialize derived element
	virtual void _dist_to_face_nodes (Str_t Key, double FaceValue, Array<Node*> const & FConn); ///< Distribute values from face/edges to nodes

}; // class ProbElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void ProbElem::Initialize(GeomElem * GE, bool IsAct)
{
	_ge      = GE;
	IsActive = IsAct;
	_initialize ();
}

inline bool ProbElem::IsEssen(Str_t Name) const
{
	for (int i=0; i<_nd; ++i) if (strcmp(Name,UD[i])==0) return true;
	return false;
}

inline void ProbElem::SetProps(Str_t Properties)
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

inline void ProbElem::SetConn(int iNod, FEM::Node * ptNode, int ID)
{
	// Check
	if (_ge->NNodes<1)      throw new Fatal("ProbElem::SetConn: There is a problem with the number of nodes: maybe derived elemet did not set _ge->NNodes");
	if (_ge->Conn.Size()<1) throw new Fatal("ProbElem::SetConn: There is a problem with connectivity array: maybe derived elemet did not allocate Connecitvity array");
	if (_nd<0 || _nl<0)     throw new Fatal("ProbElem::SetConn: There is a problem with _nd=%d or _nl=%d\n (_nd=number of degrees of freedom or _nl=number of additional labels)",_nd,_nl);

	// Connectivity
	_ge->Conn[iNod] = ptNode;

	if (IsActive)
	{
		// Add Degree of Freedom to a node (Essential, Natural)
		for (int i=0; i<_nd; ++i) _ge->Conn[iNod]->AddDOF (UD[i], FD[i]);

		// Set shared
		_ge->Conn[iNod]->SetSharedBy (ID);
	}
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
