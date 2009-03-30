/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#ifndef MECHSYS_FEM_ELEMENT_H
#define MECHSYS_FEM_ELEMENT_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "fem/probelem.h"
#include "fem/elems/lin2.h"
#include "util/string.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "linalg/laexpr.h"

#define CHECKGEPE(METHOD) if (_ge==NULL || _pe==NULL) throw new Fatal("Element::%s: Element must be Initialized before calling any of its methods",METHOD);

namespace FEM
{

typedef LinAlg::Matrix<double>  Mat_t;
typedef LinAlg::Vector<double>  Vec_t;
typedef char const            * Str_t;
typedef std::map<String,double> Prop_t;       ///< Properties type
typedef const char              ProName_t[8]; ///< Properties names. Ex: "gam", "s", "cq", ...

typedef double(*Fun_t)(double);

class Element
{
public:
	// Friends
	friend std::ostream & operator<< (std::ostream & os, Element const & E);

	// Constructor
	Element() : _id(-1), _tag(0), _type("__no_type__"), _ge(NULL), _pe(NULL) {}

	// Destructor
	virtual ~Element() { if (_ge!=NULL) delete _ge;  if (_pe!=NULL) delete _pe; }

	// Methods
	int   InitCtes   (int nDim, Str_t GeomT, Str_t ProbT);       ///< Initialize constants. Returns GI(Geometry Index): 3D=0, 2D=1, 2D-PStress=2, 2D-Axissym=3
	void  Initialize (size_t               i,                    ///< ID == index in the array of elements
	                  int                  Tag,                  ///< Tag of the element
	                  Array<Node*> const & CONN,                 ///< Connectivity
	                  Model              * Mdl,                  ///< Model. Ex: AllocModel("LinElastic")
	                  Str_t                Inis,                 ///< Initial values. Ex: "ZERO" or "Sx=0.0 Sy=0.0"
	                  Prop_t             * Prp,                  ///< Properties. Ex: "gam=20"
	                  Fun_t                SrcFun,               ///< Source function
	                  bool                 IsAct);               ///< Is the element active ?
	long  ID         () const           { return _id;          } ///< Return the ID of this element
	int   Tag        () const           { return _tag;         } ///< Return the Tag of this element
	Str_t Type       () const           { return _type.CStr(); } ///< Return the name/type of this element
	bool  Check      (String & Msg) const;                       ///< Check if everything is OK and element is ready for simulations
	void  OutNodes   (Mat_t & Vals, Array<String> & Lbls) const; ///< Output values at nodes
	void  Normal     (int iEdgeOrFace, Array<Vec_t> & X,
	                                   Array<Vec_t> & N) const;  ///< Output unit normal to edges/faces at the edge/face integration points

	// Methods related to GEOMETRY
	size_t       NNodes    () const                                     { CHECKGEPE("NNodes"   ) return _ge->NNodes;                          } ///< Return the number of nodes in this element
	Node       * Nod       (size_t i)                                   { CHECKGEPE("Nod"      ) return _ge->Conn[i];                         } ///< Return a pointer to a node in the connects list (read/write)
	Node const * Nod       (size_t i) const                             { CHECKGEPE("Nod"      ) return _ge->Conn[i];                         } ///< Return a pointer to a node in the connects list (read-only)
	double       Volume    () const                                     { CHECKGEPE("Volume"   ) return _ge->Volume();                        } ///< Return the volume/area/length of the element
	void         Extrap    (Vec_t & IPVals, Vec_t & NodVals) const      { CHECKGEPE("Extrap"   )        _ge->Extrap(IPVals,NodVals);          } ///< Extrapolate values from integration points to nodes
	void         InvMap    (double X, double Y, double Z,
	                        double & r, double & s, double & t) const   { CHECKGEPE("InvMap"   )        _ge->InvMap(X,Y,Z,r,s,t);             } ///< From "global" coordinates, compute the natural (local) coordinates
	bool         IsInside  (double X, double Y, double Z) const         { CHECKGEPE("IsInside" ) return _ge->IsInside(X,Y,Z);                 } ///< Check if a point is inside the element
	void         SetIPs    (int NIPs1D)                                 { CHECKGEPE("SetIPs"   )        _ge->SetIPs(NIPs1D);                  } ///< Set the number of integration points using 1D information. Must NOT be called after allocation of Models.
	int          VTKType   () const                                     { CHECKGEPE("VTKType"  ) return _ge->VTKType();                       } ///< Return the VTK (Visualization Tool Kit) cell type; used for generation of vtk files
	void         VTKConn   (String & Nodes) const                       { CHECKGEPE("VTKConn"  )        _ge->VTKConn(Nodes);                  } ///< Return the VTK list of connectivities with global nodes IDs
	void         GetFNodes (int FaceID, Array<Node*> & FConn) const     { CHECKGEPE("GetFNodes")        _ge->GetFNodes(FaceID,FConn);         } ///< Return the connectivity of a face, given the local face ID
	double       BoundDist (double r, double s, double t) const         { CHECKGEPE("BoundDist") return _ge->BoundDist(r,s,t);                } ///< Normalized distance from a specified point to the boundary fo the element

	// Methods related to PROBLEM
	int          NProps       () const                                  { CHECKGEPE("NProps"   ) return _pe->NProps();                        } ///< Number of properties
	ProName_t  * Props        () const                                  { CHECKGEPE("Props"    ) return _pe->Props();                         } ///< Names of properties
	Str_t        MdlName      () const                                  { CHECKGEPE("MdlName"  ) return _pe->MdlName();                       } ///< Name of the constitutive model
	void         AddVolForces ()                                        { CHECKGEPE("AddVols"  )        _pe->AddVolForces();                  } ///< Method to apply volumetric (body) forces as boundary condition
	void         ClearDisp ()                                           { CHECKGEPE("ClearDisp")        _pe->ClearDisp();                     } ///< Clear displacements and strains (for equilibrium/coupled problems)
	void         SetActive (bool Activate)                              { CHECKGEPE("SetActive")        _pe->SetActive(Activate,_id);         } ///< Activate element (construction/excavation)
	bool         IsActive  () const                                     { CHECKGEPE("IsActive" ) return _pe->IsActive;                        } ///< Check if this element is active
	void         CalcDeps  () const                                     { CHECKGEPE("CalcDeps" )        _pe->CalcDeps();                      } ///< Calculate dependent variables (to be called before Val() or OutNodes() for example). Necessary for output of principal stresses, for example.
	double       Val       (int iNod, Str_t Key) const                  { CHECKGEPE("Val"      ) return _pe->Val(iNod,Key);                   } ///< Return computed values at the Nodes of the element. Ex.: Key="ux", "fx", "Sx", "Sxy", "Ex", etc.
	double       Val       (          Str_t Key) const                  { CHECKGEPE("Val"      ) return _pe->Val(Key);                        } ///< Return computed values at the CG of the element. Ex.: Key="Sx", "Sxy", "Ex", etc.
	bool         IsEssen   (Str_t Key) const                            { CHECKGEPE("IsEssen"  ) return _pe->IsEssen(Key);                    } ///< Is the correspondent DOFKey (Degree of Freedom, such as "Dux") essential (such displacements)?
	void         EdgeBry   (Str_t Key, double Val, int iEdge)           { CHECKGEPE("EdgeBry"  )        _pe->EdgeBry(Key,Val,iEdge);          } ///< Set edge boundary conditions (Initialize MUST be called first)
	void         EdgeBry   (Str_t Key, double V0, double V1, int iEdge) { CHECKGEPE("EdgeBry"  )        _pe->EdgeBry(Key,V0,V1,iEdge);        } ///< Set edge boundary conditions (Initialize MUST be called first)
	void         FaceBry   (Str_t Key, double Val, int iFace)           { CHECKGEPE("FaceBry"  )        _pe->FaceBry(Key,Val,iFace);          } ///< Set face boundary conditions (Initialize MUST be called first)
	void         Update    (double t, double Dt, Vec_t const & dU, Vec_t & dFint)  { CHECKGEPE("Update"   ) _pe->Update(t,Dt,dU,dFint);       } ///< Update the internal state of this element for given dU and update the DOFs related to this element inside dFint (internal forces increment vector)
	void         Backup    ()                                           { CHECKGEPE("Backup"   )        _pe->Backup();                        } ///< Backup internal state
	void         Restore   ()                                           { CHECKGEPE("Restore"  )        _pe->Restore();                       } ///< Restore internal state from a previously backup state
	void         GetLbls   (Array<String> & Lbls) const                 { CHECKGEPE("GetLbls"  )        _pe->GetLbls(Lbls);                   } ///< Get the labels of all values to be output
	void         OutInfo   (std::ostream & os) const                    { CHECKGEPE("OutInfo"  )        _pe->OutInfo(os);                     } ///< Output extra info of the derived element
	bool         HasExtra  () const                                     { CHECKGEPE("HasExtra" ) return _pe->HasExtra();                      } ///< Has extra output ?
	void         OutState  (double Time, std::ostream & os, 
	                        bool OnlyCaption) const                     { CHECKGEPE("OutState")         _pe->OutState(Time,os,OnlyCaption);   } ///< Output element state. Ex.: Sx Sy Sz p q Ex Ey Ez Ev Ed
	void         OutExtra  (Mat_t & Coords, Vec_t & Norm,
	                        Mat_t & Vals, Array<String> & Lbls) const   { CHECKGEPE("OutExtra" )        _pe->OutExtra(Coords,Norm,Vals,Lbls); } ///< Output extra information
	size_t       NCMats    () const                                     { CHECKGEPE("NCMats"   ) return _pe->NCMats();                        } ///< Number of C matrices such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix
	size_t       NHMats    () const                                     { CHECKGEPE("NHMats"   ) return _pe->NHMats();                        } ///< Number of H matrices such as H:Permeability
	size_t       NUVecs    () const                                     { CHECKGEPE("NUVecs"   ) return _pe->NUVecs();                        } ///< Number of U vectors such as U:Displacements, P:Pore-pressures
	void         CMatrix   (size_t Idx, Mat_t & M) const                { CHECKGEPE("CMatrix"  )        _pe->CMatrix(Idx,M);                  } ///< C matrix such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix
	void         HMatrix   (size_t Idx, Mat_t & M) const                { CHECKGEPE("HMatrix"  )        _pe->HMatrix(Idx,M);                  } ///< H matrix such as H:Permeability
	void         UVector   (size_t Idx, Vec_t & V) const                { CHECKGEPE("UVector"  )        _pe->UVector(Idx,V);                  } ///< U vector such as U:Displacement, P:Pore-pressure
	void         CMatMap   (size_t Idx,
	                        Array<size_t> & RMap,
	                        Array<size_t> & CMap,
	                        Array<bool> & RUPresc,
	                        Array<bool> & CUPresc) const                { CHECKGEPE("CMatMap"  ) _pe->CMatMap(Idx,RMap,CMap,RUPresc,CUPresc); } ///< CMatrix map to convert local DOFs into global equation positions
	void         HMatMap   (size_t Idx,
	                        Array<size_t> & RMap,
	                        Array<size_t> & CMap,
	                        Array<bool> & RUPresc,
	                        Array<bool> & CUPresc) const                { CHECKGEPE("HMatMap"  ) _pe->HMatMap(Idx,RMap,CMap,RUPresc,CUPresc); } ///< HMatrix map to convert local DOFs into global equation positions
	void         UVecMap   (size_t Idx, Array<size_t> & RMap) const     { CHECKGEPE("UVecMap"  ) _pe->UVecMap(Idx,RMap);                      } ///< UVector map to convert local DOFs into global equation positions
	bool         HasDFs    () const                                     { CHECKGEPE("HasDFs")     return _pe->HasDFs();                       } ///< Has component (dF(s)tar) to add to external force vector?
	void         AddToDFext(double t, double Dt, Vec_t & DFext) const   { CHECKGEPE("AddToDFext")        _pe->AddToDFext(t,Dt,DFext);         } ///< Add extra component (dF(s)tar) to external force vector

#ifdef USE_BOOST_PYTHON
// {
	Node       & PyNod1     (size_t i)                                              { return (*Nod(i)); }
	Node const & PyNod2     (size_t i)                                        const { return (*Nod(i)); }
	double       PyVal1     (int iNod, BPy::str const & Key)                  const { return Val(iNod, BPy::extract<Str_t>(Key)()); }
	double       PyVal2     (          BPy::str const & Key)                  const { return Val(      BPy::extract<Str_t>(Key)()); }
	void         PyEdgeBry1 (BPy::str const & Key, double Val, int iEdge)           { EdgeBry(BPy::extract<Str_t>(Key)(), Val, iEdge); }
	void         PyEdgeBry2 (BPy::str const & Key, double V0, double V1, int iEdge) { EdgeBry(BPy::extract<Str_t>(Key)(), V0, V1, iEdge); }
	void         PyFaceBry  (BPy::str const & Key, double Val, int iFace)           { EdgeBry(BPy::extract<Str_t>(Key)(), Val, iFace); }
	void         PyOutExtra (BPy::dict & Coords, BPy::list & Norm, BPy::dict & Values) const;
// }
#endif // USE_BOOST_PYTHON

private:
	// Data
	long       _id;   ///< The ID of this element
	int        _tag;  ///< The Tag of this element
	String     _type; ///< The name of this element
	GeomElem * _ge;   ///< Geometry element
	ProbElem * _pe;   ///< Problem element

}; // class Element


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline int Element::InitCtes(int nDim, Str_t GeomT, Str_t ProbT)
{
	_type.Printf ("%s_%s_%dD",GeomT,ProbT,nDim);

	// Check
	if (strcmp(GeomT,"")==0) throw new Fatal("Element::InitCtes: GeometryType must be provided");
	if (strcmp(ProbT,"")==0) throw new Fatal("Element::InitCtes: ProblemType must be provided");
	if (strcmp(ProbT,"Rod")==0 || strcmp(ProbT,"Beam")==0 || strcmp(ProbT,"Spring")==0)
	{
		if (!(strcmp(GeomT,"Lin2")==0)) throw new Fatal("Element::InitCtes: For Rod, Beam, or Spring, GeometryType must be Lin2");
	}

	// Geometry element
	_ge = AllocGeomElem (GeomT);
	_ge->Initialize     (nDim);

	// Problem element
	_pe = AllocProbElem  (ProbT);
	return _pe->InitCtes (nDim);
}

inline void Element::Initialize(size_t i, int Tag, Array<Node*> const & CONN, Model * Mdl, Str_t Inis, Prop_t * Prp, Fun_t SrcFun, bool IsAct)
{
	_id  = i;
	_tag = Tag;
	_pe->Initialize (_ge, _id, CONN, Mdl, Inis, Prp, SrcFun, IsAct);
}

inline bool Element::Check(String & Msg) const
{
	if (_ge->CheckConn()==false) { Msg.Printf("\n  %s","CONNECTIVITY NOT SET"); return false; }
	return true;
}

inline void Element::OutNodes(Mat_t & Vals, Array<String> & Lbls) const
{
	// Get the labels of all values to output from derived elements
	_pe->GetLbls (Lbls);

	// Resize matrix with values at nodes
	size_t nlabels = Lbls.Size();
	Vals.Resize (_ge->NNodes, nlabels);

	// Fill matrix with values
	for (size_t i=0; i<_ge->NNodes; ++i)
	for (size_t j=0; j<nlabels;     ++j)
		Vals(i,j) = _pe->Val (i, Lbls[j].CStr());
}

inline void Element::Normal(int iEdgeOrFace, Array<Vec_t> & X, Array<Vec_t> & N) const
{
	CHECKGEPE("Normal")

	// Face nodes
	Array<Node*> fconn;
	_ge->GetFNodes (iEdgeOrFace, fconn);
	
	N.Resize(0);
	Mat_t J;
	Vec_t FN(_ge->NFNodes); // Face shape
	Vec_t n;                // Vector perpendicular to the face
	Vec_t x;                // Position of IP
	for (size_t i=0; i<_ge->NFIPs; i++)
	{
		// Calculate perpendicular vector
		if (_ge->NDim==3)
		{
			_ge->FaceJacob (fconn, _ge->FIPs[i].r, _ge->FIPs[i].s, J);
			Vec_t V(3); V = J(0,0), J(0,1), J(0,2);
			Vec_t W(3); W = J(1,0), J(1,1), J(1,2);
			n.Resize(3);
			x.Resize(3);
			n = V(1)*W(2) - V(2)*W(1),     
				V(2)*W(0) - V(0)*W(2),
				V(0)*W(1) - V(1)*W(0); // vectorial product
		}
		else // (_ge->NDim==2)
		{
			_ge->FaceJacob (fconn, _ge->FIPs[i].r, _ge->FIPs[i].s, J);
			n.Resize(2);
			x.Resize(2);
			n = J(0,1), -J(0,0);
		}
		n = n/LinAlg::Norm(n);
		N.Push(n);

		// Calculate position of IP
		_ge->FaceShape (_ge->FIPs[i].r, _ge->FIPs[i].s, FN);
		x.SetValues(0.0);
		for (size_t j=0; j<_ge->NFNodes; ++j)
		for (size_t k=0; k<_ge->NDim;    ++k)
			x(k) = x(k) + FN(j)*fconn[j]->Coord(k);
		X.Push(x);
	}
}


#ifdef USE_BOOST_PYTHON
// {
inline void Element::PyOutExtra (BPy::dict & Coords, BPy::list & Norm, BPy::dict & Values) const
{
	Mat_t          coords;
	Vec_t          norm;
	Mat_t          values;
	Array<String>  labels;
	OutExtra (coords, norm, values, labels);

	// Coords
	BPy::list X;
	BPy::list Y;
	for (int i=0; i<coords.Rows(); ++i)
	{
		X.append (coords(i,0));
		Y.append (coords(i,1));
	}
	Coords['X'] = X;
	Coords['Y'] = Y;

	// Norm
	Norm.append (norm(0));
	Norm.append (norm(1));

	// Values
	for (int j=0; j<values.Cols(); ++j)
	{
		BPy::list vals;
		for (int i=0; i<values.Rows(); ++i) vals.append (values(i,j));
		Values[labels[j].CStr()] = vals;
	}
}
// }
#endif // USE_BOOST_PYTHON

std::ostream & operator<< (std::ostream & os, FEM::Element const & E)
{
	os << "[" << E.ID() << "] Act=" << E.IsActive() << " Tag=" << E.Tag() << " " << E.Type() << " " << E.MdlName() << "\n";
	os << "   ";  E.OutInfo(os);  os << "\n";
	for (size_t i=0; i<E.NNodes(); ++i)
		if (E.Nod(i)!=NULL) os << "   " << (*E.Nod(i)) << "\n";
	return os;
}

}; // namespace FEM

#undef CHECKGEPE

#endif // MECHSYS_FEM_ELEMENT
