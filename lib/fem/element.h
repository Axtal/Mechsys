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

#ifndef MECHSYS_FEM_ELEMENT_H
#define MECHSYS_FEM_ELEMENT_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "fem/probelem.h"
#include "util/string.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "linalg/laexpr.h"

#define CHECKGEPE(METHOD) if (_ge==NULL || _pe==NULL) throw new Fatal("Element::%s: Element must be Initialized before calling any of its methods",METHOD);

namespace FEM
{

typedef LinAlg::Matrix<double> Mat_t;
typedef LinAlg::Vector<double> Vec_t;
typedef char const *           Str_t;

class Element
{
public:
	// Friends
	friend std::ostream & operator<< (std::ostream & os, Element const & E);

	// Constructor
	Element() : _id(-1), _tag(0), _type("__no_type__"), _ge(NULL), _pe(NULL) {}

	// Destructor
	virtual ~Element() { if (_ge!=NULL) delete _ge;  if (_pe!=NULL) delete _pe; }

	/* Initialize the element. */
	void Initialize (size_t               i,       ///< ID == index in the array of elements
	                 int                  Tag,     ///< Tag of the element
	                 int                  nDim,    ///< Space dimension
	                 Array<Node*> const & CONN,    ///< Connectivity
	                 Str_t                GeomT,   ///< Geometry type. Ex: "Hex8"
	                 Str_t                ProbT,   ///< Problem type. Ex: "Equilib"
	                 Str_t                Model,   ///< Model name. Ex: "LinElastic"
	                 Str_t                Prms,    ///< Model parameters. Ex: "E=200 nu=0.2"
	                 Str_t                Inis,    ///< Initial values. Ex: "ZERO" or "Sx=0.0 Sy=0.0"
	                 Str_t                Props,   ///< Properties. Ex: "gam=20"
	                 bool                 IsAct);  ///< Is the element active ?

	// Methods
	long  GetID    () const            { return _id;          } ///< Return the ID of this element
	int   Tag      () const            { return _tag;         } ///< Return the Tag of this element
	Str_t Type     () const            { return _type.CStr(); } ///< Return the name/type of this element
	bool  Check    (String & Msg) const;                        ///< Check if everything is OK and element is ready for simulations
	void  OutNodes (Mat_t & Vals, Array<String> & Lbls) const;  ///< Output values at nodes

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
	double       BoundDist (double r, double s, double t) const         { CHECKGEPE("BoundDist") return _ge->BoundDist(r,s,t);                } ///< TODO: Bound distance

	// Methods related to PROBLEM
	void         AddVolForces ()                                        { CHECKGEPE("AddVols"  )        _pe->AddVolForces();                  } ///< Method to apply volumetric (body) forces as boundary condition
	void         ClearDisp ()                                           { CHECKGEPE("ClearDisp")        _pe->ClearDisp();                     } ///< Clear displacements and strains (for equilibrium/coupled problems)
	void         SetActive (bool Activate)                              { CHECKGEPE("SetActive")        _pe->SetActive(Activate,_id);         } ///< Activate element (construction/excavation)
	bool         IsActive  () const                                     { CHECKGEPE("IsActive" ) return _pe->IsActive;                        } ///< Check if this element is active
	void         CalcDeps  () const                                     { CHECKGEPE("CalcDeps" )        _pe->CalcDeps();                      } ///< Calculate dependent variables (to be called before Val() or OutNodes() for example). Necessary for output of principal stresses, for example.
	Str_t        ModelName () const                                     { CHECKGEPE("ModelName") return _pe->ModelName();                     } ///< Return the name of the model of the first IP of this element
	double       Val       (int iNod, Str_t Key) const                  { CHECKGEPE("Val"      ) return _pe->Val(iNod,Key);                   } ///< Return computed values at the Nodes of the element. Ex.: Key="ux", "fx", "Sx", "Sxy", "Ex", etc.
	double       Val       (          Str_t Key) const                  { CHECKGEPE("Val"      ) return _pe->Val(Key);                        } ///< Return computed values at the CG of the element. Ex.: Key="Sx", "Sxy", "Ex", etc.
	bool         IsEssen   (Str_t Key) const                            { CHECKGEPE("IsEssen"  ) return _pe->IsEssen(Key);                    } ///< Is the correspondent DOFKey (Degree of Freedom, such as "Dux") essential (such displacements)?
	void         EdgeBry   (Str_t Key, double Val, int iEdge)           { CHECKGEPE("EdgeBry"  )        _pe->EdgeBry(Key,Val,iEdge);          } ///< Set edge boundary conditions (Initialize MUST be called first)
	void         EdgeBry   (Str_t Key, double V0, double V1, int iEdge) { CHECKGEPE("EdgeBry"  )                                              } ///< Set edge boundary conditions (Initialize MUST be called first)
	void         FaceBry   (Str_t Key, double Val, int iFace)           { CHECKGEPE("FaceBry"  )        _pe->FaceBry(Key,Val,iFace);          } ///< Set face boundary conditions (Initialize MUST be called first)
	void         Update    (double h, Vec_t const & dU, Vec_t & dFint)  { CHECKGEPE("Update"   )        _pe->Update(h,dU,dFint);              } ///< Update the internal state of this element for given dU and update the DOFs related to this element inside dFint (internal forces increment vector)
	void         Backup    ()                                           { CHECKGEPE("Backup"   )        _pe->Backup();                        } ///< Backup internal state
	void         Restore   ()                                           { CHECKGEPE("Restore"  )        _pe->Restore();                       } ///< Restore internal state from a previously backup state
	void         GetLbls   (Array<String> & Lbls) const                 { CHECKGEPE("GetLbls"  )        _pe->GetLbls(Lbls);                   } ///< Get the labels of all values to be output
	void         OutInfo   (std::ostream & os) const                    { CHECKGEPE("OutInfo"  )        _pe->OutInfo(os);                     } ///< Output extra info of the derived element
	bool         HasExtra  () const                                     { CHECKGEPE("HasExtra" ) return _pe->HasExtra();                      } ///< Has extra output ?
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

inline void Element::Initialize(size_t i, int Tag, int nDim, Array<Node*> const & CONN, Str_t GeomT, Str_t ProbT, Str_t Model, Str_t Prms, Str_t Inis, Str_t Props, bool IsAct)
{
	// Data
	_id  = i;
	_tag = Tag;
	_type.Printf ("%s_%s_%dD",GeomT,ProbT,nDim);

	// Geometry element
	if (strcmp(GeomT,"")==0) _ge = AllocGeomElem ("Lin2");
	else                     _ge = AllocGeomElem (GeomT);
	_ge->Initialize (nDim);

	// Problem element
	if (strcmp(ProbT,"")==0) throw new Fatal("Element::Initialize: ProblemType must be provided");
	_pe = AllocProbElem (ProbT);
	_pe->Initialize     (_ge, _id, CONN, Model, Prms, Inis, Props, IsAct);
}

inline bool Element::Check(String & Msg) const
{
	if (_ge->CheckConn ()==false) { Msg.Printf("\n  %s","CONNECTIVITY NOT SET");       return false; }
	if (_pe->CheckModel()==false) { Msg.Printf("\n  %s","CONSTITUTIVE MODEL NOT SET"); return false; }
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

std::ostream & operator<< (std::ostream & os, FEM::Element const & E)
{
	//os << "[" << E.GetID() << "] Act=" << E.IsActive() << " Tag=" << E.Tag() << " " << E.Name() << " " << E.ModelName() << "\n";
	//os << "   ";  E.OutInfo(os);  os << "\n";
	//for (size_t i=0; i<E.NNodes(); ++i)
		//if (E.Nod(i)!=NULL) os << "   " << (*E.Nod(i)) << "\n";
	return os;
}

}; // namespace FEM

#undef CHECKGEPE

#endif // MECHSYS_FEM_ELEMENT
