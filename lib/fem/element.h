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
//#include "fem/probelem.h"
#include "util/string.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "linalg/laexpr.h"

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
	Element() : _id(-1), _tag(0), _type("__no_type__"), _ge(NULL)/*, _pe(NULL)*/ {}

	// Methods
	void  Initialize (Str_t Type, long ID, int Tag, int NDim, bool IsActive); ///< Initialize the element
	long  GetID      () const                        { return _id;          } ///< Return the ID of this element
	int   Tag        () const                        { return _tag;         } ///< Return the Tag of this element
	Str_t Type       () const                        { return _type.CStr(); } ///< Return the name/type of this element
	bool  Check      (String & Msg) const;                                    ///< Check if everything is OK and element is ready for simulations

	// Methods related to GEOMETRY
	bool         CheckConn () const                                   { return _ge->CheckConn();             } ///< Check if connectivity is OK
	size_t       NNodes    () const                                   { return _ge->NNodes;                  } ///< Return the number of nodes in this element
	Node       * Nod       (size_t i)                                 { return _ge->Conn[i];                 } ///< Return a pointer to a node in the connects list (read/write)
	Node const * Nod       (size_t i) const                           { return _ge->Conn[i];                 } ///< Return a pointer to a node in the connects list (read-only)
	double       Volume    () const                                   { return _ge->Volume();                } ///< Return the volume/area/length of the element
	void         Extrap    (Vec_t & IPVals, Vec_t & NodVals) const    {        _ge->Extrap(IPVals,NodVals);  } ///< Extrapolate values from integration points to nodes
	void         InvMap    (double X, double Y, double Z,
	                        double & r, double & s, double & t) const {        _ge->InvMap(X,Y,Z,r,s,t);     } ///< From "global" coordinates, compute the natural (local) coordinates
	bool         IsInside  (double X, double Y, double Z) const       { return _ge->IsInside(X,Y,Z);         } ///< Check if a point is inside the element
	void         SetIPs    (int NIPs1D)                               {        _ge->SetIPs(NIPs1D);          } ///< Set the number of integration points using 1D information. Must NOT be called after allocation of Models.
	int          VTKType   () const                                   { return _ge->VTKType();               } ///< Return the VTK (Visualization Tool Kit) cell type; used for generation of vtk files
	void         VTKConn   (String & Nodes) const                     {        _ge->VTKConn(Nodes);          } ///< Return the VTK list of connectivities with global nodes IDs
	void         GetFNodes (int FaceID, Array<Node*> & FConn) const   {        _ge->GetFNodes(FaceID,FConn); } ///< Return the connectivity of a face, given the local face ID
	double       BoundDist (double r, double s, double t) const       { return _ge->BoundDist(r,s,t);        } ///< TODO: Bound distance

	// Methods related to PROBLEM
	void      AddVolForces ()                                           {/*        _pe->AddVolForces();                     */ } ///< Method to apply volumetric (body) forces as boundary condition
	void      ClearDisp    ()                                           {/*        _pe->ClearDisp();                        */ } ///< Clear displacements and strains (for equilibrium/coupled problems)
	void      SetActive    (bool Activate)                              {/*        _pe->SetActive(Activate);                */ } ///< Activate element (construction/excavation)
	bool      CheckMdl     () const                                     {return false;/* return _pe->CheckMdl(); }                       */  }///< Check if constitutive models are OK
	Element * EdgeBry      (Str_t Key, double Val, int iEdge)           {return this;/*        _pe->EdgeBry(Key,Val,iEdge); return this;*/ } ///< Set edge boundary conditions (Initialize MUST be called first)
	Element * EdgeBry      (Str_t Key, double V0, double V1, int iEdge) {return this;/*                                     return this;*/ } ///< Set edge boundary conditions (Initialize MUST be called first)
	Element * FaceBry      (Str_t Key, double Val, int iFace)           {return this;/*        _pe->FaceBry(Key,Val,iFace); return this;*/ } ///< Set face boundary conditions (Initialize MUST be called first)
	bool      IsActive     () const                                     {return false;/* return _pe->IsActive;                           */ } ///< Check if this element is active
	void      CalcDeps     () const                                     {/*        _pe->CalcDeps();                         */ } ///< Calculate dependent variables (to be called before Val() or OutNodes() for example). Necessary for output of principal stresses, for example.
	Str_t     ModelName    () const                                     {return "hi";/* return _pe->ModelName();                        */ } ///< Return the name of the model of the first IP of this element
	double    Val          (int iNod, Str_t Key) const                  {return 0;/* return _pe->Val(iNod,Key);                      */ } ///< Return computed values at the Nodes of the element. Ex.: Key="ux", "fx", "Sx", "Sxy", "Ex", etc.
	double    Val          (          Str_t Key) const                  {return 0;/* return _pe->Val(Key);                           */ } ///< Return computed values at the CG of the element. Ex.: Key="Sx", "Sxy", "Ex", etc.
	bool      IsEssen      (Str_t Key) const                            {return false;/*     return   _pe->IsEssen(name);                      */ } ///< Is the correspondent DOFKey (Degree of Freedom, such as "Dux") essential (such displacements)?
	void      SetProps     (Str_t Properties)                           {/*        _pe->SetProps(Properties);               */ } ///< Set element properties such as body forces, internal heat source, water pumping, etc.
	void      SetModel     (Str_t ModelName, Str_t Prms, Str_t Inis)    {/*        _pe->SetModel(ModelName,Prms,Inis);      */ } ///< (Re)allocate model with parameters and initial values
	Element * SetConn      (int iNod, FEM::Node * ptNode)               {return this;/*        _pe->SetConn(iNod,ptNode);   return this;*/ } ///< Set connectivity, by linking the local node ID with the pointer to the connection node
	void      Update       (double h, Vec_t const & dU, Vec_t & dFint)  {/*        _pe->Update(h,dU,dFint);                 */ } ///< Update the internal state of this element for given dU and update the DOFs related to this element inside dFint (internal forces increment vector)
	void      Backup       ()                                           {/*        _pe->Backup();                           */ } ///< Backup internal state
	void      Restore      ()                                           {/*        _pe->Restore();                          */ } ///< Restore internal state from a previously backup state
	void      GetLbls      (Array<String> & Lbls) const                 {/*        _pe->GetLbls();                          */ } ///< Get the labels of all values to be output
	void      OutNodes     (Mat_t & Vals, Array<String> & Lbls) const   {/*        _pe->OutNodes(Vals,Lbls);                */ } ///< Output values at nodes
	void      OutInfo      (std::ostream & os) const                    {/*        _pe->OutInfo();                          */ } ///< Output extra info of the derived element
	bool      HasExtra     () const                                     {return false;/*        _pe->HasExtra();                         */ } ///< Has extra output ?
	void      OutExtra     (Mat_t & Coords, Vec_t & Norm,                /*                                                 */ 
	                        Mat_t & Vals, Array<String> & Lbls) const   {/*        _pe->OutExtra(Coords,Norm,Vals,Lbls);    */ } ///< Output extra information
	size_t   NCMats        () const                                     {return 0;/* return _pe->NCMats();                           */ } ///< Number of C matrices such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix
	size_t   NHMats        () const                                     {return 0;/* return _pe->NHMats();                           */ } ///< Number of H matrices such as H:Permeability
	size_t   NUVecs        () const                                     {return 0;/* return _pe->NUVecs();                           */ } ///< Number of U vectors such as U:Displacements, P:Pore-pressures
	void     CMatrix       (size_t Idx, Mat_t & M) const                {/*        _pe->CMatrix(Idx,M);                     */ } ///< C matrix such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix
	void     HMatrix       (size_t Idx, Mat_t & M) const                {/*        _pe->HMatrix(Idx,M);                     */ } ///< H matrix such as H:Permeability
	void     UVector       (size_t Idx, Vec_t & V) const                {/*        _pe->UVector(Idx,V);                     */ } ///< U vector such as U:Displacement, P:Pore-pressure
	void     CMatMap       (size_t Idx,
	                        Array<size_t> & RMap,
	                        Array<size_t> & CMap,
	                        Array<bool> & RUPresc,
	                        Array<bool> & CUPresc) const                {/*     _pe->CMatMap(Idx,RMap,CMap,RUPresc,CUPresc);*/ } ///< CMatrix map to convert local DOFs into global equation positions
	void     HMatMap       (size_t Idx,
	                        Array<size_t> & RMap,
	                        Array<size_t> & CMap,
	                        Array<bool> & RUPresc,
	                        Array<bool> & CUPresc) const                {/*     _pe->HMatMap(Idx,RMap,CMap,RUPresc,CUPresc);*/ } ///< HMatrix map to convert local DOFs into global equation positions
	void     UVecMap       (size_t Idx, Array<size_t> & RMap) const     {/*     _pe->UVecMap(Index,RMap);                   */ } ///< UVector map to convert local DOFs into global equation positions

private:
	// Data
	long       _id;   ///< The ID of this element
	int        _tag;  ///< The Tag of this element
	String     _type; ///< The name of this element
	GeomElem * _ge;   ///< Geometry element
	//ProbElem * _pe;   ///< Problem element

}; // class Element


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline void Element::Initialize(Str_t Type, long ID, int Tag, int NDim, bool IsActive)
{
	_id   = ID;
	_tag  = Tag;
	_type = Type;
	_ge->Initialize (NDim);
	//_pe->Initialize (_ge, NDim, IsActive);
}

inline bool Element::Check(String & Msg) const
{
	/*
	if (_ge->CheckConn()==false) { Msg.Printf("\n  %s","CONNECTIVITY NOT SET");       return false; }
	if (_pe->CheckModel()==false)   { Msg.Printf("\n  %s","CONSTITUTIVE MODEL NOT SET"); return false; }
	*/
	return true;
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

#endif // MECHSYS_FEM_ELEMENT
