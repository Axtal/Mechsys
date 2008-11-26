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


/* OBS.:
 *          For 2D meshes, FACE means EDGES
 *          and EDGES related arrays are unavailable.
 *
 *          For 3D meshes, FACE and EDGES correspond to the normal
 *          meanings.
 */


// STL
#include <map>
#include <cstdarg>  // for va_list, va_start, va_end

// MechSys
#include "fem/quadrature.h"
#include "fem/node.h"
#include "util/string.h"
#include "linalg/vector.h"
#include "linalg/matrix.h"
#include "linalg/lawrap.h"
#include "linalg/laexpr.h"

namespace FEM
{

/// Elements
class Element
{
public:
	// Default constructor
	Element() : _my_id(-1), _ndim(-1), _n_nodes(0), _n_face_nodes(0), _n_int_pts(0), _n_face_int_pts(0) {}

	// Destructor
	virtual ~Element() {}

	// Set methods
	void              SetID     (long ID)            { _my_id     = ID;       }    ///< Set the ID of this element
	void              SetDim    (int nDim)           { _set_ndim(nDim);       }    ///< Set the number of dimension of the problem
	void              SetActive (bool IsActive=true) { _is_active = IsActive; }    ///< Activate/deactivate the element
	virtual Element * FaceBry   (char const * Key, double Value, int FaceLocalID); ///< Set face boundary conditions (SetDim MUST be called first)
	virtual Element * EdgeBry   (char const * Key, double Value, int EdgeLocalID); ///< Set edge boundary conditions (SetDim MUST be called first)

	// Set methods
	virtual Element * EdgeBry (char const * Key, double Value0, double Value1, int EdgeLocalID) { return this; } ///< Set edge boundary conditions (SetDim MUST be called first)
	
	// Get methods
	bool         CheckConnect ()         const;                         ///< Check if connectivity is OK
	bool         Check        (String & Message) const;                 ///< Check if everything is OK and element is ready for simulations
	long         GetID        ()         const { return _my_id;       } ///< Return the ID of this element
	bool         IsActive     ()         const { return _is_active;   } ///< Check if this element is active
	size_t       NNodes       ()         const { return _n_nodes;     } ///< Return the number of nodes in this element
	Node       * Nod          (size_t i)       { return _connects[i]; } ///< Return a pointer to a node in the connects list (read/write)
	Node const * Nod          (size_t i) const { return _connects[i]; } ///< Return a pointer to a node in the connects list (read-only)
	double       Volume       ()         const;                         ///< Return the volume/area/length of the element
	bool         IsInside     (double x, double y, double z) const;     ///< Check if a node is inside the element

	// Methods that MUST be overriden by derived classes
	virtual bool         CheckModel  ()                                  const =0; ///< Check if constitutive models are OK
	virtual void         CalcDepVars ()                                  const =0; ///< Calculate dependent variables (to be called before Val() or OutNodes() for example). Necessary for output of principal stresses, for example.
	virtual double       Val         (int iNodeLocal, char const * Name) const =0; ///< Return computed values at the Nodes of the element. Ex.: Name="ux", "fx", "Sx", "Sxy", "Ex", etc.
	virtual double       Val         (                char const * Name) const =0; ///< Return computed values at the CG of the element. Ex.: Name="Sx", "Sxy", "Ex", etc.
	virtual char const * Name        ()                                  const =0; ///< Return the name/type of this element
	virtual char const * ModelName   ()                                  const =0; ///< Return the name of the model of the first IP of this element

	// Methods related to PROBLEM (pure virtual) that MUST be overriden by derived classes
	virtual bool      IsEssential (char const * Name) const =0;                                                           ///< Is the correspondent DOFName (Degree of Freedom, such as "Dux") essential (such displacements)?
	virtual void      SetModel    (char const * ModelName, char const * Prms, char const * Inis) =0;                      ///< (Re)allocate model with parameters and initial values
	virtual void      SetProps    (char const * Properties) =0;                                                           ///< Set element properties such as body forces, internal heat source, water pumping, etc.
	virtual Element * Connect     (int iNodeLocal, FEM::Node * ptNode) =0;                                                ///< Set connectivity, by linking the local node ID with the pointer to the connection node
	virtual void      UpdateState (double TimeInc, LinAlg::Vector<double> const & dU, LinAlg::Vector<double> & dFint) =0; ///< Update the internal state of this element for given dU and update the DOFs related to this element inside dFint (internal forces increment vector)
	virtual void      GetLabels   (Array<String> & Labels) const =0;                                                      ///< Get the labels of all values to be output

	// Methods related to GEOMETRY (pure virtual) that MAY be overriden by derived classes
	virtual void SetIntPoints (int NumGaussPoints1D)                                                { throw new Fatal("FEM::Element: SetIntPoints() is not available"); } ///< Set the number of integration points using 1D information. Must NOT be called after allocation of Models.
	virtual int  VTKCellType  () const                                                              { throw new Fatal("FEM::Element: VTKCellType() is not available"); }  ///< Return the VTK (Visualization Tool Kit) cell type; used for generation of vtk files
	virtual void VTKConnect   (String & Nodes) const                                                { throw new Fatal("FEM::Element: VTKConnect() is not available"); }   ///< Return the VTK list of connectivities with global nodes IDs
	virtual void GetFaceNodes (int FaceID, Array<Node*> & FaceConnects) const                       { throw new Fatal("FEM::Element: GetFaceNodes() is not available"); } ///< Return the connectivity of a face, given the local face ID
	virtual void Shape        (double r, double s, double t, LinAlg::Vector<double> & Shape) const  { throw new Fatal("FEM::Element: Shape() is not available"); }        ///< Shape functions
	virtual void Derivs       (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const { throw new Fatal("FEM::Element: Derivs() is not available"); }       ///< Derivatives
	virtual void FaceShape    (double r, double s, LinAlg::Vector<double> & FaceShape) const        { throw new Fatal("FEM::Element: FaceShape() is not available"); }    ///< Face shape functions
	virtual void FaceDerivs   (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const       { throw new Fatal("FEM::Element: FaceDerivs() is not available"); }   ///< Face derivatives

	// Methods that MAY be overriden by derived classes
	virtual void   InverseMap    (double x, double y, double z, double & r, double & s, double & t) const;                                   ///< From "global" coordinates, compute the natural (local) coordinates
	virtual void   Jacobian      (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> & J) const;                                  ///< Jacobian matrix
	virtual void   Jacobian      (double r, double s, double t, LinAlg::Matrix<double> & J) const;                                           ///< (alternative) method to compute the Jacobian matrix
	virtual void   FaceJacobian  (Array<FEM::Node*> const & FaceConnects, double const r, double const s, LinAlg::Matrix<double> & J) const; ///< Jacobian matrix of a face
	virtual void   FaceJacobian  (Array<FEM::Node*> const & FaceConnects, double const r, LinAlg::Matrix<double> & J) const;                 ///< Jacobian matrix of a edge
	virtual void   Coords        (LinAlg::Matrix<double> & coords) const;                                                                    ///< Return the coordinates of the nodes
	virtual void   LocalCoords   (LinAlg::Matrix<double> & coords) const {};                                                                 ///< Return the local coordinates of the nodes
	virtual void   OutNodes      (LinAlg::Matrix<double> & Values, Array<String> & Labels) const;                                            ///< Output values at nodes
	virtual bool   HasExtra      () const { return false; }                                                                                                         ///< Has extra output ?
	virtual void   OutExtra      (LinAlg::Matrix<double> & Coords, LinAlg::Vector<double> & Norm, LinAlg::Matrix<double> & Values, Array<String> & Labels) const {} ///< Extra output for elements
	virtual double BoundDistance (double r, double s, double t) const { return -1; };                                                        ///< TODO
	virtual void   Extrapolate   (LinAlg::Vector<double> & IPValues, LinAlg::Vector<double> & NodalValues) const;                            ///< Extrapolate values from integration points to nodes
	virtual bool   HasVolForces  () const { return false; }                                                                                  ///< TODO
	virtual void   AddVolForces  (LinAlg::Vector<double> & FVol) const {}                                                                    ///< TODO
	virtual void   BackupState   () {}                                                                                                       ///< Backup internal state
	virtual void   RestoreState  () {}                                                                                                       ///< Restore internal state from a previously backup state

	// Methods to assemble DAS matrices; MAY be overriden by derived classes
	virtual size_t nOrder0Matrices () const { return 0; }                                                                                                                ///< Number of zero order matrices such as H:Permeability.
	virtual size_t nOrder1Matrices () const { return 0; }                                                                                                                ///< Number of first order matrices such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix.
	virtual void   Order0MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const {} ///< Order0Matrix' map to convert local DOFs into global equation positions.
	virtual void   Order0VecMap    (size_t Index, Array<size_t> & RowsMap)                                                                                      const {} ///< Order0Vector' map to convert local DOFs into global equation positions.
	virtual void   Order0Matrix    (size_t Index, LinAlg::Matrix<double> & M)                                                                                   const {} ///< Zero order matrix such as H:Permeability.
	virtual void   Order0Vector    (size_t Index, LinAlg::Vector<double> & V)                                                                                   const {} ///< Zero order vector such as [U P]^T, where U is displacement and P, porepressure.
	virtual void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const {} ///< Order0Matrix' map to convert local DOFs into global equation positions.
	virtual void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & M)                                                                                   const {} ///< First order matrix such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix.

protected:
	// Data (may be accessed by derived classes)
	long               _my_id;          ///< The ID of this element
	int                _ndim;           ///< Number of dimensions of the problem
	size_t             _n_nodes;        ///< GEOMETRY: Number of nodes in the element
	size_t             _n_face_nodes;   ///< GEOMETRY: Number of nodes in a face
	size_t             _n_int_pts;      ///< GEOMETRY: Number of integration (Gauss) points
	size_t             _n_face_int_pts; ///< GEOMETRY: Number of integration points in a face
	Array<Node*>       _connects;       ///< GEOMETRY: Connectivity (pointers to nodes in this element). size=_n_nodes
	bool               _is_active;      ///< Flag for active/inactive condition
	IntegPoint const * _a_int_pts;      ///< Array of Integration Points
	IntegPoint const * _a_face_int_pts; ///< Array of Integration Points of Faces/Edges
	
	// Methods related to GEOMETRY (pure virtual) that MUST be overriden by derived classes
	virtual void _set_ndim(int nDim) =0;
	virtual void _dist_to_face_nodes (char const * Key, double Value, Array<Node*> const & FaceConnects) const; ///< Distribute value to face nodes. FaceConnects => In: Array of ptrs to face nodes. FaceValue => In: A value applied on a face to be converted to nodes

}; // class Element


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

// Set and get methods

inline Element * Element::EdgeBry(char const * Key, double Value, int EdgeLocalID)
{
	if (_ndim<3) // For 1D/2D meshes, edges correspond to faces
	{

		// Skip if key is "Qb" TODO: Move to EquilibElem
		if (strcmp(Key,"Qb")==0) return this;

		Array<Node*> fnodes;
		GetFaceNodes        (EdgeLocalID, fnodes);
		_dist_to_face_nodes (Key, Value, fnodes);
	}
	else
	{
		throw new Fatal("Element::EdgeBry: Method not yet implemented for 3D meshes.");
	}
	return this;
}

inline Element * Element::FaceBry(char const * Key, double Value, int FaceLocalID)
{
	if (_ndim==2) throw new Fatal("Element::FaceBry: This method must be called only for 3D meshes.");
	else
	{
		Array<Node*> fnodes;
		GetFaceNodes        (FaceLocalID, fnodes);
		_dist_to_face_nodes (Key, Value, fnodes);
	}
	return this;
}

inline bool Element::CheckConnect() const
{
	if (_connects.Size()!=_n_nodes) return false;
	for (size_t i=0; i<_n_nodes; ++i) if (_connects[i]==NULL) return false;
	return true;
}

inline bool Element::Check(String & Message) const
{
	if (CheckConnect()==false) { Message.Printf("\n  %s","CONNECTIVITY NOT SET");       return false; }
	if (CheckModel()==false)   { Message.Printf("\n  %s","CONSTITUTIVE MODEL NOT SET"); return false; }
	return true;
}

inline double Element::Volume() const
{
	// Allocate entities used for every integration point
	LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
	LinAlg::Matrix<double> J;       // Jacobian matrix

	// Loop along integration points
	double vol = 0.0;
	for (size_t i=0; i<_n_int_pts; ++i)
	{
		// Temporary Integration Points
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;

		// Jacobian
		Derivs   (r,s,t, derivs); // Calculate Derivatives of Shape functions w.r.t local coordinate system
		Jacobian (derivs, J);     // Calculate J (Jacobian) matrix for i Integration Point

		// Calculate internal force vector;
		vol += det(J);
	}
	return vol;
}

inline bool Element::IsInside(double x, double y, double z) const
{
	double tiny = 1e-4;
	double huge = 1e+20;
	double max, min;
	
	//fast search in X -----------------------------------------------------------------------
	max = -huge;
	for (size_t i=0; i < NNodes(); i++) if (_connects[i]->Coord(0) > max) max=_connects[i]->Coord(0);
	if ( x > max ) return false;
	min = +huge;
	for (size_t i=0; i < NNodes(); i++) if (_connects[i]->Coord(0) < min) min=_connects[i]->Coord(0);
	if ( x < min ) return false;

	//fast search in Y -----------------------------------------------------------------------
	max = -huge;
	for (size_t i=0; i < NNodes(); i++) if (_connects[i]->Coord(1) > max) max=_connects[i]->Coord(1);
	if ( y > max ) return false;
	min = +huge;
	for (size_t i=0; i < NNodes(); i++) if (_connects[i]->Coord(1) < min) min=_connects[i]->Coord(1);
	if ( y < min ) return false;
	
	//fast search in Z -----------------------------------------------------------------------
	max = -huge;
	for (size_t i=0; i < NNodes(); i++) if (_connects[i]->Coord(2) > max) max=_connects[i]->Coord(2);
	if ( z > max ) return false;
	min = +huge;
	for (size_t i=0; i < NNodes(); i++) if (_connects[i]->Coord(2) < min) min=_connects[i]->Coord(2);
	if ( z < min ) return false;
	
	double r, s, t;
	InverseMap(x,y,z,r,s,t);
	if (BoundDistance(r,s,t)>-tiny) return true;
	else return false;
}

// Methods that MAY be overriden by derived classes

inline void Element::InverseMap(double x, double y, double z, double & r, double & s, double & t) const
{
	LinAlg::Vector<double> shape;
	LinAlg::Matrix<double> derivs;
	LinAlg::Matrix<double> J; //Jacobian matrix
	LinAlg::Vector<double> f(3);
	LinAlg::Vector<double> delta(3);
	double tx, ty, tz; //x, y, z trial
	double norm_f;
	r = s = t =0; // first suposition for natural coordinates
	int k=0;
	do
	{
		k++;
		Shape (r, s, t, shape);
		Derivs(r, s, t, derivs);
		Jacobian(derivs, J);
		tx = ty = tz = 0; 
		//calculate trial of real coordinates
		for (size_t j=0; j<_n_nodes; j++) 
		{
			tx += shape(j)*_connects[j]->Coord(0); //ok
			ty += shape(j)*_connects[j]->Coord(1); //ok
			tz += shape(j)*_connects[j]->Coord(2); //ok
		}
		
		// calculate the error
		f(0) = tx - x;
		f(1) = ty - y;
		f(2) = tz - z;
		
		delta = trn(inv(J))*f;
		
		r -= delta(0);
		s -= delta(1);
		t -= delta(2);

		norm_f = sqrt(trn(f)*f);
		if (k>25) break;
	} while(norm_f>1e-4);
}

inline void Element::Jacobian(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> & J) const
{
	// Calculate a matrix with nodal coordinates
	LinAlg::Matrix<double> cmatrix;  // size = _n_nodes x _ndim
	cmatrix.Resize(_n_nodes, _ndim);
	for (size_t i=0; i<_n_nodes; i++)
	for (int j=0; j<_ndim; ++j)
		cmatrix(i,j) = _connects[i]->Coord(j);

	// Calculate the Jacobian; 
	J = derivs * cmatrix;
}

inline void Element::Jacobian(double r, double s, double t, LinAlg::Matrix<double> & J) const
{
	LinAlg::Matrix<double> derivs; Derivs(r, s, t, derivs);
	LinAlg::Matrix<double> coords; Coords(coords);
	J = derivs*coords;
}

inline void Element::FaceJacobian(Array<FEM::Node*> const & FaceConnects, double const r, double const s, LinAlg::Matrix<double> & J) const
{
	if (_n_face_nodes>0)
	{
		// Calculate the shape function derivatives for a face
		LinAlg::Matrix<double> m_face_derivs;
		FaceDerivs(r, s, m_face_derivs);

		// Get the face coordinates
		LinAlg::Matrix<double> m_face_coords(_n_face_nodes,3);
		for (size_t i=0; i<_n_face_nodes; i++)
		for (int j=0; j<_ndim; ++j)
			m_face_coords(i,j) = FaceConnects[i]->Coord(j);

		// Determine face jacobian (2x3)
		J = m_face_derivs*m_face_coords;
	}
}

inline void Element::FaceJacobian(Array<FEM::Node*> const & FaceConnects, double const r, LinAlg::Matrix<double> & J) const
{
	if (_n_face_nodes>0)
	{
		// Calculate the shape function derivatives for a face
		LinAlg::Matrix<double> m_face_derivs;
		FaceDerivs(r, 0.0, m_face_derivs);

		// Get the face coordinates
		LinAlg::Matrix<double> m_face_coords(_n_face_nodes,2);
		for (size_t i=0; i<_n_face_nodes; i++)
		{
			m_face_coords(i,0) = FaceConnects[i]->Coord(0);
			m_face_coords(i,1) = FaceConnects[i]->Coord(1);
		}

		// Determine face jacobian (1x2)
		J = m_face_derivs*m_face_coords;
	}
}

inline void Element::Coords(LinAlg::Matrix<double> & coords) const
{
	// Calculate a matrix with nodal coordinates
	coords.Resize(_n_nodes,_ndim);

	// Loop along all nodes of this element
	for (size_t i=0; i<_n_nodes; i++)
	for (int j=0; j<_ndim; ++j)
		coords(i,j) = _connects[i]->Coord(j);
}

inline void Element::OutNodes(LinAlg::Matrix<double> & Values, Array<String> & Labels) const
{
	// Get the labels of all values to output from derived elements
	GetLabels (Labels);

	// Resize matrix with values at nodes
	size_t nlabels = Labels.Size();
	Values.Resize (_n_nodes,nlabels);

	// Fill matrix with values
	for (size_t i=0; i<_n_nodes; ++i)
	for (size_t j=0; j< nlabels; ++j)
		Values(i,j) = Val(i, Labels[j].CStr());
}

inline void Element::Extrapolate(LinAlg::Vector<double> & IPValues, LinAlg::Vector<double> & NodalValues) const 
{
	// Check
	if (IPValues.Size()!=static_cast<int>(_n_int_pts)) throw new Fatal("Element::Extrapolate: IPValues.Size()==%d must be equal to _n_int_pts==%d",IPValues.Size(),_n_int_pts);
	if (NodalValues.Size()!=static_cast<int>(_n_nodes)) throw new Fatal("Element::Extrapolate: NodalValues.Size()==%d must be equal to _n_nodes==%d",NodalValues.Size(),_n_nodes);

	// Data
	size_t m = IPValues   .Size();
	size_t n = NodalValues.Size();

	/* N: shape functions matrix:
	          _                                                 _
	         |   N11      N12      N13      ...  N1(nNode)       |
	         |   N21      N22      N23      ...  N2(nNode)       |
	     N = |   N31      N32      N33      ...  N3(nNode)       |
	         |          ...............................          |
	         |_  N(nIP)1  N(nIP)2  N(nIP)3  ...  N(nIP)(nNode)  _| [m=nIP x n=nNode]
	*/
	LinAlg::Matrix<double> N(m, n);
	LinAlg::Vector<double> shape(n);
	LinAlg::Matrix<double> IP_coord(m, _ndim+1);   // IP coordinates matrix (local)
	LinAlg::Matrix<double> node_coord;             // nodal coordinates matrix (local)

	// Mount N and IP_coord matrices
	for (size_t i=0; i<m; i++)
	{
		double r = _a_int_pts[i].r;
		double s = _a_int_pts[i].s;
		double t = _a_int_pts[i].t;
		Shape (r, s, t, shape);
		for (size_t j=0; j<n; j++) N(i,j) = shape(j);
		IP_coord(i, 0) = r;
		IP_coord(i, 1) = s;
		if (_ndim==2)
			IP_coord(i, 2) = 1.0;
		else
		{
			IP_coord(i, 2) = t;
			IP_coord(i, 3) = 1.0;
		}
	}

	// Mount node_coord matrix
	LocalCoords(node_coord);

	// Extrapolator matrix
	LinAlg::Matrix<double> E(n, m);

	//calculate extrapolator matrix
	if (n==m)     E = inv(N);
	else if (m>n) E = invg(N);
	else
	{
		LinAlg::Matrix<double> I;
		identity(m,I);
		E = invg(N)* (I - IP_coord*invg(IP_coord)) + node_coord*invg(IP_coord);
	}
	
	// Extrapolate
	NodalValues = E * IPValues;
}


/* private */

inline void Element::_dist_to_face_nodes(char const * Key, double const FaceValue, Array<Node*> const & FaceConnects) const
{
	if (IsEssential(Key)) // Assign directly
	{
		// Set nodes Brys
		for (size_t i=0; i<_n_face_nodes; ++i)
			FaceConnects[i]->Bry(Key,FaceValue);
	}
	else // Integrate along area/length
	{
		// Compute face nodal values (integration along the face)
		LinAlg::Vector<double> values;  values.Resize(_n_face_nodes);  values.SetValues(0.0);
		LinAlg::Matrix<double> J;                         // Jacobian matrix. size = [1,2] x 3
		LinAlg::Vector<double> face_shape(_n_face_nodes); // Shape functions of a face/edge. size = _n_face_nodes
		for (size_t i=0; i<_n_face_int_pts; i++)
		{
			double r = _a_face_int_pts[i].r;
			double s = _a_face_int_pts[i].s;
			double w = _a_face_int_pts[i].w;
			FaceShape    (r, s, face_shape);
			FaceJacobian (FaceConnects, r, s, J);
			values += FaceValue*face_shape*det(J)*w;
		}

		// Set nodes Brys
		for (size_t i=0; i<_n_face_nodes; ++i)
			FaceConnects[i]->Bry(Key,values(i));
	}
}


/** Outputs an element. */
std::ostream & operator<< (std::ostream & os, FEM::Element const & E)
{
	os << "[" << E.GetID() << "] " << E.Name() << " " << E.ModelName() << "\n";
	for (size_t i=0; i<E.NNodes(); ++i)
		if (E.Nod(i)!=NULL) os << "   " << (*E.Nod(i)) << "\n";
	return os;
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new element
typedef Element * (*ElementMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes elements
typedef std::map<String, ElementMakerPtr, std::less<String> > ElementFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes elements
ElementFactory_t ElementFactory;

// Allocate a new element according to a string giving the name of the element
Element * AllocElement(char const * ElementName)
{
	ElementMakerPtr ptr=NULL;
	ptr = ElementFactory[ElementName];
	if (ptr==NULL)
		throw new Fatal(_("FEM::AllocElement: There is no < %s > implemented in this library"), ElementName);

	return (*ptr)();
}

}; // namespace FEM


#ifdef USE_BOOST_PYTHON
// {

namespace BPy = boost::python;
class PyElem
{
public:
	                     PyElem   (FEM::Element * ptElem) : _elem(ptElem)                               { }
	long                 GetID    ()                                                              const { return _elem->GetID();   }
	size_t               NNodes   ()                                                              const { return _elem->NNodes();  }
	FEM::Node const    & Nod      (size_t i)                                                      const { return (*_elem->Nod(i)); }
	PyElem             & Connect  (int iNodeLocal, FEM::Node & refNode)                                 { _elem->Connect  (iNodeLocal, &refNode); return (*this); }
	PyElem             & SetModel (BPy::str const & Name, BPy::str const & Prms, BPy::str const & Inis) { _elem->SetModel (BPy::extract<char const *>(Name)(), BPy::extract<char const *>(Prms)(), BPy::extract<char const *>(Inis)()); return (*this); }
	PyElem             & EdgeBry  (BPy::str const & Key, double Value, int EdgeLocalID)                 { _elem->EdgeBry  (BPy::extract<char const *>(Key)(), Value, EdgeLocalID); return (*this); }
	PyElem             & FaceBry  (BPy::str const & Key, double Value, int FaceLocalID)                 { _elem->FaceBry  (BPy::extract<char const *>(Key)(), Value, FaceLocalID); return (*this); }
	double               Val1     (int iNodeLocal, BPy::str const & Name)                               { return _elem->Val(iNodeLocal, BPy::extract<char const *>(Name)()); }
	double               Val2     (                BPy::str const & Name)                               { return _elem->Val(            BPy::extract<char const *>(Name)()); }
	FEM::Element const * GetElem  ()                                                              const { return _elem; }
	void CalcDepsVars() const { _elem->CalcDepVars(); }
	bool HasExtra () const { return _elem->HasExtra(); }
	void OutExtra (BPy::dict & Coords, BPy::list & Norm, BPy::dict & Values) const
	{
		Matrix<double> coords;
		Vector<double> norm;
		Matrix<double> values;
		Array<String>  labels;
		_elem->OutExtra (coords, norm, values, labels);

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
	void SetProps (BPy::list & Props)
	{
		int nprops = BPy::len(Props);
		Array<double> props(nprops);
		for (int i=0; i<nprops; ++i) props[i] = BPy::extract<double>(Props[i])();
		_elem->SetProps (props);
	}
private:
	FEM::Element * _elem;
}; // class PyElement

std::ostream & operator<< (std::ostream & os, PyElem const & E) { if (E.GetElem()!=NULL) os<<(*E.GetElem()); return os; }

// }
#endif // USE_BOOST_PYTHON


#endif // MECHSYS_FEM_ELEMENT
