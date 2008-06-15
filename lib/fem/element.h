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

// STL
#include <map>

// MechSys
#include "util/string.h"
#include "fem/node.h"
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

	// Destructor
	virtual ~Element() { }

	// Auxiliar structure (integration point)
	struct IntegPoint
	{
		double r; ///< xsi coordinate
		double s; ///< eta coordinate
		double t; ///< zeta coordinate
		double w; ///< weight (coeficient) W 
	};

	// Methods
	void                 SetID       (int ID) { _my_id = ID;        }                     ///< Set the ID of this element
	int                  GetID       () const { return _my_id;      }                     ///< Return the ID of this element
	void                 Activate    ()       { _is_active = true;  }                     ///< Activate the element
	bool                 IsActive    () const { return _is_active;  }                     ///< Check if this element is active
	int                  nNodes      () const { return _n_nodes;    }                     ///< Return the number of nodes in this element
	size_t               nIntPoints  () const { return _n_int_pts;  }                     ///< Return the number of integration points in this element
	void                 IntegPoints (IntegPoint const * & IPs) const { IPs=_a_int_pts; } ///< Return a pointer to the array of integration points
	Array<Node*> const & Connects    () const;                                            ///< Return an array containing pointers to nodes
	bool                 IsInside    (double x, double y, double z) const;                ///< Check if a node is inside the element

	// Methods that MUST be overriden by derived classes
	virtual String Name() const =0;

	// Methods related to PROBLEM (pure virtual) that MUST be overriden by derived classes
	virtual bool      IsReady         () const=0;                                          ///< Check if element is ready for analysis
	virtual bool      IsEssential     (char const * DOFName) const =0;                                                                      ///< Is the correspondent DOFName (Degree of Freedom, such as "Dux") essential (such displacements)?
	virtual void      SetModel        (char const * ModelName, char const * Prms, char const * Inis) =0;                                    ///< (Re)allocate model with parameters and initial values
	virtual Element * SetNode         (int iNodeLocal, int iNodeGlobal) =0;                                                                 ///< TODO: Setup the DOFs of a node according to the DOFs needed by this element  ***** Copy a pointer of node iNode to the internal connects array
	virtual void      UpdateState     (double TimeInc, LinAlg::Vector<double> const & dU, LinAlg::Vector<double> & dFint) =0;               ///< Update the internal state of this element for given dU and update the DOFs related to this element inside dFint (internal forces increment vector)
	virtual void      BackupState     () =0;                                                                                                ///< Backup internal state
	virtual void      RestoreState    () =0;                                                                                                ///< Restore internal state from a previously backup state
	virtual void      SetGeometryType (int Geom) =0;                                                                                        ///< Set geometry type: 1:1D, 2:2D(plane-strain), 3:3D, 4:2D(axis-symmetric), 5:2D(plane-stress)
	virtual void      SetProperties   (Array<double> const & EleProps) =0;                                                                  ///< Set interal properties
	virtual String    OutCenter       (bool PrintCaptionOnly=false) const =0;                                                               ///< Output internal values computed (averaged) in the center of the element

	// Methods related to GEOMETRY (pure virtual) that MUST be overriden by derived classes
	virtual int  VTKCellType    () const =0;                                                                                                ///< Return the VTK (Visualization Tool Kit) cell type; used for generation of vtk files
	virtual void Shape          (double r, double s, double t, LinAlg::Vector<double> & Shape) const =0;                                    ///< Shape functions
	virtual void Derivs         (double r, double s, double t, LinAlg::Matrix<double> & Derivs) const =0;                                   ///< Derivatives
	virtual void FaceShape      (double r, double s, LinAlg::Vector<double> & FaceShape) const =0;                                          ///< Face shape functions
	virtual void FaceDerivs     (double r, double s, LinAlg::Matrix<double> & FaceDerivs) const =0;                                         ///< Face derivatives
	virtual void Dist2FaceNodes (Array<Node*> const & FaceConnects, double const FaceValue, LinAlg::Vector<double> & NodalValues) const =0; ///< FaceConnects => In: Array of ptrs to face nodes. FaceValue => In: A value applied on a face to be converted to nodes. NodalValues => Out:The resultant nodal values

	// Methods that MAY be overriden by derived classes
	virtual void   InverseMap    (double x, double y, double z, double & r, double & s, double & t) const;                                   ///< From "global" coordinates, compute the natural (local) coordinates
	virtual void   Jacobian      (LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> & J) const;                                  ///< Jacobian matrix
	virtual void   Jacobian      (double r, double s, double t, LinAlg::Matrix<double> & J) const;                                           ///< (alternative) method to compute the Jacobian matrix
	virtual void   FaceJacobian  (Array<FEM::Node*> const & FaceConnects, double const r, double const s, LinAlg::Matrix<double> & J) const; ///< Jacobian matrix of a face
	virtual void   FaceJacobian  (Array<FEM::Node*> const & FaceConnects, double const r, LinAlg::Matrix<double> & J) const;                 ///< Jacobian matrix of a edge
	virtual void   Coords        (LinAlg::Matrix<double> & coords) const;                                                                    ///< Return the coordinates of the nodes
	virtual void   Extrapolate   (LinAlg::Vector<double> & IPValues, LinAlg::Vector<double> & NodalValues) const;                            ///< Extrapolate values from integration points to nodes
	virtual void   OutNodes      (LinAlg::Matrix<double> & Values, Array<String> & Labels) const {};                                         ///< Output values at nodes
	virtual void   Deactivate    () { _is_active = false; }                                                                                  ///< Deactivate this element
	virtual double BoundDistance (double r, double s, double t) const { return -1; };                                                        ///< ???
	virtual void   FaceNodalVals (char const * FaceDOFName, double const FaceDOFValue, Array<FEM::Node*> const & APtrFaceNodes, String & NodalDOFName, LinAlg::Vector<double> & NodalValues) const {} ///< Return the array with nodal values in a face

	// Methods to assemble DAS matrices; MAY be overriden by derived classes
	virtual size_t nOrder0Matrices () const { return 0; }                                                                                                                ///< Number of zero order matrices such as H:Permeability.
	virtual size_t nOrder1Matrices () const { return 0; }                                                                                                                ///< Number of first order matrices such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix.
	virtual void   Order0MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const {} ///< Order0Matrix' map to convert local DOFs into global equation positions.
	virtual void   Order0VecMap    (size_t Index, Array<size_t> & RowsMap)                                                                                      const {} ///< Order0Vector' map to convert local DOFs into global equation positions.
	virtual void   Order0Matrix    (size_t Index, LinAlg::Matrix<double> & M)                                                                                   const {} ///< Zero order matrix such as H:Permeability.
	virtual void   Order0Vector    (size_t Index, LinAlg::Vector<double> & V)                                                                                   const {} ///< Zero order vector such as [U P]^T, where U is displacement and P, porepressure.
	virtual void   Order1MatMap    (size_t Index, Array<size_t> & RowsMap, Array<size_t> & ColsMap, Array<bool> & RowsEssenPresc, Array<bool> & ColsEssenPresc) const {} ///< Order0Matrix' map to convert local DOFs into global equation positions.
	virtual void   Order1Matrix    (size_t Index, LinAlg::Matrix<double> & M)                                                                                   const {} ///< First order matrix such as K:Stiffness, L1:CouplingMatrix1, L2:CouplingMatrix2 and M:MassMatrix.

	// Access methods
	virtual double Val (int iNodeLocal, char const * Name) const =0; ///< Return computed values at the CG of the element. Ex.: Name="Sx", "Sxy", "Ex", etc.
	
	// Methods
	double Volume () const
	{
		// Allocate entities used for every integration point
		LinAlg::Matrix<double> derivs;  // size = NumLocalCoords(ex.: r,s,t) x _n_nodes
		LinAlg::Matrix<double> J;       // Jacobian matrix

		// Loop along integration points
		double vol = 0.0;
		for (int i=0; i<_n_int_pts; ++i)
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

protected:
	// Data (may be accessed by derived classes)
	int                _my_id;          ///< The ID of this element
	int                _n_dim;          ///< Number of dimensions
	int                _geom;           ///< Geometry type
	int                _n_nodes;        ///< Number of nodes in the element
	int	               _n_int_pts;      ///< Number of integration (Gauss) points
	int                _n_face_nodes;   ///< Number of nodes in a face
	int	               _n_face_int_pts; ///< Number of integration points in a face
	Array<Node*>       _connects;       ///< Connectivity (pointers to nodes in this element). size=_n_nodes
	bool               _is_active;      ///< Flag for active/inactive condition
	IntegPoint const * _a_int_pts;      ///< Array of Integration Points

}; // class Element

Array<Element*> Elems; ///< Array with all elements (or only the elements in this processor)


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

inline Array<Node*> const & Element::Connects() const
{
	return _connects;
}

inline bool Element::IsInside(double x, double y, double z) const
{
	double tiny = 1e-4;
	double huge = 1e+20;
	double max, min;
	
	//fast search in X -----------------------------------------------------------------------
	max = -huge;
	for (int i=0; i < nNodes(); i++) if (_connects[i]->X() > max) max=_connects[i]->X();
	if ( x > max ) return false;
	min = +huge;
	for (int i=0; i < nNodes(); i++) if (_connects[i]->X() < min) min=_connects[i]->X();
	if ( x < min ) return false;

	//fast search in Y -----------------------------------------------------------------------
	max = -huge;
	for (int i=0; i < nNodes(); i++) if (_connects[i]->Y() > max) max=_connects[i]->Y();
	if ( y > max ) return false;
	min = +huge;
	for (int i=0; i < nNodes(); i++) if (_connects[i]->Y() < min) min=_connects[i]->Y();
	if ( y < min ) return false;
	
	//fast search in Z -----------------------------------------------------------------------
	max = -huge;
	for (int i=0; i < nNodes(); i++) if (_connects[i]->Z() > max) max=_connects[i]->Z();
	if ( z > max ) return false;
	min = +huge;
	for (int i=0; i < nNodes(); i++) if (_connects[i]->Z() < min) min=_connects[i]->Z();
	if ( z < min ) return false;
	
	double r, s, t;
	InverseMap(x,y,z,r,s,t);
	if (BoundDistance(r,s,t)>-tiny) return true;
	else return false;
}

inline void Element::InverseMap(double x, double y, double z, double & r, double & s, double & t) const
{
	LinAlg::Vector<double> shape;
	LinAlg::Matrix<double> derivs;
	LinAlg::Matrix<double> J; //Jacobian matrix
	//LinAlg::Matrix<double> inv_jacobian;
	//LinAlg::Matrix<double> trn_inv_jacobian;
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
		for (int j=0; j<_n_nodes; j++) 
		{
			tx += shape(j)*_connects[j]->X(); //ok
			ty += shape(j)*_connects[j]->Y(); //ok
			tz += shape(j)*_connects[j]->Z(); //ok
		}
		
		// calculate the error
		f(0) = tx - x;
		f(1) = ty - y;
		f(2) = tz - z;
		
		//jacobian.Inv(inv_jacobian);
		//inv_jacobian.Trn(trn_inv_jacobian);
		delta = trn(inv(J))*f;
		//Gemv(1.0, trn_inv_jacobian, f, 0.0, delta);
		
		r -= delta(0);
		s -= delta(1);
		t -= delta(2);

		//norm_f = sqrt(f(0)*f(0)+f(1)*f(1)+f(2)*f(2));
		norm_f = sqrt(trn(f)*f);
		if (k>25) break;
	} while(norm_f>1e-4);
}

inline void Element::Jacobian(LinAlg::Matrix<double> const & derivs, LinAlg::Matrix<double> & J) const
{
	// Calculate a matrix with nodal coordinates
	LinAlg::Matrix<double> cmatrix; // size = _n_nodes x 3

    if (derivs.Rows()==3)
	{
		cmatrix.Resize(_n_nodes,3);   // 3 => X,Y,Z (all nodes have X,Y and Z)

		// Loop along all nodes of this element
		for (int i=0; i<_n_nodes; i++)
		{
			cmatrix(i,0) = _connects[i]->X();
			cmatrix(i,1) = _connects[i]->Y();
			cmatrix(i,2) = _connects[i]->Z();
		}
	}
	else
	{
		cmatrix.Resize(_n_nodes,2);   // 2 => X,Y (all nodes have X and Z)

		// Loop along all nodes of this element
		for (int i=0; i<_n_nodes; i++)
		{
			cmatrix(i,0) = _connects[i]->X();
			cmatrix(i,1) = _connects[i]->Y();
		}
	}

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
		for (int i=0; i<_n_face_nodes; i++)
		{
			m_face_coords(i,0) = FaceConnects[i]->X();
			m_face_coords(i,1) = FaceConnects[i]->Y();
			m_face_coords(i,2) = FaceConnects[i]->Z();
		}

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
		for (int i=0; i<_n_face_nodes; i++)
		{
			m_face_coords(i,0) = FaceConnects[i]->X();
			m_face_coords(i,1) = FaceConnects[i]->Y();
		}

		// Determine face jacobian (1x2)
		J = m_face_derivs*m_face_coords;
	}
}


inline void Element::Coords(LinAlg::Matrix<double> & coords) const
{
	// Calculate a matrix with nodal coordinates
	coords.Resize(_n_nodes,3);   // 3 => X,Y,Z (all nodes have X,Y and Z)

	// Loop along all nodes of this element
	for (int i=0; i<_n_nodes; i++)
	{
		coords(i,0) = _connects[i]->X();
		coords(i,1) = _connects[i]->Y();
		coords(i,2) = _connects[i]->Z();
	}
}

inline void Element::Extrapolate(LinAlg::Vector<double> & IPValues, LinAlg::Vector<double> & NodalValues) const 
{
	//Extrapolation:
	//                  IPValues = E * NodalValues;
	//  where:
	//                              t           t
	//                         E = N * inv(N * N )
	//
	//  and            N = [shape functions matrix]
	//	                	  1		2		...		nNodes
	//	               1	[[N_11 N_12
	//	               2	 [N_21
	//	               :	 [
	//	              nIP	 [N_ ...					]]

	LinAlg::Matrix<double> N(_n_int_pts, _n_nodes);  // matrix of all IP shape functions
	LinAlg::Matrix<double> E(_n_nodes, _n_int_pts);  // Extrapolator matrix
	LinAlg::Vector<double> shape(_n_nodes);
	
	//filling N matrix
	for (int i_ip=0; i_ip<_n_int_pts; i_ip++)
	{
		double r = _a_int_pts[i_ip].r;
		double s = _a_int_pts[i_ip].s;
		double t = _a_int_pts[i_ip].t;
		Shape(r, s, t, shape);
		for (int j_node=0; j_node<_n_nodes; j_node++)
			N(i_ip, j_node) = shape(j_node);
	}

	//calculate extrapolator matrix
	assert(_n_nodes>=_n_int_pts); // TODO: replace with trow
	E = trn(N)*inv(N*trn(N));

	//perform extrapolation
	NodalValues = E*IPValues;
	
};


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


/////////////////////////////////////////////////////////////////////////////////////// Maps ///////////////////////


// Define a structure to contain DOF (Degree of Freedom) information from derived classes
struct DOFInfo
{
	Array<String> NodeEssential;
	Array<String> NodeNatural;
	Array<String> FaceEssential;
	Array<String> FaceNatural;
};

// Typdef of the map that contains all dof info from problem type elements
typedef std::map<String, DOFInfo, std::less<String> > DOFInfoMap_t;

// Instantiate the map that contains all dof info from problem type elements
DOFInfoMap_t DOFInfoMap;

#endif // MECHSYS_FEM_ELEMENT
