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

#ifndef MECHSYS_FEM_GEOMELEM_H
#define MECHSYS_FEM_GEOMELEM_H

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

typedef LinAlg::Matrix<double> Mat_t;
typedef LinAlg::Vector<double> Vec_t;
typedef char const *           Str_t;

class GeomElem
{
public:
	// Constructor
	GeomElem() : NDim(0) {}

	// Destructor
	virtual ~GeomElem() {}

	// Methods
	void Initialize (int nDim) { NDim=nDim; _initialize(); }

	// Methods related to GEOMETRY
	        bool         CheckConn () const;
	        Node       * Nod       (size_t i)       { return Conn[i]; }
	        Node const * Nod       (size_t i) const { return Conn[i]; }
	        double       Volume    () const;
	        void         Extrap    (Vec_t & IPVals, Vec_t & NodVals) const;
	        void         InvMap    (double X, double Y, double Z,
	                                double & r, double & s, double & t) const;
	        bool         IsInside  (double X, double Y, double Z) const;
	virtual void         SetIPs    (int NIPs1D)                             =0;
	virtual int          VTKType   () const                                 =0;
	virtual void         VTKConn   (String & Nodes) const                   =0;
	virtual void         GetFNodes (int FaceID, Array<Node*> & FConn) const =0;
	virtual double       BoundDist (double r, double s, double t) const     =0;

	// Public data (read only)
	size_t             NDim;       ///< Space dimension (2D or 3D)
	size_t             NNodes;     ///< Number of nodes
	size_t             NFaceNodes; ///< Number of face nodes
	size_t             NIPs;       ///< Number of integration points
	size_t             NFaceIPs;   ///< Number of integration points of face
	IntegPoint const * IPs;        ///< Integration points
	IntegPoint const * FaceIPs;    ///< Integration points of Faces/Edges
	Array<Node*>       Conn;       ///< Connectivity (pointers to nodes of this element) Size==NNodes

private:
	// To be derived
	virtual void _initialize () =0;
	virtual void _shape      (double r, double s, double t, Vec_t & N)  const =0;
	virtual void _face_shape (double r, double s, Vec_t & N)            const =0;
	virtual void _derivs     (double r, double s, double t, Mat_t & dN) const =0;

	// Private methods
	void _jacobian (Mat_t const & dN, Mat_t & J) const;
	//_local_coords

}; // class GeomElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


/* public */

inline bool GeomElem::CheckConn() const
{
	if (Conn.Size()!=NNodes) return false;
	for (size_t i=0; i<NNodes; ++i) if (Conn[i]==NULL) return false;
	return true;
}

inline double GeomElem::Volume() const
{
	Mat_t dN; // derivs: NDim x NNodes
	Mat_t J;  // Jacobian
	double vol=0.0;
	for (size_t i=0; i<NIPs; ++i)
	{
		Derivs    (IPs[i].r, IPs[i].s, IPs[i].t, dN);
		_jacobian (dN, J);
		vol += det(J);
	}
	return vol;
}

inline void GeomElem::Extrap(Vec_t & IPVals, Vec_t & NodVals) const 
{
	// Check
	if (IPVals.Size()!=static_cast<int>(NIPs)) throw new Fatal("Element::Extrap: IPVals.Size()==%d must be equal to NIPs==%d",IPVals.Size(),NIPs);

	// Data
	size_t m = IPVals.Size();
	size_t n = NNodes;
	NodVals.Resize(n);

	/* N: shape functions matrix:
	          _                                                 _
	         |   N11      N12      N13      ...  N1(nNode)       |
	         |   N21      N22      N23      ...  N2(nNode)       |
	     N = |   N31      N32      N33      ...  N3(nNode)       |
	         |          ...............................          |
	         |_  N(nIP)1  N(nIP)2  N(nIP)3  ...  N(nIP)(nNode)  _| [m=nIP x n=nNode]
	*/
	Mat_t N(m, n);
	Vec_t shape(n);
	Mat_t ip_coord(m, NDim+1);   // IP coordinates matrix (local)
	Mat_t node_coord;             // nodal coordinates matrix (local)

	// Mount N and ip_coord matrices
	for (size_t i=0; i<m; i++)
	{
		Shape (IPs[i].r, IPs[i].s, IPs[i].t, shape);
		for (size_t j=0; j<n; j++) N(i,j) = shape(j);
		ip_coord(i, 0) = IPs[i].r;
		ip_coord(i, 1) = IPs[i].s;
		if (NDim==2) ip_coord(i, 2) = 1.0;
		else
		{
			ip_coord(i, 2) = IPs[i].t;
			ip_coord(i, 3) = 1.0;
		}
	}

	// Mount node_coord matrix
	//_local_coords (node_coord);

	// Extrapolator matrix
	Mat_t E(n, m);

	//calculate extrapolator matrix
	if (n==m)     E = inv(N);
	else if (m>n) E = invg(N);
	else
	{
		Mat_t I;
		identity (m,I);
		E = invg(N)* (I - ip_coord*invg(ip_coord)) + node_coord*invg(ip_coord);
	}
	
	// Extrapolate
	NodVals = E * IPVals;
}

inline void GeomElem::InvMap(double x, double y, double z, double & r, double & s, double & t) const
{
	Vec_t shape;
	Mat_t derivs;
	Mat_t J; //Jacobian matrix
	Vec_t f;
	Vec_t delta;
	     if (NDim==2) f.Resize(2);
	else if (NDim==3) f.Resize(3);
	double tx, ty, tz; //x, y, z trial
	double norm_f;
	r = s = t =0; // first suposition for natural coordinates
	int max_steps= 25;
	int k=0;
	for (k=0; k<max_steps; k++)
	{
		Shape (r, s, t, shape);
		Derivs(r, s, t, derivs);
		Jacobian(derivs, J);
		tx = ty = tz = 0; 

		//calculate trial of real coordinates
		for (size_t j=0; j<NNodes; j++) 
		{
			tx += shape(j)*Conn[j]->Coord(0); //ok
			ty += shape(j)*Conn[j]->Coord(1); //ok
			if (NDim==3)
			tz += shape(j)*Conn[j]->Coord(2); //ok
		}
		
		// Calculate the error
		f(0) = tx - x;
		f(1) = ty - y;
		if (NDim==3)
		f(2) = tz - z;
		
		// Calculate the corrector
		delta = trn(inv(J))*f;
		
		r -= delta(0);
		s -= delta(1);
		if (NDim==3)
		t -= delta(2);

		norm_f = sqrt(trn(f)*f);

		if (norm_f<1.0E4) break;
	} 
	if (k==max_steps) throw new Fatal("Element::InverseMap: InverseMap did not converge after %d steps in element", max_steps);
}

inline bool GeomElem::IsInside(double X, double Y, double Z) const
{
	double tiny = 1.0e-4;  // TODO: use 0.0001*MinDiagonalOfElement
	double huge = 1.0e+20;

	// Search along X
	double max = -huge;
	for (size_t i=0; i<NNodes; i++) if (Conn[i]->Coord(0)>max) max=Conn[i]->Coord(0);
	if (X>max) return false;
	double min = +huge;
	for (size_t i=0; i<NNodes; i++) if (Conn[i]->Coord(0)<min) min=Conn[i]->Coord(0);
	if (X<min) return false;

	// Search along Y
	max = -huge;
	for (size_t i=0; i<NNodes; i++) if (Conn[i]->Coord(1)>max) max=Conn[i]->Coord(1);
	if (Y>max) return false;
	min = +huge;
	for (size_t i=0; i<NNodes; i++) if (Conn[i]->Coord(1)<min) min=Conn[i]->Coord(1);
	if (Y<min) return false;

	// Search along Z
	max = -huge;
	for (size_t i=0; i<NNodes; i++) if (Conn[i]->Coord(2)>max) max=Conn[i]->Coord(2);
	if (Z>max) return false;
	min = +huge;
	for (size_t i=0; i<NNodes; i++) if (Conn[i]->Coord(2)<min) min=Conn[i]->Coord(2);
	if (Z<min) return false;

	double r, s, t;
	InverseMap (X,Y,Z, r,s,t);
	if (BoundDist(r,s,t)>-tiny) return true;
	else return false;
}

/* private */

inline void GeomElem::_jacobian(Mat_t const & dN, Mat_t & J) const
{
	Mat_t coords(NNodes, NDim);
	for (size_t i=0; i<NNodes; ++i)
	for (size_t j=0; j<NDim;   ++j)
		coords(i,j) = Conn[i]->Coord(j);
	J = dN * coords;
}

}; // namespace FEM

#endif // MECHSYS_FEM_GEOMELEM
