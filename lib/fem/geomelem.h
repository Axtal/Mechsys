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

	// Methods related to GEOMETRY
	        double Volume    () const;
	        void   Extrap    (Vec_t & IPVals, Vec_t & NodVals) const;
	        void   InvMap    (double X, double Y, double Z, double & r, double & s, double & t) const;
	        bool   IsInside  (double X, double Y, double Z) const;
	virtual void   SetIPs    (int NIPs1D)                             =0;
	virtual int    VTKType   () const                                 =0;
	virtual void   VTKConn   (String & Nodes) const                   =0;
	virtual void   GetFNodes (int FaceID, Array<Node*> & FConn) const =0;
	virtual double BoundDist (double r, double s, double t) const { return -1; };

	// Methods
	        void Initialize (int nDim);                                                             ///< Initialize this element
	        bool CheckConn  () const;                                                               ///< Check connectivity
	        void Jacobian   (Mat_t const & dN, Mat_t & J) const;                                    ///< Jacobian matrix
	        void FaceJacob  (Array<FEM::Node*> const & FConn, double r, double s, Mat_t & J) const; ///< Face Jacobian matrix
	virtual void Shape      (double r, double s, double t, Vec_t & N)  const =0;                    ///< Shape functions
	virtual void Derivs     (double r, double s, double t, Mat_t & dN) const =0;                    ///< Derivatives of shape functions
	virtual void FaceShape  (double r, double s, Vec_t & FN)           const =0;                    ///< Face shape functions
	virtual void FaceDerivs (double r, double s, Mat_t & FdN)          const =0;                    ///< Face derivatives of shape functions

	// Public data (read only)
	size_t             NDim;    ///< Space dimension (2D or 3D)
	size_t             NNodes;  ///< Number of nodes
	size_t             NFNodes; ///< Number of face nodes
	size_t             NIPs;    ///< Number of integration points
	size_t             NFIPs;   ///< Number of integration points of face
	IntegPoint const * IPs;     ///< Integration points
	IntegPoint const * FIPs;    ///< Integration points of Faces/Edges
	Array<Node*>       Conn;    ///< Connectivity (size==NNodes). Initialized by ProbElem

private:
	virtual void _local_coords (Mat_t & C) const =0; ///< Return the local coordinates of the nodes

}; // class GeomElem


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline void GeomElem::Initialize(int nDim)
{
	if (nDim<2 || nDim>3) throw new Fatal("GeomElem::Initialize: nDim must be either 2 or 3 (%d is invalid)",nDim);
	NDim = nDim;
}

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
		Derivs   (IPs[i].r, IPs[i].s, IPs[i].t, dN);
		Jacobian (dN, J);
		vol += det(J);
	}
	return vol;
}

inline void GeomElem::Extrap(Vec_t & IPVals, Vec_t & NodVals) const 
{
	// Check
	if (IPVals.Size()!=static_cast<int>(NIPs)) throw new Fatal("GeomElem::Extrap: IPVals.Size()==%d must be equal to NIPs==%d",IPVals.Size(),NIPs);

	// Data
	size_t m = IPVals.Size();
	size_t n = NNodes;
	NodVals.Resize(n);

	/* Nmat: shape functions matrix:
	             _                                                 _
	            |   N11      N12      N13      ...  N1(nNode)       |
	            |   N21      N22      N23      ...  N2(nNode)       |
	     Nmat = |   N31      N32      N33      ...  N3(nNode)       |
	            |          ...............................          |
	            |_  N(nIP)1  N(nIP)2  N(nIP)3  ...  N(nIP)(nNode)  _| [m=nIP x n=nNode]
	*/

	// Mount Nmat and ip_coords matrices
	Mat_t Nmat      (m, n);
	Vec_t N         (n);
	Mat_t ip_coords (m, NDim+1); // IP coordinates matrix (local)
	for (size_t i=0; i<m; i++)
	{
		Shape (IPs[i].r, IPs[i].s, IPs[i].t, N);
		for (size_t j=0; j<n; j++) Nmat(i,j) = N(j);
		ip_coords(i, 0) = IPs[i].r;
		ip_coords(i, 1) = IPs[i].s;
		if (NDim==2) ip_coords(i, 2) = 1.0;
		else
		{
			ip_coords(i, 2) = IPs[i].t;
			ip_coords(i, 3) = 1.0;
		}
	}

	// Nodal coordinates matrix (local)
	Mat_t coords;
	_local_coords (coords);

	// Extrapolation matrix
	Mat_t E(n, m);
	     if (n==m) E = inv(Nmat);
	else if (m>n)  E = invg(Nmat);
	else
	{
		Mat_t I;
		identity (m,I);
		E = invg(Nmat)* (I-ip_coords*invg(ip_coords)) + coords*invg(ip_coords);
	}
	
	// Extrapolate
	NodVals = E * IPVals;
}

inline void GeomElem::InvMap(double x, double y, double z, double & r, double & s, double & t) const
{
	Vec_t N;  // Shape
	Mat_t dN; // Derivs
	Mat_t J;  // Jacobian matrix
	Vec_t f;  // Residual between real and trial global coordinates
	Vec_t delta;
	     if (NDim==2) f.Resize(2);
	else if (NDim==3) f.Resize(3);
	double tx, ty, tz;  // x, y, z trial
	r = s = t = 0;      // First suposition for natural coordinates
	int max_steps = 25; // Number of maximum steps
	int k = 0;
	for (k=0; k<max_steps; k++)
	{
		Shape    (r, s, t, N);
		Derivs   (r, s, t, dN);
		Jacobian (dN, J);
		tx = ty = tz = 0; 

		// Calculate trial of real coordinates
		for (size_t j=0; j<NNodes; j++) 
		{
			             tx += N(j)*Conn[j]->Coord(0); 
			             ty += N(j)*Conn[j]->Coord(1); 
			if (NDim==3) tz += N(j)*Conn[j]->Coord(2); 
		}
		
		// Calculate residual
		             f(0) = tx - x;
		             f(1) = ty - y;
		if (NDim==3) f(2) = tz - z;
		
		// Calculate corrector
		delta = trn(inv(J))*f;
		
		             r -= delta(0);
		             s -= delta(1);
		if (NDim==3) t -= delta(2);

		double norm_f = sqrt(trn(f)*f); // Norm of the residual

		if (norm_f<1.0e+4) break;
	}
	if (k==max_steps) throw new Fatal("GeomElem::InvMap: InvMap did not converge after %d steps", max_steps);
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
	InvMap (X,Y,Z, r,s,t);
	if (BoundDist(r,s,t)>-tiny) return true;
	else return false;
}

inline void GeomElem::Jacobian(Mat_t const & dN, Mat_t & J) const
{
	Mat_t coords(NNodes, NDim);
	for (size_t i=0; i<NNodes; ++i)
	for (size_t j=0; j<NDim;   ++j)
		coords(i,j) = Conn[i]->Coord(j);
	J = dN * coords;
}

inline void GeomElem::FaceJacob(Array<FEM::Node*> const & FConn, double r, double s, Mat_t & J) const
{
	if (NFNodes>0)
	{
		// Calculate shape function derivatives for a face
		Mat_t FdN;
		FaceDerivs(r, s, FdN);

		// Get the face coordinates
		Mat_t fcoords(NFNodes,NDim);
		for (size_t i=0; i<NFNodes; i++)
		for (size_t j=0; j<NDim;    j++)
			fcoords(i,j) = FConn[i]->Coord(j);

		// Calculate face Jacobian
		J = FdN*fcoords;
	}
	else throw new Fatal("GeomElem::FaceJacob: This element does not have any face/edge (NFNodes==%d)",NFNodes);
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Define a pointer to a function that makes (allocate) a new GeomElem
typedef GeomElem * (*GeomElemMakerPtr)();

// Typdef of the array map that contains all the pointers to the functions that makes elements
typedef std::map<String, GeomElemMakerPtr, std::less<String> > GeomElemFactory_t;

// Instantiate the array map that contains all the pointers to the functions that makes elements
GeomElemFactory_t GeomElemFactory;

// Allocate a new GeomElem according to a string giving the name of the element
GeomElem * AllocGeomElem(char const * GeomElemName)
{
	GeomElemMakerPtr ptr=NULL;
	ptr = GeomElemFactory[GeomElemName];
	if (ptr==NULL) throw new Fatal("FEM::AllocGeomElem: There is no < %s > implemented in this library", GeomElemName);
	return (*ptr)();
}

}; // namespace FEM

#endif // MECHSYS_FEM_GEOMELEM
