/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso                                *
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

/* Definitions - Copyright (C) 2007 Dorival de Moraes Pedroso */

#ifndef MPM_DEFS_H
#define MPM_DEFS_H

// Std Lib
#include <map>

// Blitz++
#include <blitz/tinyvec-et.h>
#include <blitz/tinymat.h>

// FLTK
#include <FL/Enumerations.H> // for Fl_Color

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/maps.h>

namespace MPM {

/** Curve type. */
enum CurveType { CT_POINTS, CT_LINES, CT_BOTH };

/** Color Maps. */
enum ClrMapType { CMT_BWP=0, CMT_HOT=1, CMT_JET=2 };

/** Fields. */
enum FieldType { FLD_EPSXX, FLD_EPSYY, FLD_EPSZZ, FLD_EPSXY,
                 FLD_SIGXX, FLD_SIGYY, FLD_SIGZZ, FLD_SIGXY,
                 FLD_CAMP,  FLD_CAMQ,  FLD_VELN,  FLD_TAG};

/** Nodes fixing type. */
enum FixType { FIX_X, FIX_Y, FIX_Z, FIX_XY, FIX_YZ, FIX_ZX, FIX_XYZ };

/* Curve properties. */
struct CurveProps
{
	CurveType Typ;      ///< Curve type
	Fl_Color  Clr;      ///< Color
	int       Lty;      ///< Line type
	int       Lwd;      ///< Line width
	int       Pch;      ///< Point type
	int       Psz;      ///< Point size
	char      Nam[256]; ///< Name
};

// Constants
const double SQ2    = sqrt(2.0);
const double SQ3    = sqrt(3.0);
const double SQ6    = sqrt(6.0);
const double SQ2BY3 = sqrt(2.0/3.0);
const double PI     = 4.0*atan(1.0);

// Typedefs
typedef blitz::TinyVector<size_t,4>   Connec2D; ///< for quadrilateral connectivities
typedef blitz::TinyVector<double,3>   Vector3D; ///< for position, veloc, accel, etc.
typedef blitz::TinyVector<double,6>   STensor2; ///< Symmetric 2nd order tensors: for stress, strain, ...
typedef blitz::TinyVector<double,9>   ATensor2; ///< Assymmetric 2nd order tensors: for deformation gradient, ...
typedef blitz::TinyMatrix<double,6,6> STensor4; ///< Symmetric 4th order tensors: for stiffness ...
typedef blitz::TinyMatrix<double,9,9> ATensor4; ///< Assymmetric 4th order tensors: for stiffness ...
typedef blitz::TinyVector<size_t,6>   ConnecW6; ///< Connectivities for wedge with 6 nodes
typedef blitz::TinyVector<size_t,8>   ConnecH8; ///< Connectivities for hexahedra with 8 nodes

// Structures
const size_t NRVEC = 2;
typedef std::map<int,int>::const_iterator n2i_it;
struct ShapeAndGrads2D
{
	double            N[16];     ///< Shape functions (C0)
	double            S[16];     ///< Shape functions (C1)
	Vector3D          G[16];     ///< Gradients of the shape functions
	Vector3D          R0[NRVEC]; ///< Initial vectors: centre-to-right-side
	Vector3D          R [NRVEC]; ///< Vectors: centre-to-top-side
    std::map<int,int> n2i;       ///< maps gridnode to index in N,S,G arrays
};

// Typedefs for callback functions
typedef bool      (*ptIsFixed)       (Vector3D const & N, FixType & FType);           ///< Callback for applying boundary conditions
typedef bool      (*ptIsPointInGeom) (Vector3D const & P, int & ClrIdx);              ///< Pointer to a function which checks if a mat point is inside the geometry.
typedef double    (*ptDensity)       (Vector3D const & P);                            ///< Density function
typedef void      (*ptModelData)     (Vector3D const & P, String & Name, int & Tag, int NDim, SDPair & Prms, SDPair & Inis, String & Scheme); ///< Get model data
typedef void      (*ptIniVelocity)   (Vector3D const & P, Vector3D & v);              ///< Initialize velocities
typedef bool      (*ptHasTraction)   (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & t);   ///< Apply traction to some points inside boundary cells.
typedef bool      (*ptHasAppDisp)    (Vector3D const & P, Array<Vector3D> const & N, Vector3D const & L, int nPCell, Vector3D & fu);  ///< Apply displacements to some points inside boundary cells.
typedef void      (*ptAppDisp)       (Vector3D const & P, STensor2 const & s, Vector3D const & u, Vector3D const & L, int nPCell, Vector3D & fu);
typedef void      (*ptB     )        (double t,                      Vector3D & B  ); ///< Body force
typedef void      (*ptLdM   )        (double t                     , double   & M  ); ///< Multiplier for applied external forces
typedef void      (*ptVeloc )        (double t, Vector3D const & XY, Vector3D & Vel); ///< Correct velocity at time t and position XY
typedef void      (*ptStress)        (double t, Vector3D const & XY, STensor2 & Sig); ///< Correct stress at time t and position XY

}; // namespace MPM

#endif // MPM_DEFS_H
