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

#ifndef MECHSYS_FEM_GEOMELEM_H
#define MECHSYS_FEM_GEOMELEM_H

// MechSys
#include <mechsys/geomtype.h>
#include <mechsys/fem/quadrature.h>
#include <mechsys/fem/node.h>
#include <mechsys/util/string.h>
#include <mechsys/linalg/matvec.h>

namespace FEM
{

class GeomElem
{
public:
    // Constructor
    GeomElem (int          NDim,                         ///< Space dimension
              size_t       NN,                           ///< Number of nodes
              size_t       NFN,                          ///< Number of nodes on face/edge
              double       rCt,                          ///< Natural r coordinate of centroid
              double       sCt,                          ///< Natural s coordinate of centroid
              double       tCt,                          ///< Natural t coordinate of centroid
              char const * Name="__unnamed_geomelem__"); ///< Name. Ex: Tri3, Quad4, ...

    // Methods to be overloaded
    virtual void   SetIPs     (int TotNIP)                            =0;
    virtual int    VTKType    ()                                const =0;
    virtual size_t FNode      (size_t IdxFace, size_t IdxFNode) const =0;
    virtual void   Shape      (double r, double s, double t)    const =0;
    virtual void   Derivs     (double r, double s, double t)    const =0;
    virtual void   FaceShape  (double r, double s)              const =0;
    virtual void   FaceDerivs (double r, double s)              const =0;

    // Constants
    size_t             NN;   ///< Number of nodes
    size_t             NFN;  ///< Number of face nodes
    size_t             NIP;  ///< Number of integration points
    size_t             NFIP; ///< Number of integration points of face
    IntegPoint const * IPs;  ///< Integration points
    IntegPoint const * FIPs; ///< Integration points of Faces/Edges
    IntegPoint         Rct;  ///< r,s,t coordinates of centroid
    String             Name; ///< Name. Ex: Tri3, Quad4, ...

    // Mutable data
    mutable Vec_t N;     ///< (current) (size=NN)           Shape functions
    mutable Mat_t dNdR;  ///< (current) (size=NDim x NN)    Derivative of shape functions w.r.t natural coordinates
    mutable Vec_t FN;    ///< (current) (size=NFN)          Face shape functions
    mutable Mat_t FdNdR; ///< (current) (size=NDim-1 x NFN) Derivative of face shape functions w.r.t. natural coordinates
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline GeomElem::GeomElem (int NDim, size_t TheNN, size_t TheNFN, double rCt, double sCt, double tCt, char const * TheName)
    : NN(TheNN), NFN(TheNFN), Name(TheName)
{
    Rct.r = rCt;
    Rct.s = sCt;
    Rct.t = tCt;
    N    .change_dim (NN);
    dNdR .change_dim (NDim,NN);
    FN   .change_dim (NFN);
    FdNdR.change_dim (NDim-1,NFN);
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


SDPair GEOM;

typedef GeomElem * (*GeomElemMakerPtr)(int NDim);

typedef std::map<String, GeomElemMakerPtr> GeomElemFactory_t;

GeomElemFactory_t GeomElemFactory;

GeomElem * AllocGeomElem(String const & Name, int NDim)
{
    GeomElemFactory_t::iterator it = GeomElemFactory.find(Name);
    if (it==GeomElemFactory.end()) throw new Fatal("AllocGeomElem: '%s' is not available", Name.CStr());

    GeomElem * ptr = (*it->second)(NDim);

    return ptr;
}

}; // namespace FEM

#ifdef USE_BOOST_PYTHON
double PyGEOM (BPy::str const & Key) { return FEM::GEOM(BPy::extract<char const *>(Key)()); }
#endif

#endif // MECHSYS_FEM_GEOMELEM
