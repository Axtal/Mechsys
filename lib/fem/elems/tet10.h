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

#ifndef MECHSYS_FEM_TET10_H
#define MECHSYS_FEM_TET10_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "vtkcelltype.h"

namespace FEM
{

class Tet10 : public GeomElem
{
public:
    // Auxiliar structure to map local face IDs to local node IDs
    static const int Face2Node[4][6]; // 4 faces, 6 nodes/face

    // Constructor
    Tet10 (int NDim);

    // Derived methods
    void   SetIPs     (int TotNIP);
    int    VTKType    () const { return VTK_QUADRATIC_TETRA; }
    size_t FNode      (size_t IdxFace, size_t IdxFNode) const { return Face2Node[IdxFace][IdxFNode]; }
    void   Shape      (double r, double s, double t)    const;
    void   Derivs     (double r, double s, double t)    const;
    void   FaceShape  (double r, double s)              const;
    void   FaceDerivs (double r, double s)              const;
};

const int Tet10::Face2Node[4][6] = {{ 0, 4, 7, 3 },
                                    { 1, 2, 6, 5 },
                                    { 0, 3, 2, 1 }, 
                                    { 4, 5, 6, 7 }};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Tet10::Tet10 (int NDim)
    : GeomElem(NDim, /*NN*/10, /*NFN*/6, /*rCt*/1.0/3.0, /*sCt*/1.0/3.0, /*tCt*/1.0/3.0, "Tet10")
{
    SetIPs (4);
}

inline void Tet10::SetIPs (int TotNIP)
{
         if (TotNIP==  4) IPs = TET_IP4;
    else throw new Fatal("Tet10::SetIPs: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = TRI_IP6;
    NFIP = 6; 
}

inline void  Tet10::Shape(double r, double s, double t) const
{
    /*                         t
     *                         |
     *                         |                                           
     *                         | 3                                         
     *                         @,                                          
     *                        /|`                                          
     *                        ||  `,                                       
     *                       / |    ',                                     
     *                       | |      \                                    
     *                      /  |       `.                                  
     *                      |  |         `,  9                             
     *                     /   @ 7         `@                              
     *                     |   |             \                             
     *                    /    |              `.                           
     *                    |    |                ',                         
     *                 8 @     |                  \                        
     *                   |     @.,,_       6       `.                      
     *                  |     / 0   ``'-.,,@_        `.                    
     *                  |    /              ``''-.,,_  ', 2                
     *                 |    /                        ``'@.,,,              
     *                 |   '                       ,.-``     ``''- s       
     *                |  ,@ 4                 _,-'`                        
     *                ' /                 ,.'`                             
     *               | /             _.@``                                 
     *               '/          ,-'`   5                                  
     *              |/      ,.-``                                          
     *              /  _,-``                                               
     *            .@ '`                                                    
     *           / 1                                                       
     *          /                                                          
     *         /                                                           
     *        r                                                            
     */
    double u = 1.0 - r - s - t;
    
    // corners
    N(0) = u*(2.0*u - 1.0);
    N(1) = r*(2.0*r - 1.0);
    N(2) = s*(2.0*s - 1.0);
    N(3) = t*(2.0*t - 1.0);
    
    // midedge
    N(4) = 4.0 * u * r;
    N(5) = 4.0 * r * s;
    N(6) = 4.0 * s * u;
    N(7) = 4.0 * u * t;
    N(8) = 4.0 * r * t;
    N(9) = 4.0 * s * t;
}

inline void Tet10::Derivs(double r, double s, double t) const
{
    /*           _     _ T
     *          |  dNi  |
     * Derivs = |  ---  |   , where cj = r, s, t
     *          |_ dcj _|
     *
     * Derivs(j,i), j=>local coordinate and i=>shape function
     */
    // r-derivatives: dN0/dr to dN9/dr
    dNdR(0,0) =  4.0*(r + s + t) - 3.0;
    dNdR(0,1) =  4.0*r - 1.0;
    dNdR(0,2) =  0.0;
    dNdR(0,3) =  0.0;
    dNdR(0,4) =  4.0 - 8.0*r - 4.0*s - 4.0*t;
    dNdR(0,5) =  4.0*s;
    dNdR(0,6) = -4.0*s;
    dNdR(0,7) = -4.0*t;
    dNdR(0,8) =  4.0*t;
    dNdR(0,9) =  0.0;
    
    // s-derivatives: dN0/ds to dN9/ds
    dNdR(1,0) =  4.0*(r + s + t) - 3.0;
    dNdR(1,1) =  0.0;
    dNdR(1,2) =  4.0*s - 1.0;
    dNdR(1,3) =  0.0;
    dNdR(1,4) = -4.0*r;
    dNdR(1,5) =  4.0*r;
    dNdR(1,6) =  4.0 - 4.0*r - 8.0*s - 4.0*t;
    dNdR(1,7) = -4.0*t;
    dNdR(1,8) =  0.0;
    dNdR(1,9) =  4.0*t;
    
    // t-derivatives: dN0/dt to dN9/dt
    dNdR(2,0) =  4.0*(r + s + t) - 3.0;
    dNdR(2,1) =  0.0;
    dNdR(2,2) =  0.0;
    dNdR(2,3) =  4.0*t - 1.0;
    dNdR(2,4) = -4.0*r;
    dNdR(2,5) =  0.0;
    dNdR(2,6) = -4.0*s;
    dNdR(2,7) =  4.0 - 4.0*r - 4.0*s - 8.0*t;
    dNdR(2,8) =  4.0*r;
    dNdR(2,9) =  4.0*s;
}

inline void Tet10::FaceShape(double r, double s) const
{
    /*           s
     *           ^
     *           |
     *         2 @.
     *           | '.
     *           |   '. 4
     *         5 @     @.  
     *           |       '.
     *           |         '.
     *           @-----@-----@-> r
     *           0     3     1
     */
    FN(0) = 1.0-(r+s)*(3.0-2.0*(r+s));
    FN(1) = r*(2.0*r-1.0);
    FN(2) = s*(2.0*s-1.0);
    FN(3) = 4.0*r*(1.0-(r+s));
    FN(4) = 4.0*r*s;
    FN(5) = 4.0*s*(1.0-(r+s));
}

inline void Tet10::FaceDerivs(double r, double s) const
{
    /*           _     _ T
     *          |  dNi  |
     * Derivs = |  ---  |   , where cj = r, s
     *          |_ dcj _|
     *
     * Derivs(j,i), j=>local coordinate and i=>shape function
     */
    FdNdR(0,0) = -3.0 + 4.0 * (r + s);       FdNdR(1,0) = -3.0 + 4.0*(r + s);
    FdNdR(0,1) =  4.0 * r - 1.;              FdNdR(1,1) =  0.0 ;
    FdNdR(0,2) =  0.0;                       FdNdR(1,2) =  4.0 * s - 1.0;
    FdNdR(0,3) =  4.0 - 8.0 * r - 4.0 * s;   FdNdR(1,3) = -4.0 * r;
    FdNdR(0,4) =  4.0 * s;                   FdNdR(1,4) =  4.0 * r;
    FdNdR(0,5) = -4.0 * s;                   FdNdR(1,5) =  4.0 - 4.0 * r - 8.0*s;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Tet10Maker (int NDim) { return new Tet10(NDim); }

// Register element
int Tet10Register ()
{
    GeomElemFactory["Tet10"] = Tet10Maker;
    GEOM.Set ("Tet10", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Tet10_dummy_int  = Tet10Register();

}; // namespace FEM

#endif // MECHSYS_FEM_TET10_H
