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

#ifndef MECHSYS_FEM_LIN2_H
#define MECHSYS_FEM_LIN2_H

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/geomelem.h>

namespace FEM
{

class Lin2: public GeomElem
{
public:
    // Auxiliary structure to map local face IDs to local node IDs
    static const int Face2Node[2][1]; // 2 extreme points, 1 point/extremity

    // Constructor
    Lin2 (int NDim);

    // Derived methods
    void   SetIPs     (int TotNIP);
    size_t FNode      (size_t IdxFace, size_t IdxFNode) const { return Face2Node[IdxFace][IdxFNode]; }
    void   Shape      (double r, double s, double t)    const;
    void   Derivs     (double r, double s, double t)    const;
    void   FaceShape  (double r, double s)              const {}
    void   FaceDerivs (double r, double s)              const {}
    void   NatCoords  (Mat_t & C)                       const;
};


/* Local IDs
             Nodes              Extremities
   y
   |
   +--x  @-----------@         @-----------@
        0             1        0           1
*/
const int Lin2::Face2Node[2][1] = {{ 0 },
                                   { 1 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Lin2::Lin2 (int NDim)
    : GeomElem(/*NDim*/1, /*NN*/2, /*NFN*/2, "Lin2")
{
    SetIPs (3);
}

inline void Lin2::SetIPs (int TotNIP)
{
         if (TotNIP== 2) IPs = LIN_IP2;
    else if (TotNIP== 3) IPs = LIN_IP3;
    else if (TotNIP== 4) IPs = LIN_IP4;
    else if (TotNIP== 5) IPs = LIN_IP5;
    else throw new Fatal("Lin2::SetIPs: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = NULL;
    NFIP = 0;
}

inline void Lin2::Shape (double r, double s, double t) const
{

    /*   @---------|---------@  --> r
     *  -1         0         1
     */
    N(0) = (1.0-r)/2.0;
    N(1) = (1.0+r)/2.0;
}

inline void Lin2::Derivs (double r, double s, double t) const
{
    /*           _     _ T
     *          |  dNi  |
     *   dNdR = |  ---  |   , where Rj = r
     *          |_ dRj _|
     *  
     *   dNdR(j,i), j=>local coordinate and i=>shape function
     */
    dNdR(0,0) = -1.0/2.0;
    dNdR(0,1) =  1.0/2.0;
}

inline void Lin2::NatCoords (Mat_t & C) const
{
    C.change_dim(2,3);
    C = -1.0, 0.0, 1.0,
         1.0, 0.0, 1.0;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Lin2Maker (int NDim) { return new Lin2(NDim); }

// Register element
int Lin2Register ()
{
    GeomElemFactory["Lin2"] = Lin2Maker;
    GEOM.Set ("Lin2", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Lin2_dummy_int  = Lin2Register();

}; // namespace FEM

#endif // MECHSYS_FEM_LIN2_H
