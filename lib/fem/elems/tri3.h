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

#ifndef MECHSYS_FEM_TRI3_H
#define MECHSYS_FEM_TRI3_H

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/geomelem.h>

namespace FEM
{

class Tri3: public GeomElem
{
public:
    // Auxiliar structure to map local face IDs to local node IDs
    static const int Face2Node[3][2]; // 3 edges, 2 nodes/edge

    // Constructor
    Tri3 (int NDim);

    // Derived methods
    void   SetIPs     (int TotNIP);
    size_t FNode      (size_t IdxFace, size_t IdxFNode) const { return Face2Node[IdxFace][IdxFNode]; }
    void   Shape      (double r, double s, double t)    const;
    void   Derivs     (double r, double s, double t)    const;
    void   FaceShape  (double r, double s)              const;
    void   FaceDerivs (double r, double s)              const;
    void   NatCoords  (Mat_t & C)                       const;
};


/* Local IDs
             Nodes                 Faces

   y           2
   |           @                     @
   +--x       / \                   / \
             /   \                 /   \
            /     \             2 /     \ 1
           /       \             /       \
          /         \           /         \
         @-----------@         @-----------@
        0             1              0
*/
const int Tri3::Face2Node[3][2] = {{ 0, 1 },
                                   { 1, 2 },
                                   { 2, 0 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Tri3::Tri3 (int NDim)
    : GeomElem(NDim, /*NN*/3, /*NFN*/2, "Tri3")
{
    SetIPs (3);
}

inline void Tri3::SetIPs (int TotNIP)
{
         if (TotNIP== 1) IPs = TRI_IP1;
    else if (TotNIP== 3) IPs = TRI_IP3;
    else if (TotNIP== 4) IPs = TRI_IP4;
    else if (TotNIP== 6) IPs = TRI_IP6;
    else if (TotNIP== 7) IPs = TRI_IP7;
    else if (TotNIP==13) IPs = TRI_IP13;
    else throw new Fatal("Tri3::SetIPs: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = LIN_IP2;
    NFIP = 2;
}

inline void Tri3::Shape (double r, double s, double t) const
{

    /*    s
     *    ^
     *    |
     *  2 
     *    @,(0,1)
     *    | ',
     *    |   ',
     *    |     ',
     *    |       ',    
     *    |         ',
     *    |           ', 
     *    |             ', 
     *    |               ', 
     *    |(0,0)            ', (1,0)
     *    @-------------------@  --> r
     *  0                      1
     */
    N(0) = 1.0-r-s;
    N(1) = r;
    N(2) = s;
}

inline void Tri3::Derivs (double r, double s, double t) const
{
    /*           _     _ T
     *          |  dNi  |
     *   dNdR = |  ---  |   , where Rj = r, s
     *          |_ dRj _|
     *  
     *   dNdR(j,i), j=>local coordinate and i=>shape function
     */
    dNdR(0,0) = -1.0;    dNdR(1,0) = -1.0;
    dNdR(0,1) =  1.0;    dNdR(1,1) =  0.0;
    dNdR(0,2) =  0.0;    dNdR(1,2) =  1.0;
}

inline void Tri3::FaceShape (double r, double s) const
{
    /*  
     *       0           |           1
     *       @-----------+-----------@-> r
     *      -1           |          +1
     */
    FN(0) = 0.5*(1.0-r);
    FN(1) = 0.5*(1.0+r);
}

inline void Tri3::FaceDerivs (double r, double s) const
{
    /*            _     _ T
     *           |  dNi  |
     *   FdNdR = |  ---  |   , where Rj = r, s
     *           |_ dRj _|
     *
     *   FdNdR(j,i), j=>local coordinate and i=>shape function
     */
    FdNdR(0,0) = -0.5;
    FdNdR(0,1) =  0.5;
}

inline void Tri3::NatCoords (Mat_t & C) const
{
    C.change_dim(3,3);
    C = 0.0, 0.0, 1.0,
        1.0, 0.0, 1.0,
        0.0, 1.0, 1.0;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Tri3Maker (int NDim) { return new Tri3(NDim); }

// Register element
int Tri3Register ()
{
    GeomElemFactory["Tri3"] = Tri3Maker;
    GEOM.Set ("Tri3", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Tri3_dummy_int  = Tri3Register();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI3_H
