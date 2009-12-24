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

#ifndef MECHSYS_FEM_TRI6_H
#define MECHSYS_FEM_TRI6_H

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/geomelem.h>
#include <mechsys/vtkcelltype.h>

namespace FEM
{

class Tri6: public GeomElem
{
public:
    // Auxiliar structure to map local face IDs to local node IDs
    static const int Face2Node[3][3]; // 3 edges, 3 nodes/edge

    // Constructor
    Tri6 (int NDim);

    // Derived methods
    void   SetIPs     (int TotNIP);
    int    VTKType    () const { return VTK_QUADRATIC_TRIANGLE; }
    size_t FNode      (size_t IdxFace, size_t IdxFNode) const { return Face2Node[IdxFace][IdxFNode]; }
    void   Shape      (double r, double s, double t)    const;
    void   Derivs     (double r, double s, double t)    const;
    void   FaceShape  (double r, double s)              const;
    void   FaceDerivs (double r, double s)              const;
};


/* Local IDs
             Nodes                 Faces

   y           2
   |           @                     @
   +--x       / \                   / \
           5 /   \ 4               /   \
            @     @             2 /     \ 1
           /       \             /       \
          /         \           /         \
         @-----@-----@         @-----------@
        0      3      1              0
*/
const int Tri6::Face2Node[3][3]= {{ 0, 1, 3 },
                                  { 1, 2, 4 },
                                  { 2, 0, 5 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Tri6::Tri6 (int NDim)
    : GeomElem(NDim, /*NN*/6, /*NFN*/3, /*rCt*/1.0/3.0, /*sCt*/1.0/3.0, /*tCt*/0.0, "Tri6")
{
    SetIPs (6);
}

inline void Tri6::SetIPs (int TotNIP)
{
         if (TotNIP==3)  IPs = TRI_IP3;
    else if (TotNIP==4)  IPs = TRI_IP4;
    else if (TotNIP==6)  IPs = TRI_IP6;
    else if (TotNIP==7)  IPs = TRI_IP7;
    else if (TotNIP==13) IPs = TRI_IP13;
    else throw new Fatal("tri6::SetIntPoints: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = LIN_IP3;
    NFIP = 3;
}

inline void Tri6::Shape (double r, double s, double t) const
{

    /*    s
     *    ^
     *    |
     *  2
     *    @,(0,1)
     *    | ',
     *    |   ',
     *    |     ',
     *    |       ',   4
     *  5 @          @
     *    |           ',
     *    |             ',
     *    |               ',
     *    |(0,0)            ', (1,0)
     *    @---------@---------@  --> r
     *  0           3          1
     */
    N(0) = 1.0-(r+s)*(3.0-2.0*(r+s));
    N(1) = r*(2.0*r-1.0);
    N(2) = s*(2.0*s-1.0);
    N(3) = 4.0*r*(1.0-(r+s));
    N(4) = 4.0*r*s;
    N(5) = 4.0*s*(1.0-(r+s));
}

inline void Tri6::Derivs (double r, double s, double t) const
{
    /*           _     _ T
     *          |  dNi  |
     *   dNdR = |  ---  |   , where Rj = r, s
     *          |_ dRj _|
     *  
     *   dNdR(j,i), j=>local coordinate and i=>shape function
     */
    dNdR(0,0) = -3.0 + 4.0 * (r + s);       dNdR(1,0) = -3.0 + 4.0*(r + s);
    dNdR(0,1) =  4.0 * r - 1.;              dNdR(1,1) =  0.0 ;
    dNdR(0,2) =  0.0;                       dNdR(1,2) =  4.0 * s - 1.0;
    dNdR(0,3) =  4.0 - 8.0 * r - 4.0 * s;   dNdR(1,3) = -4.0 * r;
    dNdR(0,4) =  4.0 * s;                   dNdR(1,4) =  4.0 * r;
    dNdR(0,5) = -4.0 * s;                   dNdR(1,5) =  4.0 - 4.0 * r - 8.0*s;
}

inline void Tri6::FaceShape (double r, double s) const
{
    /*
     *       @-----------@-----------@-> r
     *       0           2           1
     *       |           |           |
     *      r=-1         r=0        r=+1
     */
    FN(0) = 0.5 * (r*r-r);
    FN(1) = 0.5 * (r*r+r);
    FN(2) = 1.0 -  r*r;
}

inline void Tri6::FaceDerivs (double r, double s) const
{
    /*            _     _ T
     *           |  dNi  |
     *   FdNdR = |  ---  |   , where Rj = r, s
     *           |_ dRcj _|
     *
     *   FdNdR(j,i), j=>local coordinate and i=>shape function
     */
    FdNdR(0,0) =  r  - 0.5;
    FdNdR(0,1) =  r  + 0.5;
    FdNdR(0,2) = -2.0* r;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Tri6Maker (int NDim) { return new Tri6(NDim); }

// Register element
int Tri6Register ()
{
    GeomElemFactory["Tri6"] = Tri6Maker;
    GEOM.Set ("Tri6", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Tri6_dummy_int  = Tri6Register();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI6_H
