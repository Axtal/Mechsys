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

#ifndef MECHSYS_FEM_QUAD4_H
#define MECHSYS_FEM_QUAD4_H

// MechSys
#include "fem/node.h"
#include "fem/geomelem.h"
#include "vtkcelltype.h"

namespace FEM
{

class Quad4: public GeomElem
{
public:
    // Auxiliar structure to map local face IDs to local node IDs
    static const int Face2Node[4][2]; // 4 edges, 2 nodes/edge

    // Constructor
    Quad4 (int NDim);

    // Derived methods
    void   SetIPs     (int TotNIP);
    int    VTKType    () const { return VTK_QUAD; }
    size_t FNode      (size_t IdxFace, size_t IdxFNode) const { return Face2Node[IdxFace][IdxFNode]; }
    void   Shape      (double r, double s, double t)    const;
    void   Derivs     (double r, double s, double t)    const;
    void   FaceShape  (double r, double s)              const;
    void   FaceDerivs (double r, double s)              const;
};


/* Local IDs
                 Nodes                  Faces(edges)
   y
   |        3             2                  2
   +--x      @-----------@             +-----------+
             |           |             |           |
             |           |             |           |
             |           |            3|           |1
             |           |             |           |
             |           |             |           |
             @-----------@             +-----------+
            0             1                  0
*/
const int Quad4::Face2Node[4][2] = {{ 0, 1 },
                                    { 1, 2 },
                                    { 2, 3 },
                                    { 3, 0 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Quad4::Quad4 (int NDim)
    : GeomElem(NDim, /*NN*/4, /*NFN*/2, /*rCt*/0.0, /*sCt*/0.0, /*tCt*/0.0, "Quad4")
{
    SetIPs (4);
}

inline void Quad4::SetIPs (int TotNIP)
{
         if (TotNIP== 4) IPs = QUAD_IP2;
    else if (TotNIP== 9) IPs = QUAD_IP3;
    else if (TotNIP==16) IPs = QUAD_IP4;
    else if (TotNIP==25) IPs = QUAD_IP5;
    else throw new Fatal("Quad4::SetIPs: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = LIN_IP2;
    NFIP = 2;
}

inline void Quad4::Shape (double r, double s, double t) const
{
    /*      3                        2
     *        @--------------------@
     *        |               (1,1)|
     *        |       s ^          |
     *        |         |          |
     *        |         |          |
     *        |         +----> r   |
     *        |       (0,0)        |
     *        |                    |
     *        |                    |
     *        |(-1,-1)             |
     *        @--------------------@
     *      0                        1
     */
    N(0) = 0.25*(1.0-r-s+r*s);
    N(1) = 0.25*(1.0+r-s-r*s);
    N(2) = 0.25*(1.0+r+s+r*s);
    N(3) = 0.25*(1.0-r+s-r*s);
}

inline void Quad4::Derivs (double r, double s, double t) const
{
    /*           _     _ T
     *          |  dNi  |
     *   dNdR = |  ---  |   , where Rj = r, s
     *          |_ dRj _|
     *  
     *   dNdR(j,i), j=>local coordinate and i=>shape function
     */
    dNdR(0,0) = 0.25*(-1.0+s);   dNdR(1,0) = 0.25*(-1.0+r);
    dNdR(0,1) = 0.25*(+1.0-s);   dNdR(1,1) = 0.25*(-1.0-r);
    dNdR(0,2) = 0.25*(+1.0+s);   dNdR(1,2) = 0.25*(+1.0+r);
    dNdR(0,3) = 0.25*(-1.0-s);   dNdR(1,3) = 0.25*(+1.0-r);
}

inline void Quad4::FaceShape (double r, double s) const
{
    /*  
     *       0           |           1
     *       @-----------+-----------@-> r
     *      -1           |          +1
     */
    FN(0) = 0.5*(1.0-r);
    FN(1) = 0.5*(1.0+r);
}

inline void Quad4::FaceDerivs (double r, double s) const
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


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Quad4Maker (int NDim) { return new Quad4(NDim); }

// Register element
int Quad4Register ()
{
    GeomElemFactory["Quad4"] = Quad4Maker;
    GEOM.Set ("Quad4", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Quad4_dummy_int  = Quad4Register();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD4_H
