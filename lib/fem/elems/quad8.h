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

#ifndef MECHSYS_FEM_QUAD8_H
#define MECHSYS_FEM_QUAD8_H

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/geomelem.h>

namespace FEM
{

class Quad8: public GeomElem
{
public:
    // Auxiliar structure to map local face IDs to local node IDs
    static const int Face2Node[4][3]; // 4 edges, 3 nodes/edge

    // Constructor
    Quad8 (int NDim);

    // Derived methods
    void   SetIPs     (int TotNIP);
    size_t FNode      (size_t IdxFace, size_t IdxFNode) const { return Face2Node[IdxFace][IdxFNode]; }
    void   Shape      (double r, double s, double t)    const;
    void   Derivs     (double r, double s, double t)    const;
    void   FaceShape  (double r, double s)              const;
    void   FaceDerivs (double r, double s)              const;
};


/* Local IDs
                 Nodes                  Faces(edges)
   y
   |        3      6      2                  2
   +--x      @-----@-----@             +-----------+
             |           |             |           |
             |           |             |           |
           7 @           @ 5          3|           |1
             |           |             |           |
             |           |             |           |
             @-----@-----@             +-----------+
            0      4      1                  0
*/
const int Quad8::Face2Node[4][3] = {{ 0, 1, 4 },
                                    { 1, 2, 5 },
                                    { 2, 3, 6 },
                                    { 3, 0, 7 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Quad8::Quad8 (int NDim)
    : GeomElem(NDim, /*NN*/8, /*NFN*/3, /*rCt*/0.0, /*sCt*/0.0, /*tCt*/0.0, "Quad8")
{
    SetIPs (4);
}

inline void Quad8::SetIPs (int TotNIP)
{
         if (TotNIP== 4) IPs = QUAD_IP2;
    else if (TotNIP== 9) IPs = QUAD_IP3;
    else if (TotNIP==16) IPs = QUAD_IP4;
    else if (TotNIP==25) IPs = QUAD_IP5;
    else throw new Fatal("Quad8::SetIPs: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = LIN_IP3;
    NFIP = 3;
}

inline void Quad8::Shape (double r, double s, double t) const
{
    /*      3           6            2
     *        @---------@----------@
     *        |               (1,1)|
     *        |       s ^          |
     *        |         |          |
     *        |         |          |
     *      7 @         +----> r   @ 5
     *        |       (0,0)        |
     *        |                    |
     *        |                    |
     *        |(-1,-1)             |
     *        @---------@----------@
     *      0           4            1
     */
    double rp1=1.0+r; double rm1=1.0-r;
    double sp1=1.0+s; double sm1=1.0-s;

    N(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
    N(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
    N(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
    N(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
    N(4) = 0.50*sm1*(1.0-r*r);
    N(5) = 0.50*rp1*(1.0-s*s);
    N(6) = 0.50*sp1*(1.0-r*r);
    N(7) = 0.50*rm1*(1.0-s*s);
}

inline void Quad8::Derivs (double r, double s, double t) const
{
    /*         _     _ T
     *        |  dNi  |
     * dNdR = |  ---  |   , where Rj = r, s
     *        |_ dRj _|
     *
     * dNdR(j,i), j=>local coordinate and i=>shape function
     */
    double rp1=1.0+r; double rm1=1.0-r;
    double sp1=1.0+s; double sm1=1.0-s;

    dNdR(0,0) = - 0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0);
    dNdR(0,1) =   0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0);
    dNdR(0,2) =   0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0);
    dNdR(0,3) = - 0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0);
    dNdR(0,4) = - r * sm1;
    dNdR(0,5) =   0.50 * (1.0 - s * s);
    dNdR(0,6) = - r * sp1;
    dNdR(0,7) = - 0.5 * (1.0 - s * s);

    dNdR(1,0) = - 0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0);
    dNdR(1,1) = - 0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0);
    dNdR(1,2) =   0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0);
    dNdR(1,3) =   0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0);
    dNdR(1,4) = - 0.50 * (1.0 - r * r);
    dNdR(1,5) = - s * rp1;
    dNdR(1,6) =   0.50 * (1.0 - r * r);
    dNdR(1,7) = - s * rm1;
}

inline void Quad8::FaceShape (double r, double s) const
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

inline void Quad8::FaceDerivs (double r, double s) const
{
    /*            _     _ T
     *           |  dNi  |
     *   FdNdR = |  ---  |   , where Rj = r, s
     *           |_ dRj _|
     *
     *   FdNdR(j,i), j=>local coordinate and i=>shape function
     */
    FdNdR(0,0) =  r  - 0.5;
    FdNdR(0,1) =  r  + 0.5;
    FdNdR(0,2) = -2.0* r;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Quad8Maker (int NDim) { return new Quad8(NDim); }

// Register element
int Quad8Register ()
{
    GeomElemFactory["Quad8"] = Quad8Maker;
    GEOM.Set ("Quad8", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Quad8_dummy_int  = Quad8Register();

}; // namespace FEM

#endif // MECHSYS_FEM_QUAD8_H
