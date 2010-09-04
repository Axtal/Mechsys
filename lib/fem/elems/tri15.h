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

#ifndef MECHSYS_FEM_TRI15_H
#define MECHSYS_FEM_TRI15_H

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/geomelem.h>

namespace FEM
{

class Tri15: public GeomElem
{
public:
    // Auxiliar structure to map local face IDs to local node IDs
    static const int Face2Node[3][5]; // 3 edges, 5 nodes/edge

    // Constructor
    Tri15 (int NDim);

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
                   Nodes                                  Faces

   y                 2
   |                 @                                      @
   |                / \                                    / \
   +--x            /   \                                  /   \
               10 @     @ 9                              /     \
                 /       \                              /       \
                /   14    \                            /         \
             5 @     @     @ 4                      2 /           \ 1
              /             \                        /             \
             /               \                      /               \
         11 @     @     @     @ 8                  /                 \
           /     12     13     \                  /                   \
          /                     \                /                     \
         @-----@-----@-----@-----@              @-----------------------@
         0     6     3     7     1                          0
*/
const int Tri15::Face2Node[3][5]= {{ 0, 1, 3,  6,  7 },
                                   { 1, 2, 4,  8,  9 },
                                   { 2, 0, 5, 10, 11 }}; // order of nodes is important


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Tri15::Tri15 (int NDim)
    : GeomElem(NDim, /*NN*/15, /*NFN*/5, "Tri15")
{
    SetIPs (16);
}

inline void Tri15::SetIPs (int TotNIP)
{
         if (TotNIP==12)  IPs = TRI_IP12;
    else if (TotNIP==13)  IPs = TRI_IP12;
    else if (TotNIP==16)  IPs = TRI_IP16;
    else if (TotNIP==25)  IPs = TRI_IP25;
    else throw new Fatal("tri6::SetIntPoints: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = LIN_IP5;
    NFIP = 5;
}

inline void Tri15::Shape (double r, double s, double t) const
{
    /*    s
     *    ^
     *    |
     *  2
     *    @,(0,1)
     *    | ',
     *    |   ', 9
     * 10 @     @,
     *    |  14   ',   4
     *  5 @    @     @
     *    |           ',  8
     * 11 @  12@   @    '@
     *    |       13      ',
     *    |(0,0)            ', (1,0)
     *    @----@----@----@----@  --> r
     *  0      6    3    7     1
     */
    double pt1 = 42.666666666666667;
    double pt2 = 10.666666666666667;
    double cc  = 1.0 - r - s;
    double t1  = r  - 0.25;
    double t2  = r  - 0.5;
    double t3  = r -  0.75;
    double t4  = s  - 0.25;
    double t5  = s  - 0.5;
    double t6  = s -  0.75;
    double t7  = cc - 0.25;
    double t8  = cc - 0.5;
    double t9  = cc - 0.75;
    N(0)  = pt2*cc*t7*t8*t9;
    N(1)  = pt2*r*t1*t2*t3;
    N(2)  = pt2*s*t4*t5*t6;
    N(3)  = 64.0*cc*r*t1*t7;
    N(4)  = 64.0*r*s*t1*t4;
    N(5)  = 64.0*s*cc*t4*t7;
    N(6)  = pt1*cc*r*t7*t8;
    N(7)  = pt1*cc*r*t1*t2;
    N(8)  = pt1*r*s*t1*t2;
    N(9)  = pt1*r*s*t4*t5;
    N(10) = pt1*s*cc*t4*t5;
    N(11) = pt1*s*cc*t7*t8;
    N(12) = 128.0*r*s*cc*t7;
    N(13) = 128.0*r*s*t1*cc;
    N(14) = 128.0*r*s*cc*t4;
}

inline void Tri15::Derivs (double r, double s, double t) const
{
    /*           _     _ T
     *          |  dNi  |
     *   dNdR = |  ---  |   , where Rj = r, s
     *          |_ dRj _|
     *
     *   dNdR(j,i), j=>local coordinate and i=>shape function
     */
    double pt1 = 42.666666666666667;
    double pt2 = 10.666666666666667;
    double cc  = 1.0 - r - s;
    double t1  = r  - 0.25;
    double t2  = r  - 0.5;
    double t3  = r  - 0.75;
    double t4  = s  - 0.25;
    double t5  = s  - 0.5;
    double t6  = s  - 0.75;
    double t7  = cc - 0.25;
    double t8  = cc - 0.5;
    double t9  = cc - 0.75;

    dNdR(0, 0)  = - pt2 * (t8 * t9 * (t7 + cc) + cc * t7 * (t8 + t9) );
    dNdR(0, 1)  = pt2 * (t2 * t3 * (t1 + r) + r * t1 * (t3 + t2) );
    dNdR(0, 2)  = 0.0;
    dNdR(0, 3)  = 64.0 * (cc * t7 * (t1 + r) - r * t1 * (t7 + cc) );
    dNdR(0, 4)  = 64.0 * s * t4 * (t1 + r);
    dNdR(0, 5)  = - 64.0 * s * t4 * (t7 + cc);
    dNdR(0, 6)  = pt1 * (cc * t7 * t8 - r * (t8 * (t7 + cc) + cc * t7) );
    dNdR(0, 7)  = pt1 * (cc * (t2 * (t1 + r) + r * t1) - r * t1 * t2);
    dNdR(0, 8)  = pt1 * s * (t2 * (t1 + r) + r * t1);
    dNdR(0, 9)  = pt1 * s * t4 * t5;
    dNdR(0, 10) = - pt1 * s * t4 * t5;
    dNdR(0, 11) = - pt1 * s * (t8 * (t7 + cc) + cc * t7);
    dNdR(0, 12) = 128.0 * s * (cc * t7 - r * (t7 + cc) );
    dNdR(0, 13) = 128.0 * s * (cc * (t1 + r) - r * t1);
    dNdR(0, 14) = 128.0 * s * t4 * (cc - r);

    dNdR(1, 0)  = - pt2 * (t8 * t9 * (t7 + cc) + cc * t7 * (t8 + t9) );
    dNdR(1, 1)  = 0.0;
    dNdR(1, 2)  = pt2 * (t5 * t6 * (t4 + s) + s * t4 * (t6 + t5) );
    dNdR(1, 3)  = - 64.0 * r * t1 * (t7 + cc);
    dNdR(1, 4)  = 64.0 * r * t1 * (t4 + s);
    dNdR(1, 5)  = 64.0 * (cc * t7 * (t4 + s) - s * t4 * (t7 + cc) );
    dNdR(1, 6)  = - pt1 * r * (t8 * (t7 + cc) + cc * t7);
    dNdR(1, 7)  = - pt1 * r * t1 * t2;
    dNdR(1, 8)  = pt1 * r * t1 * t2;
    dNdR(1, 9)  = pt1 * r * (t5 * (t4 + s) + s * t4);
    dNdR(1, 10) = pt1 * ((cc * (t5 * (t4 + s) + s * t4)) - s * t4 * t5);
    dNdR(1, 11) = pt1 * (cc * t7 * t8 - s * (t8 * (t7 + cc) + cc * t7));
    dNdR(1, 12) = 128.0 * r * (cc * t7 - s * (cc + t7) );
    dNdR(1, 13) = 128.0 * r * t1 * (cc - s);
    dNdR(1, 14) = 128.0 * r * (cc * (t4 + s) - s * t4);
}

inline void Tri15::FaceShape (double r, double s) const
{
    /*
     *       @-----@-----@-----@-----@-> r
     *       0     3     2     4     1
     *       |           |           |
     *      r=-1  -1/2   r=0  1/2   r=+1
     */
    FN(0) = (r-1.)*(1.-2.*r)*r*(-1.-2.*r)/6.;
    FN(1) = (1.-2.*r)*r*(-1.-2.*r)*(1.+r)/6.;
    FN(2) = (1.-r)*(1.-2.*r)*(-1.-2.*r)*(-1.-r);
    FN(3) = 4.*(1.-r)*(1.-2.*r)*r*(-1.-r)/3.;
    FN(4) = 4.*(1.-r)*r*(-1.-2.*r)*(-1.-r)/3.;
}

inline void Tri15::FaceDerivs (double r, double s) const
{
    /*            _     _ T
     *           |  dNi  |
     *   FdNdR = |  ---  |   , where Rj = r, s
     *           |_ dRj _|
     *
     *   FdNdR(j,i), j=>local coordinate and i=>shape function
     */
    FdNdR(0,0) = -((1.-2.*r)*(r-1.)*r)/3.-((-2.*r-1.)*(r-1.)*r)/3.+((-2.*r-1.)*(1.-2.*r)*r)/6.+((-2.*r-1.)*(1.-2.*r)*(r-1.))/6.;
    FdNdR(0,1) = -((1.-2.*r)*r*(r+1.))/3.-((-2.*r-1.)*r*(r+1.))/3.+((-2.*r-1.)*(1.-2.*r)*(r+1.))/6.+((-2.*r-1.)*(1.-2.*r)*r)/6.;
    FdNdR(0,2) = -2.*(1.-2.*r)*(-r-1.)*(1.-r)-2.*(-2.*r-1.)*(-r-1.)*(1.-r)-(-2.*r-1.)*(1.-2.*r)*(1.-r)-(-2.*r-1.)*(1.-2.*r)*(-r-1.);
    FdNdR(0,3) = -(8.*(-r-1.)*(1.-r)*r)/3.-(4.*(1.-2.*r)*(1.-r)*r)/3.-(4.*(1.-2.*r)*(-r-1.)*r)/3.+(4.*(1.-2.*r)*(-r-1.)*(1.-r))/3.;
    FdNdR(0,4) = -(8.*(-r-1.)*(1.-r)*r)/3.-(4.*(-2.*r-1.)*(1.-r)*r)/3.-(4.*(-2.*r-1.)*(-r-1.)*r)/3.+(4.*(-2.*r-1.)*(-r-1.)*(1.-r))/3.;
}

inline void Tri15::NatCoords (Mat_t & C) const
{
    throw new Fatal("Tri15::NatCoords: method is not available yet");
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Tri15Maker (int NDim) { return new Tri15(NDim); }

// Register element
int Tri15Register ()
{
    GeomElemFactory["Tri15"] = Tri15Maker;
    GEOM.Set ("Tri15", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Tri15_dummy_int  = Tri15Register();

}; // namespace FEM

#endif // MECHSYS_FEM_TRI15_H
