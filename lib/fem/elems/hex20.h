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

#ifndef MECHSYS_FEM_HEX20_H
#define MECHSYS_FEM_HEX20_H

// STL
#include <algorithm> 

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/geomelem.h>
#include <mechsys/vtkcelltype.h>

namespace FEM
{

class Hex20 : public GeomElem
{
public:
    // Auxiliar structure to map local face IDs to local node IDs
    static const int Face2Node[6][8]; // 6 faces, 8 nodes/face

    // Constructor
    Hex20 (int NDim);

    // Derived methods
    void   SetIPs     (int TotNIP);
    int    VTKType    () const { return VTK_QUADRATIC_HEXAHEDRON; }
    size_t FNode      (size_t IdxFace, size_t IdxFNode) const { return Face2Node[IdxFace][IdxFNode]; }
    void   Shape      (double r, double s, double t)    const;
    void   Derivs     (double r, double s, double t)    const;
    void   FaceShape  (double r, double s)              const;
    void   FaceDerivs (double r, double s)              const;
};


/* Local IDs
                  Vertices                               Faces
    t
    |           4        15        7    
   ,+--s         @-------@--------@                   +----------------+ 
 r'            ,'|              ,'|                 ,'|              ,'| 
          12 @'  |         14 ,'  |               ,'  |  ___       ,'  | 
           ,'    |16        ,@    |19           ,'    |,'5,'  [0],'    | 
     5   ,'      @      6 ,'      @           ,'      |~~~     ,'      | 
       @'=======@=======@'        |         +'===============+'  ,'|   | 
       |      13 |      |         |         |   ,'|   |      |   |3|   | 
       |         |      |  11     |         |   |2|   |      |   |,'   | 
    17 |       0 @- - - | @- - - -@         |   |,'   +- - - | +- - - -+ 
       @       ,'       @       ,' 3        |       ,'       |       ,'  
       |   8 @'      18 |     ,'            |     ,' [1]  ___|     ,'    
       |   ,'           |   ,@ 10           |   ,'      ,'4,'|   ,'      
       | ,'             | ,'                | ,'        ~~~  | ,'        
       @-------@--------@'                  +----------------+'          
     1         9         2
*/
const int Hex20::Face2Node[6][8] = {{ 0, 4, 7, 3, 16, 15, 19, 11 },
                                    { 1, 2, 6, 5,  9, 18, 13, 17 },
                                    { 0, 1, 5, 4,  8, 17, 12, 16 },
                                    { 2, 3, 7, 6, 10, 19, 14, 18 },
                                    { 0, 3, 2, 1, 11, 10,  9,  8 }, 
                                    { 4, 5, 6, 7, 12, 13, 14, 15 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Hex20::Hex20 (int NDim)
    : GeomElem(NDim, /*NN*/20, /*NFN*/8, /*rCt*/0.0, /*sCt*/0.0, /*tCt*/0.0, "Hex20")
{
    SetIPs (27);
}

inline void Hex20::SetIPs (int TotNIP)
{
         if (TotNIP==  8) IPs = HEX_IP2;
    else if (TotNIP== 27) IPs = HEX_IP3;
    else if (TotNIP== 64) IPs = HEX_IP4;
    else if (TotNIP==125) IPs = HEX_IP5;
    else throw new Fatal("Hex20::SetIPs: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = QUAD_IP2;
    NFIP = 4; 
}

inline void Hex20::Shape (double r, double s, double t) const
{
    double rp1=1.0+r; double rm1=1.0-r;
    double sp1=1.0+s; double sm1=1.0-s;
    double tp1=1.0+t; double tm1=1.0-t;

    N( 0) = 0.125*rm1*sm1*tm1*(-r-s-t-2);
    N( 1) = 0.125*rp1*sm1*tm1*( r-s-t-2);
    N( 2) = 0.125*rp1*sp1*tm1*( r+s-t-2);
    N( 3) = 0.125*rm1*sp1*tm1*(-r+s-t-2);
    N( 4) = 0.125*rm1*sm1*tp1*(-r-s+t-2);
    N( 5) = 0.125*rp1*sm1*tp1*( r-s+t-2);
    N( 6) = 0.125*rp1*sp1*tp1*( r+s+t-2);
    N( 7) = 0.125*rm1*sp1*tp1*(-r+s+t-2);
    N( 8) = 0.25*(1-r*r)*sm1*tm1;
    N( 9) = 0.25*rp1*(1-s*s)*tm1;
    N(10) = 0.25*(1-r*r)*sp1*tm1;
    N(11) = 0.25*rm1*(1-s*s)*tm1;
    N(12) = 0.25*(1-r*r)*sm1*tp1;
    N(13) = 0.25*rp1*(1-s*s)*tp1;
    N(14) = 0.25*(1-r*r)*sp1*tp1;
    N(15) = 0.25*rm1*(1-s*s)*tp1;
    N(16) = 0.25*rm1*sm1*(1-t*t);
    N(17) = 0.25*rp1*sm1*(1-t*t);
    N(18) = 0.25*rp1*sp1*(1-t*t);
    N(19) = 0.25*rm1*sp1*(1-t*t);
}

inline void Hex20::Derivs(double r, double s, double t) const
{
    /*       _     _ T
     *      |  dNi  |
     * dN = |  ---  |   , where cj = r, s
     *      |_ dcj _|
     *
     * dN(j,i), j=>local coordinate and i=>shape function
     */
    double rp1=1.0+r; double rm1=1.0-r;
    double sp1=1.0+s; double sm1=1.0-s;
    double tp1=1.0+t; double tm1=1.0-t;
    
    //Derivatives with respect to r
    dNdR(0, 0)= -0.125*sm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
    dNdR(0, 1)=  0.125*sm1*tm1*( r-s-t-2)+0.125*rp1*sm1*tm1;
    dNdR(0, 2)=  0.125*sp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1;
    dNdR(0, 3)= -0.125*sp1*tm1*(-r+s-t-2)-0.125*rm1*sp1*tm1;
    dNdR(0, 4)= -0.125*sm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1;
    dNdR(0, 5)=  0.125*sm1*tp1*( r-s+t-2)+0.125*rp1*sm1*tp1;
    dNdR(0, 6)=  0.125*sp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
    dNdR(0, 7)= -0.125*sp1*tp1*(-r+s+t-2)-0.125*rm1*sp1*tp1;
    dNdR(0, 8)= -0.5*r*sm1*tm1;
    dNdR(0, 9)=  0.25*(1-s*s)*tm1;
    dNdR(0,10)= -0.5*r*sp1*tm1;
    dNdR(0,11)= -0.25*(1-s*s)*tm1;
    dNdR(0,12)= -0.5*r*sm1*tp1;
    dNdR(0,13)=  0.25*(1-s*s)*tp1;
    dNdR(0,14)= -0.5*r*sp1  *tp1;
    dNdR(0,15)= -0.25*(1-s*s)*tp1;
    dNdR(0,16)= -0.25*sm1*(1-t*t);
    dNdR(0,17)=  0.25*sm1*(1-t*t);
    dNdR(0,18)=  0.25*sp1*(1-t*t);
    dNdR(0,19)= -0.25*sp1*(1-t*t);

    //Derivatives with respect to s
    dNdR(1, 0)= -0.125*rm1*tm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
    dNdR(1, 1)= -0.125*rp1*tm1*( r-s-t-2)-0.125*rp1*sm1*tm1;
    dNdR(1, 2)=  0.125*rp1*tm1*( r+s-t-2)+0.125*rp1*sp1*tm1;
    dNdR(1, 3)=  0.125*rm1*tm1*(-r+s-t-2)+0.125*rm1*sp1*tm1;
    dNdR(1, 4)= -0.125*rm1*tp1*(-r-s+t-2)-0.125*rm1*sm1*tp1;
    dNdR(1, 5)= -0.125*rp1*tp1*( r-s+t-2)-0.125*rp1*sm1*tp1;
    dNdR(1, 6)=  0.125*rp1*tp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
    dNdR(1, 7)=  0.125*rm1*tp1*(-r+s+t-2)+0.125*rm1*sp1*tp1;
    dNdR(1, 8)= -0.25*(1-r*r)*tm1;
    dNdR(1, 9)= -0.5*s*rp1*tm1;
    dNdR(1,10)=  0.25*(1-r*r)*tm1;
    dNdR(1,11)= -0.5*s*rm1*tm1;
    dNdR(1,12)= -0.25*(1-r*r)*tp1;
    dNdR(1,13)= -0.5*s*rp1*tp1;
    dNdR(1,14)=  0.25*(1-r*r)*tp1;
    dNdR(1,15)= -0.5*s*rm1*tp1;
    dNdR(1,16)= -0.25*rm1*(1-t*t);
    dNdR(1,17)= -0.25*rp1*(1-t*t);
    dNdR(1,18)=  0.25*rp1*(1-t*t);
    dNdR(1,19)=  0.25*rm1*(1-t*t);

    //Derivatives with respect to t
    dNdR(2, 0)= -0.125*rm1*sm1*(-r-s-t-2)-0.125*rm1*sm1*tm1;
    dNdR(2, 1)= -0.125*rp1*sm1*( r-s-t-2)-0.125*rp1*sm1*tm1;
    dNdR(2, 2)= -0.125*rp1*sp1*( r+s-t-2)-0.125*rp1*sp1*tm1;
    dNdR(2, 3)= -0.125*rm1*sp1*(-r+s-t-2)-0.125*rm1*sp1*tm1;
    dNdR(2, 4)=  0.125*rm1*sm1*(-r-s+t-2)+0.125*rm1*sm1*tp1;
    dNdR(2, 5)=  0.125*rp1*sm1*( r-s+t-2)+0.125*rp1*sm1*tp1;
    dNdR(2, 6)=  0.125*rp1*sp1*( r+s+t-2)+0.125*rp1*sp1*tp1;
    dNdR(2, 7)=  0.125*rm1*sp1*(-r+s+t-2)+0.125*rm1*sp1*tp1;
    dNdR(2, 8)= -0.25*(1-r*r)*sm1;
    dNdR(2, 9)= -0.25*rp1*(1-s*s);
    dNdR(2,10)= -0.25*(1-r*r)*sp1;
    dNdR(2,11)= -0.25*rm1*(1-s*s);
    dNdR(2,12)=  0.25*(1-r*r)*sm1;
    dNdR(2,13)=  0.25*rp1*(1-s*s);
    dNdR(2,14)=  0.25*(1-r*r)*sp1;
    dNdR(2,15)=  0.25*rm1*(1-s*s);
    dNdR(2,16)= -0.5*t*rm1*sm1;
    dNdR(2,17)= -0.5*t*rp1*sm1;
    dNdR(2,18)= -0.5*t*rp1*sp1;
    dNdR(2,19)= -0.5*t*rm1*sp1;
}

inline void Hex20::FaceShape(double r, double s) const
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

    FN(0) = 0.25*rm1*sm1*(rm1+sm1-3.0);
    FN(1) = 0.25*rp1*sm1*(rp1+sm1-3.0);
    FN(2) = 0.25*rp1*sp1*(rp1+sp1-3.0);
    FN(3) = 0.25*rm1*sp1*(rm1+sp1-3.0);
    FN(4) = 0.50*sm1*(1.0-r*r);
    FN(5) = 0.50*rp1*(1.0-s*s);
    FN(6) = 0.50*sp1*(1.0-r*r);
    FN(7) = 0.50*rm1*(1.0-s*s);
}

inline void Hex20::FaceDerivs(double r, double s) const
{
    /*          _     _ T
     *         |  dNi  |
     *   FdN = |  ---  |   , where cj = r, s
     *         |_ dcj _|
     *
     *   FdNdR(j,i), j=>local coordinate and i=>shape function
     */
    double rp1=1.0+r; double rm1=1.0-r;
    double sp1=1.0+s; double sm1=1.0-s;

    FdNdR(0,0) = - 0.25 * sm1 * (rm1 + rm1 + sm1 - 3.0);
    FdNdR(0,1) =   0.25 * sm1 * (rp1 + rp1 + sm1 - 3.0);
    FdNdR(0,2) =   0.25 * sp1 * (rp1 + rp1 + sp1 - 3.0);
    FdNdR(0,3) = - 0.25 * sp1 * (rm1 + rm1 + sp1 - 3.0);
    FdNdR(0,4) = - r * sm1;
    FdNdR(0,5) =   0.50 * (1.0 - s * s);
    FdNdR(0,6) = - r * sp1;
    FdNdR(0,7) = - 0.5 * (1.0 - s * s);

    FdNdR(1,0) = - 0.25 * rm1 * (sm1 + rm1 + sm1 - 3.0);
    FdNdR(1,1) = - 0.25 * rp1 * (sm1 + rp1 + sm1 - 3.0);
    FdNdR(1,2) =   0.25 * rp1 * (sp1 + rp1 + sp1 - 3.0);
    FdNdR(1,3) =   0.25 * rm1 * (sp1 + rm1 + sp1 - 3.0);
    FdNdR(1,4) = - 0.50 * (1.0 - r * r);
    FdNdR(1,5) = - s * rp1;
    FdNdR(1,6) =   0.50 * (1.0 - r * r);
    FdNdR(1,7) = - s * rm1;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Hex20Maker (int NDim) { return new Hex20(NDim); }

// Register element
int Hex20Register ()
{
    GeomElemFactory["Hex20"] = Hex20Maker;
    GEOM.Set ("Hex20", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Hex20_dummy_int  = Hex20Register();

}; // namespace FEM

#endif // MECHSYS_FEM_HEX20_H
