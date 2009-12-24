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

#ifndef MECHSYS_FEM_HEX8_H
#define MECHSYS_FEM_HEX8_H

// MechSys
#include <mechsys/fem/node.h>
#include <mechsys/fem/geomelem.h>
#include <mechsys/vtkcelltype.h>

namespace FEM
{

class Hex8 : public GeomElem
{
public:
    // Auxiliar structure to map local face IDs to local node IDs
    static const int Face2Node[6][4]; // 6 faces, 4 nodes/face

    // Constructor
    Hex8 (int NDim);

    // Derived methods
    void   SetIPs     (int TotNIP);
    int    VTKType    () const { return VTK_HEXAHEDRON; }
    size_t FNode      (size_t IdxFace, size_t IdxFNode) const { return Face2Node[IdxFace][IdxFNode]; }
    void   Shape      (double r, double s, double t)    const;
    void   Derivs     (double r, double s, double t)    const;
    void   FaceShape  (double r, double s)              const;
    void   FaceDerivs (double r, double s)              const;
};


/* Local IDs
                 Nodes                                   Faces
    z
    |           4                  7
   ,+--y         @________________@                    +________________+ 
 x'            ,'|              ,'|                  ,'|              ,'| 
             ,'  |            ,'  |                ,'  |  ___       ,'  | 
           ,'    |          ,'    |              ,'    |,'5,'  [0],'    | 
     5   ,'      |      6 ,'      |            ,'      |~~~     ,'      | 
       @'===============@'        |          +'===============+'  ,'|   | 
       |         |      |         |          |   ,'|   |      |   |3|   | 
       |         |      |         |          |   |2|   |      |   |,'   | 
       |       0 @______|_________@          |   |,'   +______|_________+ 
       |       ,'       |       ,' 3         |       ,'       |       ,'  
       |     ,'         |     ,'             |     ,' [1]  ___|     ,'    
       |   ,'           |   ,'               |   ,'      ,'4,'|   ,'      
       | ,'             | ,'                 | ,'        ~~~  | ,'        
       @________________@'                   +________________+'          
     1                   2
*/
const int Hex8::Face2Node[6][4] = {{ 0, 4, 7, 3 },
                                   { 1, 2, 6, 5 },
                                   { 0, 1, 5, 4 },
                                   { 2, 3, 7, 6 },
                                   { 0, 3, 2, 1 }, 
                                   { 4, 5, 6, 7 }};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Hex8::Hex8 (int NDim)
    : GeomElem(NDim, /*NN*/8, /*NFN*/4, /*rCt*/0.0, /*sCt*/0.0, /*tCt*/0.0, "Hex8")
{
    SetIPs (8);
}

inline void Hex8::SetIPs (int TotNIP)
{
         if (TotNIP==  8) IPs = HEX_IP2;
    else if (TotNIP== 27) IPs = HEX_IP3;
    else if (TotNIP== 64) IPs = HEX_IP4;
    else if (TotNIP==125) IPs = HEX_IP5;
    else throw new Fatal("Hex8::SetIPs: Total number of integration points = %d is invalid",TotNIP);

    NIP  = TotNIP;
    FIPs = QUAD_IP2;
    NFIP = 4; 
}

inline void Hex8::Shape (double r, double s, double t) const
{
    /*                    t
     *                    ^
     *                    |     
     *                   4                7
     *                    @________________@
     *                  ,'|              ,'|
     *                ,'  |            ,'  |
     *              ,'    |          ,'    |
     *        5   ,'      |     6  ,'      |
     *          @'===============@'        |
     *          |         |      |         |
     *          |         |      |         |            
     *          |       0 @_____ | ________@ 3 --> s
     *          |       ,'       |       ,' 
     *          |     ,'         |     ,' 
     *          |   ,'           |   ,' 
     *          | ,'             | ,' 
     *          @________________@'
     *         1                  2 
     *      ,'
     *    |_
     *   r
     */
    N(0) = 0.125*(1.0-r-s+r*s-t+s*t+r*t-r*s*t);
    N(1) = 0.125*(1.0+r-s-r*s-t+s*t-r*t+r*s*t);
    N(2) = 0.125*(1.0+r+s+r*s-t-s*t-r*t-r*s*t);
    N(3) = 0.125*(1.0-r+s-r*s-t-s*t+r*t+r*s*t);
    N(4) = 0.125*(1.0-r-s+r*s+t-s*t-r*t+r*s*t);
    N(5) = 0.125*(1.0+r-s-r*s+t-s*t+r*t-r*s*t);
    N(6) = 0.125*(1.0+r+s+r*s+t+s*t+r*t+r*s*t);
    N(7) = 0.125*(1.0-r+s-r*s+t+s*t-r*t-r*s*t);
}

inline void Hex8::Derivs (double r, double s, double t) const
{
    /*       _     _ T
     *      |  dNi  |
     * dN = |  ---  |   , where cj = r, s
     *      |_ dcj _|
     *
     * dN(j,i), j=>local coordinate and i=>shape function
     */
    dNdR(0,0) = 0.125*(-1.0+s+t-s*t);   dNdR(1,0)=0.125*(-1.0+r+t-r*t);   dNdR(2,0)=0.125*(-1.0+r+s-r*s);
    dNdR(0,1) = 0.125*(+1.0-s-t+s*t);   dNdR(1,1)=0.125*(-1.0-r+t+r*t);   dNdR(2,1)=0.125*(-1.0-r+s+r*s);
    dNdR(0,2) = 0.125*(+1.0+s-t-s*t);   dNdR(1,2)=0.125*(+1.0+r-t-r*t);   dNdR(2,2)=0.125*(-1.0-r-s-r*s);
    dNdR(0,3) = 0.125*(-1.0-s+t+s*t);   dNdR(1,3)=0.125*(+1.0-r-t+r*t);   dNdR(2,3)=0.125*(-1.0+r-s+r*s);
    dNdR(0,4) = 0.125*(-1.0+s-t+s*t);   dNdR(1,4)=0.125*(-1.0+r-t+r*t);   dNdR(2,4)=0.125*(+1.0-r-s+r*s);
    dNdR(0,5) = 0.125*(+1.0-s+t-s*t);   dNdR(1,5)=0.125*(-1.0-r-t-r*t);   dNdR(2,5)=0.125*(+1.0+r-s-r*s);
    dNdR(0,6) = 0.125*(+1.0+s+t+s*t);   dNdR(1,6)=0.125*(+1.0+r+t+r*t);   dNdR(2,6)=0.125*(+1.0+r+s+r*s);
    dNdR(0,7) = 0.125*(-1.0-s-t-s*t);   dNdR(1,7)=0.125*(+1.0-r+t-r*t);   dNdR(2,7)=0.125*(+1.0-r+s-r*s);
}

inline void Hex8::FaceShape (double r, double s) const
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
    FN(0) = 0.25*(1.0-r-s+r*s);
    FN(1) = 0.25*(1.0+r-s-r*s);
    FN(2) = 0.25*(1.0+r+s+r*s);
    FN(3) = 0.25*(1.0-r+s-r*s);
}

inline void Hex8::FaceDerivs (double r, double s) const
{
    /*            _     _ T
     *           |  dNi  |
     *   FdNdR = |  ---  |   , where Rj = r, s
     *           |_ dRj _|
     *
     *   FdNdR(j,i), j=>local coordinate and i=>shape function
     */

    FdNdR(0,0) = 0.25*(-1.0+s);   FdNdR(1,0) = 0.25*(-1.0+r);
    FdNdR(0,1) = 0.25*(+1.0-s);   FdNdR(1,1) = 0.25*(-1.0-r);
    FdNdR(0,2) = 0.25*(+1.0+s);   FdNdR(1,2) = 0.25*(+1.0+r);
    FdNdR(0,3) = 0.25*(-1.0-s);   FdNdR(1,3) = 0.25*(+1.0-r);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


// Allocate a new element
GeomElem * Hex8Maker (int NDim) { return new Hex8(NDim); }

// Register element
int Hex8Register ()
{
    GeomElemFactory["Hex8"] = Hex8Maker;
    GEOM.Set ("Hex8", (double)GEOM.Keys.Size());
    return 0;
}

// Call register
int __Hex8_dummy_int  = Hex8Register();

}; // namespace FEM

#endif // MECHSYS_FEM_HEX8_H
