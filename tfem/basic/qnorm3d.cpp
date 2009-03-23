/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/hex20.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/mesh.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_8s;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// Description:
	// Test of normal boundary conditions (traction) applied normal to a face.

	double H1 = 0.5;    // height
	double H2 = 1.5;    // height
	double L  = 1.0;    // length
	double th = 1.0;    // thickness
	double E  = 2000.0; // Young
	double nu = 0.25;   // Poisson
	double q  = -1.0;   // Load

	/*     
                                 6
                                ,@
                              ,' |',
                            ,'      '@13
                          ,'     |    ',
                     14 ,'              ',
                      ,@         |       `@ 5    ---
                    ,'         18@       ,'|      |
                  ,'                   ,'  |      |
                ,'    Normal     |   ,'    |      |
            7 ,'     pressure      ,'      |      | 
      ---    @`,     on face     @         |      | 
       |     |  `,             ,'12        |      | 
      H1   19@    @,15       ,'            @17    H2
       |     |      `,     ,'    |         |      |
      ---   3@,- - -  '4 ,'@ - - -@2       |      | 
          ,    ',      @   10      ',  9   |      |          z
        ,',    11'@    |             '@    |      |    y,    |
           ',      ',  @16             ',  |      |      ',  |
          th ',      ',|         8       ',|      |        ',|
               ', ,  0 @---------@---------@ 1   ---         '-------> x
                ,'    /_\       /_\       /_\
                      ///       o o       o o
               
                       |-------- L --------|
	 */

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	Mesh::Generic mesh(/*Is3D*/true);
	mesh.SetNVerts  (20);
	mesh.SetNElems  (1);
	mesh.SetVert    ( 0, true,   0.0 ,      0.0,         0.0);
	mesh.SetVert    ( 1, true,     L ,      0.0,         0.0);
	mesh.SetVert    ( 2, true,     L ,       th,         0.0);
	mesh.SetVert    ( 3, true,   0.0 ,       th,         0.0);
	mesh.SetVert    ( 4, true,   0.0 ,      0.0,          H1);
	mesh.SetVert    ( 5, true,     L ,      0.0,          H2);
	mesh.SetVert    ( 6, true,     L ,       th,          H2);
	mesh.SetVert    ( 7, true,   0.0 ,       th,          H1);
	mesh.SetVert    ( 8, true, L/2.0 ,      0.0,         0.0);
	mesh.SetVert    ( 9, true,     L ,   th/2.0,         0.0);
	mesh.SetVert    (10, true, L/2.0 ,       th,         0.0);
	mesh.SetVert    (11, true,   0.0 ,   th/2.0,         0.0);
	mesh.SetVert    (12, true, L/2.0 ,      0.0, (H1+H2)/2.0);
	mesh.SetVert    (13, true,     L ,   th/2.0,          H2);
	mesh.SetVert    (14, true, L/2.0 ,       th, (H1+H2)/2.0);
	mesh.SetVert    (15, true,   0.0 ,   th/2.0,          H1);
	mesh.SetVert    (16, true,   0.0 ,      0.0,      H1/2.0);
	mesh.SetVert    (17, true,     L ,      0.0,      H2/2.0);
	mesh.SetVert    (18, true,     L ,       th,      H2/2.0);
	mesh.SetVert    (19, true,   0.0 ,       th,      H1/2.0);
	mesh.SetElem    (0, -1, true, VTK_QUADRATIC_HEXAHEDRON);
	mesh.SetElemCon (0,  0,  0);  mesh.SetElemCon(0,  1,  1);  mesh.SetElemCon(0,  2,  2);  mesh.SetElemCon(0,  3,  3);
	mesh.SetElemCon (0,  4,  4);  mesh.SetElemCon(0,  5,  5);  mesh.SetElemCon(0,  6,  6);  mesh.SetElemCon(0,  7,  7);
	mesh.SetElemCon (0,  8,  8);  mesh.SetElemCon(0,  9,  9);  mesh.SetElemCon(0, 10, 10);  mesh.SetElemCon(0, 11, 11);
	mesh.SetElemCon (0, 12, 12);  mesh.SetElemCon(0, 13, 13);  mesh.SetElemCon(0, 14, 14);  mesh.SetElemCon(0, 15, 15);
	mesh.SetElemCon (0, 16, 16);  mesh.SetElemCon(0, 17, 17);  mesh.SetElemCon(0, 18, 18);  mesh.SetElemCon(0, 19, 19);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and solver
	FEM::Data   dat (3);
	FEM::Solver sol (dat,"tpstrain04");

	// Elements attributes
	FEM::EAtts_T eatts(1);
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	eatts = T(-1, "Hex20", "Equilib", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 -----------------------------------------------------------
	dat.Nod(0)->Bry    ("uy",0.0)->Bry("ux",0.0);
	dat.Nod(3)->Bry    ("uy",0.0)->Bry("ux",0.0);
	dat.Nod(11)->Bry   ("uy",0.0)->Bry("ux",0.0);
	dat.Ele(0)->FaceBry("uy",0.0,2);
	dat.Ele(0)->FaceBry("uy",0.0,3);
	dat.Ele(0)->FaceBry("uz",0.0,4);
	dat.Ele(0)->FaceBry( "Q",  q,5); // 5 => top face
	sol.SolveWithInfo  (/*NDiv*/1, /*DTime*/0.0);

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Error summary
	double err_ux = 0.0;
	double err_uy = 0.0;
	double tol    = 1E-7;

	err_ux += fabs( dat.Ele(0)->Val( 4,"ux") - ( 2.4422E-3) );
	err_ux += fabs( dat.Ele(0)->Val(12,"ux") - ( 4.7371E-3) );
	err_ux += fabs( dat.Ele(0)->Val( 5,"ux") - ( 9.2043E-3) );

	err_uy += fabs( dat.Ele(0)->Val( 4,"uz") - ( 2.8685E-4) );
	err_uy += fabs( dat.Ele(0)->Val(12,"uz") - (-1.0436E-4) );
	err_uy += fabs( dat.Ele(0)->Val( 5,"uz") - (-3.1506E-3) );

	cout << _8s << "Err x" << _8s << "Err z" << endl;
	cout << _8s << (err_ux>tol?"[1;31m":"[1;32m") << _8s <<err_ux << "[0m";
	cout << _8s << (err_uy>tol?"[1;31m":"[1;32m") << _8s <<err_uy << "[0m" << endl;
	cout << endl;

	// Return error flag
	if (err_ux>tol || err_uy>tol) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
