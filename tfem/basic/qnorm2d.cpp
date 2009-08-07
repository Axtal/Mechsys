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
#include "fem/elems/quad8.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "util/fatal.h"
#include "util/numstreams.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_8s;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	double H1 = 0.5;    // height
	double H2 = 1.5;    // height
	double L  = 1.0;    // length
	double E  = 2000.0; // Young
	double nu = 0.25;   // Poisson
	double q  = -1.0;    // Load

	/*     
			                      @ 2  ---
                                ,'|     |
                              ,'  |     |
	        Diagonal        ,'    |     |
	          load        ,'      |     | 
	      on the face   @         |     | 
	                  ,'6         |     | 
	                ,'            @ 5   H2  
	              ,'              |     |
	            ,'     e[0]       |     | 
	  ---   3 @                   |     | 
	   |      |                   |     | 
	  H1    7 @                   |     | 
	   |      |         4         |     | 
	  ---   0 @---------@---------@ 1  ---
	         /_\       /_\       /_\
	         ///       o o       o o
	
	          |-------- L --------|
	 */

	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	Mesh::Generic mesh(/*Is3D*/false);
	mesh.SetNVerts  (8);
	mesh.SetNElems  (1);
	mesh.SetVert    (0, true,   0.0 ,         0.0);
	mesh.SetVert    (1, true,     L ,         0.0);
	mesh.SetVert    (2, true,     L ,          H2);
	mesh.SetVert    (3, true,   0.0 ,          H1);
	mesh.SetVert    (4, true, L/2.0 ,         0.0);
	mesh.SetVert    (5, true,     L ,      H2/2.0);
	mesh.SetVert    (6, true, L/2.0 , (H2+H1)/2.0);
	mesh.SetVert    (7, true,   0.0 ,      H1/2.0);
	mesh.SetElem    (0, -1, true, VTK_QUADRATIC_QUAD);
	mesh.SetElemCon (0, 0, 0);  mesh.SetElemCon(0, 1, 1);  mesh.SetElemCon(0, 2, 2);  mesh.SetElemCon(0, 3, 3);
	mesh.SetElemCon (0, 4, 4);  mesh.SetElemCon(0, 5, 5);  mesh.SetElemCon(0, 6, 6);  mesh.SetElemCon(0, 7, 7);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and solver
	FEM::Data   dat (2); // 2D
	FEM::Solver sol (dat,"tpstrain03");

	// Elements attributes
	FEM::EAtts_T eatts(1);
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	eatts = T(-1, "Quad8", "PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 -----------------------------------------------------------
	dat.Nod(0)->Bry     ("uy",0.0)->Bry("ux",0.0);
	dat.Nod(4)->Bry     ("uy",0.0);
	dat.Nod(1)->Bry     ("uy",0.0);
	dat.Ele(0)->EdgeBry ("Q",q, 3); // 3 => top edge
	sol.SolveWithInfo   (/*NDiv*/1, /*DTime*/0.0);

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Error summary
	double err_ux = 0.0;
	double err_uy = 0.0;
	double tol    = 1E-7;

	err_ux += fabs( dat.Ele(0)->Val(3,"ux") - ( 2.4422E-3) );
	err_ux += fabs( dat.Ele(0)->Val(6,"ux") - ( 4.7371E-3) );
	err_ux += fabs( dat.Ele(0)->Val(2,"ux") - ( 9.2043E-3) );

	err_uy += fabs( dat.Ele(0)->Val(3,"uy") - ( 2.8685E-4) );
	err_uy += fabs( dat.Ele(0)->Val(6,"uy") - (-1.0436E-4) );
	err_uy += fabs( dat.Ele(0)->Val(2,"uy") - (-3.1506E-3) );

	cout << _8s << "Err x" << _8s << "Err y" << endl;
	cout << _8s << (err_ux>tol?"[1;31m":"[1;32m") << _8s <<err_ux << "[0m";
	cout << _8s << (err_uy>tol?"[1;31m":"[1;32m") << _8s <<err_uy << "[0m" << endl;
	cout << endl;

	// Return error flag
	if (err_ux>tol || err_uy>tol) return 1;
	else return 0;
}
MECHSYS_CATCH
