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
#include "fem/elems/quad8pstrain.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"
#include "util/numstreams.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_8s;


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

	// 0) Problem dimension
	FEM::Data dat(2); // 2D

	// 1) Nodes
	dat.SetNNodes (8);
	dat.SetNode   (0,   0.0 ,         0.0);
	dat.SetNode   (1,     L ,         0.0);
	dat.SetNode   (2,     L ,          H2);
	dat.SetNode   (3,   0.0 ,          H1);
	dat.SetNode   (4, L/2.0 ,         0.0);
	dat.SetNode   (5,     L ,      H2/2.0);
	dat.SetNode   (6, L/2.0 , (H2+H1)/2.0);
	dat.SetNode   (7,   0.0 ,      H1/2.0);

	// 2) Elements
	dat.SetNElems (1);
	dat.SetElem   (0, "Quad8PStrain", /*IsActive*/true);

	// 3) Set connectivity (list of nodes must be LOCAL)
	dat.Ele(0)->Connect(0, dat.Nod(0))
	        ->Connect(1, dat.Nod(1))
	        ->Connect(2, dat.Nod(2))
	        ->Connect(3, dat.Nod(3))
	        ->Connect(4, dat.Nod(4))
	        ->Connect(5, dat.Nod(5))
	        ->Connect(6, dat.Nod(6))
	        ->Connect(7, dat.Nod(7));

	// 4) Boundary conditions (must be after connectivity)
	dat.Nod(0)->Bry     ("uy",0.0)->Bry("ux",0.0);
	dat.Nod(4)->Bry     ("uy",0.0);
	dat.Nod(1)->Bry     ("uy",0.0);
	dat.Ele(0)->EdgeBry ("Q",q, 3); // 3 => top edge

	// 5) Parameters and initial values
	String prms; prms.Printf("E=%f  nu=%f",E,nu);
	dat.Ele(0)->SetModel("LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;

	// 6) Solve
	FEM::Solver sol(dat,"tpstrain03");
	sol.SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);

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

	Output out;
	out.VTU(&dat, "out.vtu");

	// Return error flag
	if (err_ux>tol || err_uy>tol) return 1;
	else return 0;
}
catch (Exception * e) 
{
	e->Cout();
	if (e->IsFatal()) {delete e; exit(1);}
	delete e;
}
catch (char const * m)
{
	std::cout << "Fatal: " << m << std::endl;
	exit (1);
}
catch (...)
{
	std::cout << "Some exception (...) ocurred\n";
} 
