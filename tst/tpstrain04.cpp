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
#include "fem/elems/hex20equilib.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "fem/output.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_8s;


int main(int argc, char **argv) try
{
	// Description:
	// Test of normal boundary conditions (traction) applied normal to a face.


	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);

	double H1 = 0.5;    // height
	double H2 = 1.5;    // height
	double L  = 1.0;    // length
	double T  = 1.0;    // thickness
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
           T ',      ',|         8       ',|      |        ',|
               ', ,  0 @---------@---------@ 1   ---         '-------> x
                ,'    /_\       /_\       /_\
                      ///       o o       o o
               
                       |-------- L --------|
	 */

	// Problem dimension
	FEM::Data dat(3); // 3D

	// Nodes
	dat.SetNNodes (20);
	dat.SetNode   ( 0,   0.0 ,      0.0,         0.0);
	dat.SetNode   ( 1,     L ,      0.0,         0.0);
	dat.SetNode   ( 2,     L ,        T,         0.0);
	dat.SetNode   ( 3,   0.0 ,        T,         0.0);
	dat.SetNode   ( 4,   0.0 ,      0.0,          H1);
	dat.SetNode   ( 5,     L ,      0.0,          H2);
	dat.SetNode   ( 6,     L ,        T,          H2);
	dat.SetNode   ( 7,   0.0 ,        T,          H1);
	dat.SetNode   ( 8, L/2.0 ,      0.0,         0.0);
	dat.SetNode   ( 9,     L ,    T/2.0,         0.0);
	dat.SetNode   (10, L/2.0 ,        T,         0.0);
	dat.SetNode   (11,   0.0 ,    T/2.0,         0.0);
	dat.SetNode   (12, L/2.0 ,      0.0, (H1+H2)/2.0);
	dat.SetNode   (13,     L ,    T/2.0,          H2);
	dat.SetNode   (14, L/2.0 ,        T, (H1+H2)/2.0);
	dat.SetNode   (15,   0.0 ,    T/2.0,          H1);
	dat.SetNode   (16,   0.0 ,      0.0,      H1/2.0);
	dat.SetNode   (17,     L ,      0.0,      H2/2.0);
	dat.SetNode   (18,     L ,        T,      H2/2.0);
	dat.SetNode   (19,   0.0 ,        T,      H1/2.0);

	// Elements
	dat.SetNElems (1);
	dat.SetElem   (0, "Hex20Equilib", /*IsActive*/true, /*Tag*/-1);

	// Set connectivity (list of nodes must be LOCAL)
	dat.Ele(0)->Connect(0, dat.Nod( 0))->Connect(1, dat.Nod( 1))->Connect(2, dat.Nod( 2))->Connect(3, dat.Nod( 3))
	        ->Connect(4, dat.Nod( 4))->Connect(5, dat.Nod( 5))->Connect(6, dat.Nod( 6))->Connect(7, dat.Nod( 7))
	        ->Connect(8, dat.Nod( 8))->Connect(9, dat.Nod( 9))->Connect(10,dat.Nod(10))->Connect(11,dat.Nod(11))
	        ->Connect(12,dat.Nod(12))->Connect(13,dat.Nod(13))->Connect(14,dat.Nod(14))->Connect(15,dat.Nod(15))
	        ->Connect(16,dat.Nod(16))->Connect(17,dat.Nod(17))->Connect(18,dat.Nod(18))->Connect(19,dat.Nod(19));

	// Parameters and initial values
	String prms; prms.Printf("E=%f  nu=%f",E,nu);
	dat.Ele(0)->SetModel("LinElastic", prms.CStr(), "ZERO");

	// 4) Boundary conditions (must be after connectivity)
	dat.Nod(0)->Bry    ("uy",0.0)->Bry("ux",0.0);
	dat.Nod(3)->Bry    ("uy",0.0)->Bry("ux",0.0);
	dat.Nod(11)->Bry   ("uy",0.0)->Bry("ux",0.0);

	dat.Ele(0)->FaceBry("uy",0.0,2);
	dat.Ele(0)->FaceBry("uy",0.0,3);
	dat.Ele(0)->FaceBry("uz",0.0,4);
	dat.Ele(0)->FaceBry( "Q",  q,5); // 5 => top face

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&dat)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);
	delete sol;

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
