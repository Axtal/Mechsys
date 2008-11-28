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

/*       F1        F2        F3    F1= 0.0 F2=F3 = 1.0
          ^         ^         ^
          |         |         |

        2 @---------@---------@ 6
          |',       7         |
          |  ',        e[1]   |
          |    ',             |
          |      ',  4        |
        5 @        '@         @ 8
          |          ',       |
          |   e[0]     ',     |
          |              ',   |
          |         3      ', |
        0 @---------@---------@ 1
         /_\       /_\       /_\
         ///       o o       o o
 */


// STL
#include <iostream>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/tri6pstrain.h"
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
	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Nodes
	g.SetNNodes (9);
	g.SetNode   (0, 0.0, 0.0);
	g.SetNode   (1, 1.0, 0.0);
	g.SetNode   (2, 0.0, 1.0);
	g.SetNode   (3, 0.5, 0.0);
	g.SetNode   (4, 0.5, 0.5);
	g.SetNode   (5, 0.0, 0.5);
	g.SetNode   (6, 1.0, 1.0);
	g.SetNode   (7, 0.5, 1.0);
	g.SetNode   (8, 1.0, 0.5);

	// Elements
	g.SetNElems (2);
	g.SetElem   (0, "Tri6PStrain", /*IsActive*/true, /*Tag*/-1);
	g.SetElem   (1, "Tri6PStrain", /*IsActive*/true, /*Tag*/-1);

	// Set connectivity
	g.Ele(0)->Connect(0, g.Nod(0))
	        ->Connect(1, g.Nod(1))
			->Connect(2, g.Nod(2))
			->Connect(3, g.Nod(3))
			->Connect(4, g.Nod(4))
			->Connect(5, g.Nod(5))->SetIntPoints(3);
	
	g.Ele(1)->Connect(0, g.Nod(6))
	        ->Connect(1, g.Nod(2))
	        ->Connect(2, g.Nod(1))
	        ->Connect(3, g.Nod(7))
	        ->Connect(4, g.Nod(4))
	        ->Connect(5, g.Nod(8))->SetIntPoints(3);

	// Boundary conditions (must be after connectivity)
	g.Nod(0)->Bry("ux",0.0)->Bry("uy",0.0);
	g.Nod(1)->Bry("uy",0.0);
	g.Nod(3)->Bry("uy",0.0);
	g.Nod(2)->Bry("fy",0.0);
	g.Nod(7)->Bry("fy",1.0);
	g.Nod(6)->Bry("fy",1.0);

	// Parameters and initial values
	g.Ele(0)->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");
	g.Ele(1)->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);
	delete sol;

	Output out;
	out.VTU(&g, "out.vtu");

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Check
    Array<double> err_eps;
    Array<double> err_sig;
    Array<double> err_dis;

	Vector<double> sx0 (6);
	Vector<double> sy0 (6);
	Vector<double> sxy0(6);
	Vector<double> sx1 (6);
	Vector<double> sy1 (6);
	Vector<double> sxy1(6);

	for (int i=0; i<6; i++)
	{
		sx0 (i) = g.Ele(0)->Val(i, "Sx" );
		sy0 (i) = g.Ele(0)->Val(i, "Sy" );
		sxy0(i) = g.Ele(0)->Val(i, "Sxy");
		sx1 (i) = g.Ele(1)->Val(i, "Sx" );
		sy1 (i) = g.Ele(1)->Val(i, "Sy" );
		sxy1(i) = g.Ele(1)->Val(i, "Sxy");
	}

	Vector<double> sx (9);
	Vector<double> sy (9);
	Vector<double> sxy(9);

	sx(0) = sx0(0);                 sy(0) = sy0(0);                 sxy(0) = sxy0(0);                 
	sx(1) = (sx0(1) + sx1(2))*0.5;  sy(1) = (sy0(1) + sy1(2))*0.5;  sxy(1) = (sxy0(1) + sxy1(2))*0.5;
	sx(2) = (sx0(2) + sx1(1))*0.5;  sy(2) = (sy0(2) + sy1(1))*0.5;  sxy(2) = (sxy0(2) + sxy1(1))*0.5;
	sx(3) = sx0(3);                 sy(3) = sy0(3);                 sxy(3) = sxy0(3);
	sx(4) = (sx0(4) + sx1(4))*0.5;  sy(4) = (sy0(4) + sy1(4))*0.5;  sxy(4) = (sxy0(4) + sxy1(4))*0.5;
	sx(5) = sx0(5);                 sy(5) = sy0(5);                 sxy(5) = sxy0(5);
	sx(6) = sx1(0);                 sy(6) = sy1(0);                 sxy(6) = sxy1(0);
	sx(7) = sx1(3);                 sy(7) = sy1(3);                 sxy(7) = sxy1(3);
	sx(8) = sx1(5);                 sy(8) = sy1(5);                 sxy(8) = sxy1(5);                 

	err_sig.Push( fabs(sx(0) - (-1.1241E-2)) );
	err_sig.Push( fabs(sx(1) - ( 3.3784E-1)) );
	err_sig.Push( fabs(sx(2) - (-8.1807E-2)) );
	err_sig.Push( fabs(sx(3) - ( 9.0341E-2)) );
	err_sig.Push( fabs(sx(4) - ( 1.2802E-1)) );
	err_sig.Push( fabs(sx(5) - (-8.5346E-3)) );
	err_sig.Push( fabs(sx(6) - (-5.0083E-1)) );
	err_sig.Push( fabs(sx(7) - (-3.2931E-1)) );
	err_sig.Push( fabs(sx(8) - (-8.5345E-3)) );

	err_sig.Push( fabs(sy(0) - ( -8.5054E-1 )) );
	err_sig.Push( fabs(sy(1) - (  4.8378)) );
	err_sig.Push( fabs(sy(2) - ( -1.1622)) );
	err_sig.Push( fabs(sy(3) - (  1.8464)) );
	err_sig.Push( fabs(sy(4) - (  1.8378)) );
	err_sig.Push( fabs(sy(5) - ( -6.8422E-1)) );
	err_sig.Push( fabs(sy(6) - (  5.4992)) );
	err_sig.Push( fabs(sy(7) - (  1.8464)) );
	err_sig.Push( fabs(sy(8) - (  5.3158)) );

	err_sig.Push( fabs(sxy(0) - (  1.1241E-2 )) );
	err_sig.Push( fabs(sxy(1) - (  8.1807E-2)) );
	err_sig.Push( fabs(sxy(2) - ( -3.3784E-1)) );
	err_sig.Push( fabs(sxy(3) - ( -8.5345E-3)) );
	err_sig.Push( fabs(sxy(4) - ( -1.2802E-1)) );
	err_sig.Push( fabs(sxy(5) - ( -7.3272E-2)) );
	err_sig.Push( fabs(sxy(6) - (  5.0083E-1)) );
	err_sig.Push( fabs(sxy(7) - ( -8.5346E-3)) );
	err_sig.Push( fabs(sxy(8) - (  3.4638E-1)) );


	// Error summary
	double tol_sig     = 1.0e-4;
	double min_err_sig = err_sig[err_sig.Min()];
	double max_err_sig = err_sig[err_sig.Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Sig" << _8s<<min_err_sig << _8s<<err_sig.Mean() << (max_err_sig>tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig << "[0m" << _8s<<err_sig.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_sig>tol_sig) return 1;
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
