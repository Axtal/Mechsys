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

/*       F1        F2        F3    F1=F2=F3 = 1.0
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
	g.SetElem   (0, "Tri6PStrain", /*IsActive*/true);
	g.SetElem   (1, "Tri6PStrain", /*IsActive*/true);

	// Set connectivity
	g.Ele(0)->Connect(0, g.Nod(0))
	        ->Connect(1, g.Nod(1))
			->Connect(2, g.Nod(2))
			->Connect(3, g.Nod(3))
			->Connect(4, g.Nod(4))
			->Connect(5, g.Nod(5));
	g.Ele(1)->Connect(0, g.Nod(6))
	        ->Connect(1, g.Nod(2))
	        ->Connect(2, g.Nod(1))
	        ->Connect(3, g.Nod(7))
	        ->Connect(4, g.Nod(4))
	        ->Connect(5, g.Nod(8));

	// Boundary conditions (must be after connectivity)
	g.Nod(0)->Bry("ux",0.0)->Bry("uy",0.0);
	g.Nod(1)->Bry("uy",0.0);
	g.Nod(3)->Bry("uy",0.0);
	g.Nod(2)->Bry("fy",1.0);
	g.Nod(7)->Bry("fy",1.0);
	g.Nod(6)->Bry("fy",1.0);

	// Parameters and initial values
	g.Ele(0)->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");
	g.Ele(1)->SetModel("LinElastic", "E=10000.0 nu=0.25", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Check
    Array<double> err_eps;
    Array<double> err_sig;
    Array<double> err_dis;

	// Element 0
	err_sig.Push( fabs(g.Ele(0)->Val(0, "Sx") - ( 1.56432140e-01)) );
	/*
	err_sig.Push( fabs(g.Ele(0)->Val(1, "Sx") - (-3.00686928e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(2, "Sx") - ( 1.44254788e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(3, "Sx") - (-3.19109076e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(4, "Sx") - (-3.31286428e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(5, "Sx") - ( 1.25832639e-01)) );

	err_sig.Push( fabs(g.Ele(0)->Val(0, "Sy") - (-2.05141549e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(1, "Sy") - ( 1.15872190e+00)) );
	err_sig.Push( fabs(g.Ele(0)->Val(2, "Sy") - (-9.53580350e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(3, "Sy") - (-2.22127394e+00)) );
	err_sig.Push( fabs(g.Ele(0)->Val(4, "Sy") - (-2.96971274e+00)) );
	err_sig.Push( fabs(g.Ele(0)->Val(5, "Sy") - (-4.33357619e+00)) );

	err_sig.Push( fabs(g.Ele(0)->Val(0, "Sxy") - (-1.56432140e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(1, "Sxy") - (-6.74437968e-02)) );
	err_sig.Push( fabs(g.Ele(0)->Val(2, "Sxy") - ( 2.23875937e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(3, "Sxy") - (-4.90216486e-02)) );
	err_sig.Push( fabs(g.Ele(0)->Val(4, "Sxy") - ( 3.31286428e-01)) );
	err_sig.Push( fabs(g.Ele(0)->Val(5, "Sxy") - ( 2.42298085e-01)) );

	// Element 1
	err_sig.Push( fabs(g.Ele(1)->Val(0, "Sx")  - ( 9.95732723e-01)) );
	err_sig.Push( fabs(g.Ele(1)->Val(1, "Sx")  - ( 2.23875937e-01)) );
	err_sig.Push( fabs(g.Ele(1)->Val(2, "Sx")  - (-1.21960866e+00)) );
	err_sig.Push( fabs(g.Ele(1)->Val(3, "Sx")  - ( 1.39446295e+00)) );
	err_sig.Push( fabs(g.Ele(1)->Val(4, "Sx")  - (-8.20878435e-01)) );
	err_sig.Push( fabs(g.Ele(1)->Val(5, "Sx")  - (-4.90216486e-02)) );

	err_sig.Push( fabs(g.Ele(1)->Val(0, "Sy")  - (-1.25426728e+00)) );
	err_sig.Push( fabs(g.Ele(1)->Val(1, "Sy")  - ( 1.68328476e+00)) );
	err_sig.Push( fabs(g.Ele(1)->Val(2, "Sy")  - (-4.29017485e-01)) );
	err_sig.Push( fabs(g.Ele(1)->Val(3, "Sy")  - (-2.39612823e+00)) );
	err_sig.Push( fabs(g.Ele(1)->Val(4, "Sy")  - (-1.57087843e+00)) );
	err_sig.Push( fabs(g.Ele(1)->Val(5, "Sy")  - (-4.50843047e+00)) );

	err_sig.Push( fabs(g.Ele(1)->Val(0, "Sxy") - (-9.95732723e-01)) );
	err_sig.Push( fabs(g.Ele(1)->Val(1, "Sxy") - ( 1.29641965e+00)) );
	err_sig.Push( fabs(g.Ele(1)->Val(2, "Sxy") - (-3.00686928e-01)) );
	err_sig.Push( fabs(g.Ele(1)->Val(3, "Sxy") - ( 1.25832639e-01)) );
	err_sig.Push( fabs(g.Ele(1)->Val(4, "Sxy") - ( 8.20878435e-01)) );
	err_sig.Push( fabs(g.Ele(1)->Val(5, "Sxy") - (-1.47127394e+00)) );
	*/

	// Error summary
	//double tol_eps     = 1.0e-16;
	double tol_sig     = 1.0e-14;
	//double tol_dis     = 1.0e-16;
	//double min_err_eps = err_eps[err_eps.Min()];
	double min_err_sig = err_sig[err_sig.Min()];
	//double min_err_dis = err_dis[err_dis.Min()];
	//double max_err_eps = err_eps[err_eps.Max()];
	double max_err_sig = err_sig[err_sig.Max()];
	//double max_err_dis = err_dis[err_dis.Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	//cout << _4<< "Eps" << _8s<<min_err_eps << _8s<<err_eps.Mean() << (max_err_eps>tol_eps?"[1;31m":"[1;32m") << _8s<<max_err_eps << "[0m" << _8s<<err_eps.Norm() << endl;
	cout << _4<< "Sig" << _8s<<min_err_sig << _8s<<err_sig.Mean() << (max_err_sig>tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig << "[0m" << _8s<<err_sig.Norm() << endl;
	//cout << _4<< "Dis" << _8s<<min_err_dis << _8s<<err_dis.Mean() << (max_err_dis>tol_dis?"[1;31m":"[1;32m") << _8s<<max_err_dis << "[0m" << _8s<<err_dis.Norm() << endl;
	cout << endl;

	// Return error flag
	//if (max_err_eps>tol_eps || max_err_sig>tol_sig || max_err_dis>tol_dis) return 1;
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
