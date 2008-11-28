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
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/quad8pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_8s;

int main(int argc, char **argv) try
{
	double H  = 2.0;   // height
	double L  = 2.0;   // length
	double E  = 207.0; // Young
	double nu = 0.3;   // Poisson
	double q  = -1.0;  // Load

	/*        | | | | | | | | | | |  q=1
	          V V V V V V V V V V V 
	  ---   3 @---------@---------@ 2
	   |      |         6         |
	   |      |                   |
	   |      |                   |
	          |                   |
	   H    7 @        e[0]       @ 5
	          |                   |
	   |      |                   |
	   |      |                   |
	   |      |         4         |
	  ---   0 @---------@---------@ 1
	         /_\       /_\       /_\
	         o o       ///       o o
	
	          |-------- L --------|
	 */

	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);

	// 0) Problem dimension
	FEM::Geom g(2); // 2D

	// 1) Nodes
	g.SetNNodes (8);
	g.SetNode   (0,   0.0 ,   0.0);
	g.SetNode   (1,     L ,   0.0);
	g.SetNode   (2,     L ,     H);
	g.SetNode   (3,   0.0 ,     H);
	g.SetNode   (4, L/2.0 ,   0.0);
	g.SetNode   (5,     L , H/2.0);
	g.SetNode   (6, L/2.0 ,     H);
	g.SetNode   (7,   0.0 , H/2.0);

	// 2) Elements
	g.SetNElems (1);
	g.SetElem   (0, "Quad8PStrain", /*IsActive*/true, /*Tag*/-1);

	// 3) Set connectivity (list of nodes must be LOCAL)
	g.Ele(0)->Connect(0, g.Nod(0))
	        ->Connect(1, g.Nod(1))
	        ->Connect(2, g.Nod(2))
	        ->Connect(3, g.Nod(3))
	        ->Connect(4, g.Nod(4))
	        ->Connect(5, g.Nod(5))
	        ->Connect(6, g.Nod(6))
	        ->Connect(7, g.Nod(7));

	g.Ele(0)->SetIntPoints(4);

	// 4) Boundary conditions (must be after connectivity)
	g.Nod(0)->Bry     ("uy",0.0);
	g.Nod(4)->Bry     ("uy",0.0)->Bry("ux",0.0);
	g.Nod(1)->Bry     ("uy",0.0);
	g.Ele(0)->EdgeBry ("fy",q,3); // 3 => top edge

	// 5) Parameters and initial values
	String prms; prms.Printf("E=%f  nu=%f",E,nu);
	g.Ele(0)->SetModel("LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;
	g.Ele(0)->Order1Matrix(0,Ke0);
	//cout << "Ke0=\n" << Ke0 << endl;

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);
	delete sol;

	// Output
	cout << "Node 3: ux = " << g.Nod(3)->Val("ux") << " : uy = " << g.Nod(3)->Val("uy") << " : fy = "  << g.Nod(3)->Val("fy")  << endl;
	cout << "Node 6: ux = " << g.Nod(6)->Val("ux") << " : uy = " << g.Nod(6)->Val("uy") << " : fy = "  << g.Nod(6)->Val("fy")  << endl;
	cout << "Node 2: ux = " << g.Nod(2)->Val("ux") << " : uy = " << g.Nod(2)->Val("uy") << " : fy = "  << g.Nod(2)->Val("fy")  << endl << endl;;
	cout << "Node 0: ux = " << g.Nod(0)->Val("ux") << " : uy = " << g.Nod(0)->Val("uy") << " : fy = "  << g.Nod(0)->Val("fy")  << endl;
	cout << "Node 4: ux = " << g.Nod(4)->Val("ux") << " : uy = " << g.Nod(4)->Val("uy") << " : fy = "  << g.Nod(4)->Val("fy")  << endl;
	cout << "Node 1: ux = " << g.Nod(1)->Val("ux") << " : uy = " << g.Nod(1)->Val("uy") << " : fy = "  << g.Nod(1)->Val("fy")  << endl << endl;;
	cout << "Elem 0: Sx = " << g.Ele(0)->Val("Sx") << " : Sy = " << g.Ele(0)->Val("Sy") << " : Sxy = " << g.Ele(0)->Val("Sxy") << endl;
	cout << "Elem 0: Ex = " << g.Ele(0)->Val("Ex") << " : Ey = " << g.Ele(0)->Val("Ey") << " : Exy = " << g.Ele(0)->Val("Exy") << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
    Array<double> err_eps;
    Array<double> err_sig;
    Array<double> err_dis;

	double Sx  = 0.0;
	double Sy  = q;
	double Ex  = -nu*(1.0+nu)*Sy/E;
	double Ey  =  (1.0-nu*nu)*Sy/E;
	double Ez  = 0.0;
	double Exy = 0.0;
	double Sz  = (E/(1.0+nu))*(nu/(1.0-2.0*nu))*(Ex+Ey);
	double Sxy = 0.0;

	// Stress and strain
	for (size_t i=0; i<g.NElems(); ++i)
	{
		for (size_t j=0; j<g.Ele(i)->NNodes(); ++j)
		{
			err_eps.Push ( fabs(g.Ele(i)->Val(j,"Ex" ) - Ex ) / (1.0+fabs(Ex )) );
			err_eps.Push ( fabs(g.Ele(i)->Val(j,"Ey" ) - Ey ) / (1.0+fabs(Ey )) );
			err_eps.Push ( fabs(g.Ele(i)->Val(j,"Ez" ) - Ez ) / (1.0+fabs(Ez )) );
			err_eps.Push ( fabs(g.Ele(i)->Val(j,"Exy") - Exy) / (1.0+fabs(Exy)) );
			err_sig.Push ( fabs(g.Ele(i)->Val(j,"Sx" ) - Sx ) / (1.0+fabs(Sx )) );
			err_sig.Push ( fabs(g.Ele(i)->Val(j,"Sy" ) - Sy ) / (1.0+fabs(Sy )) );
			err_sig.Push ( fabs(g.Ele(i)->Val(j,"Sz" ) - Sz ) / (1.0+fabs(Sz )) );
			err_sig.Push ( fabs(g.Ele(i)->Val(j,"Sxy") - Sxy) / (1.0+fabs(Sxy)) );
		}
	}

	// Displacements
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		double ux_correct = Ex*(g.Nod(i)->X()-L/2.0);
		double uy_correct = Ey* g.Nod(i)->Y();
		err_dis.Push ( fabs(g.Nod(i)->Val("ux") - ux_correct) / (1.0+fabs(ux_correct)) );
		err_dis.Push ( fabs(g.Nod(i)->Val("uy") - uy_correct) / (1.0+fabs(uy_correct)) );
	}

	// Error summary
	double tol_eps     = 1.0e-16;
	double tol_sig     = 1.0e-14;
	double tol_dis     = 1.0e-16;
	double min_err_eps = err_eps[err_eps.Min()];
	double min_err_sig = err_sig[err_sig.Min()];
	double min_err_dis = err_dis[err_dis.Min()];
	double max_err_eps = err_eps[err_eps.Max()];
	double max_err_sig = err_sig[err_sig.Max()];
	double max_err_dis = err_dis[err_dis.Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Eps" << _8s<<min_err_eps << _8s<<err_eps.Mean() << (max_err_eps>tol_eps?"[1;31m":"[1;32m") << _8s<<max_err_eps << "[0m" << _8s<<err_eps.Norm() << endl;
	cout << _4<< "Sig" << _8s<<min_err_sig << _8s<<err_sig.Mean() << (max_err_sig>tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig << "[0m" << _8s<<err_sig.Norm() << endl;
	cout << _4<< "Dis" << _8s<<min_err_dis << _8s<<err_dis.Mean() << (max_err_dis>tol_dis?"[1;31m":"[1;32m") << _8s<<max_err_dis << "[0m" << _8s<<err_dis.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_eps>tol_eps || max_err_sig>tol_sig || max_err_dis>tol_dis) return 1;
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
