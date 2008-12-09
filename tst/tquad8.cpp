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
	FEM::Data dat(2); // 2D

	// 1) Nodes
	dat.SetNNodes (8);
	dat.SetNode   (0,   0.0 ,   0.0);
	dat.SetNode   (1,     L ,   0.0);
	dat.SetNode   (2,     L ,     H);
	dat.SetNode   (3,   0.0 ,     H);
	dat.SetNode   (4, L/2.0 ,   0.0);
	dat.SetNode   (5,     L , H/2.0);
	dat.SetNode   (6, L/2.0 ,     H);
	dat.SetNode   (7,   0.0 , H/2.0);

	// 2) Elements
	dat.SetNElems (1);
	dat.SetElem   (0, "Quad8PStrain", /*IsActive*/true, /*Tag*/-1);

	// 3) Set connectivity (list of nodes must be LOCAL)
	dat.Ele(0)->Connect(0, dat.Nod(0))
	        ->Connect(1, dat.Nod(1))
	        ->Connect(2, dat.Nod(2))
	        ->Connect(3, dat.Nod(3))
	        ->Connect(4, dat.Nod(4))
	        ->Connect(5, dat.Nod(5))
	        ->Connect(6, dat.Nod(6))
	        ->Connect(7, dat.Nod(7));

	dat.Ele(0)->SetIntPoints(4);

	// 4) Boundary conditions (must be after connectivity)
	dat.Nod(0)->Bry     ("uy",0.0);
	dat.Nod(4)->Bry     ("uy",0.0)->Bry("ux",0.0);
	dat.Nod(1)->Bry     ("uy",0.0);
	dat.Ele(0)->EdgeBry ("fy",q,3); // 3 => top edge

	// 5) Parameters and initial values
	String prms; prms.Printf("E=%f  nu=%f",E,nu);
	dat.Ele(0)->SetModel("LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0");

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;
	dat.Ele(0)->Order1Matrix(0,Ke0);
	//cout << "Ke0=\n" << Ke0 << endl;

	// 6) Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&dat)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);
	delete sol;

	// Output
	cout << "Node 3: ux = " << dat.Nod(3)->Val("ux") << " : uy = " << dat.Nod(3)->Val("uy") << " : fy = "  << dat.Nod(3)->Val("fy")  << endl;
	cout << "Node 6: ux = " << dat.Nod(6)->Val("ux") << " : uy = " << dat.Nod(6)->Val("uy") << " : fy = "  << dat.Nod(6)->Val("fy")  << endl;
	cout << "Node 2: ux = " << dat.Nod(2)->Val("ux") << " : uy = " << dat.Nod(2)->Val("uy") << " : fy = "  << dat.Nod(2)->Val("fy")  << endl << endl;;
	cout << "Node 0: ux = " << dat.Nod(0)->Val("ux") << " : uy = " << dat.Nod(0)->Val("uy") << " : fy = "  << dat.Nod(0)->Val("fy")  << endl;
	cout << "Node 4: ux = " << dat.Nod(4)->Val("ux") << " : uy = " << dat.Nod(4)->Val("uy") << " : fy = "  << dat.Nod(4)->Val("fy")  << endl;
	cout << "Node 1: ux = " << dat.Nod(1)->Val("ux") << " : uy = " << dat.Nod(1)->Val("uy") << " : fy = "  << dat.Nod(1)->Val("fy")  << endl << endl;;
	cout << "Elem 0: Sx = " << dat.Ele(0)->Val("Sx") << " : Sy = " << dat.Ele(0)->Val("Sy") << " : Sxy = " << dat.Ele(0)->Val("Sxy") << endl;
	cout << "Elem 0: Ex = " << dat.Ele(0)->Val("Ex") << " : Ey = " << dat.Ele(0)->Val("Ey") << " : Exy = " << dat.Ele(0)->Val("Exy") << endl;

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
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j)
		{
			err_eps.Push ( fabs(dat.Ele(i)->Val(j,"Ex" ) - Ex ) / (1.0+fabs(Ex )) );
			err_eps.Push ( fabs(dat.Ele(i)->Val(j,"Ey" ) - Ey ) / (1.0+fabs(Ey )) );
			err_eps.Push ( fabs(dat.Ele(i)->Val(j,"Ez" ) - Ez ) / (1.0+fabs(Ez )) );
			err_eps.Push ( fabs(dat.Ele(i)->Val(j,"Exy") - Exy) / (1.0+fabs(Exy)) );
			err_sig.Push ( fabs(dat.Ele(i)->Val(j,"Sx" ) - Sx ) / (1.0+fabs(Sx )) );
			err_sig.Push ( fabs(dat.Ele(i)->Val(j,"Sy" ) - Sy ) / (1.0+fabs(Sy )) );
			err_sig.Push ( fabs(dat.Ele(i)->Val(j,"Sz" ) - Sz ) / (1.0+fabs(Sz )) );
			err_sig.Push ( fabs(dat.Ele(i)->Val(j,"Sxy") - Sxy) / (1.0+fabs(Sxy)) );
		}
	}

	// Displacements
	for (size_t i=0; i<dat.NNodes(); ++i)
	{
		double ux_correct = Ex*(dat.Nod(i)->X()-L/2.0);
		double uy_correct = Ey* dat.Nod(i)->Y();
		err_dis.Push ( fabs(dat.Nod(i)->Val("ux") - ux_correct) / (1.0+fabs(ux_correct)) );
		err_dis.Push ( fabs(dat.Nod(i)->Val("uy") - uy_correct) / (1.0+fabs(uy_correct)) );
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
