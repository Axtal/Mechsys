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
#include "fem/elems/rod.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/output.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_6;
using Util::_8s;

int main(int argc, char **argv) try
{
	/*  Small truss
	                                2 |
	                                  |
	                                  V 
	                                2 o----> 1
	                                ,'|
	                              ,'  |
	                   E=1.0    ,'    |
	                   A=1.0  ,'      | E=1.0
	                    [2] ,'        | A=1.0
	                      ,'          | [1]
	                    ,'            |
	                  ,'              |
	   y            ,'                |
	   |         0,'        [0]       |1
	   |         o--------------------o
	  (z)___x   /_\        E=1.0     /_\
	           ////        A=1.0     ///  
	*/

	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);

	// Geometry
	FEM::Geom g(2); // 2D

	// Nodes
	g.SetNNodes (3);
	g.SetNode   (0, 0.0, 0.0);
	g.SetNode   (1, 1.0, 0.0);
	g.SetNode   (2, 1.0, 1.0);

	// Elements
	g.SetNElems (3);
	g.SetElem   (0, "Rod", /*IsActive*/true);
	g.SetElem   (1, "Rod", /*IsActive*/true);
	g.SetElem   (2, "Rod", /*IsActive*/true);

	// Number of integration points
	g.Ele(0)->SetIntPoints(2);
	g.Ele(1)->SetIntPoints(2);
	g.Ele(2)->SetIntPoints(2);

	// Set connectivity
	g.Ele(0)->Connect(0, g.Nod(0))->Connect(1, g.Nod(1));
	g.Ele(1)->Connect(0, g.Nod(1))->Connect(1, g.Nod(2));
	g.Ele(2)->Connect(0, g.Nod(0))->Connect(1, g.Nod(2));

	// Parameters and initial value
	g.Ele(0)->SetModel("LinElastic", "E=1.0 A=1.0", "ZERO");
	g.Ele(1)->SetModel("LinElastic", "E=1.0 A=1.0", "ZERO");
	g.Ele(2)->SetModel("LinElastic", "E=1.0 A=1.0", "ZERO");

	// Boundary conditions (must be after set connectivity)
	g.Nod(0)->Bry("ux", 0.0)->Bry("uy",  0.0);
	g.Nod(1)->Bry("ux", 0.0)->Bry("uy",  0.0);
	g.Nod(2)->Bry("fx", 1.0)->Bry("fy", -2.0);

	// Check stiffness matrices
	double err_ke = 0.0;
	LinAlg::Matrix<double> Ke0,  Ke1,  Ke2;
	LinAlg::Matrix<double> Ke0c, Ke1c, Ke2c; // correct matrices
	double c = 0.5/sqrt(2.0);
	Ke0c.Resize(4,4);
	Ke1c.Resize(4,4);
	Ke2c.Resize(4,4);
	g.Ele(0)->Order1Matrix(0,Ke0);
	g.Ele(1)->Order1Matrix(0,Ke1);
	g.Ele(2)->Order1Matrix(0,Ke2);
	Ke0c =  1.0,  0.0, -1.0,  0.0,
	        0.0,  0.0,  0.0,  0.0,
	       -1.0,  0.0,  1.0,  0.0,
	        0.0,  0.0,  0.0,  0.0;
	Ke1c =  0.0,  0.0,  0.0,  0.0,
	        0.0,  1.0,  0.0, -1.0,
	        0.0,  0.0,  0.0,  0.0,
	        0.0, -1.0,  0.0,  1.0;
	Ke2c =    c,    c,   -c,   -c,
		      c,    c,   -c,   -c,
		     -c,   -c,    c,    c,
		     -c,   -c,    c,    c;
	for (int i=0; i<4; ++i)
	for (int j=0; j<4; ++j)
	{
		err_ke += fabs(Ke0(i,j)-Ke0c(i,j));
		err_ke += fabs(Ke1(i,j)-Ke1c(i,j));
		err_ke += fabs(Ke2(i,j)-Ke2c(i,j));
	}
	if (err_ke>8.9e-16) throw new Fatal("tex831: err_ke=%e is bigger than %e.",err_ke,8.9e-16);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();
	double norm_resid = LinAlg::Norm(sol->Resid());
	cout << "\n[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
	cout << "[1;32mNumber of DOFs          = " << sol->nDOF() << "[0m\n";
	if (norm_resid>1.1e-15) throw new Fatal("tex831: norm_resid=%e is bigger than %e.",norm_resid,1.1e-15);

	// Output: VTU
	Output o; o.VTU (&g, "tex461.vtu");
	cout << "[1;34mFile <tex461.vtu> saved.[0m\n\n";

	// Output: Nodes
	cout << _6<<"Node #" << _8s<<"ux" << _8s<<"uy" << _8s<<"fx"<< _8s<<"fy" << endl;
	for (size_t i=0; i<g.NNodes(); ++i)
		cout << _6<<i << _8s<<g.Nod(i)->Val("ux") <<  _8s<<g.Nod(i)->Val("uy") << _8s<<g.Nod(i)->Val("fx") << _8s<<g.Nod(i)->Val("fy") << endl;
	cout << endl;

	// Output: Elements
	cout << _6<<"Elem #" << _8s<<"Sa" << _8s<<"Ea" << _8s<<"Sa(left)" << _8s<<"Sa(right)" << _8s<<"Ea(left)" << _8s<<"Ea(right)" << endl;
	for (size_t i=0; i<g.NElems(); ++i)
	{
		cout << _6<<i << _8s<<g.Ele(i)->Val("Sa") <<  _8s<<g.Ele(i)->Val("Ea");
		for (size_t j=0; j<g.Ele(i)->NNodes(); ++j) cout << _8s<<g.Ele(i)->Val(j, "Sa");
		for (size_t j=0; j<g.Ele(i)->NNodes(); ++j) cout << _8s<<g.Ele(i)->Val(j, "Ea");
		cout << endl;
	}
	cout << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////
	
	// Displacements
	Array<double> err_u(6);
	err_u[0] = fabs(g.Nod(0)->Val("ux")-0.0);
	err_u[1] = fabs(g.Nod(0)->Val("uy")-0.0);
	err_u[2] = fabs(g.Nod(1)->Val("ux")-0.0);
	err_u[3] = fabs(g.Nod(1)->Val("uy")-0.0);
	err_u[4] = fabs(g.Nod(2)->Val("ux")-(3.0+2.0*sqrt(2.0)));
	err_u[5] = fabs(g.Nod(2)->Val("uy")-(-3.0));

	// Forces
	Array<double> err_f(6);
	err_f[0] = fabs(g.Nod(0)->Val("fx")-(-1.0));
	err_f[1] = fabs(g.Nod(0)->Val("fy")-(-1.0));
	err_f[2] = fabs(g.Nod(1)->Val("fx")-  0.0 );
	err_f[3] = fabs(g.Nod(1)->Val("fy")-  3.0 );
	err_f[4] = fabs(g.Nod(2)->Val("fx")-  1.0 );
	err_f[5] = fabs(g.Nod(2)->Val("fy")-(-2.0));

	// Stresses
	Array<double> err_s(3);
	err_s[0] = fabs(g.Ele(0)->Val("Sa")- 0.0);
	err_s[1] = fabs(g.Ele(1)->Val("Sa")-(3.0));
	err_s[2] = fabs(g.Ele(2)->Val("Sa")-(-sqrt(2.0)));

	// Strains
	Array<double> err_e(3);
	err_e[0] = fabs(g.Ele(0)->Val("Ea")- 0.0);
	err_e[1] = fabs(g.Ele(1)->Val("Sa")-(3.0));
	err_e[2] = fabs(g.Ele(2)->Val("Sa")-(-sqrt(2.0)));

	// Error summary
	double tol_u     = DBL_EPSILON;
	double tol_f     = 4.45e-16;
	double tol_s     = 4.45e-16;
	double tol_e     = 4.45e-16;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	double min_err_f = err_f[err_f.Min()];
	double max_err_f = err_f[err_f.Max()];
	double min_err_s = err_s[err_s.Min()];
	double max_err_s = err_s[err_s.Max()];
	double min_err_e = err_e[err_e.Min()];
	double max_err_e = err_e[err_e.Max()];
	cout << _4<< ""    << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u"   << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << _4<< "f"   << _8s<<min_err_f << _8s<<err_f.Mean() << (max_err_f>tol_f?"[1;31m":"[1;32m") << _8s<<max_err_f << "[0m" << _8s<<err_f.Norm() << endl;
	cout << _4<< "Sig" << _8s<<min_err_s << _8s<<err_s.Mean() << (max_err_s>tol_s?"[1;31m":"[1;32m") << _8s<<max_err_s << "[0m" << _8s<<err_s.Norm() << endl;
	cout << _4<< "Eps" << _8s<<min_err_e << _8s<<err_e.Mean() << (max_err_e>tol_e?"[1;31m":"[1;32m") << _8s<<max_err_e << "[0m" << _8s<<err_e.Norm() << endl;
	cout << endl;

	return 0;

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
