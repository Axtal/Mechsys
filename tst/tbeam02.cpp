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
#include "fem/elems/beam.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
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
	/*  _  0             1     2     3     4      5
	 *  _|\@-------------@-----@-----@-----@------@
	 *  _|/       0      |  1     2     3,-|   4
	 *                   |            ,-'  |
	 *                  5|         ,-'     |
	 *                   |      ,-'6       |7
	 *                   |   ,-'           |
	 *                   |,-'              |
	 *                 6 @                 |
	 *                  ###                @ 7
	 *                                    ###
	 */
	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);

	// Geometry
	FEM::Geom g(2); // 2D

	// Nodes
	g.SetNNodes (8);
	g.SetNode   (0,  0.0, 5.0);
	g.SetNode   (1,  6.0, 5.0);
	g.SetNode   (2,  8.0, 5.0);
	g.SetNode   (3, 10.0, 5.0);
	g.SetNode   (4, 12.0, 5.0);
	g.SetNode   (5, 14.0, 5.0);
	g.SetNode   (6,  6.0, 1.0);
	g.SetNode   (7, 12.0, 0.0);

	// Elements
	g.SetNElems (8);
	g.SetElem   (0, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, g.Nod(0))->Connect(1, g.Nod(1));
	g.SetElem   (1, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, g.Nod(1))->Connect(1, g.Nod(2));
	g.SetElem   (2, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, g.Nod(2))->Connect(1, g.Nod(3));
	g.SetElem   (3, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, g.Nod(3))->Connect(1, g.Nod(4));
	g.SetElem   (4, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, g.Nod(4))->Connect(1, g.Nod(5));
	g.SetElem   (5, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, g.Nod(6))->Connect(1, g.Nod(1));
	g.SetElem   (6, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, g.Nod(6))->Connect(1, g.Nod(4));
	g.SetElem   (7, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, g.Nod(7))->Connect(1, g.Nod(4));

	// Parameters and initial value
	g.Ele(0)->SetModel("LinElastic", "E=1.0 A=5e+9 Izz=6e+4", "ZERO");
	g.Ele(1)->SetModel("LinElastic", "E=1.0 A=5e+9 Izz=6e+4", "ZERO");
	g.Ele(2)->SetModel("LinElastic", "E=1.0 A=5e+9 Izz=6e+4", "ZERO");
	g.Ele(3)->SetModel("LinElastic", "E=1.0 A=5e+9 Izz=6e+4", "ZERO");
	g.Ele(4)->SetModel("LinElastic", "E=1.0 A=5e+9 Izz=6e+4", "ZERO");
	g.Ele(5)->SetModel("LinElastic", "E=1.0 A=1e+9 Izz=2e+4", "ZERO");
	g.Ele(6)->SetModel("LinElastic", "E=1.0 A=1e+9 Izz=2e+4", "ZERO");
	g.Ele(7)->SetModel("LinElastic", "E=1.0 A=1e+9 Izz=2e+4", "ZERO");

	// Boundary conditions (must be after set connectivity)
	g.Ele(0)->EdgeBry("q", -20.0, -20.0, 0);
	g.Ele(1)->EdgeBry("q", -20.0, -20.0, 0);
	g.Ele(2)->EdgeBry("q", -20.0, -20.0, 0);
	g.Ele(3)->EdgeBry("q", -20.0, -20.0, 0);
	g.Ele(4)->EdgeBry("q", -20.0, -20.0, 0);
	g.Nod(2)->Bry("fy", -60.0);
	g.Nod(3)->Bry("fy", -60.0);
	g.Nod(0)->Bry("ux", 0.0)->Bry("uy", 0.0);
	g.Nod(6)->Bry("ux", 0.0)->Bry("uy", 0.0)->Bry("wz", 0.0);
	g.Nod(7)->Bry("ux", 0.0)->Bry("uy", 0.0)->Bry("wz", 0.0);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve(); ///*tIni*/0.0, /*tFin*/1.0, /*hIni*/0.001, /*DtOut*/0.01);
	double norm_resid = LinAlg::Norm(sol->Resid());
	cout << "\n[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
	cout << "[1;32mNumber of DOFs          = " << sol->nDOF() << "[0m\n";
	if (norm_resid>1.0e-12) throw new Fatal("tex831: norm_resid=%e is bigger than %e.",norm_resid,1.0e-12);
	cout << endl;
	delete sol;

	// Output: Nodes
	cout << _6<<"Node #" << _8s<<"ux" << _8s<<"uy" << _8s<<"wz" << _8s<<"fx"<< _8s<<"fy" << _8s<<"mz" << endl;
	for (size_t i=0; i<g.NNodes(); ++i)
		cout << _6<<i << _8s<<g.Nod(i)->Val("ux") <<  _8s<<g.Nod(i)->Val("uy") << _8s<<g.Nod(i)->Val("wz") << _8s<<g.Nod(i)->Val("fx") << _8s<<g.Nod(i)->Val("fy") << _8s<<g.Nod(i)->Val("mz") << endl;
	cout << endl;

	// Output: Elements
	cout << _6<<"Elem #" << _8s<<"N0" << _8s<<"M0" << _8s<<"V0" << _8s<<"N1" << _8s<<"M1" << _8s<<"V1" << endl;
	for (size_t i=0; i<g.NElems(); ++i)
	{
		g.Ele(i)->CalcDepVars();
		cout << _6<<i << _8s<<g.Ele(i)->Val(0, "N") << _8s<<g.Ele(i)->Val(0, "M") << _8s<<g.Ele(i)->Val(0, "V");
		cout <<          _8s<<g.Ele(i)->Val(1, "N") << _8s<<g.Ele(i)->Val(1, "M") << _8s<<g.Ele(i)->Val(1, "V") << endl;
	}
	cout << endl;

	// Output: VTU
	Output o; o.VTU (&g, "tbeam02.vtu");
	cout << "[1;34mFile <tbeam02.vtu> saved.[0m\n\n";

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	/*
	// Displacements
	Array<double> err_u(9);
	err_u[ 0] = fabs(g.Nod(0)->Val("ux") - (0.0));

	// Forces
	Array<double> err_f(9);
	err_f[ 0] = fabs(g.Nod(0)->Val("fx") - (0.0));

	// Stresses
	Array<double> err_s(12);
	err_s[ 0] = fabs(g.Ele(0)->Val(0,"N") - (  0.0));   err_s[ 6] = fabs(g.Ele(0)->Val(1,"N") - (  0.0));

	// Error summary
	double tol_u     = 1.0e-12;
	double tol_f     = 1.0e-13;
	double tol_s     = 1.0e-13;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	double min_err_f = err_f[err_f.Min()];
	double max_err_f = err_f[err_f.Max()];
	double min_err_s = err_s[err_s.Min()];
	double max_err_s = err_s[err_s.Max()];
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << _4<< "f" << _8s<<min_err_f << _8s<<err_f.Mean() << (max_err_f>tol_f?"[1;31m":"[1;32m") << _8s<<max_err_f << "[0m" << _8s<<err_f.Norm() << endl;
	cout << _4<< "s" << _8s<<min_err_s << _8s<<err_s.Mean() << (max_err_s>tol_s?"[1;31m":"[1;32m") << _8s<<max_err_s << "[0m" << _8s<<err_s.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_u>tol_u || max_err_f>tol_f || max_err_s>tol_s) return 1;
	else return 0;
	*/

	return 1;
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
