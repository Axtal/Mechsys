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

/* J. N. Reddy's Finite Element Method:   Example 8.3.1   */

// STL
#include <iostream>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "fem/elems/tri3diffusion.h"
#include "models/diffusions/lindiffusion.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"

using std::cout;
using std::endl;
using LinAlg::Vector;

int main(int argc, char **argv) try
{
	// Input
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou can call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	// Geometry
	FEM::Geom g(2); // 2D

	// Nodes
	g.SetNNodes (6);
	g.SetNode   (0, 0.0, 0.0);
	g.SetNode   (1, 0.5, 0.0);
	g.SetNode   (2, 0.5, 0.5);
	g.SetNode   (3, 1.0, 0.0);
	g.SetNode   (4, 1.0, 0.5);
	g.SetNode   (5, 1.0, 1.0);

	// Elements
	g.SetNElems (4);
	g.SetElem   (0, "Tri3Diffusion");
	g.SetElem   (1, "Tri3Diffusion");
	g.SetElem   (2, "Tri3Diffusion");
	g.SetElem   (3, "Tri3Diffusion");

	// Set connectivity
	g.Ele(0)->Connect(0, g.Nod(0))->Connect(1, g.Nod(1))->Connect(2, g.Nod(2));
	g.Ele(1)->Connect(0, g.Nod(4))->Connect(1, g.Nod(2))->Connect(2, g.Nod(1));
	g.Ele(2)->Connect(0, g.Nod(1))->Connect(1, g.Nod(3))->Connect(2, g.Nod(4));
	g.Ele(3)->Connect(0, g.Nod(2))->Connect(1, g.Nod(4))->Connect(2, g.Nod(5));

	// Parameters and initial values
	g.Ele(0)->SetModel("LinDiffusion", "k=1.0", "");
	g.Ele(1)->SetModel("LinDiffusion", "k=1.0", "");
	g.Ele(2)->SetModel("LinDiffusion", "k=1.0", "");
	g.Ele(3)->SetModel("LinDiffusion", "k=1.0", "");
	
	// Properties (heat source)
	Array<double> source(1);
	//source.SetValues(6.0/0.25);
	source.SetValues(1.0);
	g.Ele(0)->SetProps(source);
	g.Ele(1)->SetProps(source);
	g.Ele(2)->SetProps(source);
	g.Ele(3)->SetProps(source);

	// Boundary conditions (must be after connectivity)
	g.Ele(0)->EdgeBry("q", 0.0, 0)->EdgeBry("q", 0.0, 2);
	g.Ele(2)->EdgeBry("q", 0.0, 0)->EdgeBry("u", 0.0, 1);
	g.Ele(3)->EdgeBry("u", 0.0, 1)->EdgeBry("q", 0.0, 2);
	g.Nod(3)->Bry("u",0.0);
	g.Nod(5)->Bry("u",0.0);

	// Stiffness
	LinAlg::Matrix<double> Ke0, Ke1, Ke2, Ke3;
	LinAlg::Matrix<double> Ke_correct;  Ke_correct.Resize(3,3);
	g.Ele(0)->Order1Matrix(0,Ke0);
	g.Ele(1)->Order1Matrix(0,Ke1);
	g.Ele(2)->Order1Matrix(0,Ke2);
	g.Ele(3)->Order1Matrix(0,Ke3);
	Ke_correct =  0.5, -0.5,  0.0,
	             -0.5,  1.0, -0.5,
	              0.0, -0.5,  0.5;

	// Output
	Ke0.SetNS(Util::_6_3);
	Ke1.SetNS(Util::_6_3);
	Ke2.SetNS(Util::_6_3);
	Ke3.SetNS(Util::_6_3);
	cout << Ke0 << endl;
	cout << Ke1 << endl;
	cout << Ke2 << endl;
	cout << Ke3 << endl;

	// Check
	double errors = 0.0;
	for (int i=0; i<3; ++i)
	for (int j=0; j<3; ++j)
	{
		errors += fabs(Ke0(i,j)-Ke_correct(i,j));
		errors += fabs(Ke1(i,j)-Ke_correct(i,j));
		errors += fabs(Ke2(i,j)-Ke_correct(i,j));
		errors += fabs(Ke3(i,j)-Ke_correct(i,j));
	}

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve();
	cout << "DFext = \n"               << sol->DFext() << endl;
	cout << "DFint = \n"               << sol->DFint() << endl;
	cout << "Resid = DFext-DFint = \n" << sol->Resid() << endl;

	if (fabs(errors)>1.0e-14) cout << "[1;31mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                      cout << "[1;32mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	// Output
	cout << "  [ID](X,Y,Z) : {EssentialBryName,NaturalBryName,EssentialBry,NaturalBry,IsEssenPresc,EqID,EssentialVal,NaturalVal}" << endl;;
	cout << g << endl;
	Output o; o.VTU (&g, "tex831.vtu");

	// Output
	LinAlg::Matrix<double> vals0,vals1,vals2,vals3;
	Array<String>          labs0,labs1,labs2,labs3;
	g.Ele(0)->OutNodes (vals0, labs0);
	g.Ele(1)->OutNodes (vals1, labs1);
	g.Ele(2)->OutNodes (vals2, labs2);
	g.Ele(3)->OutNodes (vals3, labs3);
	std::cout << "Element # 0\n" << labs0 << "\n" << vals0 << std::endl;
	std::cout << "Element # 1\n" << labs1 << "\n" << vals1 << std::endl;
	std::cout << "Element # 2\n" << labs2 << "\n" << vals2 << std::endl;
	std::cout << "Element # 3\n" << labs3 << "\n" << vals3 << std::endl;

	for (int i=0; i<6; ++i)
		cout << "Node # " << i << ":  u=" << g.Nod(i)->Val("u") << ",  q=" << g.Nod(i)->Val("q") << endl;

	return 1;

	// Element 0
	errors += fabs(g.Ele(0)->Val(0, "Sx") - ( 1.56432140e-01));
	errors += fabs(g.Ele(0)->Val(1, "Sx") - (-3.00686928e-01));
	errors += fabs(g.Ele(0)->Val(2, "Sx") - ( 1.44254788e-01));
	errors += fabs(g.Ele(0)->Val(3, "Sx") - (-3.19109076e-01));
	errors += fabs(g.Ele(0)->Val(4, "Sx") - (-3.31286428e-01));
	errors += fabs(g.Ele(0)->Val(5, "Sx") - ( 1.25832639e-01));

	errors += fabs(g.Ele(0)->Val(0, "Sy") - (-2.05141549e-01));
	errors += fabs(g.Ele(0)->Val(1, "Sy") - ( 1.15872190e+00));
	errors += fabs(g.Ele(0)->Val(2, "Sy") - (-9.53580350e-01));
	errors += fabs(g.Ele(0)->Val(3, "Sy") - (-2.22127394e+00));
	errors += fabs(g.Ele(0)->Val(4, "Sy") - (-2.96971274e+00));
	errors += fabs(g.Ele(0)->Val(5, "Sy") - (-4.33357619e+00));

	errors += fabs(g.Ele(0)->Val(0, "Sxy") - (-1.56432140e-01));
	errors += fabs(g.Ele(0)->Val(1, "Sxy") - (-6.74437968e-02));
	errors += fabs(g.Ele(0)->Val(2, "Sxy") - ( 2.23875937e-01));
	errors += fabs(g.Ele(0)->Val(3, "Sxy") - (-4.90216486e-02));
	errors += fabs(g.Ele(0)->Val(4, "Sxy") - ( 3.31286428e-01));
	errors += fabs(g.Ele(0)->Val(5, "Sxy") - ( 2.42298085e-01));

	// Element 1
	errors += fabs(g.Ele(1)->Val(0, "Sx")  - ( 9.95732723e-01));
	errors += fabs(g.Ele(1)->Val(1, "Sx")  - ( 2.23875937e-01));
	errors += fabs(g.Ele(1)->Val(2, "Sx")  - (-1.21960866e+00));
	errors += fabs(g.Ele(1)->Val(3, "Sx")  - ( 1.39446295e+00));
	errors += fabs(g.Ele(1)->Val(4, "Sx")  - (-8.20878435e-01));
	errors += fabs(g.Ele(1)->Val(5, "Sx")  - (-4.90216486e-02));

	errors += fabs(g.Ele(1)->Val(0, "Sy")  - (-1.25426728e+00));
	errors += fabs(g.Ele(1)->Val(1, "Sy")  - ( 1.68328476e+00));
	errors += fabs(g.Ele(1)->Val(2, "Sy")  - (-4.29017485e-01));
	errors += fabs(g.Ele(1)->Val(3, "Sy")  - (-2.39612823e+00));
	errors += fabs(g.Ele(1)->Val(4, "Sy")  - (-1.57087843e+00));
	errors += fabs(g.Ele(1)->Val(5, "Sy")  - (-4.50843047e+00));

	errors += fabs(g.Ele(1)->Val(0, "Sxy") - (-9.95732723e-01));
	errors += fabs(g.Ele(1)->Val(1, "Sxy") - ( 1.29641965e+00));
	errors += fabs(g.Ele(1)->Val(2, "Sxy") - (-3.00686928e-01));
	errors += fabs(g.Ele(1)->Val(3, "Sxy") - ( 1.25832639e-01));
	errors += fabs(g.Ele(1)->Val(4, "Sxy") - ( 8.20878435e-01));
	errors += fabs(g.Ele(1)->Val(5, "Sxy") - (-1.47127394e+00));

	if (fabs(errors)>1.0e-7) cout << "[1;31mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;
	else                     cout << "[1;32mErrors(" << linsol << ") = " << errors << "[0m\n" << endl;

	// Return error flag
	if (fabs(errors)>1.0e-7) return 1;
	else                     return 0;
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
