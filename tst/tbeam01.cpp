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
#include "fem/elems/beam.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/output.h"
#include "fem/solvers/autome.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_6;
using Util::_8s;

inline double fy (double t)
{
	return 0.0;
}

int main(int argc, char **argv) try
{
	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);

	// Constants
	double M   = -20.0;   // kN*m
	double P   = -10.0;   // kN
	double L   =  1.0;    // m
	double E   =  2.1e+8; // kPa
	double A   =  4.0e-2; // m^2
	double Izz =  4.0e-4; // m^4

	// Geometry
	FEM::Data dat(2); // 2D

	// Nodes
	dat.SetNNodes (4);
	dat.SetNode   (0, 0.0,   L);
	dat.SetNode   (1,   L,   L);
	dat.SetNode   (2,   L, 0.0);
	dat.SetNode   (3, L+L, 0.0);

	// Elements
	dat.SetNElems (3);
	dat.SetElem   (0, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, dat.Nod(0))->Connect(1, dat.Nod(1));
	dat.SetElem   (1, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, dat.Nod(1))->Connect(1, dat.Nod(2));
	dat.SetElem   (2, "Beam", /*Active*/true, /*Tag*/-5)->Connect(0, dat.Nod(2))->Connect(1, dat.Nod(3));

	// Parameters and initial value
	String prms; prms.Printf("E=%f A=%f Izz=%f", E, A, Izz);
	dat.Ele(0)->SetModel("LinElastic", prms.CStr(), "ZERO");
	dat.Ele(1)->SetModel("LinElastic", prms.CStr(), "ZERO");
	dat.Ele(2)->SetModel("LinElastic", prms.CStr(), "ZERO");

	// Boundary conditions (must be after set connectivity)
	dat.Nod(0)->Bry("ux", 0.0)->Bry("uy", 0.0);
	dat.Nod(1)->Bry("fy", P);
	dat.Nod(3)->Bry("uy", 0.0)->Bry("mz", M);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&dat)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo();
	delete sol;

	// Output: VTU
	Output o; o.VTU (&dat, "tbeam01.vtu");
	cout << "\n[1;34mFile <tbeam01.vtu> saved.[0m\n\n";

	// Output: Nodes
	cout << "" << _6<<"Node #" << _8s<<"ux" << _8s<<"uy" << _8s<<"wz" << _8s<<"fx"<< _8s<<"fy" << _8s<<"mz" << endl;
	for (size_t i=0; i<dat.NNodes(); ++i)
	{
		if (dat.Nod(i)->HasVar("wz"))
			cout << _6<<i << _8s<<dat.Nod(i)->Val("ux") <<  _8s<<dat.Nod(i)->Val("uy") << _8s<<dat.Nod(i)->Val("wz") << _8s<<dat.Nod(i)->Val("fx") << _8s<<dat.Nod(i)->Val("fy") << _8s<<dat.Nod(i)->Val("mz") << endl;
	}
	cout << endl;

	// Output: Elements
	cout << _6<<"Elem #" << _8s<<"N0" << _8s<<"M0" << _8s<<"V0" << _8s<<"N1" << _8s<<"M1" << _8s<<"V1" << endl;
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		dat.Ele(i)->CalcDepVars();
		cout << _6<<i << _8s<<dat.Ele(i)->Val(0, "N") << _8s<<dat.Ele(i)->Val(0, "M") << _8s<<dat.Ele(i)->Val(0, "V");
		cout <<          _8s<<dat.Ele(i)->Val(1, "N") << _8s<<dat.Ele(i)->Val(1, "M") << _8s<<dat.Ele(i)->Val(1, "V") << endl;
	}
	cout << endl;

	//ExtraOut test
	Matrix<double> Coords;
	Vector<double> Norm;
	Matrix<double> Values;
	Array<String>  Labels;

	dat.Ele(2)->OutExtra(Coords, Norm, Values, Labels);
	cout << "Coords: \n"       << Coords << endl;
	cout << "Normal vector: "  << Norm   << endl;
	cout << "Values:\n"        << Values << endl;
	cout << "Labels: "         << Labels << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Displacements
	Array<double> err_u(12);
	err_u[ 0] = fabs(dat.Nod(0)->Val("ux") - (0.0));
	err_u[ 1] = fabs(dat.Nod(0)->Val("uy") - (0.0));
	err_u[ 2] = fabs(dat.Nod(0)->Val("wz") - (7.84722222e-5));
	err_u[ 3] = fabs(dat.Nod(1)->Val("ux") - (0.0));
	err_u[ 4] = fabs(dat.Nod(1)->Val("uy") - (6.85515873e-05));
	err_u[ 5] = fabs(dat.Nod(1)->Val("wz") - (4.87103175e-05));
	err_u[ 6] = fabs(dat.Nod(2)->Val("ux") - (1.89484127e-05));
	err_u[ 7] = fabs(dat.Nod(2)->Val("uy") - (7.03373016e-05));
	err_u[ 8] = fabs(dat.Nod(2)->Val("wz") - (-1.08134921e-05));
	err_u[ 9] = fabs(dat.Nod(3)->Val("ux") - (1.89484127e-05));
	err_u[10] = fabs(dat.Nod(3)->Val("uy") - (0.0));
	err_u[11] = fabs(dat.Nod(3)->Val("wz") - (-1.59623016e-04));

	// Forces
	Array<double> err_f(12);
	err_f[ 0] = fabs(dat.Nod(0)->Val("fx") - (0.0));
	err_f[ 1] = fabs(dat.Nod(0)->Val("fy") - (-5.0));
	err_f[ 2] = fabs(dat.Nod(0)->Val("mz") - (0.0));
	err_f[ 3] = fabs(dat.Nod(1)->Val("fx") - (0.0));
	err_f[ 4] = fabs(dat.Nod(1)->Val("fy") - (-10.0));
	err_f[ 5] = fabs(dat.Nod(1)->Val("mz") - (0.0));
	err_f[ 6] = fabs(dat.Nod(2)->Val("fx") - (0.0));
	err_f[ 7] = fabs(dat.Nod(2)->Val("fy") - (0.0));
	err_f[ 8] = fabs(dat.Nod(2)->Val("mz") - (0.0));
	err_f[ 9] = fabs(dat.Nod(3)->Val("fx") - (0.0));
	err_f[10] = fabs(dat.Nod(3)->Val("fy") - (15.0));
	err_f[11] = fabs(dat.Nod(3)->Val("mz") - (-20.0));

	// Stresses
	Array<double> err_s(18);
	err_s[ 0] = fabs(dat.Ele(0)->Val(0,"N") -     (  0.0));   err_s[ 9] = fabs(dat.Ele(0)->Val(1,"N") -     (  0.0));
	err_s[ 1] = fabs(dat.Ele(0)->Val(0,"M") - fabs(  0.0));   err_s[10] = fabs(dat.Ele(0)->Val(1,"M") - fabs( -5.0));
	err_s[ 2] = fabs(dat.Ele(0)->Val(0,"V") - fabs( -5.0));   err_s[11] = fabs(dat.Ele(0)->Val(1,"V") - fabs( -5.0));
	err_s[ 3] = fabs(dat.Ele(1)->Val(0,"N") -     (-15.0));   err_s[12] = fabs(dat.Ele(1)->Val(1,"N") -     (-15.0));
	err_s[ 4] = fabs(dat.Ele(1)->Val(0,"M") - fabs( -5.0));   err_s[13] = fabs(dat.Ele(1)->Val(1,"M") - fabs( -5.0));
	err_s[ 5] = fabs(dat.Ele(1)->Val(0,"V") - fabs(  0.0));   err_s[14] = fabs(dat.Ele(1)->Val(1,"V") - fabs(  0.0));
	err_s[ 6] = fabs(dat.Ele(2)->Val(0,"N") -     (  0.0));   err_s[15] = fabs(dat.Ele(2)->Val(1,"N") -     (  0.0));
	err_s[ 7] = fabs(dat.Ele(2)->Val(0,"M") - fabs( -5.0));   err_s[16] = fabs(dat.Ele(2)->Val(1,"M") - fabs(-20.0));
	err_s[ 8] = fabs(dat.Ele(2)->Val(0,"V") - fabs(-15.0));   err_s[17] = fabs(dat.Ele(2)->Val(1,"V") - fabs(-15.0));

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
