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

/*  Small truss
                              1.0 ^
                                  |
                                  |2
                                  o----> 2.0
                                ,'|
                              ,'  |
                   E=200    ,'    |
                   A=SQ2  ,'      | E=50
                    [2] ,'        | A=1
                      ,'          | [1]
                    ,'            |
                  ,'              |
   y            ,'                |
   |         0,'        [0]       |1
   |         o--------------------o
  (z)___x   /_\        E=100     /_\
           ////        A=1       ___  
*/

// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/elems/rod.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "util/exception.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_6;
using Util::_8s;

int main(int argc, char **argv) try
{
	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);

	// Geometry
	FEM::Data dat(2); // 2D

	// Nodes
	dat.SetNNodes (3);
	dat.SetNode   (0,  0.0,  0.0);
	dat.SetNode   (1, 10.0,  0.0);
	dat.SetNode   (2, 10.0, 10.0);

	// Elements
	dat.SetNElems (3);
	dat.SetElem   (0, "Rod", /*IsActive*/true, /*Tag*/-1);
	dat.SetElem   (1, "Rod", /*IsActive*/true, /*Tag*/-1);
	dat.SetElem   (2, "Rod", /*IsActive*/true, /*Tag*/-1);

	// Set connectivity
	dat.Ele(0)->Connect(0, dat.Nod(0))->Connect(1, dat.Nod(1));
	dat.Ele(1)->Connect(0, dat.Nod(1))->Connect(1, dat.Nod(2));
	dat.Ele(2)->Connect(0, dat.Nod(0))->Connect(1, dat.Nod(2));

	// Parameters and initial value
	dat.Ele(0)->SetModel("", "E=100.0 A=1.0"              , "Sx=0.0");
	dat.Ele(1)->SetModel("", "E= 50.0 A=1.0"              , "Sx=0.0");
	dat.Ele(2)->SetModel("", "E=200.0 A=1.414213562373095", "Sx=0.0");

	// Boundary conditions (must be after set connectivity)
	dat.Nod(0)->Bry("ux", 0.0)->Bry("uy", -0.5); // Essential
	dat.Nod(1)->                Bry("uy",  0.4); // Essential
	dat.Nod(2)->Bry("fx", 2.0)->Bry("fy",  1.0); // Natural

	// Check stiffness matrices
	double err_ke = 0.0;
	Array<size_t>  map;
	Array<bool>    pre;
	Matrix<double> Ke0,  Ke1,  Ke2;
	Matrix<double> Ke0c, Ke1c, Ke2c; // correct matrices
	Ke0c.Resize(4,4);
	Ke1c.Resize(4,4);
	Ke2c.Resize(4,4);
	dat.Ele(0)->Order1Matrix(0,Ke0);
	dat.Ele(1)->Order1Matrix(0,Ke1);
	dat.Ele(2)->Order1Matrix(0,Ke2);
	Ke0c =  10.0,   0.0, -10.0,   0.0,
	         0.0,   0.0,   0.0,   0.0,
	       -10.0,   0.0,  10.0,   0.0,
	         0.0,   0.0,   0.0,   0.0;
	Ke1c =   0.0,   0.0,   0.0,   0.0,
	         0.0,   5.0,   0.0,  -5.0,
	         0.0,   0.0,   0.0,   0.0,
	         0.0,  -5.0,   0.0,   5.0;
	Ke2c =  10.0,  10.0, -10.0, -10.0,
	        10.0,  10.0, -10.0, -10.0,
	       -10.0, -10.0,  10.0,  10.0,
	       -10.0, -10.0,  10.0,  10.0;
	for (int i=0; i<4; ++i)
	for (int j=0; j<4; ++j)
	{
		err_ke += fabs(Ke0(i,j)-Ke0c(i,j));
		err_ke += fabs(Ke1(i,j)-Ke1c(i,j));
		err_ke += fabs(Ke2(i,j)-Ke2c(i,j));
	}
	if (err_ke>8.55e-14) throw new Fatal("ttruss01: err_ke=%e is bigger than %e.",err_ke,8.55e-14);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&dat)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);
	delete sol;

	// Output: Nodes
	cout << _6<<"Node #" << _8s<<"ux" << _8s<<"uy" << _8s<<"fx"<< _8s<<"fy" << endl;
	for (size_t i=0; i<dat.NNodes(); ++i)
		cout << _6<<i << _8s<<dat.Nod(i)->Val("ux") <<  _8s<<dat.Nod(i)->Val("uy") << _8s<<dat.Nod(i)->Val("fx") << _8s<<dat.Nod(i)->Val("fy") << endl;
	cout << endl;

	// Output: Elements
	cout << _6<<"Elem #" << _8s<<"Sa(left)" << _8s<<"Sa(right)" << _8s<<"Ea(left)" << _8s<<"Ea(right)" << endl;
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		dat.Ele(i)->CalcDepVars();
		cout << _6<<i;
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j) cout << _8s<<dat.Ele(i)->Val(j, "Sa");
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j) cout << _8s<<dat.Ele(i)->Val(j, "Ea");
		cout << endl;
	}
	cout << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Displacements
	Array<double> err_u(6);
	err_u[0] =  fabs(dat.Nod(0)->Val("ux") - ( 0.0));
	err_u[1] =  fabs(dat.Nod(0)->Val("uy") - (-0.5));
	err_u[2] =  fabs(dat.Nod(1)->Val("ux") - ( 0.0));
	err_u[3] =  fabs(dat.Nod(1)->Val("uy") - ( 0.4));
	err_u[4] =  fabs(dat.Nod(2)->Val("ux") - (-0.5));
	err_u[5] =  fabs(dat.Nod(2)->Val("uy") - ( 0.2));

	// Forces
	Array<double> err_f(6);
	err_f[0] = fabs(dat.Nod(0)->Val("fx") - (-2.0));
	err_f[1] = fabs(dat.Nod(0)->Val("fy") - (-2.0));
	err_f[2] = fabs(dat.Nod(1)->Val("fx") - ( 0.0));
	err_f[3] = fabs(dat.Nod(1)->Val("fy") - ( 1.0));
	err_f[4] = fabs(dat.Nod(2)->Val("fx") - ( 2.0));
	err_f[5] = fabs(dat.Nod(2)->Val("fy") - ( 1.0));

	// Correct axial normal stresses
	Array<double> err_s(dat.NElems());
	err_s.SetValues(0.0);
	err_s[1] = fabs(dat.Ele(1)->Val(0, "N") - (-1.0));

	// Error summary
	double tol_u     = DBL_EPSILON;
	double tol_f     = 8.89e-16;
	double tol_s     = 4.4e-4;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	double min_err_f = err_f[err_f.Min()];
	double max_err_f = err_f[err_f.Max()];
	double min_err_s = err_s[err_s.Min()];
	double max_err_s = err_s[err_s.Max()];
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << _4<< "f" << _8s<<min_err_f << _8s<<err_f.Mean() << (max_err_f>tol_f?"[1;31m":"[1;32m") << _8s<<max_err_f << "[0m" << _8s<<err_f.Norm() << endl;
	cout << _4<< "N" << _8s<<min_err_s << _8s<<err_s.Mean() << (max_err_s>tol_s?"[1;31m":"[1;32m") << _8s<<max_err_s << "[0m" << _8s<<err_s.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_u>tol_u || max_err_f>tol_f || max_err_s>tol_s) return 1;
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
