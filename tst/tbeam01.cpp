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
	double M   = 20.0;   // kN*m
	double P   = 10.0;   // kN
	double L   = 1.0;    // m
	double E   = 2.1e+8; // kPa
	double A   = 4.0e-2; // m^2
	double Izz = 4.0e-4; // m^4

	// Geometry
	FEM::Geom g(2); // 2D

	// Nodes
	g.SetNNodes (4);
	g.SetNode   (0, 0.0,   L);
	g.SetNode   (1,   L,   L);
	g.SetNode   (2,   L, 0.0);
	g.SetNode   (3, L+L, 0.0);

	// Elements
	g.SetNElems (3);
	g.SetElem   (0, "Beam")->Connect(0, g.Nod(0))->Connect(1, g.Nod(1));
	g.SetElem   (1, "Beam")->Connect(0, g.Nod(1))->Connect(1, g.Nod(2));
	g.SetElem   (2, "Beam")->Connect(0, g.Nod(2))->Connect(1, g.Nod(3));

	// Parameters and initial value
	String prms; prms.Printf("E=%f A=%f Izz=%f", E, A, Izz);
	g.Ele(0)->SetModel("LinElastic", prms.CStr(), "ZERO");
	g.Ele(1)->SetModel("LinElastic", prms.CStr(), "ZERO");
	g.Ele(2)->SetModel("LinElastic", prms.CStr(), "ZERO");

	// Boundary conditions (must be after set connectivity)
	g.Nod(0)->Bry("ux", 0.0)->Bry("uy", 0.0);
	g.Nod(1)->Bry("fy", -P);
	g.Nod(3)->Bry("uy", 0.0)->Bry("mz", -M);

	cout << g << endl;

	/*
	Matrix<double> Ke;  Ke.SetNS(Util::_6_3);
	g.Ele(2)->Order1Matrix(0, Ke);
	Ke = Ke/1.0e+6;
	cout << Ke << endl;
	*/

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol(linsol.CStr()) -> SetNumDiv(1) -> SetDeltaTime(0.0);
	sol -> Solve(); ///*tIni*/0.0, /*tFin*/1.0, /*hIni*/0.001, /*DtOut*/0.01);
	double norm_resid = LinAlg::Norm(sol->Resid());
	cout << "\n[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
	cout << "[1;32mNumber of DOFs          = " << sol->nDOF() << "[0m\n";
	if (norm_resid>1.0e-13) throw new Fatal("tex831: norm_resid=%e is bigger than %e.",norm_resid,1.3e-15);
	cout << endl;

	// Output: Nodes
	cout << _6<<"Node #" << _8s<<"ux" << _8s<<"uy" << _8s<<"wz" << _8s<<"fx"<< _8s<<"fy" << _8s<<"mz" << endl;
	for (size_t i=0; i<g.NNodes(); ++i)
		cout << _6<<i << _8s<<g.Nod(i)->Val("ux") <<  _8s<<g.Nod(i)->Val("uy") << _8s<<g.Nod(i)->Val("wz") << _8s<<g.Nod(i)->Val("fx") << _8s<<g.Nod(i)->Val("fy") << _8s<<g.Nod(i)->Val("mz") << endl;
	cout << endl;

	/*
	// Set number of prescribed and unknown DOFs
	size_t ndofs  = 0; // total
	size_t npdofs = 0; // prescribed
	size_t nudofs = 0; // unknown
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		if (g.Nod(i)->nSharedBy()>0)
		{
			for (size_t j=0; j<g.Nod(i)->nDOF(); ++j)
			{
				if (g.Nod(i)->DOFVar(j).IsEssenPresc) npdofs++;
				else                                  nudofs++;
				g.Nod(i)->DOFVar(j).EqID = ndofs;
				ndofs++;
			}
		}
	}

	Vector<double> u(ndofs);
	Vector<double> u1(ndofs);

	// Assemble matrices
	Matrix<double> Me;
	Array<size_t>  rmap, cmap;
	Array<bool>    rpre, cpre;
	Matrix<double> K(ndofs,ndofs);
	Matrix<double> M(ndofs,ndofs);
	K.SetValues(0.0);
	M.SetValues(0.0);
	for (size_t i=0; i<g.NElems(); ++i)
	{
		g.Ele(i)->Order1MatMap (0, rmap, cmap, rpre, cpre);
		g.Ele(i)->Order1Matrix (0, Ke);
		g.Ele(i)->Order2MatMap (0, rmap, cmap, rpre, cpre);
		g.Ele(i)->Order2Matrix (0, Me);
		for (int p=0; p<Ke.Rows(); ++p)
		for (int q=0; q<Ke.Cols(); ++q)
		{
			size_t I = rmap[p];
			size_t J = cmap[q];
			K(I,J) += Ke(p,q);
			M(I,J) += Me(p,q);
		}
	}

	// Assemble
	Vector<double> f(ndofs);

	Matrix<double> A(ndofs,ndofs);
	Matrix<double> C(ndofs,ndofs);
	Vector<double> v(ndofs);
	Vector<double> a(ndofs);
	Vector<double> r(ndofs);
	v.SetValues(0.0);
	C.SetValues(0.0);
	double Dt  = 0.05;
	double bet = 0.25;
	double gam = 0.5;
	double c1  = gam/(bet*Dt);
	double c2  = 1.0/(bet*Dt*Dt);
	A = K + c1*C + c2*M;
	r = f - C*v - K*u;

	a = inv(M)*dif;
	*/


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
