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
#include <fstream>
#include <cmath>
#include <cfloat> // for DBL_EPSILON

// MechSys
#include "fem/data.h"
#include "fem/elems/tri3.h"
#include "fem/elems/quad4.h"
#include "fem/diffusionelem.h"
#include "models/diffusions/lindiffusion.h"
#include "fem/solver.h"
#include "fem/output.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Vector;
using boost::make_tuple;
using Util::_4;
using Util::_6;
using Util::_8s;

double u_correct(double s, double k, double x, double y)
{
	// s: source, k:conductivity
	// x,y: position
	size_t nterms = 200;
	double sum    = 0.0;
	for (size_t i=1; i<=nterms; ++i)
	{
		double alp = 0.5*(2.0*i-1.0)*Util::PI;
		sum += pow(-1.0,i)*cos(alp*y)*cosh(alp*x)/(pow(alp,3.0)*cosh(alp));
	}
	return 0.5*s*(1.0-y*y+4.0*sum)/k;
}

int main(int argc, char **argv) try
{
	// Constants
	size_t ndiv       = 2;     // number of divisions along x and y
	bool   check_conv = false; // check convergence ?
	String linsol("UM");       // linear solver: UMFPACK
	
	// Input
	cout << "Input: " << argv[0] << "  ndiv  check_conv  linsol(LA,UM,SLU)\n";
	if (argc>=2) ndiv        = atoi(argv[1]);
	if (argc>=3) check_conv  = atoi(argv[2]);
	if (argc>=4) linsol.Printf("%s",argv[3]);

	double err_1 = 0.0;
	if (check_conv==false)
	{
		cout << "\n[0;1m================================================================== Coarse triangular mesh[0m\n";

		////////////////////////////////////////////////////////////////////////////////////////// FEM /////

		// Data and Solver
		FEM::Data   dat (2);
		FEM::Solver sol (dat);

		// Nodes
		dat.SetNNodes (6);
		dat.SetNode   (0, 0.0, 0.0);
		dat.SetNode   (1, 0.5, 0.0);
		dat.SetNode   (2, 0.5, 0.5);
		dat.SetNode   (3, 1.0, 0.0);
		dat.SetNode   (4, 1.0, 0.5);
		dat.SetNode   (5, 1.0, 1.0);

		// Elements
		dat.SetNElems (4);
		dat.SetElem   (0, "Tri3", "Diffusion", true, -1);
		dat.SetElem   (1, "Tri3", "Diffusion", true, -1);
		dat.SetElem   (2, "Tri3", "Diffusion", true, -1);
		dat.SetElem   (3, "Tri3", "Diffusion", true, -1);

		// Set connectivity
		dat.Ele(0)->Connect(0, dat.Nod(0))->Connect(1, dat.Nod(1))->Connect(2, dat.Nod(2));
		dat.Ele(1)->Connect(0, dat.Nod(4))->Connect(1, dat.Nod(2))->Connect(2, dat.Nod(1));
		dat.Ele(2)->Connect(0, dat.Nod(1))->Connect(1, dat.Nod(3))->Connect(2, dat.Nod(4));
		dat.Ele(3)->Connect(0, dat.Nod(2))->Connect(1, dat.Nod(4))->Connect(2, dat.Nod(5));

		// Parameters and initial values
		dat.Ele(0)->SetModel("LinDiffusion", "k=1.0", "");
		dat.Ele(1)->SetModel("LinDiffusion", "k=1.0", "");
		dat.Ele(2)->SetModel("LinDiffusion", "k=1.0", "");
		dat.Ele(3)->SetModel("LinDiffusion", "k=1.0", "");
		
		// Properties (heat source)
		dat.Ele(0)->SetProps("s=1.0");
		dat.Ele(1)->SetProps("s=1.0");
		dat.Ele(2)->SetProps("s=1.0");
		dat.Ele(3)->SetProps("s=1.0");

		// Boundary conditions (must be after connectivity)
		dat.Ele(0)->EdgeBry("q", 0.0, 0)->EdgeBry("q", 0.0, 2);
		dat.Ele(2)->EdgeBry("q", 0.0, 0)->EdgeBry("u", 0.0, 1);
		dat.Ele(3)->EdgeBry("u", 0.0, 1)->EdgeBry("q", 0.0, 2);
		dat.Nod(3)->Bry("u",0.0);
		dat.Nod(5)->Bry("u",0.0);

		// Check conductivity matrices
		double err_ke = 0.0;
		LinAlg::Matrix<double> Ke0, Ke1, Ke2, Ke3;
		LinAlg::Matrix<double> Ke_correct;  Ke_correct.Resize(3,3);
		dat.Ele(0)->Order1Matrix(0,Ke0);
		dat.Ele(1)->Order1Matrix(0,Ke1);
		dat.Ele(2)->Order1Matrix(0,Ke2);
		dat.Ele(3)->Order1Matrix(0,Ke3);
		Ke_correct =  0.5, -0.5,  0.0,
					 -0.5,  1.0, -0.5,
					  0.0, -0.5,  0.5;
		for (int i=0; i<3; ++i)
		for (int j=0; j<3; ++j)
		{
			err_ke += fabs(Ke0(i,j)-Ke_correct(i,j));
			err_ke += fabs(Ke1(i,j)-Ke_correct(i,j));
			err_ke += fabs(Ke2(i,j)-Ke_correct(i,j));
			err_ke += fabs(Ke3(i,j)-Ke_correct(i,j));
		}
		if (err_ke>DBL_EPSILON) throw new Fatal("tex831: err_ke=%e for coarse triangular mesh is bigger than %e.",err_ke,DBL_EPSILON);

		// Solve
		sol.SetLinSol(linsol.CStr());
		sol.SolveWithInfo();

		// Output: Nodes
		cout << _6<<"Node #" << _8s<<"u" << _8s<<"q" << endl;
		for (size_t i=0; i<dat.NNodes(); ++i)
			cout << _6<<i << _8s<<dat.Nod(i)->Val("u") << _8s<<dat.Nod(i)->Val("q") << endl;

		//////////////////////////////////////////////////////////////////////////////////////// Check /////

		std::ofstream of("tex831_1.cal", std::ios::out);
		of << _8s<<"x" << _8s<<"u" << _8s<<"ucorr" << endl;
		Array<double> err_u;
		for (size_t i=0; i<dat.NNodes(); ++i)	
		{
			double x     = dat.Nod(i)->X();
			double y     = dat.Nod(i)->Y();
			double u     = dat.Nod(i)->Val("u");
			double ucorr = u_correct(1.0,1.0,x,y);
			err_u.Push ( fabs(u-ucorr) / (1.0+fabs(ucorr)) );
			if (fabs(y)<=1e-5) of << _8s<<x << _8s<<u << _8s<<ucorr << endl;
		}
		of.close();

		// Error summary
		double tol_u     = 1.0e-2;
		double min_err_u = err_u[err_u.Min()];
		double max_err_u = err_u[err_u.Max()];
		cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
		cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
		cout << endl;
		err_1 = max_err_u;
	}

	cout << endl;
	double err_2  = 0.0;
	size_t ntests = (check_conv ? 4 : 1);
	Array<size_t> ndivs;   ndivs  .Resize(ntests);
	Array<size_t> ndofs;   ndofs  .Resize(ntests);
	Array<double> max_err; max_err.Resize(ntests);
	if (check_conv)
	{
		ndivs[0] = 4;
		if (ntests>=2) ndivs[1] =  38;
		if (ntests>=3) ndivs[2] =  78;
		if (ntests>=4) ndivs[3] = 100;
		if (ntests>=5) ndivs[4] = 220;
	}
	else ndivs.SetValues(ndiv);
	for (size_t k=0; k<ntests; ++k)
	{
		cout << "[0;1m================================================================== Fine quadrangular mesh[0m\n";

		///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

		// Blocks
		Mesh::Block b;
		b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
		b.SetCoords (false, 4,           // Is3D, NNodes
					 0., 1., 1., 0.,     // x coordinates
					 0., 0., 1., 1.);    // y coordinates
		b.SetNx     (ndivs[k]);          // x weights and num of divisions along x
		b.SetNy     (ndivs[k]);          // y weights and num of divisions along y
		b.SetETags  (4,-10,-20,-30,-40); // edge tags
		Array<Mesh::Block*> blocks;
		blocks.Push (&b);

		// Generate
		Mesh::Structured mesh(/*Is3D*/false);
		mesh.SetBlocks (blocks);
		mesh.Generate  (true);

		////////////////////////////////////////////////////////////////////////////////////////// FEM /////

		// Data and Solver
		FEM::Data   dat (2);
		FEM::Solver sol (dat);

		// Elements attributes
		FEM::EAtts_T eatts;
		eatts.Push (make_tuple(-1, "Quad4", "Diffusion", "LinDiffusion", "k=1.0", "", "s=1.0", true));

		// Set geometry: nodes and elements
		dat.SetNodesElems (&mesh, &eatts);

		// Edges brys (the order matters!)
		FEM::EBrys_T ebrys;
		ebrys.Push  (make_tuple(-10, "q", 0.0));
		ebrys.Push  (make_tuple(-30, "q", 0.0));
		ebrys.Push  (make_tuple(-20, "u", 0.0));
		ebrys.Push  (make_tuple(-40, "u", 0.0));
		dat.SetBrys (&mesh, NULL, &ebrys, NULL);

		// Check conductivity matrices
		double max_err_ke = 0.0;
		LinAlg::Matrix<double> Ke_correct;  Ke_correct.Resize(4,4);
		Ke_correct =  4.0, -1.0, -2.0, -1.0,
		             -1.0,  4.0, -1.0, -2.0,
		             -2.0, -1.0,  4.0, -1.0,
		             -1.0, -2.0, -1.0,  4.0;
		Ke_correct = (1.0/6.0)*Ke_correct;
		for (size_t i=0; i<dat.NElems(); ++i)
		{
			LinAlg::Matrix<double> Ke;
			dat.Ele(i)->Order1Matrix(0,Ke);
			double err_ke = 0.0;
			for (int i=0; i<4; ++i)
			for (int j=0; j<4; ++j)
				err_ke += fabs(Ke(i,j)-Ke_correct(i,j));
			if (err_ke>max_err_ke) max_err_ke = err_ke;
		}
		if (max_err_ke>1.0e-12) throw new Fatal("tex831: max_err_ke==%e for quadrangular mesh is bigger than %e.",max_err_ke,1.0e-12);

		// Solve
		sol.SetLinSol(linsol.CStr());
		sol.SolveWithInfo();
		ndofs[k] = sol.nDOF();

		// Output: Nodes
		cout << endl;
		if (ndivs[k]<3)
		{
			cout << _6<<"Node #" << _8s<<"u" << _8s<<"q" << endl;
			for (size_t i=0; i<dat.NNodes(); ++i)
				cout << _6<<i << _8s<<dat.Nod(i)->Val("u") << _8s<<dat.Nod(i)->Val("q") << endl;
			cout << endl;
		}

		//////////////////////////////////////////////////////////////////////////////////////// Check /////

		// Check
		std::ofstream of("tex831_2.cal", std::ios::out);
		of << _8s<<"x" << _8s<<"u" << _8s<<"ucorr" << endl;
		Array<double> err_u;
		for (size_t i=0; i<dat.NNodes(); ++i)	
		{
			double x     = dat.Nod(i)->X();
			double y     = dat.Nod(i)->Y();
			double u     = dat.Nod(i)->Val("u");
			double ucorr = u_correct(1.0,1.0,x,y);
			err_u.Push ( fabs(u-ucorr) / (1.0+fabs(ucorr)) );
			if (fabs(y)<=1e-5) of << _8s<<x << _8s<<u << _8s<<ucorr << endl;
		}
		of.close();

		// Error summary
		double tol_u     = 1.0e-2;
		double min_err_u = err_u[err_u.Min()];
		double max_err_u = err_u[err_u.Max()];
		cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
		cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
		cout << endl;

		// Error -- convergence
		max_err[k] = max_err_u;
	}
	err_2 = max_err[0];

	// Check convergence
	if (check_conv)
	{
		// Janicke, L. & Kost, A. (1999). Convergence properties of the finite element solutions. IEEE Transactions on Magnetics, 35(3), 1414-1417.
		cout << _6<<"nDOF1" << _6 <<"nDOF2" << _8s<<"conv" << endl;
		for (size_t k=1; k<ntests; ++k)
		{
			double conv = -2.0*(log(max_err[k-1])-log(max_err[k]))/(log(ndofs[k-1])-log(ndofs[k]));
			cout << _6<<ndofs[k-1] << _6<<ndofs[k] << "[1;35m" << _8s<<conv << "[0m\n";
		}
		// First -> Last tests
		double conv = -2.0*(log(max_err[0])-log(max_err[ntests-1]))/(log(ndofs[0])-log(ndofs[ntests-1]));
		cout << _6<<ndofs[0] << _6<<ndofs[ntests-1] << "[1;35m" << _8s<<conv << "[0m\n";
	}

	// Return error flag
	if (err_1>1.38e-2 || err_2>1.24e-2)
	{
		cout << "[1;31mERROR too big: err_1 = " << _8s<<err_1 << ",    err_2 = " << _8s<<err_2 << "[0m" << endl;
		return 1;
	}
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
