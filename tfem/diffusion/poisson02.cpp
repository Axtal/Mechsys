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
using Util::_4;
using Util::_6;
using Util::_8s;

#define T boost::make_tuple

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

	double error  = 0.0;
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
		Array<Mesh::Block> bks(1);

		// Block # 0 --------------------------------
		Mesh::Verts_T ve0(4);
		Mesh::Edges_T ed0(4);
		Mesh::ETags_T et0(4);
		ve0 = T(0,0.0,0.0,0.0), T(1,1.0,0.0,0.0), T(2,1.0,1.0,0.0), T(3,0.0,1.0,0.0);
		ed0 = T(0,1), T(1,2), T(2,3), T(0,3);
		et0 = T(0,3,-10), T(1,2,-20), T(0,1,-30), T(2,3,-40);
		bks[0].Set   (-1, ve0, ed0, &et0, NULL, /*orig*/0, /*xplus*/1, /*yplus*/3);
		bks[0].SetNx (ndivs[k]);
		bks[0].SetNy (ndivs[k]);

		// Generate
		Mesh::Structured mesh(/*Is3D*/false);
		mesh.SetBlocks (bks);
		mesh.Generate  (true);

		////////////////////////////////////////////////////////////////////////////////////////// FEM /////

		// Data and Solver
		FEM::Data   dat (2);
		FEM::Solver sol (dat);
		sol.SetLinSol   (linsol.CStr());

		// Elements attributes
		FEM::EAtts_T eatts(1);
		eatts = T(-1, "Quad4", "Diffusion", "LinDiffusion", "k=1.0", "", "s=1.0", FNULL, true);

		// Set geometry: nodes and elements
		dat.SetNodesElems (&mesh, &eatts);

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
			dat.Ele(i)->CMatrix(0,Ke);
			double err_ke = 0.0;
			for (int i=0; i<4; ++i)
			for (int j=0; j<4; ++j)
				err_ke += fabs(Ke(i,j)-Ke_correct(i,j));
			if (err_ke>max_err_ke) max_err_ke = err_ke;
		}
		if (max_err_ke>1.0e-12) throw new Fatal("tex831: max_err_ke==%e for quadrangular mesh is bigger than %e.",max_err_ke,1.0e-12);

		// Stage # 1 -----------------------------------------------------------
		FEM::EBrys_T ebrys;
		ebrys.Push        (T(-10, "f", 0.0));
		ebrys.Push        (T(-30, "f", 0.0));
		ebrys.Push        (T(-20, "u", 0.0));
		ebrys.Push        (T(-40, "u", 0.0));
		dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
		dat.AddVolForces  ();
		sol.SolveWithInfo ();
		ndofs[k] = sol.nDOF();

		// Output: Nodes
		cout << endl;
		if (ndivs[k]<3)
		{
			cout << _6<<"Node #" << _8s<<"u" << _8s<<"f" << endl;
			for (size_t i=0; i<dat.NNodes(); ++i)
				cout << _6<<i << _8s<<dat.Nod(i)->Val("u") << _8s<<dat.Nod(i)->Val("f") << endl;
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
	error = max_err[0];

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
	if (error>1.24e-2)
	{
		cout << "[1;31mERROR too big:   error = " << _8s<<error << "[0m" << endl;
		return 1;
	}
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
