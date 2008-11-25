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
#include "fem/elems/quad4biot.h" // << plane strain
//#include "fem/elems/quad8pstrain.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::_8_4;
using boost::make_tuple;

void CallSolve (int iStage, FEM::Solver * Sol)
{
	double start = std::clock();
	Sol -> Solve();
	double total = std::clock() - start;
	double norm_resid = LinAlg::Norm(Sol->Resid());
	cout << "[1;31mStage # " << iStage << "[0m" << endl;
	cout << "\tTime elapsed (solution) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";
	cout << "\t[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
	cout << "\t[1;32mNumber of DOFs          = " << Sol->nDOF() << "[0m\n";
	// Output: VTU
	//Output o; o.VTU (&g, "tconsolid01.vtu");
	//cout << "[1;34mFile <tconsolid01.vtu> saved.[0m\n\n";
}

double Terz(double Z, double T) // Pore-pressure excess for one dimensional consolidation
{
	double ue=0.0;
	for (int m=0; m<100; m++)
	{
		double M = 3.14159*(2.0*m+1.0)/2.0;
		ue+=2.0/M*sin(M*Z)*exp(-M*M*T);
	}
	return ue;
}

int main(int argc, char **argv) try
{
	// Constants
	double W     = 1.0;    // Width
	double H     = 10.0;   // Height
	double E     = 5000.0; // Young
	double nu    = 0.25;   // Poisson
	double gw    = 10.0;   // GammaW
	double k     = 1.0e-5; // Isotropic permeability
	int    ndivy = 10;     // number of divisions along x and y
	bool   is_o2 = false;  // use high order elements?

	// More constants related with the one-dimensional consolidation
	double load  = -10;
	double mv    = (1+nu)*(1-2*nu)/(E*(1-nu));
	double cv    = k/(mv*gw);
	Vector<int>    SampleNodes(11);  // Nodes where pwp is evaluated
	SampleNodes = 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21;
	Vector<double> NormTimes(9); 
	NormTimes   = 0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0;  // List of normalized times

	// Calculate time increments
	Vector<double> TimeIncs(8);  // Time increments
	for (int i=0; i<TimeIncs.Size(); i++) 
		TimeIncs(i) = (NormTimes(i+1)-NormTimes(i))*pow(H,2)/cv;

	// Output information
	Matrix<double> OutPwp(SampleNodes.Size(), TimeIncs.Size());  OutPwp.SetValues(0.0);

	// Normalized depths of the sample nodes
	Vector<double> NDepth(SampleNodes.Size()); 
	for (int i=0; i<NDepth.Size(); i++) 
		NDepth(i) = (10.0 - i)/10.0;

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ndivy\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndivy =  atof(argv[2]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Blocks
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (false, 4,                // Is3D, NNodes
	             0.0,   W,  W, 0.0,       // x coordinates
	             0.0, 0.0,  H,   H);      // y coordinates
	b.SetNx     (1);                      // x weights and num of divisions along x
	b.SetNy     (ndivy);                  // y weights and num of divisions along y
	b.SetETags  (4,  -10, -20, -30, -40); // edge tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();                // Non-linear elements
	clock_t start = std::clock();           // Initial time
	size_t  ne    = mesh.Generate (blocks); // Discretize domain
	clock_t total = std::clock() - start;   // Time elapsed
	if (is_o2) cout << "\nNum of quadrangles (o2) = " << ne << endl;
	else       cout << "\nNumber of quadrangles   = " << ne << endl;
	cout << "Time elapsed (mesh)     = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Elements attributes
	String prms; prms.Printf("gw=%f E=%f nu=%f k=%f",gw,E,nu,k);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8Biot", "", prms.CStr(), "ZERO"));
	else       eatts.Push (make_tuple(-1, "Quad4Biot", "", prms.CStr(), "ZERO"));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);

	// Solver
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g)->SetLinSol("UM");

	// Edges boundaries
	FEM::EBrys_T ebrys;

	// Stage # 1
	sol->SetNumDiv(4)->SetDeltaTime(1.0);
	ebrys.Resize (0);
	ebrys.Push   (make_tuple(-10, "ux",    0.0));
	ebrys.Push   (make_tuple(-20, "ux",    0.0));
	ebrys.Push   (make_tuple(-30, "uy",    0.0));
	ebrys.Push   (make_tuple(-40, "fy",  -10.0));
	ebrys.Push   (make_tuple(-40, "pwp",   0.0));
	FEM::SetBrys (&mesh, NULL, &ebrys, NULL, &g);
	CallSolve    (1, sol);

	// Output: VTU 
	Output o; o.VTU (&g, "tbiot01_01.vtu");
	cout << "[1;34mFile <tbiot01_01.vtu> saved.[0m\n\n";

	// Stage # 2..
	for (int i=0; i<TimeIncs.Size(); i++)
	{
		cout << "Stage " << i+2 << endl;
		sol->SetNumDiv(10)->SetDeltaTime(TimeIncs(i));
		ebrys.Resize (0);
		ebrys.Push   (make_tuple(-10, "ux",    0.0));
		ebrys.Push   (make_tuple(-20, "ux",    0.0));
		ebrys.Push   (make_tuple(-30, "uy",    0.0));
		ebrys.Push   (make_tuple(-40, "pwp",   0.0));
		FEM::SetBrys (&mesh, NULL, &ebrys, NULL, &g);
		CallSolve    (2, sol);
		for (int j=0; j<SampleNodes.Size(); j++)
			OutPwp(j,i) = g.Nod(SampleNodes(j))->Val("pwp")/load; // Saving normalized pore-pressure
	}

	OutPwp.SetNS(Util::_8_4);
	cout << "OutPwp :" << endl << OutPwp << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Analytical values
	Matrix<double> AValues(SampleNodes.Size(), TimeIncs.Size());  
	for (int i=0; i<NDepth.Size(); i++)
		for (int j=0; j<TimeIncs.Size(); j++)
			AValues(i,j) = Terz(NDepth(i), NormTimes(j+1));
	
	AValues.SetNS(Util::_8_4);
	cout << "AValues:" << endl << AValues<< endl;

	// Test
	Array<double> err;

	for (int i=0; i<OutPwp.Rows(); i++)
		for (int j=0; j<OutPwp.Cols(); j++)
			err.Push(fabs(OutPwp(i,j)-AValues(i,j)));

	// Error summary
	double tol     = 1.0e-1;
	double min_err = err[err.Min()];
	double max_err = err[err.Max()];
	cout << _4<< ""    << _8s<<"Min"   << _8s<<"Mean"                                            << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Eps" << _8s<<min_err << _8s<<err.Mean() << (max_err>tol?"[1;31m":"[1;32m") << _8s<<max_err << "[0m" << _8s<<err.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err>tol) return 1;
	else return 0;


	// Output: VTU
	//o.VTU (&g, "tbiot01_02.vtu");
	//cout << "[1;34mFile <tbiot01_02.vtu> saved.[0m\n\n";

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
