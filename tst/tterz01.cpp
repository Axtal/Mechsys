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
#include "fem/solver.h"
#include "fem/elems/quad4.h"
#include "fem/elems/quad8.h"
#include "fem/biotelem.h"
#include "models/equilibs/biotelastic.h"
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
		TimeIncs(i) = (NormTimes(i+1)-NormTimes(i))*pow(H,2.0)/cv;

	// Output information
	Matrix<double> OutPwp(SampleNodes.Size(), TimeIncs.Size());  // Pore-pressure values at the sample points obtained from de analises
	Vector<double> Pwp0  (SampleNodes.Size());  // Pore-pressure values before the load application

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
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (blocks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and solver
	FEM::Data   dat (2); // 2D
	FEM::Solver sol (dat,"tterz01");

	// Elements attributes
	FEM::EAtts_T eatts;
	String prms; prms.Printf("E=%f nu=%f k=%f", E,nu,k);
	String prps; prps.Printf("gam=20 gw=%f",    gw);
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8", "Biot", "BiotElastic", prms.CStr(), "ZERO", prps.CStr(), true));
	else       eatts.Push (make_tuple(-1, "Quad4", "Biot", "BiotElastic", prms.CStr(), "ZERO", prps.CStr(), true));

	// Set geometry: nodes, elements, attributes, and boundaries
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 0: Stage performed to approach a stationary condition ------------------------------
	FEM::EBrys_T ebrys;
	ebrys.Resize (0);
	ebrys.Push   (make_tuple(-10, "ux",    0.0));
	ebrys.Push   (make_tuple(-20, "ux",    0.0));
	ebrys.Push   (make_tuple(-30, "uy",    0.0));
	ebrys.Push   (make_tuple(-40, "pwp",   0.0));
	dat.SetBrys (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo(/*NDiv*/4, /*DTime*/1000000, /*iStage*/0);

	for (int i=0; i<SampleNodes.Size(); i++) 
		Pwp0(i) = dat.Nod(SampleNodes(i))->Val("pwp");

	// Stage # 1: Load application ------------------------------
	ebrys.Resize (0);
	ebrys.Push   (make_tuple(-10, "ux",    0.0));
	ebrys.Push   (make_tuple(-20, "ux",    0.0));
	ebrys.Push   (make_tuple(-30, "uy",    0.0));
	ebrys.Push   (make_tuple(-40, "fy",   load));
	ebrys.Push   (make_tuple(-40, "pwp",   0.0));
	dat.SetBrys (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo(/*NDiv*/4, /*DTime*/1.0, /*iStage*/1);

	// Stage # 2.. : Consolidation ------------------------------
	for (int i=0; i<TimeIncs.Size(); i++)
	{
		ebrys.Resize (0);
		ebrys.Push   (make_tuple(-10, "ux",    0.0));
		ebrys.Push   (make_tuple(-20, "ux",    0.0));
		ebrys.Push   (make_tuple(-30, "uy",    0.0));
		ebrys.Push   (make_tuple(-40, "pwp",   0.0));
		dat.SetBrys (&mesh, NULL, &ebrys, NULL);
		sol.SolveWithInfo(/*NDiv*/4, /*DTime*/TimeIncs(i), /*iStage*/i+2);
		for (int j=0; j<SampleNodes.Size(); j++)
			OutPwp(j,i) = (dat.Nod(SampleNodes(j))->Val("pwp")-Pwp0(j))/load; // Saving normalized excess of pore-pressure
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
