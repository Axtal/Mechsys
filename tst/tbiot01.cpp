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

/*                            | | | | p
                              V V V V    
      8> @---------------------------@ <8     -
         |                    .      |        |
         |                    .      |        |
         |                    .      |        |
         |                    .      |        |
         |                    .      |        |
      8> |                    .      | <8     H
         |                    .      |        |
         |                    .      |        |
         |                    .      |        |
         |                    .      |        |
         |                    .      |        |
      8>@---------------------------@  <8     -
        /_\      /_\       /_\      /_\ 
         o        o         o        o  

         <---------- W ------------->

*/


// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/quad4.h"
#include "fem/elems/quad8.h"
#include "fem/biotelem.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::_8_4;
using Util::PI;
using boost::make_tuple;

// Biot f function
double f(double e)
{
	e = e + 0.000000001; // to avoid zero when e==0
	double  ee = e*e;
	double r_pi = sqrt(PI);
	return (1.0/(4.0*r_pi))*e*log(1.0 + 4.0/(PI*ee)) + 1.0/PI*atan(r_pi*e/2.0) + e/(2.0*r_pi*(3.24+ee));
}

// Normalized vertical displacement function
double Biot(double X, double T)
{
	double r_pi = sqrt(PI);
	// input values
	double b  = 1.0, p=1.0, v=0.0, E=1.0, cv=1.0;
	double l  = 2.0*b;
	double mv = (1.0+v)*(1.0-2.0*v)/(E*(1.0-v));
	double w_inf = mv*p*l/(4.0*r_pi); // quantity used to normalize
	double t = pow(T*l,2)/cv;
	double x = X*b; //
	double r_cvt = sqrt(cv*t);
	double ws = 2.0*mv*p*(r_cvt/r_pi)*(f((x+b)/r_cvt)-f((x-b)/r_cvt));
	return ws/w_inf;
}

int main(int argc, char **argv) try
{
	// Description:
	// Two dimensional analysis of the settlement of a footing
	// compared with the analytical solution obtained from Biot.

	// Constants
	double W     = 12.0;    // Width
	double H     = 12.0;    // Height
	double b     =  2.0;    // Load application length
	double E     = 10000.0; // Young
	double nu    = 0.0;     // Poisson
	double gw    = 10.0;    // GammaW
	double k     = 1.0e-6;  // Isotropic permeability
	int    ndivy = 12;      // number of divisions along x and y
	bool   is_o2 = false;   // use high order elements?

	// More constants related with the one-dimensional consolidation
	double load  = -100.0;
	double mv    = (1+nu)*(1-2*nu)/(E*(1-nu));
	double cv    = k/(mv*gw);
	double l     = 2*b;
	double r_pi  = sqrt(PI);
	double winf  = mv*load*l/(4*r_pi); // Quantity used to normalize the settlement
	Vector<int>    SampleNodes(19);  // Nodes where pwp is evaluated
	SampleNodes = 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21;
	SampleNodes = 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 239, 240, 241, 242, 243, 244, 245, 246;

	Vector<double> NormTimes(6); 
	NormTimes   = 0.0, 1.0/8.0, 2.0/8.0, 3.0/8.0, 4.0/8.0, 5.0/8.0; // List of normalized times

	// Calculate time increments
	Vector<double> TimeIncs(NormTimes.Size()-1);  // Time increments
	for (int i=0; i<TimeIncs.Size(); i++) 
		TimeIncs(i) = l*l/cv* (pow(NormTimes(i+1),2.0)-pow(NormTimes(i),2.0));

	// Output information
	Matrix<double> OutUy(SampleNodes.Size(), TimeIncs.Size());  
	Vector<double> Uy0  (SampleNodes.Size());  

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ndivy\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndivy =  atof(argv[2]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Block # 1
	Mesh::Block b1;
	b1.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b1.SetCoords (false, 4,               // Is3D, NNodes
	             0.0, W-b,  W-b, 0.0,     // x coordinates
	             0.0, 0.0,  H,   H);      // y coordinates
	b1.SetNx     (10);                    // x weights and num of divisions along x
	b1.SetNy     (ndivy);                 // y weights and num of divisions along y
	b1.SetETags  (4,  -10, 0, -30, -50);  // edge tags

	// Block # 2
	Mesh::Block b2;
	b2.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b2.SetCoords (false, 4,               // Is3D, NNodes
	             W-b,   W,  W, W-b,       // x coordinates
	             0.0, 0.0,  H,   H);      // y coordinates
	b2.SetNx     (8);                     // x weights and num of divisions along x
	b2.SetNy     (ndivy);                 // y weights and num of divisions along y
	b2.SetETags  (4,  0, -20, -30, -40);  // edge tags

	// Blocks
	Array<Mesh::Block*> blocks;
	blocks.Push (&b1);
	blocks.Push (&b2);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (blocks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Data dat(2); // 2D

	// Elements attributes
	String prms; prms.Printf("gw=%f E=%f nu=%f k=%f",gw,E,nu,k);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8", "Biot", "", prms.CStr(), "ZERO", "gam=20", true));
	else       eatts.Push (make_tuple(-1, "Quad4", "Biot", "", prms.CStr(), "ZERO", "gam=20", true));

	// Set geometry: nodes, elements, attributes, and boundaries
	dat.SetNodesElems (&mesh, &eatts);

	// Solver
	FEM::Solver sol(dat, "tbiot01");

	// Edges boundaries
	FEM::EBrys_T ebrys;

	// Stage # 0 --------------------------------------------------------------
	ebrys.Resize      (0);
	ebrys.Push        (make_tuple(-10, "ux",    0.0));
	ebrys.Push        (make_tuple(-20, "ux",    0.0));
	ebrys.Push        (make_tuple(-30, "uy",    0.0));
	ebrys.Push        (make_tuple(-40, "pwp",   0.0));
	ebrys.Push        (make_tuple(-50, "pwp",   0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	dat.AddVolForces  ();
	sol.SolveWithInfo (4, 1e+6, 0, "  Initial stress state due to self weight (zero displacements)\n");
	dat.ClearDisp     ();

	// Stage # 1 --------------------------------------------------------------
	ebrys.Resize      (0);
	ebrys.Push        (make_tuple(-10, "ux",    0.0));
	ebrys.Push        (make_tuple(-20, "ux",    0.0));
	ebrys.Push        (make_tuple(-30, "uy",    0.0));
	ebrys.Push        (make_tuple(-40, "fy",   load));
	ebrys.Push        (make_tuple(-40, "pwp",   0.0));
	ebrys.Push        (make_tuple(-50, "pwp",   0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo (4, 0.0001, 1, "  Apply surface (footing) loading\n");

	// Calculate displacements after first stage
	for (int i=0; i<SampleNodes.Size(); i++) 
		Uy0(i) = dat.Nod(SampleNodes(i))->Val("uy");

	// Stage # 2+ -------------------------------------------------------------
	for (int i=0; i<TimeIncs.Size(); i++)
	{
		ebrys.Resize      (0);
		ebrys.Push        (make_tuple(-10, "ux",    0.0));
		ebrys.Push        (make_tuple(-20, "ux",    0.0));
		ebrys.Push        (make_tuple(-30, "uy",    0.0));
		ebrys.Push        (make_tuple(-40, "pwp",   0.0));
		ebrys.Push        (make_tuple(-50, "pwp",   0.0));
		dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
		sol.SolveWithInfo (10, TimeIncs(i), i+2, "  Consolidation\n");
		for (int j=0; j<SampleNodes.Size(); j++)
			OutUy(j,i) = (dat.Nod(SampleNodes(j))->Val("uy") - Uy0(j))/(-winf); // Saving normalized vertical displacement
	}

	OutUy.SetNS(Util::_8_4);
	cout << "\nOutUy :" << endl << OutUy << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Normalized X coordinate of the sample nodes
	Vector<double> XNorm(SampleNodes.Size()); 
	for (int i=0; i<XNorm.Size(); i++) 
		XNorm(i) = (W-dat.Nod(SampleNodes(i))->X())/b;

	// Analytical values
	Matrix<double> AValues(SampleNodes.Size(), TimeIncs.Size());  
	for (int i=0; i<XNorm.Size(); i++)
		for (int j=0; j<TimeIncs.Size(); j++)
			AValues(i,j) = -Biot(XNorm(i), NormTimes(j+1));
	
	AValues.SetNS(Util::_8_4);
	cout << "AValues:" << endl << AValues<< endl;

	// Test
	Array<double> err;

	for (int i=0; i<OutUy.Rows(); i++)
		for (int j=0; j<OutUy.Cols(); j++)
			err.Push(fabs(OutUy(i,j)-AValues(i,j)));

	// Error summary
	double tol     = 1.0e-0;
	double min_err = err[err.Min()];
	double max_err = err[err.Max()];
	cout << _4<< ""    << _8s<<"Min"   << _8s<<"Mean"                                            << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Eps" << _8s<<min_err << _8s<<err.Mean() << (max_err>tol?"[1;31m":"[1;32m") << _8s<<max_err << "[0m" << _8s<<err.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err>tol) return 1;
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
