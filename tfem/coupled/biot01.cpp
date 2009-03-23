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
using Util::PI;

#define T boost::make_tuple

// Biot f function
double f(double e)
{
	e = e + 0.000000001; // to avoid zero when e==0
	double  ee = e*e;
	double r_pi = sqrt(PI);
	return (1.0/(4.0*r_pi))*e*log(1.0 + 4.0/(PI*ee)) + 1.0/PI*atan(r_pi*e/2.0) + e/(2.0*r_pi*(3.24+ee));
}

// Normalized vertical displacement function
double Biot(double X, double Time)
{
	double r_pi = sqrt(PI);
	// input values
	double b  = 1.0, p=1.0, v=0.0, E=1.0, cv=1.0;
	double l  = 2.0*b;
	double mv = (1.0+v)*(1.0-2.0*v)/(E*(1.0-v));
	double w_inf = mv*p*l/(4.0*r_pi); // quantity used to normalize
	double t = pow(Time*l,2)/cv;
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
	double gam   = 20.0;    // Specific weight
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

	// Blocks
	Array<Mesh::Block> bks(2);

	// Block # 0 --------------------------------
    Mesh::Verts_T ve0(4);
    Mesh::Edges_T ed0(4);
    Mesh::ETags_T et0(3);
    ve0 = T(0,0.0,0.0,0.0), T(1,W-b,0.0,0.0), T(2,W-b,H,0.0), T(3,0.0,H,0.0);
    ed0 = T(0,1), T(1,2), T(2,3), T(0,3);
    et0 = T(0,3,-10), T(0,1,-30), T(2,3,-50);
    bks[0].Set   (-1, ve0, ed0, &et0, NULL, /*orig*/0, /*xplus*/1, /*yplus*/3);
	bks[0].SetNx (10);
	bks[0].SetNy (ndivy);

	// Block # 1
    Mesh::Verts_T ve1(4);
    Mesh::Edges_T ed1(4);
    Mesh::ETags_T et1(3);
    ve1 = T(0,W-b,0.0,0.0), T(1,W,0.0,0.0), T(2,W,H,0.0), T(3,W-b,H,0.0);
    ed1 = T(0,1), T(1,2), T(2,3), T(0,3);
    et1 = T(1,2,-20), T(0,1,-30), T(2,3,-40);
    bks[1].Set   (-1, ve1, ed1, &et1, NULL, /*orig*/0, /*xplus*/1, /*yplus*/3);
	bks[1].SetNx (8);
	bks[1].SetNy (ndivy);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (bks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Data dat(2); // 2D

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f k=%f",E,nu,k);
	String prps; prps.Printf("gam=%f gw=%f",gam,gw);
	FEM::EAtts_T eatts(1);
	if (is_o2) eatts = T(-1, "Quad8", "Biot", "BiotElastic", prms.CStr(), "ZERO", prps.CStr(), FNULL, true);
	else       eatts = T(-1, "Quad4", "Biot", "BiotElastic", prms.CStr(), "ZERO", prps.CStr(), FNULL, true);

	// Set geometry: nodes, elements, attributes, and boundaries
	dat.SetNodesElems (&mesh, &eatts);

	// Solver
	FEM::Solver sol(dat, "tbiot01");

	// Edges boundaries
	FEM::EBrys_T ebrys;

	// Stage # 0 --------------------------------------------------------------
	ebrys.Resize      (0);
	ebrys.Push        (T(-10, "ux",    0.0));
	ebrys.Push        (T(-20, "ux",    0.0));
	ebrys.Push        (T(-30, "uy",    0.0));
	ebrys.Push        (T(-40, "pwp",   0.0));
	ebrys.Push        (T(-50, "pwp",   0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	dat.AddVolForces  ();
	sol.SolveWithInfo (4, 1e+6, 0, "  Initial stress state due to self weight (zero displacements)\n", /*ClearDisp*/true);

	// Stage # 1 --------------------------------------------------------------
	ebrys.Resize      (0);
	ebrys.Push        (T(-10, "ux",    0.0));
	ebrys.Push        (T(-20, "ux",    0.0));
	ebrys.Push        (T(-30, "uy",    0.0));
	ebrys.Push        (T(-40, "fy",   load));
	ebrys.Push        (T(-40, "pwp",   0.0));
	ebrys.Push        (T(-50, "pwp",   0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo (4, 0.0001, 1, "  Apply surface (footing) loading\n");

	// Calculate displacements after first stage
	for (int i=0; i<SampleNodes.Size(); i++) 
		Uy0(i) = dat.Nod(SampleNodes(i))->Val("uy");

	// Stage # 2+ -------------------------------------------------------------
	for (int i=0; i<TimeIncs.Size(); i++)
	{
		ebrys.Resize      (0);
		ebrys.Push        (T(-10, "ux",    0.0));
		ebrys.Push        (T(-20, "ux",    0.0));
		ebrys.Push        (T(-30, "uy",    0.0));
		ebrys.Push        (T(-40, "pwp",   0.0));
		ebrys.Push        (T(-50, "pwp",   0.0));
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
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
