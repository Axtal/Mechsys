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

// Python
#include <Python.h>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/hex8.h"
#include "fem/elems/hex20.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "models/equilibs/pyequilib.h"
#include "models/equilibs/camclay.h"
#include "models/equilibs/toto.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// Input
	int  mdl   = 1;     // Model
	bool is_o2 = false; // use high order elements?
	cout << "Input: " << argv[0] << "  model(1:LinElastic, 2:CamClay, 3:PyEquilib, 4:Toto)  is_o2\n";
	if (argc>=2) mdl   =  atoi(argv[1]);
	if (argc>=3) is_o2 = (atoi(argv[2])>0 ? true : false);

	// Init Python
    Py_Initialize();

	// States
	LinAlg::Matrix<double> states(2,3); // Stress states
	states =  -2.0,  -2.0,  -2.0,
	          -1.08, -1.08, -3.84;

	// Number of divisions
	int gndiv = 20; // global number of divisions for states
	int sndiv = 20; // solver number of divisions

	// Parameters
	double E    = 1300.0;
	double lam  = 0.0778;
	double kap  = 0.00824;
	double phi  = 35.0;
	double nu   = 0.0;
	double vini = 1.699;

	double k1   = 10.0;
	double l1   = 0.1;
	double b1   = 1.0;
	double psi1 = 2.0;

	double k2   = 10.0;
	double l2   = 0.1;
	double b2   = 1.0;
	double psi2 = 2.0;

	double k3   = 10.0;
	double l3   = 0.1;
	double b3   = 0.1;
	double ev3  = -2.0;

	double k4   = 10.0;
	double l4   = 0.1;
	double b4   = 1.0;
	double ev4  = -2.0;

	// Mesh
	//int nxyz = 3;  // Divisions along x
	//int iele = 13; // Index of element for output (Centre)
	int nxyz = 1;
	int iele = 0;

	//############################################################################## Mesh

	/*
						4----------------7
					  ,'|              ,'|
					,'  |            ,'  |
				  ,'    | -6    -1 ,'    |
				,'      |        ,'      |
			  5'===============6'        |
			  |         |      |    -4   |
			  |    -3   |      |         |
			  |         0- - - | -  - - -3
			  |       ,'       |       ,'
			  |     ,' -2      |     ,'
			  |   ,'        -5 |   ,'
			  | ,'             | ,'
			  1----------------2'                */

	// Generate
	Mesh::Structured mesh(/*Is3D*/true);
	if (is_o2) mesh.SetO2();
	mesh.GenBox (nxyz,nxyz,nxyz);

	//############################################################################### FEM

	// Data and Solver
	FEM::Data   dat (3); // 3D
	FEM::Solver sol (dat,"ttriax");

	// Elements attributes
	double Sx = states(0,0);
	double Sy = states(0,1);
	double Sz = states(0,2);
	String geom; geom = (is_o2 ? "Hex20" : "Hex8");
	String stat; stat.Printf("Sx=%f Sy=%f Sz=%f v=%f",Sx,Sy,Sz,vini);
	String prms;
	FEM::EAtts_T eatts(1);
	switch (mdl)
	{
		case 1:
		{
			cout << "[1;34m ########################################## LinElastic ########################################## [0m" << endl;
			prms.Printf("E=%f nu=%f",E,nu);
			eatts = T(-1, geom.CStr(), "Equilib", "LinElastic", prms.CStr(), stat.CStr(), "gam=20", true);
			break;
		}
		case 2:
		{
			cout << "[1;34m ########################################## CamClay ########################################## [0m" << endl;
			prms.Printf("lam=%f kap=%f phics=%f nu=%f",lam,kap,phi,nu);
			eatts = T(-1, geom.CStr(), "Equilib", "CamClay", prms.CStr(), stat.CStr(), "gam=20", true);
			break;
		}
		case 3:
		{
			cout << "[1;34m ########################################## PyEquilib ########################################## [0m" << endl;
			prms.Printf("a0=%f a1=%f a2=%f a3=%f",lam,kap,phi,nu);
			eatts = T(-1, geom.CStr(), "Equilib", "PyEquilib,camclay.py", prms.CStr(), stat.CStr(), "gam=20", true);
			break;
		}
		case 4:
		{
			cout << "[1;34m ########################################## Toto ########################################## [0m" << endl;
			prms.Printf("k1=%f l1=%f b1=%f psi1=%f  k2=%f l2=%f b2=%f psi2=%f  k3=%f l3=%f b3=%f ev3=%f  k4=%f l4=%f b4=%f ev4=%f",k1,l1,b1,psi1, k2,l2,b2,psi2, k3,l3,b3,ev3, k4,l4,b4,ev4);
			eatts = T(-1, geom.CStr(), "Equilib", "Toto", prms.CStr(), stat.CStr(), "gam=20", true);
			break;
		}
		default: throw new Fatal("main: Model # %d is invalid. Please use 1:LinElastic, 2:CamClay, 3:PyEquilib, 4:Toto",mdl);
	}

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Output
	Array<size_t> eles;
	eles.Push      (iele);
	dat.SetOutEles (eles, "ttriax");
	dat.OutEles    (true); // true => only caption
	dat.OutEles    ();     // initial state

	// Solve
	int k = 0;
	for (int i=1; i<states.Rows(); ++i)
	{
		// increments
		double dfx = (states(i,0)-Sx)/gndiv;
		double dfy = (states(i,1)-Sy)/gndiv;
		double dfz = (states(i,2)-Sz)/gndiv;
		//cout << dfx << ", " << dfy << ", " << dfz << endl;
		for (int j=0; j<gndiv; ++j)
		{
			// solve
			FEM::FBrys_T fbrys; fbrys.Resize(6);
			fbrys = T(-1, "ux", 0.0), T(-3, "uy", 0.0), T(-5, "uz", 0.0),
			        T(-2, "fx", dfx), T(-4, "fy", dfy), T(-6, "fz", dfz);
			dat.SetBrys (&mesh, NULL, NULL, &fbrys);
			sol.SolveWithInfo (sndiv,0.0,k,"");
			// results
			dat.OutEles();
			// update state
			Sx += dfx;
			Sy += dfy;
			Sz += dfz;
			k++;
		}
	}

	// End
    Py_Finalize();
	return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }