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

/*
                        | P                                                   | P = -1.0
         3            2 V              5                                      V
      8> @-----------------------------@       -                              .
         |            .'|'.            |       |                            .' '.
         |          .'  |  '.          |       |                          .'     '.
         |        .'    |    '.        |       |                        .'         '.
         |      .'      |      '.      |       |  H                   .'             '.
         |    .'        |        '.    |       |                    .'                 '.
         |  .'         1|          '.  |       |                  .'                     '.
        0|.'............|............'.|4      |                .'.........................'.
      8>@------------------------------@       -               /_\                         /_\
        /_\                           /_\                      ///                          o
        ///                            o

         <------------ W ------------->

*/


// STL
#include <iostream>

// MechSys
#include "fem/geometry.h"
#include "fem/functions.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"
#include "fem/embedded.h"
#include "fem/elems/hex8equilib.h"
#include "fem/elems/rod3.h"
#include "fem/elems/embspring.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::_8_4;
using Util::PI;
using boost::make_tuple;

int main(int argc, char **argv) try
{
	int ndivz = 1;
	int ndivx = 2*ndivz;
	int ndivy = 2*ndivz;
	double H  = 1.0;
	double W  = 2.0*H;
	bool is_o2 = false;

	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Block # 1
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (true, 8,                 // Is3D, NNodes
	             0.0,   W,   W, 0.0, 0.0,   W, W, 0.0,  // x coordinates
	             0.0, 0.0,   W,   W, 0.0, 0.0, W,   W,  // y coordinates
	             0.0, 0.0, 0.0, 0.0,   H,   H, H,   H); // y coordinates
	b.SetNx     (ndivx);                  // x weights and num of divisions along x
	b.SetNy     (ndivy);                  // y weights and num of divisions along y
	b.SetNz     (ndivz);                  // z weights and num of divisions along z
	b.SetFTags  (6, 0, 0, 0, 0, -10, 0);  // face tags

	// Blocks
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	Mesh::Structured mesh(/*Is3D*/true);
	mesh.SetBlocks (blocks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////
	
	// weight
	double P = 1.0;

	// Geometry
	FEM::Geom g(3); // 2D

	// Nodes brys
	FEM::NBrys_T nbrys;
	nbrys.Push (make_tuple(2.0, 0.0, 0.0, "ux", 0.0)); // x,y,z, key, val
	nbrys.Push (make_tuple(2.0, 0.0, 0.0, "uy", 0.0)); // x,y,z, key, val
	//nbrys.Push (make_tuple(  0.0, 0.0, 0.0, "uy", 0.0)); // x,y,z, key, val
	//nbrys.Push (make_tuple(    W, 0.0, 0.0, "uy", 0.0)); // x,y,z, key, val
	nbrys.Push (make_tuple(W/2.0,   W/2.0, H, "fz",  -P)); // x,y,z, key, val

	// Faces brys
	FEM::FBrys_T fbrys;
	fbrys.Push (make_tuple(-10, "uz", 0.0)); // tag, key, val

	// Elements attributes
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Hex8Equilib", "LinElastic", "E=1.0 nu=0.0", "", "", true)); // tag, type, model, prms, inis, props
	else       eatts.Push (make_tuple(-1, "Hex8Equilib", "LinElastic", "E=1.0 nu=0.0", "", "", true)); // tag, type, model, prms, inis, props

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);

	// Add reinforcements
	AddReinf (2.0, 0.0, 0.0, 1.0, 1.0, 1.0, "E=1.0E8 A=0.1 K=1E12", true, -10, &g);
	AddReinf (1.0, 1.0, 1.0, 0.0, 2.0, 0.0, "E=1.0E8 A=0.1 K=1E12", true, -20, &g);
	AddReinf (2.0, 0.0, 0.0, 0.0, 2.0, 0.0, "E=1.0E8 A=0.1 K=1E12", true, -30, &g);

	// Boundary conditions
	FEM::SetBrys (&mesh, &nbrys, NULL, &fbrys, &g);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol->SetGeom(&g)->SetLinSol(linsol.CStr());
	sol->SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);
	delete sol;

	Output out;
	out.VTU(&g, "out.vtu");

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Test
	Array<double> err;

	for (size_t i=0; i<g.NElems(); i++)
	{
		int tag = g.Ele(i)->Tag();
		int  nn = g.Ele(i)->NNodes();

		     if (tag==-10 && nn<=3) { err.Push(fabs(g.Ele(i)->Val(0, "Sa") - (-0.866025404*P*10))); }
		else if (tag==-20 && nn<=3) { err.Push(fabs(g.Ele(i)->Val(0, "Sa") - (-0.866025404*P*10))); }
		else if (tag==-30 && nn<=3) { err.Push(fabs(g.Ele(i)->Val(0, "Sa") - ( 0.707106781*P*10))); }
	}

	// Error summary
	double tol     = 1.0e-2;
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

