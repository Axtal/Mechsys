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
#include "fem/elems/quad4pstrain.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/solvers/autome.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

#include "fem/embedded.h"        // << embedded
#include "fem/elems/rod3.h"        // << embedded
#include "fem/elems/embspring.h"        // << embedded

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
	int ndivx = 1;
	int ndivy = 1;
	double H  = 1.0;
	double W  = 2.0*H;
	bool is_o2 = false;

	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Block # 1
	Mesh::Block b1;
	b1.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b1.SetCoords (false, 4,                // Is3D, NNodes
	             0.0, W/2.0,  W/2.0, 0.0,  // x coordinates
	             0.0,   0.0,      H,   H); // y coordinates
	b1.SetNx     (ndivx);                  // x weights and num of divisions along x
	b1.SetNy     (ndivy);                  // y weights and num of divisions along y
	b1.SetETags  (4, 0, 0, -30, 0);   // edge tags

	// Block # 2
	Mesh::Block b2;
	b2.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b2.SetCoords (false, 4,               // Is3D, NNodes
	             W/2.0,   W,  W, W/2.0,   // x coordinates
	             0.0,   0.0,  H,    H);   // y coordinates
	b2.SetNx     (ndivx);                 // x weights and num of divisions along x
	b2.SetNy     (ndivy);                 // y weights and num of divisions along y
	b2.SetETags  (4, 0, 0, -30, 0);      // edge tags

	// Blocks
	Array<Mesh::Block*> blocks;
	blocks.Push (&b1);
	blocks.Push (&b2);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();                // Non-linear elements
	mesh.SetBlocks (blocks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////
	
	// weight
	double P = 10.0;

	// Geometry
	FEM::Geom g(2); // 2D

	// Nodes brys
	FEM::NBrys_T nbrys;
	nbrys.Push (make_tuple(  0.0, 0.0, 0.0, "ux", 0.0)); // x,y,z, key, val
	nbrys.Push (make_tuple(  0.0, 0.0, 0.0, "uy", 0.0)); // x,y,z, key, val
	nbrys.Push (make_tuple(    W, 0.0, 0.0, "uy", 0.0)); // x,y,z, key, val
	nbrys.Push (make_tuple(W/2.0,   H, 0.0, "fy",  -P)); // x,y,z, key, val

	// Edges brys
	FEM::EBrys_T ebrys;
	//ebrys.Push (make_tuple(-30, "uy", 0.0)); // tag, key, val

	// Elements attributes
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", "E=1.0 nu=0.0", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0", "", true)); // tag, type, model, prms, inis, props
	else       eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", "E=1.0 nu=0.0", "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0", "", true)); // tag, type, model, prms, inis, props

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g);
	FEM::SetBrys       (&mesh, &nbrys, NULL, NULL, &g);

	// AddReinforcement
	LinAlg::Vector<double> P0(3); P0=0.0,0.0,0.0;
	LinAlg::Vector<double> P1(3); P1=1.0,1.0,0.0;
	LinAlg::Vector<double> P2(3); P2=2.0,0.0,0.0;

	AddReinf(0.0, 0.0, 0.0, 1.0, 1.0, 0.0, "E=1.0E6 A=1 K=1E12", true, -10, &g);
	AddReinf(1.0, 1.0, 0.0, 2.0, 0.0, 0.0, "E=1.0E6 A=1 K=1E12", true, -20, &g);
	AddReinf(0.0, 0.0, 0.0, 2.0, 0.0, 0.0, "E=1.0E6 A=1 K=1E12", true, -30, &g);

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

		     if (tag==-10 && nn<=3) err.Push(fabs(g.Ele(i)->Val(0, "Sa") - (-0.707106781*P)));
		else if (tag==-20 && nn<=3) err.Push(fabs(g.Ele(i)->Val(0, "Sa") - (-0.707106781*P)));
		else if (tag==-30 && nn<=3) err.Push(fabs(g.Ele(i)->Val(0, "Sa") - ( 0.500000000*P)));
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
