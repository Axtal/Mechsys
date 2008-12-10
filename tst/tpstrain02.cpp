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

/*       | | | | | | | q
         V V V V V V V    
         @-----------@
         |           |
         |           |
         |           | H
         |           |
         |     L     |
         @-----@-----@
        /_\   / \   /_\
         o    ///    o
*/

// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
//#include "fem/elems/quad4pstrain.h"
//#include "fem/elems/quad8pstrain.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using boost::make_tuple;

int main(int argc, char **argv) try
{
	// Constants
	double E     = 207.0; // Young
	double nu    = 0.3;   // Poisson
	double q     = -1.0;  // Load
	int    ndiv  = 10;    // number of divisions along x and y
	bool   is_o2 = false; // use high order elements?
	String linsol("UM");  // UMFPACK

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ndiv  linsol(LA,UM,SLU)\n";
	if (argc>=2) is_o2      = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndiv       =  atof(argv[2]);
	if (argc>=4) linsol.Printf("%s",argv[3]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Blocks
	Mesh::Block b;
	b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b.SetCoords (false, 4,            // Is3D, NNodes
	             0.0, 1.0, 1.0, 0.0,  // x coordinates
	             0.0, 0.0, 1.0, 1.0); // y coordinates
	b.SetNx     (ndiv);               // x weights and num of divisions along x
	b.SetNy     (ndiv);               // y weights and num of divisions along y
	b.SetETags  (4,  0, 0, -10, -20); // edge tags
	Array<Mesh::Block*> blocks;
	blocks.Push (&b);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	mesh.SetBlocks (blocks);           // Set Blocks
	if (is_o2) mesh.SetO2();           // Non-linear elements
	mesh.Generate (/*WithInfo*/ true); // Discretize domain

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Data dat(2); // 2D

	// Nodes brys
	FEM::NBrys_T nbrys;
	nbrys.Push (make_tuple(0.5, 0.0, 0.0, "ux", 0.0)); // x,y,z, key, val

	// Edges brys
	FEM::EBrys_T ebrys;
	ebrys.Push (make_tuple(-10, "uy", 0.0)); // tag, key, val
	ebrys.Push (make_tuple(-20, "fy",   q)); // tag, key, val

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", true));
	else       eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", true));

	// Set geometry: nodes, elements, attributes, and boundaries
	dat.SetNodesElems (&mesh, &eatts);
	dat.SetBrys       (&mesh, &nbrys, &ebrys, NULL);

	// Solve
	FEM::Solver sol(dat,"tpstrain02");
	sol.SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
    Array<double> err_eps;
    Array<double> err_sig;
    Array<double> err_dis;

	double Sx  = 0.0;
	double Sy  = q;
	double Ex  = -nu*(1.0+nu)*Sy/E;
	double Ey  =  (1.0-nu*nu)*Sy/E;
	double Ez  = 0.0;
	double Exy = 0.0;
	double Sz  = (E/(1.0+nu))*(nu/(1.0-2.0*nu))*(Ex+Ey);
	double Sxy = 0.0;

	// Stress and strains
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j)
		{
			err_eps.Push ( fabs(dat.Ele(i)->Val(j,"Ex" ) - Ex ) / (1.0+fabs(Ex )) );
			err_eps.Push ( fabs(dat.Ele(i)->Val(j,"Ey" ) - Ey ) / (1.0+fabs(Ey )) );
			err_eps.Push ( fabs(dat.Ele(i)->Val(j,"Ez" ) - Ez ) / (1.0+fabs(Ez )) );
			err_eps.Push ( fabs(dat.Ele(i)->Val(j,"Exy") - Exy) / (1.0+fabs(Exy)) );
			err_sig.Push ( fabs(dat.Ele(i)->Val(j,"Sx" ) - Sx ) / (1.0+fabs(Sx )) );
			err_sig.Push ( fabs(dat.Ele(i)->Val(j,"Sy" ) - Sy ) / (1.0+fabs(Sy )) );
			err_sig.Push ( fabs(dat.Ele(i)->Val(j,"Sz" ) - Sz ) / (1.0+fabs(Sz )) );
			err_sig.Push ( fabs(dat.Ele(i)->Val(j,"Sxy") - Sxy) / (1.0+fabs(Sxy)) );
		}
	}

	// Displacements
	for (size_t i=0; i<dat.NNodes(); ++i)
	{
		double ux_correct = Ex*(dat.Nod(i)->X()-0.5);
		double uy_correct = Ey* dat.Nod(i)->Y();
		err_dis.Push ( fabs(dat.Nod(i)->Val("ux") - ux_correct) / (1.0+fabs(ux_correct)) );
		err_dis.Push ( fabs(dat.Nod(i)->Val("uy") - uy_correct) / (1.0+fabs(uy_correct)) );
	}

	// Error summary
	double tol_eps     = 1.0e-16;
	double tol_sig     = 1.0e-14;
	double tol_dis     = 1.0e-16;
	double min_err_eps = err_eps[err_eps.Min()];
	double min_err_sig = err_sig[err_sig.Min()];
	double min_err_dis = err_dis[err_dis.Min()];
	double max_err_eps = err_eps[err_eps.Max()];
	double max_err_sig = err_sig[err_sig.Max()];
	double max_err_dis = err_dis[err_dis.Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Eps" << _8s<<min_err_eps << _8s<<err_eps.Mean() << (max_err_eps>tol_eps?"[1;31m":"[1;32m") << _8s<<max_err_eps << "[0m" << _8s<<err_eps.Norm() << endl;
	cout << _4<< "Sig" << _8s<<min_err_sig << _8s<<err_sig.Mean() << (max_err_sig>tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig << "[0m" << _8s<<err_sig.Norm() << endl;
	cout << _4<< "Dis" << _8s<<min_err_dis << _8s<<err_dis.Mean() << (max_err_dis>tol_dis?"[1;31m":"[1;32m") << _8s<<max_err_dis << "[0m" << _8s<<err_dis.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_eps>tol_eps || max_err_sig>tol_sig || max_err_dis>tol_dis) return 1;
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
