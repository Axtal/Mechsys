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
#include "fem/elems/quad4pstrain.h"
#include "fem/elems/quad8pstrain.h"
#include "fem/elems/beam.h"
#include "models/equilibs/linelastic.h"
#include "fem/solvers/forwardeuler.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "linalg/laexpr.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::PI;
using boost::make_tuple;

void Kirsch_stress(double p1, double p2, double r, double R, double th, Vector<double> & sig_k) // Calculate the Kirsch solution for a cylindrical hole
{
	sig_k.Resize(3);
	sig_k = (p1+p2)/2.0*(1.0-r*r/R/R) + (p1-p2)/2.0*(1.0 - 4.0*r*r/R/R + 3.0*pow(r,4.0)/pow(R,4.0))*cos(2.0*th), // sig_r
			(p1+p2)/2.0*(1.0+r*r/R/R) - (p1-p2)/2.0*(1.0 + 3.0*pow(r,4.0)/pow(R,4.0))*cos(2.0*th),               // sig_th
		   -(p1-p2)/2.0*(1.0+2.0*r*r/R/R - 3.0*pow(r,4.0)/pow(R,4.0))*sin(2.0*th);                               // sig_r_th
}

void Kirsch_disp(double p1, double p2, double r, double R, double th, double E, double nu, double ux_correct, double uy_correct) // Calculate the Kirsch solution for a cylindrical hole
{
	double G = E/2.0/(1+nu);
	ux_correct =  (p1+p2)/(4.0*G)*r*r/R + (p1-p2)/(4.0*G)*r*r/R*(4.0*(1.0-nu)-r*r/R/R)*cos(2.0*th), // disp_r
	uy_correct = -(p1-p2)/(4.0*G)*r*r/R   *(2.0*(1.0-2.0*nu)+r*r/R/R)*sin(2.0*th);                  // disp_theta
}
void RotateStress(Vector<double> & S, double th, Vector<double> & Q) // Calculate the stress components in a system rotated an angle equal to th
{
	double l1, l2, m1, m2;
	l1 =  cos(th); l2 = sin(th);
	m1 = -sin(th); m2 = cos(th);
	Matrix<double> T(3,3);
	T = l1*l1, m1*m1,     2*l1*m1,
		l2*l2, m2*m2,     2*l2*m2,
		l1*l2, m1*m2, l1*m2+m1*l2;
	Q = T*S;
}

int main(int argc, char **argv) try
{
	// Constants
	double E_soil   = 6000.0;   // Young [MPa]
	double nu_soil  = 0.2;      // Poisson
	double r        = 1.0;      // radius
	double L        = 10.0;     // length
	double H        = 10.0;     // height
	bool   is_o2    = false;    // use high order elements?

	// Input
	cout << "Input: " << argv[0] << "  is_o2\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	/*                |---------- L ----------|
	 *   -+- -+- -+-  o___________o_-30_______o
	 *    |   |   |   |                      ,|
	 *    |   |   |   |-40    [b1]         ,' |
	 *    |   |   d   o                  ,'   |
	 *    |   f   |   |    y       x   ,'     |
	 *    |   |   |   |     ',   ,'  ,o       |
	 *    |   |  -+-  o-,_    '+'  ,'         |
	 *    H   |           o-,    ,'           o -10
	 *    |  -+- . . . . . . 'o '      [b0]   |
	 *    |   |               .',             |
	 *    |   e               .  o  y^        |
	 *    |   |               .   \  |        |
	 *    |   |               .   |  +-->x    |
	 *   -+- -+-      +----r----->o-----o-----o
	 *                        .       -20
	 *                        .   |---- a ----|
	 *                |-- b --|------ c ------|
	 */

	// Geometry
	double a = L-r;
	double b = r*cos(2.*PI/8.);
	double c = L-b;
	double d = H-r;
	double e = r*sin(2.*PI/8.);
	double f = H-e;
	
	double K = 1;

	// Lower block -- coordinates
	Mesh::Block b0;
	b0.SetTag    (-1);
	b0.SetCoords (false, 8, // Is3D, NNodes
	               r,  L, L, b,    r+a/2.,    L, b+c/2., r*cos(PI/8.),
	              0., 0., H, e,        0., H/2., e+f/2., r*sin(PI/8.));
	b0.SetNx     (30*K,/*Ax*/2.0, /*Nonlinear*/true);
	b0.SetNy     (15*K);
	b0.SetETags  (4, 0, -10, -20, 0);

	// Upper block -- coordinates
	Mesh::Block b1;
	b1.SetTag    (-1);
	b1.SetCoords (false, 8,
	              b, L, 0., 0.,   b+c/2., L/2.,     0., r*cos(3.*PI/8.),
	              e, H,  H,  r,   e+f/2.,    H, r+d/2., r*sin(3.*PI/8.));
	b1.SetNx     (30*K,2.0, true);
	b1.SetNy     (15*K);
	b1.SetETags  (4, 0, -30,  0, -40);

	// Blocks
	Array<Mesh::Block*> blocks;  blocks.Resize(2);
	blocks[0] = &b0;
	blocks[1] = &b1;

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();               // Non-linear elements
	clock_t start = std::clock();          // Initial time
	size_t  ne    = mesh.Generate(blocks); // Discretize domain
	clock_t total = std::clock() - start;  // Time elapsed
	cout << "\nNumber of quads     = " << ne << endl;
	cout << "Time elapsed (mesh) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Geom g(2); // 2D

	// Nodes brys
	FEM::NBrys_T nbrys;

	// Edges brys
	FEM::EBrys_T ebrys;
	double p1 = -30.0;
	double p2 = -30.0;
	ebrys.Push (make_tuple(-10, "fx", p1));
	ebrys.Push (make_tuple(-20, "uy", 0.0));
	ebrys.Push (make_tuple(-30, "fy", p2));
	ebrys.Push (make_tuple(-40, "ux", 0.0));

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E_soil,nu_soil);
	FEM::EAtts_T eatts;
	if (is_o2) eatts.Push (make_tuple(-1, "Quad8PStrain", "LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0"));
	else       eatts.Push (make_tuple(-1, "Quad4PStrain", "LinElastic", prms.CStr(), "Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0"));

	// Set geometry: nodes, elements, attributes, and boundaries
	FEM::SetNodesElems (&mesh, &eatts, &g, 1.0e-5);
	FEM::SetBrys       (&mesh, &nbrys, &ebrys, NULL, &g);

	// Solve
	FEM::Solver * sol = FEM::AllocSolver("ForwardEuler");
	sol -> SetGeom(&g) -> SetLinSol("UM") -> SetNumDiv(1) -> SetDeltaTime(0.0);
	start = std::clock();
	sol -> Solve();
	total = std::clock() - start;
	double norm_resid = LinAlg::Norm(sol->Resid());
	cout << "Time elapsed (solution) = "<<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds\n";
	cout << "[1;35mNorm(Resid=DFext-DFint) = " << norm_resid << "[0m\n";
	cout << "[1;32mNumber of DOFs          = " << sol->nDOF() << "[0m\n";

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

    Vector<double> sig(3);
    Vector<double> sig_k(3);
    Array <double> err_sig_x;
    Array <double> err_sig_y;
    Array <double> err_sig_xy;
    Array <double> err_disp;

	// Stress 
	for (size_t i=0; i<g.NElems(); ++i)
	{
		for (size_t j=0; j<g.Ele(i)->NNodes(); ++j)
		{
			// Analytical
			double x  = g.Ele(i)->Nod(j)->X();
			double y  = g.Ele(i)->Nod(j)->Y();
			double th = atan(y/x);
			double R  = sqrt(x*x + y*y);
			Kirsch_stress(p1, p2, r, R, th, sig_k); // Calculate the Kirsch solution for a cylindrical hole

			// Analysis
			sig = g.Ele(i)->Val(j,"Sx"), g.Ele(i)->Val(j,"Sy"), g.Ele(i)->Val(j,"Sxy");
			Vector<double> sig_polar(3); // Vector to transform to polar coordinates
			sig_polar = sig(0)*pow(cos(th),2.0) + sig(1)*pow(sin(th),2.0) + 2.0*sig(2)*sin(th)*cos(th),
			            sig(0)*pow(sin(th),2.0) + sig(1)*pow(cos(th),2.0) - 2.0*sig(2)*sin(th)*cos(th),
			            (sig(1)-sig(0))*sin(th)*cos(th) + sig(2)*(pow(cos(th),2.0) - pow(sin(th),2.0));

			err_sig_x .Push ( fabs(sig_k(0) - sig_polar(0)) );
			err_sig_y .Push ( fabs(sig_k(1) - sig_polar(1)) );
			err_sig_xy.Push ( fabs(sig_k(2) - sig_polar(2)) );
		}
	}

	// Displacements
	for (size_t i=0; i<g.NNodes(); ++i)
	{
		// Analytical
		double x  = g.Nod(i)->X();
		double y  = g.Nod(i)->Y();
		double th = atan(y/x);
		double R  = sqrt(x*x + y*y);
		double ux_correct;
		double uy_correct;
		Kirsch_disp(p1, p2, r, R, th, E_soil, nu_soil, ux_correct, uy_correct);// Calculate the Kirsch solution for a cylindrical hole

		// Analysis
		err_disp.Push ( fabs(g.Nod(i)->Val("ux") - ux_correct) );
		err_disp.Push ( fabs(g.Nod(i)->Val("uy") - uy_correct) );
	}

	// Output: VTU
	Output o; o.VTU (&g, "texam1.vtu");
	cout << "[1;34mFile <texam1.vtu> saved.[0m\n\n";

	// Error summary
	double tol_sig        = 3.0e0;
	double tol_disp       = 1.0e-1;
	double min_err_sig_x  = err_sig_x [err_sig_x .Min()];
	double min_err_sig_y  = err_sig_y [err_sig_y .Min()];
	double min_err_sig_xy = err_sig_xy[err_sig_xy.Min()];
	double min_err_disp   = err_disp  [err_disp.Min()];
	double max_err_sig_x  = err_sig_x [err_sig_x .Max()];
	double max_err_sig_y  = err_sig_y [err_sig_y .Max()];
	double max_err_sig_xy = err_sig_xy[err_sig_xy.Max()];
	double max_err_disp   = err_disp  [err_disp.Max()];
	cout << _4 << ""      << _8s<<"Min"       << _8s << "Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4 << "SigX"  << _8s<<min_err_sig_x  << _8s<<err_sig_x .Mean() << (max_err_sig_x >tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig_x  << "[0m" << _8s<<err_sig_x.Norm() << endl;
	cout << _4 << "SigY"  << _8s<<min_err_sig_y  << _8s<<err_sig_y .Mean() << (max_err_sig_y >tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig_y  << "[0m" << _8s<<err_sig_y.Norm() << endl;
	cout << _4 << "SigXY" << _8s<<min_err_sig_xy << _8s<<err_sig_xy.Mean() << (max_err_sig_xy>tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig_xy << "[0m" << _8s<<err_sig_xy.Norm() << endl;
	cout << _4 << "Disp"  << _8s<<min_err_disp   << _8s<<err_disp  .Mean() << (max_err_disp  >tol_disp?"[1;31m":"[1;32m") << _8s<<max_err_disp   << "[0m" << _8s<<err_disp  .Norm() << endl;
	cout << endl;

	// Return error flag
	if ( max_err_sig_x>tol_sig || max_err_sig_y>tol_sig || max_err_sig_xy>tol_sig || max_err_disp>tol_disp ) return 1;
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
