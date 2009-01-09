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
#include "fem/elems/quad4.h"
#include "fem/elems/quad8.h"
#include "fem/equilibelem.h"
#include "fem/elems/beam.h"
#include "models/equilibs/linelastic.h"
#include "models/equilibs/beamelastic.h"
#include "fem/solver.h"
#include "fem/output.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;
using Util::PI;

#define T boost::make_tuple

void Einsten_Schwartz(
		double p0,   // initial stress load
		double K,    // ratio beetween horizontal and vertical initial stress load
		double R,    // radius: polar coordinate
		double th,   // angle:  polar coordinate
		double r,    // radius of the hole
		double E,    // rock Young modulus
		double nu,   // rock Poisson ratio
		double E_s,  // liner Young modulus
		double nu_s, // liner Poisson ratio
		double I,    // liner moment of inertia
		double A,    // cross sectional area of the liner for a 1m long section
		double & M,  // bending moment in the liner
		double & N,  // axial force in the liner
		double & uT, // tangential displacement
		double & uR  // radial displacement
		)
{
	double Nu   = 1.0-nu;
	double Nu2  = 1.0-nu*nu;
	double Nus2 = 1.0-nu_s*nu_s;
	double C    = E*R*Nus2 / (E_s*A*Nu2);
	double F    = E*pow(R,3)*Nus2 / ( E_s*I*(1.0-nu*nu) );
	double beta = ((6.0+F)*C*Nu + 2*F*nu) / (3*F + 3*C + 2*C*F*Nu);
	double b2   = C*Nu * 0.5 / (C*Nu + 4*nu - 6*beta - 3*beta*C*Nu);
	double a2   = beta*b2;
	double a0   = C*F*Nu / (C+F+C*F*Nu);
	double G    = E/2.0/(1+nu);

	M = p0*r*r*0.25*(1.0-K)*(1.0 - 2.0*a2 + 2.0*b2)*cos(2.0*th);
	N = p0*r*(0.5*(1.0+K)*(1.0 - a0) + 0.5*(1.0-K)*(1.0+2*a2)*cos(2.0*th) );
	uT = 0.5*p0*r/G*(1.0-K)*(a2 + (1.0 - 2*nu)*b2)*sin(2.0*th);
	uR = 0.5*p0*r/G*(-0.5*(1.0+K)*a0 + 0.5*(1.0-K)*( 4.0*(1.0-nu)*b2 - 2.0*a2 )*cos(2*th) );
}

int main(int argc, char **argv) try
{
	// Constants
	double E_soil   = 6000.0;   // Young [MPa]
	double nu_soil  = 0.2;      // Poisson
	double E_beam   = 20000.0;  // Young [MPa]
	double Izz_beam = 0.01042;  // Inertia m^4
	double A_beam   = 0.5;      // Area m^2
	double r        = 1.25;     // radius
	double L        = 50.0;     // length
	double H        = 50.0;     // height
	bool   is_o2    = false;    // use high order elements?
	int    ny       = 15;       // num divisions along Y
	double Ax       = 2.0;      // rate of increase of X divisions
	double NonLinX  = false;    // nonlinear divisions along X?

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ny(nx=2*ny)\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ny    =  atoi(argv[2]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	/*                |---------- L ----------|
	 *   -+- -+- -+-  o___________o_-30_______o
	 *    |   |   |   |                      ,|
	 *    |   |   |   |-40    [b1]         ,' |
	 *    |   |   d   o                  ,'   |
	 *    |   f   |   |    y       x   ,'     |
	 *    |   |   |   |     ',   ,'  ,o       |
	 *    |   |  -+-  o-,_    '+'  ,'         |
	 *    H   |       -55 o-,    ,'           o -20
	 *    |  -+- . . . . . . 'o '      [b0]   |
	 *    |   |               .',             |
	 *    |   e               .  o  y^        |
	 *    |   |               .-55\  |        |
	 *    |   |               .   |  +-->x    |
	 *   -+- -+-      +----r----->o-----o-----o
	 *                        .       -10
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

	// Blocks
	Array<Mesh::Block> bks(2);

	// Block # 0 --------------------------------
    Mesh::Verts_T ve0(8);
    Mesh::Edges_T ed0(8);
    Mesh::ETags_T et0(6);
	ve0 = T(0 , r            , 0.           , 0.), 
	      T(1 , L            , 0.           , 0.), 
	      T(2 , L            , H            , 0.), 
	      T(3 , b            , e            , 0.), 
	      T(4 , r+a/2.       , 0.           , 0.), 
	      T(5 , L            , H/2.         , 0.), 
	      T(6 , b+c/2.       , e+f/2.       , 0.), 
	      T(7 , r*cos(PI/8.) , r*sin(PI/8.) , 0.); 
    ed0 = T(0,4), T(4,1), T(1,5), T(5,2), T(2,6), T(6,3), T(3,7), T(7,0);
    et0 = T(0,4,-10), T(4,1,-10), T(1,5,-20), T(5,2,-20), T(0,7,-55), T(3,7,-55);
    bks[0].Set   (-1, ve0, ed0, &et0, NULL, /*orig*/0, /*xplus*/4, /*yplus*/7);
	bks[0].SetNx (2*ny, Ax, NonLinX);
	bks[0].SetNy (ny);

	// Block # 1 --------------------------------
    Mesh::Verts_T ve1(8);
    Mesh::Edges_T ed1(8);
    Mesh::ETags_T et1(6);
	ve1 = T(0 , b               , e               , 0.), 
	      T(1 , L               , H               , 0.), 
	      T(2 , 0.              , H               , 0.), 
	      T(3 , 0.              , r               , 0.), 
	      T(4 , b+c/2.          , e+f/2.          , 0.), 
	      T(5 , L/2.            , H               , 0.), 
	      T(6 , 0.              , r+d/2.          , 0.), 
	      T(7 , r*cos(3.*PI/8.) , r*sin(3.*PI/8.) , 0.); 
    ed1 = T(0,4), T(4,1), T(1,5), T(5,2), T(2,6), T(6,3), T(3,7), T(7,0);
    et1 = T(1,5,-30), T(5,2,-30), T(2,6,-40), T(6,3,-40), T(0,7,-55), T(3,7,-55);
    bks[1].Set   (-1, ve1, ed1, &et1, NULL, /*orig*/0, /*xplus*/4, /*yplus*/7);
	bks[1].SetNx (2*ny, Ax, NonLinX);
	bks[1].SetNy (ny);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (bks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and Solver
	FEM::Data   dat (2);
	FEM::Solver sol (dat, "texam5");

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E_soil,nu_soil);
	String geom; geom = (is_o2 ? "Quad8" : "Quad4");
	FEM::EAtts_T eatts(1);
	eatts = T(-1, geom.CStr(), "PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", true);

	// Beam attributes
	String beamprms;
	beamprms.Printf("E=%f A=%f Izz=%f",E_beam,A_beam,Izz_beam);
	eatts.Push (T(-55, "Lin2", "Beam", "BeamElastic", beamprms.CStr(), "ZERO", "gam=20 cq=1", true));

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Nodes brys
	FEM::NBrys_T nbrys;
	nbrys.Push (T(r,0.0,0.0,"wz",0.0));
	nbrys.Push (T(0.0,r,0.0,"wz",0.0));

	// Edges brys
	FEM::EBrys_T ebrys;
	double p0 = -15.0;
	double  K =  2.0;
	ebrys.Push (T(-10, "uy",  0.0));
	ebrys.Push (T(-20, "fx", p0*K));
	ebrys.Push (T(-30, "fy",   p0));
	ebrys.Push (T(-40, "ux",  0.0));
	ebrys.Push (T(-55, "Qb",  -1));

	// Stage # 1
	dat.SetBrys (&mesh, &nbrys, &ebrys, NULL);
	sol.SolveWithInfo();

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

    Vector<double> sig(3);
    Vector<double> sig_k(3);
    Array <double> err_M;
    Array <double> err_N;
    Array <double> err_uR;
    Array <double> err_uT;

	// Stress 
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		String ele_name;
		if (dat.Ele(i)->Tag()==-55) 
			for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j)
			{
				// Analytical
				double x  = dat.Ele(i)->Nod(j)->X();
				double y  = dat.Ele(i)->Nod(j)->Y();
				double th = atan(y/x);
				double R  = sqrt(x*x + y*y);
				double M_correct; // Bending momentun
				double N_correct; // Axial force
				double uT_correct, uR_correct;
				double nu_beam = 0.;
				Einsten_Schwartz( p0, K, R, th, r, E_soil, nu_soil, E_beam, nu_beam, Izz_beam, A_beam, M_correct, N_correct, uT_correct, uR_correct );

				// Analysis
				double M;
				double N;
				dat.Ele(i)->CalcDeps();
				M = dat.Ele(i)->Val(j,"M");
				N = dat.Ele(i)->Val(j,"N");

				err_M .Push ( fabs(M_correct - M) );
				err_N .Push ( fabs(N_correct - N) );

				int node_id = dat.Ele(i)->Nod(j)->ID();

				double ux = dat.Nod(node_id)->Val("ux");
				double uy = dat.Nod(node_id)->Val("uy");
				double uR = -(ux*cos(th)+uy*sin(th));
				double uT = -(uy*cos(th)-ux*sin(th));

				err_uR.Push (fabs(uR - uR_correct)) ;
				err_uT.Push (fabs(uT - uT_correct)) ;

				//cout << "node_id  = " << node_id << endl;

				if (node_id==4949)
				{
					cout << "uR = " << uR << "  uRc = " << uR_correct << endl;
					cout << "uT = " << uT << "  uTc = " << uT_correct << endl;
				}

			}
	}

	//// Displacements
	//for (size_t i=0; i<dat.NNodes(); ++i)
	//{
	//	// Analytical
	//	double x  = dat.Nod(i)->X();
	//	double y  = dat.Nod(i)->Y();
	//	double th = atan(y/x);
	//	double R  = sqrt(x*x + y*y);
	//	double uR_correct;
	//	double uT_correct;
	//	double nu_beam = 0.25;
	//	double M, N;

	//	Einsten_Schwartz( p0, K, R, th, r, E_soil, nu_soil, E_beam, nu_beam, Izz_beam, A_beam, M, N, uT_correct, uR_correct);

	//	// Analysis
	//	double ux = dat.Nod(i)->Val("ux");
	//	double uy = dat.Nod(i)->Val("uy");
	//	double uR = ux*cos(th)+uy*sin(th);
	//	double uT = uy*cos(th)-ux*sin(th);

	//	//if (i==50

	//	err_uR.Push (fabs(uR - uR_correct)) ;
	//	err_uT.Push (fabs(uT - uT_correct)) ;
	//}

	// Error summary
	double tol_uR       = 1.0e-1;
	double tol_uT       = 1.0e-1;
	double tol_M        = 1.0e-1;
	double tol_N        = 1.0e-1;
	double min_err_uR   = err_uR[err_uR.Min()];
	double max_err_uR   = err_uR[err_uR.Max()];
	double min_err_uT   = err_uT[err_uT.Min()];
	double max_err_uT   = err_uT[err_uT.Max()];
	double min_err_M    = err_M [err_M .Min()];
	double max_err_M    = err_M [err_M .Max()];
	double min_err_N    = err_N [err_N .Min()];
	double max_err_N    = err_N [err_N .Max()];

	cout << "\n";
	cout << _4 << ""   << _8s <<"Min"      << _8s << "Mean"                                                     << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4 << "uR" << _8s <<min_err_uR << _8s<<err_uR.Mean() << (max_err_uR  >tol_uR?"[1;31m":"[1;32m") << _8s<<max_err_uR  << "[0m" << _8s<<err_uR .Norm() << endl;
	cout << _4 << "uT" << _8s <<min_err_uT << _8s<<err_uT.Mean() << (max_err_uT  >tol_uT?"[1;31m":"[1;32m") << _8s<<max_err_uT  << "[0m" << _8s<<err_uT .Norm() << endl;
	cout << _4 << "M"  << _8s <<min_err_M  << _8s<<err_M .Mean() << (max_err_M   >tol_M?"[1;31m":"[1;32m")  << _8s<<max_err_M   << "[0m" << _8s<<err_M  .Norm() << endl;
	cout << _4 << "N"  << _8s <<min_err_N  << _8s<<err_N .Mean() << (max_err_N   >tol_N?"[1;31m":"[1;32m")  << _8s<<max_err_N   << "[0m" << _8s<<err_N  .Norm() << endl;
	cout << endl;

	return 1;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
