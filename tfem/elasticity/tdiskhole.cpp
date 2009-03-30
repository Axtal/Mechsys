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
#include "models/equilibs/linelastic.h"
#include "fem/solver.h"
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

#define T boost::make_tuple

// Calculate Kirsch's solution for a cylindrical hole (stresses)
inline void KirschStress(double ph, double pv, double r, double R, double TT, double & SigR, double & SigT, double & SigRT)
{
	double pm = (ph+pv)/2.0;
	double pd = (ph-pv)/2.0;
	double c1 = r*r/(R*R);
	SigR  =  pm*(1.0-c1) + pd*(1.0-4.0*c1+3.0*c1*c1)*cos(2.0*TT);
	SigT  =  pm*(1.0+c1) - pd*(1.0+3.0*c1*c1)*cos(2.0*TT);
	SigRT = -pd*(1.0+2.0*c1-3.0*c1*c1)*sin(2.0*TT);
}

// Calculate Kirsch's solution for a cylindrical hole (displacements)
void KirschDisp(double ph, double pv, double r, double R, double TT, double E, double nu, double & uR, double & uT)
{
	double G  = E/(2.0*(1.0+nu)); // Shear modulus
	double c1 = r*r/R;
	double qm = (ph+pv)/(4.0*G);
	double qd = (ph-pv)/(4.0*G);
	uR =  qm*c1 + qd*c1*(4.0*(1.0-nu)-c1/R)*cos(2.0*TT);
	uT = -qd*c1*(2.0*(1.0-2.0*nu)+c1/R)*sin(2.0*TT);
}

int main(int argc, char **argv) try
{
	// Constants
	double ph      = -30.0;    // horizontal pressure
	double pv      = -30.0;    // vertical pressure
	double E_soil  = 6000.0;   // Young [MPa]
	double nu_soil = 0.2;      // Poisson
	double r       = 1.0;      // radius
	double L       = 10.0;     // length
	double H       = 10.0;     // height
	bool   is_o2   = false;    // use high order elements?
	int    ny      = 15;       // num divisions along Y
	double Ax      = 2.0;      // rate of increase of X divisions
	double NonLinX = false;    // nonlinear divisions along X?

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ny(nx=2*ny)\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ny    =  atoi(argv[2]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	/*                     | | | pv | | |
	 *                     V V V    V V V
	 *
	 *                |---------- L ----------|
	 *   -+- -+- -+-  o___________o_-30_______o
	 *    |   |   |   |                      ,|
	 *    |   |   |   |-40    [b1]         ,' |
	 *    |   |   d   o                  ,'   |
	 *    |   f   |   |    y       x   ,'     |      <--
	 *    |   |   |   |     ',   ,'  ,o       |      <--
	 *    |   |  -+-  o-,_    '+'  ,'         |      <--
	 *    H   |           o-,    ,'           o -20  <-- ph
	 *    |  -+- . . . . . . 'o '      [b0]   |      <--
	 *    |   |               .',             |      <--
	 *    |   e               .  o  y^        |      <--
	 *    |   |               .   \  |        |
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
    Mesh::ETags_T et0(4);
	ve0 = T(0 , r            , 0.           , 0.), 
	      T(1 , L            , 0.           , 0.), 
	      T(2 , L            , H            , 0.), 
	      T(3 , b            , e            , 0.), 
	      T(4 , r+a/2.       , 0.           , 0.), 
	      T(5 , L            , H/2.         , 0.), 
	      T(6 , b+c/2.       , e+f/2.       , 0.), 
	      T(7 , r*cos(PI/8.) , r*sin(PI/8.) , 0.); 
    ed0 = T(0,4), T(4,1), T(1,5), T(5,2), T(2,6), T(6,3), T(3,7), T(7,0);
    et0 = T(0,4,-10), T(4,1,-10), T(1,5,-20), T(5,2,-20);
    bks[0].Set   (-1, ve0, ed0, &et0, NULL, /*orig*/0, /*xplus*/4, /*yplus*/7);
	bks[0].SetNx (2*ny, Ax, NonLinX);
	bks[0].SetNy (ny);

	// Block # 1 --------------------------------
    Mesh::Verts_T ve1(8);
    Mesh::Edges_T ed1(8);
    Mesh::ETags_T et1(4);
	ve1 = T(0 , b               , e               , 0.), 
	      T(1 , L               , H               , 0.), 
	      T(2 , 0.              , H               , 0.), 
	      T(3 , 0.              , r               , 0.), 
	      T(4 , b+c/2.          , e+f/2.          , 0.), 
	      T(5 , L/2.            , H               , 0.), 
	      T(6 , 0.              , r+d/2.          , 0.), 
	      T(7 , r*cos(3.*PI/8.) , r*sin(3.*PI/8.) , 0.); 
    ed1 = T(0,4), T(4,1), T(1,5), T(5,2), T(2,6), T(6,3), T(3,7), T(7,0);
    et1 = T(1,5,-30), T(5,2,-30), T(2,6,-40), T(6,3,-40);
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
	FEM::Solver sol (dat, "texam1");

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E_soil,nu_soil);
	String geom; geom = (is_o2 ? "Quad8" : "Quad4");
	FEM::EAtts_T eatts(1);
	eatts = T(-1, geom.CStr(), "PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1
	FEM::EBrys_T ebrys;
	ebrys.Push  (T(-10, "uy", 0.));
	ebrys.Push  (T(-20, "fx", ph));
	ebrys.Push  (T(-30, "fy", pv));
	ebrys.Push  (T(-40, "ux", 0.));
	dat.SetBrys (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo();

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Stress
	Array <double> err_sR;
	Array <double> err_sT;
	Array <double> err_sRT;
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		for (size_t j=0; j<dat.Ele(i)->NNodes(); ++j)
		{
			// Analytical
			double x = dat.Ele(i)->Nod(j)->X();
			double y = dat.Ele(i)->Nod(j)->Y();
			double t = atan(y/x);
			double R = sqrt(x*x+y*y);
			double sigRc, sigTc, sigRTc; // correct stress components
			KirschStress (ph,pv,r,R,t, sigRc,sigTc,sigRTc);

			// Numerical
			double c     = x/R;
			double s     = y/R;
			double cc    = c*c;
			double ss    = s*s;
			double sc    = s*c;
			double Sx    = dat.Ele(i)->Val(j,"Sx");
			double Sy    = dat.Ele(i)->Val(j,"Sy");
			double Sxy   = dat.Ele(i)->Val(j,"Sxy");
			double sigR  = Sx*cc + Sy*ss + 2.0*Sxy*sc;
			double sigT  = Sx*ss + Sy*cc - 2.0*Sxy*sc;
			double sigRT = (Sy-Sx)*sc + (cc-ss)*Sxy;

			// Error
			err_sR .Push (fabs(sigR  - sigRc ));
			err_sT .Push (fabs(sigT  - sigTc ));
			err_sRT.Push (fabs(sigRT - sigRTc));
		}
	}

	// Displacements
	Array <double> err_uR;
	Array <double> err_uT;
	for (size_t i=0; i<dat.NNodes(); ++i)
	{
		// Analytical
		double x  = dat.Nod(i)->X();
		double y  = dat.Nod(i)->Y();
		double t  = atan(y/x);
		double R  = sqrt(x*x+y*y);
		double uRc, uTc; // correct displacement components
		KirschDisp (ph,pv,r,R,t, E_soil,nu_soil, uRc,uTc);

		// Numerical
		double c  = x/R;
		double s  = y/R;
		double ux = dat.Nod(i)->Val("ux");
		double uy = dat.Nod(i)->Val("uy");
		double uR = ux*c + uy*s;
		double uT = uy*c - ux*s;

		// Error
		err_uR.Push (fabs(uR-uRc));
		err_uT.Push (fabs(uT-uTc));
	}

	// Error summary
	double tol_sR      = 3.0e0;
	double tol_sT      = 3.0e0;
	double tol_sRT     = 3.0e0;
	double tol_uR      = 1.0e-1;
	double tol_uT      = 1.0e-3;
	double min_err_sR  = err_sR [err_sR .Min()];   double max_err_sR  = err_sR [err_sR .Max()];
	double min_err_sT  = err_sT [err_sT .Min()];   double max_err_sT  = err_sT [err_sT .Max()];
	double min_err_sRT = err_sRT[err_sRT.Min()];   double max_err_sRT = err_sRT[err_sRT.Max()];
	double min_err_uR  = err_uR [err_uR .Min()];   double max_err_uR  = err_uR [err_uR .Max()];
	double min_err_uT  = err_uT [err_uT .Min()];   double max_err_uT  = err_uT [err_uT .Max()];
	cout << endl;
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "sR"  << _8s<<min_err_sR  << _8s<<err_sR .Mean() << (max_err_sR >tol_sR ?"[1;31m":"[1;32m") << _8s<<max_err_sR  << "[0m" << _8s<<err_sR.Norm()  << endl;
	cout << _4<< "sT"  << _8s<<min_err_sT  << _8s<<err_sT .Mean() << (max_err_sT >tol_sT ?"[1;31m":"[1;32m") << _8s<<max_err_sT  << "[0m" << _8s<<err_sT.Norm()  << endl;
	cout << _4<< "sRT" << _8s<<min_err_sRT << _8s<<err_sRT.Mean() << (max_err_sRT>tol_sRT?"[1;31m":"[1;32m") << _8s<<max_err_sRT << "[0m" << _8s<<err_sRT.Norm() << endl;
	cout << _4<< "uR"  << _8s<<min_err_uR  << _8s<<err_uR .Mean() << (max_err_uR >tol_uR ?"[1;31m":"[1;32m") << _8s<<max_err_uR  << "[0m" << _8s<<err_uR.Norm()  << endl;
	cout << _4<< "uT"  << _8s<<min_err_uT  << _8s<<err_uT .Mean() << (max_err_uT >tol_uT ?"[1;31m":"[1;32m") << _8s<<max_err_uT  << "[0m" << _8s<<err_uT.Norm()  << endl;

	// Return error flag
	if (max_err_sR>tol_sR || max_err_sT>tol_sT || max_err_sRT>tol_sRT || max_err_uR>tol_uR || max_err_uT>tol_uT) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }