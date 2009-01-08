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

/*       F1        F2        F3    F1= 0.0 F2=F3 = 1.0
          ^         ^         ^
          |         |         |

        2 @---------@---------@ 6
          |',       7         |
          |  ',        e[1]   |
          |    ',             |
          |      ',  4        |
        5 @        '@         @ 8
          |          ',       |
          |   e[0]     ',     |
          |              ',   |
          |         3      ', |
        0 @---------@---------@ 1
         /_\       /_\       /_\
         ///       o o       o o
 */


// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/tri6.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "mesh/mesh.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_8s;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	Mesh::Generic mesh(/*Is3D*/false);
	mesh.SetNVerts  (9);
	mesh.SetNElems  (2);
	mesh.SetVert    (0, true, 0.0, 0.0);
	mesh.SetVert    (1, true, 1.0, 0.0);
	mesh.SetVert    (2, true, 0.0, 1.0);
	mesh.SetVert    (3, true, 0.5, 0.0);
	mesh.SetVert    (4, true, 0.5, 0.5);
	mesh.SetVert    (5, true, 0.0, 0.5);
	mesh.SetVert    (6, true, 1.0, 1.0);
	mesh.SetVert    (7, true, 0.5, 1.0);
	mesh.SetVert    (8, true, 1.0, 0.5);
	mesh.SetElem    (0, -1, true, VTK_QUADRATIC_TRIANGLE);
	mesh.SetElem    (1, -1, true, VTK_QUADRATIC_TRIANGLE);
	mesh.SetElemCon (0, 0, 0);
	mesh.SetElemCon (0, 1, 1);
	mesh.SetElemCon (0, 2, 2);
	mesh.SetElemCon (0, 3, 3);
	mesh.SetElemCon (0, 4, 4);
	mesh.SetElemCon (0, 5, 5);
	mesh.SetElemCon (1, 0, 6);
	mesh.SetElemCon (1, 1, 2);
	mesh.SetElemCon (1, 2, 1);
	mesh.SetElemCon (1, 3, 7);
	mesh.SetElemCon (1, 4, 4);
	mesh.SetElemCon (1, 5, 8);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and solver
	FEM::Data   dat (2); // 2D
	FEM::Solver sol (dat,"ttri6");

	// Elements attributes
	FEM::EAtts_T eatts(1);
	String prms; prms.Printf("E=%f nu=%f", 10000.0, 0.25);
	eatts = T(-1, "Tri6", "PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 --------------------------------------
	dat.Nod(0)->Bry("ux",0.0)->Bry("uy",0.0);
	dat.Nod(1)->Bry("uy",0.0);
	dat.Nod(3)->Bry("uy",0.0);
	dat.Nod(2)->Bry("fy",0.0);
	dat.Nod(7)->Bry("fy",1.0);
	dat.Nod(6)->Bry("fy",1.0);
	sol.SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);

	///////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
    Array<double> err_eps;
    Array<double> err_sig;
    Array<double> err_dis;

	Vector<double> sx0 (6);
	Vector<double> sy0 (6);
	Vector<double> sxy0(6);
	Vector<double> sx1 (6);
	Vector<double> sy1 (6);
	Vector<double> sxy1(6);

	for (int i=0; i<6; i++)
	{
		sx0 (i) = dat.Ele(0)->Val(i, "Sx" );
		sy0 (i) = dat.Ele(0)->Val(i, "Sy" );
		sxy0(i) = dat.Ele(0)->Val(i, "Sxy");
		sx1 (i) = dat.Ele(1)->Val(i, "Sx" );
		sy1 (i) = dat.Ele(1)->Val(i, "Sy" );
		sxy1(i) = dat.Ele(1)->Val(i, "Sxy");
	}

	Vector<double> sx (9);
	Vector<double> sy (9);
	Vector<double> sxy(9);

	sx(0) = sx0(0);                 sy(0) = sy0(0);                 sxy(0) = sxy0(0);                 
	sx(1) = (sx0(1) + sx1(2))*0.5;  sy(1) = (sy0(1) + sy1(2))*0.5;  sxy(1) = (sxy0(1) + sxy1(2))*0.5;
	sx(2) = (sx0(2) + sx1(1))*0.5;  sy(2) = (sy0(2) + sy1(1))*0.5;  sxy(2) = (sxy0(2) + sxy1(1))*0.5;
	sx(3) = sx0(3);                 sy(3) = sy0(3);                 sxy(3) = sxy0(3);
	sx(4) = (sx0(4) + sx1(4))*0.5;  sy(4) = (sy0(4) + sy1(4))*0.5;  sxy(4) = (sxy0(4) + sxy1(4))*0.5;
	sx(5) = sx0(5);                 sy(5) = sy0(5);                 sxy(5) = sxy0(5);
	sx(6) = sx1(0);                 sy(6) = sy1(0);                 sxy(6) = sxy1(0);
	sx(7) = sx1(3);                 sy(7) = sy1(3);                 sxy(7) = sxy1(3);
	sx(8) = sx1(5);                 sy(8) = sy1(5);                 sxy(8) = sxy1(5);                 

	err_sig.Push( fabs(sx(0) - (-1.1241E-2)) );
	err_sig.Push( fabs(sx(1) - ( 3.3784E-1)) );
	err_sig.Push( fabs(sx(2) - (-8.1807E-2)) );
	err_sig.Push( fabs(sx(3) - ( 9.0341E-2)) );
	err_sig.Push( fabs(sx(4) - ( 1.2802E-1)) );
	err_sig.Push( fabs(sx(5) - (-8.5346E-3)) );
	err_sig.Push( fabs(sx(6) - (-5.0083E-1)) );
	err_sig.Push( fabs(sx(7) - (-3.2931E-1)) );
	err_sig.Push( fabs(sx(8) - (-8.5345E-3)) );

	err_sig.Push( fabs(sy(0) - ( -8.5054E-1 )) );
	err_sig.Push( fabs(sy(1) - (  4.8378)) );
	err_sig.Push( fabs(sy(2) - ( -1.1622)) );
	err_sig.Push( fabs(sy(3) - (  1.8464)) );
	err_sig.Push( fabs(sy(4) - (  1.8378)) );
	err_sig.Push( fabs(sy(5) - ( -6.8422E-1)) );
	err_sig.Push( fabs(sy(6) - (  5.4992)) );
	err_sig.Push( fabs(sy(7) - (  1.8464)) );
	err_sig.Push( fabs(sy(8) - (  5.3158)) );

	err_sig.Push( fabs(sxy(0) - (  1.1241E-2 )) );
	err_sig.Push( fabs(sxy(1) - (  8.1807E-2)) );
	err_sig.Push( fabs(sxy(2) - ( -3.3784E-1)) );
	err_sig.Push( fabs(sxy(3) - ( -8.5345E-3)) );
	err_sig.Push( fabs(sxy(4) - ( -1.2802E-1)) );
	err_sig.Push( fabs(sxy(5) - ( -7.3272E-2)) );
	err_sig.Push( fabs(sxy(6) - (  5.0083E-1)) );
	err_sig.Push( fabs(sxy(7) - ( -8.5346E-3)) );
	err_sig.Push( fabs(sxy(8) - (  3.4638E-1)) );


	// Error summary
	double tol_sig     = 1.0e-4;
	double min_err_sig = err_sig[err_sig.Min()];
	double max_err_sig = err_sig[err_sig.Max()];
	cout << _4<< ""    << _8s<<"Min"       << _8s<<"Mean"                                                        << _8s<<"Max"                  << _8s<<"Norm"         << endl;
	cout << _4<< "Sig" << _8s<<min_err_sig << _8s<<err_sig.Mean() << (max_err_sig>tol_sig?"[1;31m":"[1;32m") << _8s<<max_err_sig << "[0m" << _8s<<err_sig.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_sig>tol_sig) return 1;
	else return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
