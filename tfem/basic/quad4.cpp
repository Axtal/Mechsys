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
#include "fem/solver.h"
#include "fem/elems/quad4.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "util/fatal.h"
#include "util/numstreams.h"
#include "mesh/mesh.h"

using std::cout;
using std::endl;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	/*       F1                 F3    F1=F2=F3 = 1.0
	          ^                   ^
	          |                   |
	
	        3 @-------------------@ 2    Prof. Carlos Felippa
	          |                   |      IFEM.Ch23.pdf
	          |                   |      Pag. 23-7
	   L/2    |                   |
	          |                   |
	        0 @-------------------@ 1
	         /_\                 /_\
	         ///       L         o o
	 */

	// Input
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);
	else cout << "[1;32mYou may call this program as in:\t " << argv[0] << " LinSol\n  where LinSol:\n \tLA  => LAPACK_T  : DENSE\n \tUM  => UMFPACK_T : SPARSE\n \tSLU => SuperLU_T : SPARSE\n [0m[1;34m Now using LA (LAPACK)\n[0m" << endl;

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	Mesh::Generic mesh(/*Is3D*/false);
	mesh.SetNVerts  (4);
	mesh.SetVert    (0, true, 0.0, 0.0);
	mesh.SetVert    (1, true, 1.0, 0.0);
	mesh.SetVert    (2, true, 1.0, 0.5);
	mesh.SetVert    (3, true, 0.0, 0.5);
	mesh.SetNElems  (1);
	mesh.SetElem    (0, -1, true, VTK_QUAD);
	mesh.SetElemCon (0, 0, 0);
	mesh.SetElemCon (0, 1, 1);
	mesh.SetElemCon (0, 2, 2);
	mesh.SetElemCon (0, 3, 3);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////
	
	// Data and solver
	FEM::Data   dat (2); // 2D
	FEM::Solver sol (dat,"tquad4");

	// Elements attributes
	FEM::EAtts_T eatts(1);
	String prms; prms.Printf("E=%f nu=%f", 96.0, 1.0/3.0);
	eatts = T(-1, "Quad4", "PStress", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, true);

	// Set geometry: nodes and elements
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 --------------------------------------
	dat.Nod(0)->Bry("ux",0.0)->Bry("uy",0.0);
	dat.Nod(1)->Bry("uy",0.0);
	sol.SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);

	// Stiffness
	Array<size_t>          map;
	Array<bool>            pre;
	LinAlg::Matrix<double> Ke0;  Ke0.SetNS(Util::_6_3);
	dat.Ele(0)->CMatrix(0,Ke0);
	cout << "Ke0=\n" << Ke0 << endl;

	///////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Check
	LinAlg::Matrix<double> Ke_correct; Ke_correct.Resize(8,8);
	Ke_correct =  42.0,  18.0,  -6.0,   0.0, -21.0, -18.0, -15.0,   0.0,
	              18.0,  78.0,   0.0,  30.0, -18.0, -39.0,   0.0, -69.0,
	              -6.0,   0.0,  42.0, -18.0, -15.0,   0.0, -21.0,  18.0,
	               0.0,  30.0, -18.0,  78.0,   0.0, -69.0,  18.0, -39.0,
	             -21.0, -18.0, -15.0,   0.0,  42.0,  18.0,  -6.0,   0.0,
	             -18.0, -39.0,   0.0, -69.0,  18.0,  78.0,   0.0,  30.0,
	             -15.0,   0.0, -21.0,  18.0,  -6.0,   0.0,  42.0, -18.0,
	               0.0, -69.0,  18.0, -39.0,   0.0,  30.0, -18.0,  78.0;

	// Check
    double errors = 0.0;
	for (int i=0; i<8; ++i)
	for (int j=0; j<8; ++j)
		errors += fabs(Ke0(i,j)-Ke_correct(i,j));

	if (fabs(errors)>1.0e-3) cout << "[1;31mErrors(" << linsol << ") = " << Util::_8s << errors << "[0m\n" << endl;
	else                      cout << "[1;32mErrors(" << linsol << ") = " << Util::_8s << errors << "[0m\n" << endl;

	// Return error flag
	if (fabs(errors)>1.0e-3) return 1;
	else                     return 0;
}
MECHSYS_CATCH
