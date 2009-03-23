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

/*            12
     o|\ +-----------+ /|o
     o|/ |  embank 2 | \|o
         |-----------|
         |  embank 1 |  12
         |-----------|
         |           |
         +-----------+
        /_\         /_\
        o o         o o
*/

// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/quad4.h"
#include "fem/elems/quad8.h"
#include "fem/equilibelem.h"
#include "models/equilibs/linelastic.h"
#include "util/exception.h"
#include "linalg/matrix.h"
#include "mesh/structured.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_8s;

#define T boost::make_tuple

int main(int argc, char **argv) try
{
	// Constants
	double E     = 5000.0; // Young
	double nu    = 0.3;    // Poisson
	int    ndivx = 48;     // number of divisions along x and y (for each block)
	int    ndivy = 16;     // number of divisions along x and y (for each block)
	bool   is_o2 = false;  // use high order elements?

	// Input
	cout << "Input: " << argv[0] << "  is_o2  ndivx  ndivy\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndivx =  atof(argv[2]);
	if (argc>=4) ndivy =  atof(argv[3]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Blocks
	Array<Mesh::Block> bks(3);

	// Block # 0 --------------------------------
    Mesh::Verts_T ve0(4);
    Mesh::Edges_T ed0(4);
    Mesh::ETags_T et0(3);
	ve0 = T(0,0.0,0.0,0.), T(1,12.0,0.0,0.), T(2,12.0,4.0,0.), T(3,0.0,4.0,0.);
	ed0 = T(0,1), T(1,2), T(2,3), T(3,0);
	et0 = T(3,0,-10), T(1,2,-10), T(0,1,-11);
    bks[0].Set (-1, ve0, ed0, &et0, NULL, /*orig*/0, /*xplus*/1, /*yplus*/3);
	bks[0].SetNx (ndivx);
	bks[0].SetNy (ndivy);

	// Block # 1 --------------------------------
    Mesh::Verts_T ve1(4);
    Mesh::Edges_T ed1(4);
    Mesh::ETags_T et1(2);
	ve1 = T(0,0.0,4.0,0.), T(1,12.0,4.0,0.), T(2,12.0,8.0,0.), T(3,0.0,8.0,0.);
	ed1 = T(0,1), T(1,2), T(2,3), T(3,0);
	et1 = T(3,0,-10), T(1,2,-10);
    bks[1].Set (-2, ve1, ed1, &et1, NULL, /*orig*/0, /*xplus*/1, /*yplus*/3);
	bks[1].SetNx (ndivx);
	bks[1].SetNy (ndivy);

	// Block # 2 --------------------------------
    Mesh::Verts_T ve2(4);
    Mesh::Edges_T ed2(4);
    Mesh::ETags_T et2(3);
	ve2 = T(0,0.0,8.0,0.), T(1,12.0,8.0,0.), T(2,12.0,12.0,0.), T(3,0.0,12.0,0.);
	ed2 = T(0,1), T(1,2), T(2,3), T(3,0);
	et2 = T(3,0,-10), T(1,2,-10), T(3,2,-12);
    bks[2].Set (-3, ve2, ed2, &et2, NULL, /*orig*/0, /*xplus*/1, /*yplus*/3);
	bks[2].SetNx (ndivx);
	bks[2].SetNy (ndivy);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (bks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Data dat(2); // 2D

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f",E,nu);
	String geom; geom = (is_o2 ? "Quad8" : "Quad4");
	FEM::EAtts_T eatts(3);
	eatts = T(-1, geom.CStr(), "PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, true ),
	        T(-2, geom.CStr(), "PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, false),
	        T(-3, geom.CStr(), "PStrain", "LinElastic", prms.CStr(), "ZERO", "gam=20", FNULL, false);

	// Set geometry: nodes, elements, attributes, and boundaries
	dat.SetNodesElems (&mesh, &eatts);

	// Solver
	FEM::Solver sol(dat, "tembank01");

	// Stage # -1 --------------------------------------------------------------
	FEM::EBrys_T ebrys;
	ebrys.Push        (T(-10, "ux", 0.0));
	ebrys.Push        (T(-11, "uy", 0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	dat.AddVolForces  ();
	sol.SolveWithInfo (/*NDiv*/1, /*DTime*/1.0, /*iStage*/-1, "  Initial stress state due to self weight (zero displacements)\n", /*ClearDisp*/true);

	// Stage # 0 ---------------------------------------------------------------
	dat.Activate      (/*Tag*/-2);
    ebrys.Resize      (0);
	ebrys.Push        (T(-10, "ux", 0.0));
	ebrys.Push        (T(-11, "uy", 0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo (1, 2.0, 0, "  Construction of first layer\n");

	// Stage # 1 ---------------------------------------------------------------
	dat.Activate      (/*Tag*/-3);
    ebrys.Resize      (0);
	ebrys.Push        (T(-10, "ux", 0.0));
	ebrys.Push        (T(-11, "uy", 0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo (1, 3.0, 0, "  Construction of second layer\n");

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	return 0;
}
catch (Exception  * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout << "Fatal: "<<m<<std::endl;  exit(1); }
catch (...)            { std::cout << "Some exception (...) ocurred\n"; }
