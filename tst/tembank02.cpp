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
#include "fem/biot.h"
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
	double E     = 5000.0; // Young
	double nu    = 0.3;    // Poisson
	double gw    = 10.0;   // gamma w
	double k     = 1.0e-2; // permeability
	int    ndivx = 10;     // number of divisions along x and y (for each block)
	int    ndivy =  4;     // number of divisions along x and y (for each block)
	//int    ndivx = 48;   // number of divisions along x and y (for each block)
	//int    ndivy = 16;   // number of divisions along x and y (for each block)
	bool   is_o2 = false;  // use high order elements?

	// Input
	//cout << "Input: " << argv[0] << "  is_o2  ndivx  ndivy\n";
	if (argc>=2) is_o2 = (atoi(argv[1])>0 ? true : false);
	if (argc>=3) ndivx =  atof(argv[2]);
	if (argc>=4) ndivy =  atof(argv[3]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	// Block # 0
	Mesh::Block b0;
	b0.SetTag    (-1); // tag to be replicated to all generated elements inside this block
	b0.SetCoords (false, 4,               // Is3D, NNodes
	              0.0, 12.0, 12.0, 0.0,   // x coordinates
	              0.0,  0.0,  4.0, 4.0);  // y coordinates
	b0.SetNx     (ndivx);                 // x weights and num of divisions along x
	b0.SetNy     (ndivy);                 // y weights and num of divisions along y
	b0.SetETags  (4, -10, -10, -11, -12); // edge tags

	// Block # 1
	Mesh::Block b1;
	b1.SetTag    (-2); // tag to be replicated to all generated elements inside this block
	b1.SetCoords (false, 4,               // Is3D, NNodes
	              0.0, 12.0, 12.0, 0.0,   // x coordinates
	              4.0,  4.0,  8.0, 8.0);  // y coordinates
	b1.SetNx     (ndivx);                 // x weights and num of divisions along x
	b1.SetNy     (ndivy);                 // y weights and num of divisions along y
	b1.SetETags  (4, -10, -10,   0, -13); // edge tags

	// Block # 2
	Mesh::Block b2;
	b2.SetTag    (-3); // tag to be replicated to all generated elements inside this block
	b2.SetCoords (false, 4,               // Is3D, NNodes
	              0.0, 12.0, 12.0,  0.0,  // x coordinates
	              8.0,  8.0, 12.0, 12.0); // y coordinates
	b2.SetNx     (ndivx);                 // x weights and num of divisions along x
	b2.SetNy     (ndivy);                 // y weights and num of divisions along y
	b2.SetETags  (4, -10, -10,   0, -14); // edge tags

	// Blocks
	Array<Mesh::Block*> blocks;
	blocks.Push (&b0);
	blocks.Push (&b1);
	blocks.Push (&b2);

	// Generate
	Mesh::Structured mesh(/*Is3D*/false);
	if (is_o2) mesh.SetO2();
	mesh.SetBlocks (blocks);
	mesh.Generate  (true);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Geometry
	FEM::Data dat(2); // 2D

	// Elements attributes
	String prms; prms.Printf("E=%f nu=%f gw=%f k=%f",E,nu,gw,k);
	FEM::EAtts_T eatts;
	if (is_o2)
	{
		eatts.Push (make_tuple(-1, "Quad8", "Biot", "LinElastic", prms.CStr(), "ZERO", "gam=20", true ));
		eatts.Push (make_tuple(-2, "Quad8", "Biot", "LinElastic", prms.CStr(), "ZERO", "gam=20", false));
		eatts.Push (make_tuple(-3, "Quad8", "Biot", "LinElastic", prms.CStr(), "ZERO", "gam=20", false));
	}
	else
	{
		eatts.Push (make_tuple(-1, "Quad4", "Biot", "LinElastic", prms.CStr(), "ZERO", "gam=20", true ));
		eatts.Push (make_tuple(-2, "Quad4", "Biot", "LinElastic", prms.CStr(), "ZERO", "gam=20", false));
		eatts.Push (make_tuple(-3, "Quad4", "Biot", "LinElastic", prms.CStr(), "ZERO", "gam=20", false));
	}

	// Set geometry: nodes, elements, attributes, and boundaries
	dat.SetNodesElems (&mesh, &eatts);

	// Solver
	FEM::Solver sol(dat,"tembed02");

	// Stage # -1 --------------------------------------------------------------
	FEM::EBrys_T ebrys;
	ebrys.Push           (make_tuple(-10, "ux",  0.0));
	ebrys.Push           (make_tuple(-11, "uy",  0.0));
	ebrys.Push           (make_tuple(-12, "pwp", 0.0));
	dat.SetBrys         (&mesh, NULL, &ebrys, NULL);
	dat.ApplyBodyForces    ();
	sol.SolveWithInfo   (10, 1e+2, /*iStage*/-1, "  Initial stress state due to self weight (zero displacements)\n");
	dat.ClearDisplacements ();

	// Stage # 0 ---------------------------------------------------------------
	dat.Activate         (/*Tag*/-2);
    ebrys.Resize       (0);
	ebrys.Push         (make_tuple(-10, "ux",  0.0));
	ebrys.Push         (make_tuple(-11, "uy",  0.0));
	ebrys.Push         (make_tuple(-13, "pwp", 0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo (10, 1e+2, 0, "  Construction of first layer\n");

	// Stage # 1 ---------------------------------------------------------------
	dat.Activate         (/*Tag*/-3);
    ebrys.Resize       (0);
	ebrys.Push         (make_tuple(-10, "ux",  0.0));
	ebrys.Push         (make_tuple(-11, "uy",  0.0));
	ebrys.Push         (make_tuple(-14, "pwp", 0.0));
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo (10, 1e+2, 0, "  Construction of second layer\n");

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	return 0;
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
