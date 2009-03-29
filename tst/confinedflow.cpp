// Std Lib
#include <iostream>

// MechSys
#include "mechsys.h"
#include "util/exception.h"

#define T boost::make_tuple

int main(int argc, char **argv) try
{

	////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	Array<Mesh::Block> bks(3);

	// Block # 0 ----------------------------------------------------------------
	Mesh::Verts_T ve0(4);
	Mesh::Edges_T ed0(4);
	Mesh::ETags_T et0(3);
	ve0 = T(0,0,0,0.0), T(3,0,10,0.0), T(4,20,10,0.0), T(6,20,0,0.0);
	ed0 = T(0,3), T(3,4), T(4,6), T(0,6);
	et0 = T(0,6,-50), T(0,3,-51), T(3,4,-52);
	bks[0].Set (-1, ve0, ed0, &et0, NULL, 0, 6, 3, -1);
	bks[0].SetNx (20,0,0);
	bks[0].SetNy (10,0,0);

	// Block # 1 ----------------------------------------------------------------
	Mesh::Verts_T ve1(4);
	Mesh::Edges_T ed1(4);
	Mesh::ETags_T et1(2);
	ve1 = T(4,20,10,0.0), T(5,30,10,0.0), T(6,20,0,0.0), T(7,30,0,0.0);
	ed1 = T(4,5), T(4,6), T(5,7), T(6,7);
	et1 = T(6,7,-50), T(4,5,-50);
	bks[1].Set (-1, ve1, ed1, &et1, NULL, 6, 7, 4, -1);
	bks[1].SetNx (10,0,0);
	bks[1].SetNy (10,0,0);

	// Block # 2 ----------------------------------------------------------------
	Mesh::Verts_T ve2(4);
	Mesh::Edges_T ed2(4);
	Mesh::ETags_T et2(3);
	ve2 = T(1,50,0,0.0), T(2,50,10,0.0), T(5,30,10,0.0), T(7,30,0,0.0);
	ed2 = T(1,2), T(2,5), T(5,7), T(1,7);
	et2 = T(1,7,-50), T(1,2,-51), T(2,5,-53);
	bks[2].Set (-1, ve2, ed2, &et2, NULL, 7, 1, 5, -1);
	bks[2].SetNx (20,0,0);
	bks[2].SetNy (10,0,0);

	// Generate
	Mesh::Structured mesh(false); // false=>2D
	mesh.SetBlocks (bks);
	mesh.Generate  (true); // true=>WithInfo

	////////////////////////////////////////////////////////////////////////////////////// FEM //////

	// Data and Solver
	FEM::Data   dat (2);
	FEM::Solver sol (dat, "confinedflow");

	// Element attributes
	FEM::EAtts_T eatts(1);
	eatts = T(-1, "Quad4", "Diffusion", "LinDiffusion", "k=1e-05", "ZERO", "s=0.0", FNULL, true );  // material

	// Set nodes and elements (geometry)
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 --------------------------------------------------------------
	FEM::EBrys_T ebrys(4);
	ebrys = T(-50,"f",0), T(-51,"f",0), T(-52,"u",10), T(-53,"u",0);
	dat.SetBrys       (&mesh, NULL, &ebrys, NULL);
	sol.SolveWithInfo (1, 1, 1, "simulation stage\n");

}
catch (Exception * e) { e->Cout();  if (e->IsFatal()) {delete e; exit(1);}  delete e; }
catch (char const * m) { std::cout<<"Fatal: "<<m<<std::endl;  exit(1); }
catch (...) { std::cout << "Some exception (...) ocurred\n"; } 
