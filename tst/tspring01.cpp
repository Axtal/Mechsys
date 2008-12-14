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

/* Bhatti, Example 4.12, pag. 278 */

// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/beam.h"
#include "fem/elems/spring.h"
#include "models/equilibs/beamelastic.h"
#include "models/equilibs/springelastic.h"
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
	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

	Mesh::Generic mesh(false);
	mesh.SetNVerts   (3); // set number of vertices
	mesh.SetVert     (0, true,             0.0,            0.0); // VertIdx, OnBry==true, X, Y, Z
	mesh.SetVert     (1, true, 180.0/sqrt(2.0), 180.0/sqrt(2.0));
	mesh.SetVert     (2, true, 360.0/sqrt(2.0),            0.0);
	mesh.SetNElems   (2); // set number of elements
	mesh.SetElem     (0,-50,true,3); // ElemIdx, ETag, OnBry==true, VTK_LINE==3
	mesh.SetElemCon  (0,0,0); // ElemIdx, Vert1, Vert2
	mesh.SetElemCon  (0,1,1); // ElemIdx, Vert1, Vert2
	mesh.SetElemETag (0,0,-50); // ElemIdx, LocalEdgeID==0, Tag
	mesh.SetElem     (1,-50,true,3);
	mesh.SetElemCon  (1,0,1);
	mesh.SetElemCon  (1,1,2);
	mesh.SetElemETag (1,0,-50);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////

	// Data and Solver
	FEM::Data   dat (2);
	FEM::Solver sol (dat, "tspring01");

	// Element attributes
	FEM::EAtts_T eatts;
	eatts.Push (make_tuple(  -50, "", "Beam",   "BeamElastic",   "E=30000 A=0.5 Izz=1000", "ZERO", "gam=20 cq=1", true));
	eatts.Push (make_tuple( -200, "", "Spring", "SpringElastic", "ks=100",                 "ZERO", "",            true));

	// Set nodes and elements (geometry)
	dat.SetOnlyFrame  (); // frame (beam/truss) mesh only
	dat.SetNodesElems (&mesh, &eatts);

	// Add linear elements
	//Array<int> conn(2);  conn = 0, 2;
	//dat.PushElem (-200, "Spring", "", "ks=100", "ZERO", "gam=20", true, conn);

	// Print stiffness
	/*
	Matrix<double> Ke;  Ke.SetNS(Util::_12_4);
	dat.Ele(2)->Order1Matrix(0,Ke);
	cout << "Ke =\n" << Ke << endl;
	*/

	// Stage # 1 --------------------------------------------------------------
	dat.Nod(0)->Bry ("ux",0.0);
	dat.Nod(0)->Bry ("uy",0.0);
	dat.Nod(2)->Bry ("uy",0.0);
	dat.Nod(2)->Bry ("fx",50.0);
	dat.Nod(0)->Bry ("wz",0.0);
	dat.Nod(2)->Bry ("wz",0.0);
	sol.SolveWithInfo ();

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Displacements
	Array<double> err_u(9);
	err_u[ 0] = fabs(dat.Nod(0)->Val("ux") - (0.0));
	err_u[ 1] = fabs(dat.Nod(0)->Val("uy") - (0.0));
	err_u[ 2] = fabs(dat.Nod(0)->Val("wz") - (0.0));
	err_u[ 3] = fabs(dat.Nod(1)->Val("ux") - (0.184555));
	err_u[ 4] = fabs(dat.Nod(1)->Val("uy") - (-0.0274869));
	err_u[ 5] = fabs(dat.Nod(1)->Val("wz") - (0.0));
	err_u[ 6] = fabs(dat.Nod(2)->Val("ux") - (0.36911));
	err_u[ 7] = fabs(dat.Nod(2)->Val("uy") - (0.0));
	err_u[ 8] = fabs(dat.Nod(2)->Val("wz") - (0.0));

	Array<double> err_s(2);
	err_s[ 0] = fabs(dat.Ele(2)->Val(0,"N") - (36.911));
	err_s[ 1] = fabs(dat.Ele(2)->Val(1,"N") - (36.911));

	// Error summary
	double tol_u     = 1.0e-7;
	double tol_s     = 1.0e-5;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	double min_err_s = err_s[err_s.Min()];
	double max_err_s = err_s[err_s.Max()];
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << _4<< "s" << _8s<<min_err_s << _8s<<err_s.Mean() << (max_err_s>tol_s?"[1;31m":"[1;32m") << _8s<<max_err_s << "[0m" << _8s<<err_s.Norm() << endl;

	if (max_err_u>tol_u || max_err_s>tol_s) return 1;
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
