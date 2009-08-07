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
#include "fem/elems/beam.h"
#include "models/equilibs/beamelastic.h"
#include "util/fatal.h"
#include "mesh/mesh.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::_4;
using Util::_6;
using Util::_8s;
using boost::make_tuple;

int main(int argc, char **argv) try
{
	/*  _  0             1     2     3     4      5
	 *  _|\@-------------@-----@-----@-----@------@
	 *  _|/       0      |  1     2     3,-|   4
	 *                   |            ,-'  |
	 *                  5|         ,-'     |
	 *                   |      ,-'6       |7
	 *                   |   ,-'           |
	 *                   |,-'              |
	 *                 6 @                 |
	 *                  ###                @ 7
	 *                                    ###
	 */
	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("LA");
	if (argc==2) linsol.Printf("%s",argv[1]);

	///////////////////////////////////////////////////////////////////////////////////////// Mesh /////
	
	Mesh::Generic mesh(/*Is3D*/false);
	mesh.SetNVerts (8);
	mesh.SetNElems (8);
	mesh.SetVert   (0, true,  0.0, 5.0);
	mesh.SetVert   (1, true,  6.0, 5.0);
	mesh.SetVert   (2, true,  8.0, 5.0);
	mesh.SetVert   (3, true, 10.0, 5.0);
	mesh.SetVert   (4, true, 12.0, 5.0);
	mesh.SetVert   (5, true, 14.0, 5.0);
	mesh.SetVert   (6, true,  6.0, 1.0);
	mesh.SetVert   (7, true, 12.0, 0.0);
	mesh.SetElem   (0, -5, true, VTK_LINE);
	mesh.SetElem   (1, -5, true, VTK_LINE);
	mesh.SetElem   (2, -5, true, VTK_LINE);
	mesh.SetElem   (3, -5, true, VTK_LINE);
	mesh.SetElem   (4, -5, true, VTK_LINE);
	mesh.SetElem   (5, -6, true, VTK_LINE);
	mesh.SetElem   (6, -6, true, VTK_LINE);
	mesh.SetElem   (7, -6, true, VTK_LINE);
	mesh.SetElemCon(0, 0, 0);  mesh.SetElemCon(0, 1, 1);
	mesh.SetElemCon(1, 0, 1);  mesh.SetElemCon(1, 1, 2);
	mesh.SetElemCon(2, 0, 2);  mesh.SetElemCon(2, 1, 3);
	mesh.SetElemCon(3, 0, 3);  mesh.SetElemCon(3, 1, 4);
	mesh.SetElemCon(4, 0, 4);  mesh.SetElemCon(4, 1, 5);
	mesh.SetElemCon(5, 0, 6);  mesh.SetElemCon(5, 1, 1);
	mesh.SetElemCon(6, 0, 6);  mesh.SetElemCon(6, 1, 4);
	mesh.SetElemCon(7, 0, 7);  mesh.SetElemCon(7, 1, 4);

	////////////////////////////////////////////////////////////////////////////////////////// FEM /////
	
	// Data and Solver
	FEM::Data   dat (2); // 2D
	FEM::Solver sol (dat, "tbeam02");

	// Elements attributes
	double E1   = 1.0;
	double A1   = 5.0e+9;
	double Izz1 = 6.0e+4;
	double E2   = 1.0;
	double A2   = 1.0e+9;
	double Izz2 = 2.0e+4;
	String prms1; prms1.Printf("E=%f A=%f Izz=%f",E1,A1,Izz1);
	String prms2; prms2.Printf("E=%f A=%f Izz=%f",E2,A2,Izz2);
	FEM::EAtts_T eatts;
	eatts.Push (make_tuple(-5, "Lin2", "Beam", "BeamElastic", prms1.CStr(), "ZERO", "gam=20 cq=1", FNULL, true));
	eatts.Push (make_tuple(-6, "Lin2", "Beam", "BeamElastic", prms2.CStr(), "ZERO", "gam=20 cq=1", FNULL, true));

	// Set geometry: nodes and elements
	dat.SetOnlyFrame  (true);
	dat.SetNodesElems (&mesh, &eatts);

	// Stage # 1 -----------------------------------------------------------
	dat.Ele(0)->EdgeBry("Qb", -20.0, -20.0, 0);
	dat.Ele(1)->EdgeBry("Qb", -20.0, -20.0, 0);
	dat.Ele(2)->EdgeBry("Qb", -20.0, -20.0, 0);
	dat.Ele(3)->EdgeBry("Qb", -20.0, -20.0, 0);
	dat.Ele(4)->EdgeBry("Qb", -20.0, -20.0, 0);
	dat.Nod(2)->Bry("fy", -60.0);
	dat.Nod(3)->Bry("fy", -60.0);
	dat.Nod(0)->Bry("ux", 0.0)->Bry("uy", 0.0);
	dat.Nod(6)->Bry("ux", 0.0)->Bry("uy", 0.0)->Bry("wz", 0.0);
	dat.Nod(7)->Bry("ux", 0.0)->Bry("uy", 0.0)->Bry("wz", 0.0);
	sol.SolveWithInfo();

	//////////////////////////////////////////////////////////////////////////////////////// Output ////
	
	// Output: Nodes
	cout << _6<<"Node #" << _8s<<"ux" << _8s<<"uy" << _8s<<"wz" << _8s<<"fx"<< _8s<<"fy" << _8s<<"mz" << endl;
	for (size_t i=0; i<dat.NNodes(); ++i)
		cout << _6<<i << _8s<<dat.Nod(i)->Val("ux") <<  _8s<<dat.Nod(i)->Val("uy") << _8s<<dat.Nod(i)->Val("wz") << _8s<<dat.Nod(i)->Val("fx") << _8s<<dat.Nod(i)->Val("fy") << _8s<<dat.Nod(i)->Val("mz") << endl;
	cout << endl;

	// Output: Elements
	cout << _6<<"Elem #" << _8s<<"N0" << _8s<<"M0" << _8s<<"V0" << _8s<<"N1" << _8s<<"M1" << _8s<<"V1" << endl;
	for (size_t i=0; i<dat.NElems(); ++i)
	{
		dat.Ele(i)->CalcDeps();
		cout << _6<<i << _8s<<dat.Ele(i)->Val(0, "N") << _8s<<dat.Ele(i)->Val(0, "M") << _8s<<dat.Ele(i)->Val(0, "V");
		cout <<          _8s<<dat.Ele(i)->Val(1, "N") << _8s<<dat.Ele(i)->Val(1, "M") << _8s<<dat.Ele(i)->Val(1, "V") << endl;
	}
	cout << endl;

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Displacements
	Array<double> err_u;
    err_u.Push(fabs(dat.Nod(0)->Val("ux") - 0.000000000000000E+00));
    err_u.Push(fabs(dat.Nod(1)->Val("ux") - 0.364510505250599E-07));
    err_u.Push(fabs(dat.Nod(4)->Val("ux") - 0.643549267921723E-07));
    err_u.Push(fabs(dat.Nod(5)->Val("ux") - 0.643549267921723E-07));
    err_u.Push(fabs(dat.Nod(6)->Val("ux") - 0.000000000000000E+00));
    err_u.Push(fabs(dat.Nod(7)->Val("ux") - 0.000000000000000E+00));
    err_u.Push(fabs(dat.Nod(0)->Val("uy") -  0.000000000000000E+00));
    err_u.Push(fabs(dat.Nod(1)->Val("uy") - -0.831942398779212E-06));
    err_u.Push(fabs(dat.Nod(4)->Val("uy") - -0.628326321776739E-06));
    err_u.Push(fabs(dat.Nod(5)->Val("uy") -  0.287991471238717E-02));
    err_u.Push(fabs(dat.Nod(6)->Val("uy") -  0.000000000000000E+00));
    err_u.Push(fabs(dat.Nod(7)->Val("uy") -  0.000000000000000E+00));
    err_u.Push(fabs(dat.Nod(0)->Val("wz") - -0.102535585850665E-02));
    err_u.Push(fabs(dat.Nod(1)->Val("wz") - -0.949704254186087E-03));
    err_u.Push(fabs(dat.Nod(4)->Val("wz") -  0.177354929713225E-02));
    err_u.Push(fabs(dat.Nod(5)->Val("wz") -  0.132921596379892E-02));
    err_u.Push(fabs(dat.Nod(6)->Val("wz") -  0.000000000000000E+00));
    err_u.Push(fabs(dat.Nod(7)->Val("wz") -  0.000000000000000E+00));

	// Error summary
	double tol_u     = 1.0e-6;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	cout << endl;
	cout << _4<< ""  << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u" << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;

	// Return error flag
	if (max_err_u>tol_u) return 1;
	else return 0;
}
MECHSYS_CATCH
