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

/*  Six-bay bridge truss
    (from Prof. Carlos Felippa's FEM notes)
                                                                                                                                           
    E=1000                                                     A=10        (5)        A=10                                                 
                                                                [8] ___,----o----,___ [9]                                                  
                                                      (3)   ___,---'        |        '---,___   (7)                                        
                                         A=10      ____o---'                |                '---o____      A=10                           
                                          [7]__,--'    |',                  |                  ,'|    '--,__[10]                           
                                  (1) ___,--'          |  '-,               |               ,-'  |          '--,___ (9)                    
                                 _,o,'                 |     ',  A=1        |[14]    A=1  ,'     |                 ',o,_                   
                    A=10     _,-'  | '-,_     A=1      |       ',[18]       |A=3    [19],'       |[15]   A=1    _,-' |  '-,_     A=10      
                     [6] _,-'      |     '-,_ [17]     |[13]     ',         |         ,'         |A=3   [20]_,-'     |      '-,_ [11]      
                     _,-'          |[12]     '-,_      |A=3        '-,      |      ,-'           |      _,-'         |[16]      '-,_       
   y             _,-'              |A=3          '-,_  |              ',    |    ,'              |  _,-'             |A=3           '-,_   
   |       (0),-'        [0]       |       [1]       '-|         [2]    '-,_|_,-'     [3]        |-'      [4]        |        [5]       '-,(11)
   |         o---------------------o-------------------o--------------------o--------------------o-------------------o--------------------o
  (z)___x   /_\        A=2        (2)      A=2        (4)       A=2        (6)       A=2        (8)       A=2      (10)       A=2        /_\
           ////                    |                   |                    |                    |                   |                   ___
                                   V                   V                    V                    V                   V
                                  10                  10                   16                   10                  10                             */

// STL
#include <iostream>

// MechSys
#include "fem/data.h"
#include "fem/solver.h"
#include "fem/elems/rod2.h"
#include "fem/output.h"
#include "util/exception.h"

using std::cout;
using std::endl;
using Util::_4;
using Util::_6;
using Util::_8s;

int main(int argc, char **argv) try
{
	// Input
	cout << "Input: " << argv[0] << "  linsol(LA,UM,SLU)\n";
	String linsol("UM");
	if (argc==2) linsol.Printf("%s",argv[1]);

	// Geometry
	FEM::Data dat(2); // 2D

	// Nodes
	dat.SetNNodes (12);
	dat.SetNode   ( 0,  0.0, 0.0);
	dat.SetNode   ( 1, 10.0, 5.0);
	dat.SetNode   ( 2, 10.0, 0.0);
	dat.SetNode   ( 3, 20.0, 8.0);
	dat.SetNode   ( 4, 20.0, 0.0);
	dat.SetNode   ( 5, 30.0, 9.0);
	dat.SetNode   ( 6, 30.0, 0.0);
	dat.SetNode   ( 7, 40.0, 8.0);
	dat.SetNode   ( 8, 40.0, 0.0);
	dat.SetNode   ( 9, 50.0, 5.0);
	dat.SetNode   (10, 50.0, 0.0);
	dat.SetNode   (11, 60.0, 0.0);

	// Elements
	dat.SetNElems (21);
	for (int i=0; i<21; ++i) dat.SetElem(i, "Rod2", true, -1);

	// Set connectivity
	dat.Ele( 0)->Connect(0, dat.Nod( 0))->Connect(1, dat.Nod( 2));
	dat.Ele( 1)->Connect(0, dat.Nod( 2))->Connect(1, dat.Nod( 4));
	dat.Ele( 2)->Connect(0, dat.Nod( 4))->Connect(1, dat.Nod( 6));
	dat.Ele( 3)->Connect(0, dat.Nod( 6))->Connect(1, dat.Nod( 8));
	dat.Ele( 4)->Connect(0, dat.Nod( 8))->Connect(1, dat.Nod(10));
	dat.Ele( 5)->Connect(0, dat.Nod(10))->Connect(1, dat.Nod(11));
	dat.Ele( 6)->Connect(0, dat.Nod( 0))->Connect(1, dat.Nod( 1));
	dat.Ele( 7)->Connect(0, dat.Nod( 1))->Connect(1, dat.Nod( 3));
	dat.Ele( 8)->Connect(0, dat.Nod( 3))->Connect(1, dat.Nod( 5));
	dat.Ele( 9)->Connect(0, dat.Nod( 5))->Connect(1, dat.Nod( 7));
	dat.Ele(10)->Connect(0, dat.Nod( 7))->Connect(1, dat.Nod( 9));
	dat.Ele(11)->Connect(0, dat.Nod( 9))->Connect(1, dat.Nod(11));
	dat.Ele(12)->Connect(0, dat.Nod( 1))->Connect(1, dat.Nod( 2));
	dat.Ele(13)->Connect(0, dat.Nod( 3))->Connect(1, dat.Nod( 4));
	dat.Ele(14)->Connect(0, dat.Nod( 5))->Connect(1, dat.Nod( 6));
	dat.Ele(15)->Connect(0, dat.Nod( 7))->Connect(1, dat.Nod( 8));
	dat.Ele(16)->Connect(0, dat.Nod( 9))->Connect(1, dat.Nod(10));
	dat.Ele(17)->Connect(0, dat.Nod( 1))->Connect(1, dat.Nod( 4));
	dat.Ele(18)->Connect(0, dat.Nod( 3))->Connect(1, dat.Nod( 6));
	dat.Ele(19)->Connect(0, dat.Nod( 6))->Connect(1, dat.Nod( 7));
	dat.Ele(20)->Connect(0, dat.Nod( 8))->Connect(1, dat.Nod( 9));

	// Parameters and initial value
	dat.Ele( 0)->SetModel("", "E=1000.0  A= 2.0", "ZERO");
	dat.Ele( 1)->SetModel("", "E=1000.0  A= 2.0", "ZERO");
	dat.Ele( 2)->SetModel("", "E=1000.0  A= 2.0", "ZERO");
	dat.Ele( 3)->SetModel("", "E=1000.0  A= 2.0", "ZERO");
	dat.Ele( 4)->SetModel("", "E=1000.0  A= 2.0", "ZERO");
	dat.Ele( 5)->SetModel("", "E=1000.0  A= 2.0", "ZERO");
	dat.Ele( 6)->SetModel("", "E=1000.0  A=10.0", "ZERO");
	dat.Ele( 7)->SetModel("", "E=1000.0  A=10.0", "ZERO");
	dat.Ele( 8)->SetModel("", "E=1000.0  A=10.0", "ZERO");
	dat.Ele( 9)->SetModel("", "E=1000.0  A=10.0", "ZERO");
	dat.Ele(10)->SetModel("", "E=1000.0  A=10.0", "ZERO");
	dat.Ele(11)->SetModel("", "E=1000.0  A=10.0", "ZERO");
	dat.Ele(12)->SetModel("", "E=1000.0  A= 3.0", "ZERO");
	dat.Ele(13)->SetModel("", "E=1000.0  A= 3.0", "ZERO");
	dat.Ele(14)->SetModel("", "E=1000.0  A= 3.0", "ZERO");
	dat.Ele(15)->SetModel("", "E=1000.0  A= 3.0", "ZERO");
	dat.Ele(16)->SetModel("", "E=1000.0  A= 3.0", "ZERO");
	dat.Ele(17)->SetModel("", "E=1000.0  A= 1.0", "ZERO");
	dat.Ele(18)->SetModel("", "E=1000.0  A= 1.0", "ZERO");
	dat.Ele(19)->SetModel("", "E=1000.0  A= 1.0", "ZERO");
	dat.Ele(20)->SetModel("", "E=1000.0  A= 1.0", "ZERO");

	// Boundary conditions (must be after set connectivity)
	dat.Nod( 0)->Bry("ux",   0.0)->Bry("uy", 0.0);
	dat.Nod(11)->Bry("uy",   0.0);
	dat.Nod( 2)->Bry("fy", -10.0);
	dat.Nod( 4)->Bry("fy", -10.0);
	dat.Nod( 6)->Bry("fy", -16.0);
	dat.Nod( 8)->Bry("fy", -10.0);
	dat.Nod(10)->Bry("fy", -10.0);

	// Solve
	FEM::Solver sol(dat,"ttruss02");
	sol.SolveWithInfo(/*NDiv*/1, /*DTime*/0.0);

	//////////////////////////////////////////////////////////////////////////////////////// Check /////

	// Correct node displacements
	const double u_corr[][2] = {{ 0.000000,  0.000000 },
	                            { 0.809536, -1.775600 },
	                            { 0.280000, -1.792260 },
	                            { 0.899001, -2.291930 },
	                            { 0.560000, -2.316600 },
	                            { 0.847500, -2.385940 },
	                            { 0.847500, -2.421940 },
	                            { 0.795999, -2.291930 },
	                            { 1.135000, -2.316600 },
	                            { 0.885464, -1.775600 },
	                            { 1.415000, -1.792260 },
	                            { 1.695000,  0.000000 }};

	// Correct node forces including reactions
	const double f_corr[][2] = {{ 0.0,  28.0 },
	                            { 0.0,   0.0 },
	                            { 0.0, -10.0 },
	                            { 0.0,   0.0 },
	                            { 0.0, -10.0 },
	                            { 0.0,   0.0 },
	                            { 0.0, -16.0 },
	                            { 0.0,   0.0 },
	                            { 0.0, -10.0 },
	                            { 0.0,   0.0 },
	                            { 0.0, -10.0 },
	                            { 0.0,  28.0 }};

	Array<double> err_u;
	Array<double> err_f;
	for (size_t i=0; i<dat.NNodes(); ++i)
	{
		// Displacements
		err_u.Push ( fabs(dat.Nod(i)->Val("ux")-u_corr[i][0]) );
		err_u.Push ( fabs(dat.Nod(i)->Val("uy")-u_corr[i][1]) );

		// Forces
		err_f.Push ( fabs(dat.Nod(i)->Val("fx")-f_corr[i][0]) );
		err_f.Push ( fabs(dat.Nod(i)->Val("fy")-f_corr[i][1]) );
	}
	
	// Correct axial normal stresses
	Array<double> err_s(dat.NElems());
	err_s[ 0] = fabs(dat.Ele( 0)->Val(0, "Sa") - (28.0000));
	err_s[ 1] = fabs(dat.Ele( 1)->Val(0, "Sa") - (28.0000));
	err_s[ 2] = fabs(dat.Ele( 2)->Val(0, "Sa") - (28.7500));
	err_s[ 3] = fabs(dat.Ele( 3)->Val(0, "Sa") - (28.7500));
	err_s[ 4] = fabs(dat.Ele( 4)->Val(0, "Sa") - (28.0000));
	err_s[ 5] = fabs(dat.Ele( 5)->Val(0, "Sa") - (28.0000));
	err_s[ 6] = fabs(dat.Ele( 6)->Val(0, "Sa") - (-6.2610));
	err_s[ 7] = fabs(dat.Ele( 7)->Val(0, "Sa") - (-6.0030));
	err_s[ 8] = fabs(dat.Ele( 8)->Val(0, "Sa") - (-6.0300));
	err_s[ 9] = fabs(dat.Ele( 9)->Val(0, "Sa") - (-6.0300));
	err_s[10] = fabs(dat.Ele(10)->Val(0, "Sa") - (-6.0030));
	err_s[11] = fabs(dat.Ele(11)->Val(0, "Sa") - (-6.2610));
	err_s[12] = fabs(dat.Ele(12)->Val(0, "Sa") - ( 3.3330));
	err_s[13] = fabs(dat.Ele(13)->Val(0, "Sa") - ( 3.0830));
	err_s[14] = fabs(dat.Ele(14)->Val(0, "Sa") - ( 4.0000));
	err_s[15] = fabs(dat.Ele(15)->Val(0, "Sa") - ( 3.0830));
	err_s[16] = fabs(dat.Ele(16)->Val(0, "Sa") - ( 3.3330));
	err_s[17] = fabs(dat.Ele(17)->Val(0, "Sa") - ( 1.6770));
	err_s[18] = fabs(dat.Ele(18)->Val(0, "Sa") - ( 3.2020));
	err_s[19] = fabs(dat.Ele(19)->Val(0, "Sa") - ( 3.2020));
	err_s[20] = fabs(dat.Ele(20)->Val(0, "Sa") - ( 1.6770));

	// Error summary
	double tol_u     = 5.0e-6;
	double tol_f     = 7.0e-13;
	double tol_s     = 4.4e-4;
	double min_err_u = err_u[err_u.Min()];
	double max_err_u = err_u[err_u.Max()];
	double min_err_f = err_f[err_f.Min()];
	double max_err_f = err_f[err_f.Max()];
	double min_err_s = err_s[err_s.Min()];
	double max_err_s = err_s[err_s.Max()];
	cout << _4<< ""   << _8s<<"Min"     << _8s<<"Mean"                                                  << _8s<<"Max"                << _8s<<"Norm"       << endl;
	cout << _4<< "u"  << _8s<<min_err_u << _8s<<err_u.Mean() << (max_err_u>tol_u?"[1;31m":"[1;32m") << _8s<<max_err_u << "[0m" << _8s<<err_u.Norm() << endl;
	cout << _4<< "f"  << _8s<<min_err_f << _8s<<err_f.Mean() << (max_err_f>tol_f?"[1;31m":"[1;32m") << _8s<<max_err_f << "[0m" << _8s<<err_f.Norm() << endl;
	cout << _4<< "Sa" << _8s<<min_err_s << _8s<<err_s.Mean() << (max_err_s>tol_s?"[1;31m":"[1;32m") << _8s<<max_err_s << "[0m" << _8s<<err_s.Norm() << endl;
	cout << endl;

	// Return error flag
	if (max_err_u>tol_u || max_err_f>tol_f || max_err_s>tol_s) return 1;
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
