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
#include <cmath>

// MechSys
#include <mechsys/numerical/quadrature.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;
using Numerical::Quadrature;

class Test01
{ public: double Fun(double x) const { return cos(100.0*sin(x)); } };

class Test02
{ public: double Fun(double x) const { return sqrt(x)*log(x); } };

class Test03
{
public:
	double Fun(double x) const
	{
		double gam = 0.0001;
		double f   = pow(x,3.0)-8.0*pow(x,2.0)+3.0*x-2.0;
		double g   = 3*pow(x,2.0)-16.0*x+3.0;
		double h   = 6.0*x-16.0;
		return (f*h-g*g)/(f*f+gam*gam*g*g);
	}
};

int main(int argc, char **argv) try
{
	cout << "======================================  Test 01" << endl;
	Test01             t01;
	Quadrature<Test01> qu1(&t01, &Test01::Fun, /*method*/Numerical::QAGS_T);
	double             r01 = qu1.Integrate(0.0,Util::PI);
	cout << "Int{0,Pi}(cos(100*sin(x))) = " << r01 << endl;
	cout << "Error = " << fabs(0.062787400491491 - r01) << endl;
	cout << "Error estimative = " << qu1.LastError() << endl << endl;

	cout << "======================================  Test 02" << endl;
	Test02             t02;
	Quadrature<Test02> qu2(&t02, &Test02::Fun, /*method*/Numerical::QAGS_T);
	double             r02 = qu2.Integrate(0.0,1.0);
	cout << "Int{0,1}(sqrt(x)*log(x)) = " << r02 << endl;
	cout << "Error = " << fabs(-4.0/9.0 - r02) << endl;
	cout << "Error estimative = " << qu2.LastError() << endl << endl;

	cout << "======================================  Test 03" << endl;
	Test03             t03;
	Quadrature<Test03> qu3(&t03, &Test03::Fun, /*method*/Numerical::QAGS_T);
	double             r03 = qu3.Integrate(0.0,10.0);
	cout << "Int{0,10}(m(x)) = " << r03 << endl;
	cout << "Error estimative = " << qu3.LastError() << endl;

	return 0;
}
MECHSYS_CATCH
