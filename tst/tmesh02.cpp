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
#include <ctime>  // for std::clock()

// MechSys
#include "mesh/structured.h"
#include "util/util.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/matrix.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::PI;

inline int Max2(int A, int B)        { return (A>B ? A       : B      ); }
inline int Max3(int A, int B, int C) { return (A>B ? A>C?A:C : B>C?B:C); }

int main(int argc, char **argv) try
{
	
	/////////////////////////////////////////////////////////////////////////////////////////// 2D /////
	
	{
		/*                |---------- L ----------|
		 *   -+- -+- -+-  o___________o___________o
		 *    |   |   |   |                      ,|
		 *    |   |   |   |       [b1]         ,' |
		 *    |   |   d   o                  ,'   |
		 *    |   f   |   |    y       x   ,'     |
		 *    |   |   |   |     ',   ,'  ,o       |
		 *    |   |  -+-  o-,_    '+'  ,'         |
		 *    H   |           o-,    ,'           o
		 *    |  -+- . . . . . . 'o '      [b0]   |
		 *    |   |               .',             |
		 *    |   e               .  o  y^        |
		 *    |   |               .   \  |        |
		 *    |   |               .   |  +-->x    |
		 *   -+- -+-      +----r----->o-----o-----o
		 *                        .
		 *                        .   |---- a ----|
		 *                |-- b --|------ c ------|
		 */

		// Constants
		double r = 0.5; // radius
		double L = 1.0; // length
		double H = 1.0; // height
		double a = L-r;
		double b = r*cos(2.*PI/8.);
		double c = L-b;
		double d = H-r;
		double e = r*sin(2.*PI/8.);
		double f = H-e;

		// Lower block -- coordinates
		Matrix<double> c0(2,8);
		c0 =  r,  L, L, b,    r+a/2.,    L, b+c/2., r*cos(PI/8.),
		     0., 0., H, e,        0., H/2., e+f/2., r*sin(PI/8.);

		// Lower block -- weights
		int ndivx0 = 20;
		int ndivy0 = 30;
		Matrix<double> w0(2,Max2(ndivx0,ndivy0)+1);
		for (int i=0; i<ndivx0; ++i) w0(0,i) = 1.0;
		for (int i=0; i<ndivy0; ++i) w0(1,i) = 1.0;

		// Upper block -- coordinates
		Matrix<double> c1(2,8);
		c1 =  b, L, 0., 0.,   b+c/2., L/2.,     0., r*cos(3.*PI/8.),
		      e, H,  H,  r,   e+f/2.,    H, r+d/2., r*sin(3.*PI/8.);

		//for (int i=0; i<8; ++i) c1(1,i) += 1.;

		// Upper block -- weights
		int ndivx1 = ndivx0;
		int ndivy1 = 40;
		Matrix<double> w1(2,Max2(ndivx1,ndivy1)+1);
		for (int i=0; i<ndivx1; ++i) w1(0,i) = 1.0;
		for (int i=0; i<ndivy1; ++i) w1(1,i) = 1.0;

		// Blocks
		Array<Mesh::Block> blocks;  blocks.Resize(2);
		blocks[0].Set (c0, w0, ndivx0, ndivy0);
		blocks[1].Set (c1, w1, ndivx1, ndivy1);

		// Generate
		std::clock_t start = std::clock(); // Initial time
		Mesh::Structured m;
		size_t ne = m.Generate (blocks);
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "2D:("<<ne<<" elements) Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		// Output
		m.WriteVTU ("tmesh02_2D.vtu");
		cout << "[1;34mFile <tmesh02_2D.vtu> created[0m" << endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////// 3D /////

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
