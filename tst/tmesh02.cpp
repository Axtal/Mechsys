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
#include "util/array.h"
#include "util/util.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/matrix.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;
using Util::PI;

int main(int argc, char **argv) try
{
	double errors = 0.0;

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
		Mesh::Block b0;
		b0.SetCoords(false, 8, // Is3D, NNodes
		              r,  L, L, b,    r+a/2.,    L, b+c/2., r*cos(PI/8.),
		             0., 0., H, e,        0., H/2., e+f/2., r*sin(PI/8.));
		b0.SetNx (10);
		b0.SetNy (10);

		// Upper block -- coordinates
		Mesh::Block b1;
		b1.SetCoords(false, 8,
		              b, L, 0., 0.,   b+c/2., L/2.,     0., r*cos(3.*PI/8.),
		              e, H,  H,  r,   e+f/2.,    H, r+d/2., r*sin(3.*PI/8.));
		b1.SetNx (10);
		b1.SetNy (10);

		// Blocks
		Array<Mesh::Block*> blocks;  blocks.Resize(2);
		blocks[0] = &b0;
		blocks[1] = &b1;

		// Generate
		std::clock_t start = std::clock(); // Initial time
		Mesh::Structured m(/*Is3D*/false);
		size_t ne = m.Generate (blocks);
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "2D:("<<ne<<" elements) Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		// Output
		m.WriteVTU ("tmesh02_2D.vtu");
		cout << "[1;34mFile <tmesh02_2D.vtu> created[0m" << endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////// 3D /////

	// Return error flag
	if (fabs(errors)>1.0e-13) return 1;
	else                      return 0;
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
