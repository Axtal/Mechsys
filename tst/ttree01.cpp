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
#include <cmath>    // for fabs

// MechSys
#include "util/tree.h"
#include "util/exception.h"
#include "util/numstreams.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	double errors = 0.0;

	{
		cout << "------------------------------------------------ Test 1 ---\n";

		// Edges
		Array<long> edges(24);
		edges = 0, 3,  // edge  0
				1, 2,  // edge  1
				0, 1,  // edge  2
				2, 3,  // edge  3
				4, 7,  // edge  4
				5, 6,  // edge  5
				4, 5,  // edge  6
				6, 7,  // edge  7
				0, 4,  // edge  8
				1, 5,  // edge  9
				2, 6,  // edge 10
				3, 7;  // edge 11

		// Tree
		Util::Tree t1(edges);
		Util::Tree t2(edges);

		// Remove edge
		t1.DelEdge (0,1);
		t2 = t1;
		t1.Reset (edges);
		cout << t1 << endl;
		cout << t2 << endl;

		// Shortest path
		Array<long> path; path.SetNS(Util::_4);
		t2.ShortPath (1,3, path);
		cout << path << endl;
	}

	{
		cout << "------------------------------------------------ Test 2 ---\n";

		// Edges
		Array<long> edges(48);
		edges = 0, 11,  //  0
		        3, 11,
		        1,  9,  //  1
		        2,  9,
		        0,  8,  //  2
		        1,  8,
		        2, 10,  //  3
		        3, 10,
		        4, 15,  //  4
		        7, 15,
		        5, 13,  //  5
		        6, 13,
		        4, 12,  //  6
		        5, 12,
		        6, 14,  //  7
		        7, 14,
		        0, 16,  //  8
		        4, 16,
		        1, 17,  //  9
		        5, 17,
		        2, 18,  // 10
		        6, 18,
		        3, 19,  // 11
		        7, 19;

		// Path arrays
		Array<long> path; path.SetNS(Util::_4);
		Array<long> cor1; cor1.Resize(7);
		Array<long> cor2; cor2.Resize(7);
		cor1 =  8,  1,  9,  2, 10,  3, 11;
		cor2 = 12,  5, 13,  6, 14,  7, 15;

		// Path 1
		Util::Tree tree(edges);
		tree.DelEdge   (0,8);
		tree.ShortPath (8,11, path);
		cout << path << endl;
		for (size_t i=0; i<path.Size(); ++i) errors += fabs(path[i]-cor1[i]);

		// Path 2
		tree.Reset     (edges);
		tree.DelEdge   (4,12);
		tree.ShortPath (12,15, path);
		cout << path << endl;
		for (size_t i=0; i<path.Size(); ++i) errors += fabs(path[i]-cor2[i]);

		// Errors
		if (fabs(errors)>1.0e-14) cout << "[1;31m\nErrors = " << errors << "[0m\n" << endl;
		else                      cout << "[1;32m\nErrors = " << errors << "[0m\n" << endl;
	}

	// Return error flag
	if (fabs(errors)>1.0e-14) return 1;
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
