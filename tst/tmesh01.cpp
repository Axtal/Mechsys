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
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/matrix.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;

int main(int argc, char **argv) try
{
	double errors = 0.0;
	
	/////////////////////////////////////////////////////////////////////////////////////////// 2D /////
	
	{
		// Constants
		double L  = 1.0;   // length
		double H  = 1.0;   // height
		double hL = L/2.0; // half-length
		double hH = H/2.0; // half-height

		// Blocks
		Mesh::Block b;
		b.SetTag    (-1); // tag to be replicated to all generated elements inside this block
		b.SetCoords (false, 8,                          // Is3D, NNodes
		             0.,  L, L, 0.,    hL,  L, hL, 0.,  // x coordinates
		             0., 0., H,  H,    0., hH,  H, hH); // y coordinates
		b.SetETags  (4,  -10, -11, -12, -13);           // edge tags
		b.SetWx     (2, 1.0,1.0);                       // x weights and num of divisions along x
		b.SetWy     (2, 1.0,1.0);                       // y weights and num of divisions along y
		Array<Mesh::Block*> blocks;
		blocks.Push (&b);

		// Generate
		std::clock_t start = std::clock(); // Initial time
		Mesh::Structured ms(/*Is3D*/false);
		size_t ne = ms.Generate (blocks);
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "2D("<<ne<<" elements): Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		// Check
		//cout << "ETags:" << endl;
		//for (size_t i=0; i<ms.ElemsBry().Size(); ++i)
		//{
			//ms.ElemsBry()[i]->ETags.SetNS(Util::_6_3);
			//cout << ms.ElemsBry()[i]->ETags;
		//}

		// Output
		ms.WriteVTU ("tmesh01_2D.vtu");
		cout << "[1;34mFile <tmesh01_2D.vtu> created[0m" << endl << endl;
	}

	/////////////////////////////////////////////////////////////////////////////////////////// 3D /////
	
	{
		// Constants
		double L  = 1.0;   // length
		double H  = 1.0;   // height
		double D  = 1.0;   // depth
		double hL = L/2.0; // half-length
		double hH = H/2.0; // half-height
		double hD = D/2.0; // half-depth

		// Coordinates
		Mesh::Block b;
		b.SetTag    (-2); // tag to be replicated to all generated elements inside this block
		b.SetCoords (true, 20, // Is3D, NNodes
		             0.,  L,  L, 0.,    0.,  L, L, 0.,    hL,  L, hL, 0.,    hL,  L, hL, 0.,    0.,  L,  L, 0.,  // x coordinates
		             0., 0.,  H,  H,    0., 0., H,  H,    0., hH,  H, hH,    0., hH,  H, hH,    0., 0.,  H,  H,  // y coordinates
		             0., 0., 0., 0.,     D,  D, D,  D,    0., 0., 0., 0.,     D,  D,  D,  D,    hD, hD, hD, hD); // z coordinates
		b.SetWx     (2, 1.0,1.0); // x weights and num of divisions along x
		b.SetWy     (2, 1.0,1.0); // y weights and num of divisions along y
		b.SetWz     (2, 1.0,1.0); // z weights and num of divisions along z
		b.SetETags  (12, -100, -101, -102, -103, -104, -105, -106, -107, -108, -109, -110, -111);
		b.SetFTags  (6,  -10, -11, -12, -13, -14, -15);
		Array<Mesh::Block*> blocks;
		blocks.Push(&b);

		// Generate
		std::clock_t start = std::clock(); // Initial time
		Mesh::Structured ms(/*Is3D*/true);
		size_t ne = ms.Generate (blocks);
		std::clock_t total = std::clock() - start; // Time elapsed
		std::cout << "3D:("<<ne<<" elements) Time elapsed = [1;31m" << static_cast<double>(total)/CLOCKS_PER_SEC << "[0m [1;32mseconds[0m" << std::endl;

		// Check
		cout << "FTags:" << endl;
		for (size_t eb=0; eb<ms.NElemsBry(); ++eb)
		{
			size_t i = ms.ElemBry(eb);
			for (size_t j=0; j<ms.ElemNFTags(i); ++j)
			{
				char buf[12];
				sprintf(buf,"%d",ms.ElemFTag(i,j));
				cout << "Elem # " << i << "\tFace # " << j << "\t FTag = " << buf << endl;
			}
		}

		// Output
		ms.WriteVTU ("tmesh01_3D.vtu");
		cout << "[1;34mFile <tmesh01_3D.vtu> created[0m" << endl;
	}

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
