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
#include "mesh/alphashape.h"
#include "util/array.h"
#include "util/exception.h"
#include "util/numstreams.h"
#include "linalg/matrix.h"

using std::cout;
using std::endl;
using LinAlg::Matrix;

int main(int argc, char **argv) try
{
	Mesh::AlphaShape ma;
	ma.AddCloudPoint(0.084, 0.611); //  0, 
	ma.AddCloudPoint(0.086, 0.444); //  1, 
	ma.AddCloudPoint(0.133, 0.701); //  2, 
	ma.AddCloudPoint(0.126, 0.341); //  3, 
	ma.AddCloudPoint(0.198, 0.780); //  4, 
	ma.AddCloudPoint(0.183, 0.196); //  5, 
	ma.AddCloudPoint(0.250, 0.100); //  6, 
	ma.AddCloudPoint(0.232, 0.335); //  7, 
	ma.AddCloudPoint(0.232, 0.487); //  8, 
	ma.AddCloudPoint(0.214, 0.630); //  9, 
	ma.AddCloudPoint(0.270, 0.859); // 10, 
	ma.AddCloudPoint(0.290, 0.721); // 11, 
	ma.AddCloudPoint(0.308, 0.630); // 12, 
	ma.AddCloudPoint(0.320, 0.222); // 13, 
	ma.AddCloudPoint(0.346, 0.379); // 14, 
	ma.AddCloudPoint(0.370, 0.049); // 15, 
	ma.AddCloudPoint(0.451, 0.112); // 16, 
	ma.AddCloudPoint(0.404, 0.525); // 17, 
	ma.AddCloudPoint(0.412, 0.750); // 18, 
	ma.AddCloudPoint(0.406, 0.906); // 19, 
	ma.AddCloudPoint(0.485, 0.815); // 20, 
	ma.AddCloudPoint(0.501, 0.685); // 21, 
	ma.AddCloudPoint(0.498, 0.589); // 22, 
	ma.AddCloudPoint(0.488, 0.440); // 23, 
	ma.AddCloudPoint(0.459, 0.329); // 24, 
	ma.AddCloudPoint(0.498, 0.223); // 25, 
	ma.AddCloudPoint(0.579, 0.269); // 26, 
	ma.AddCloudPoint(0.596, 0.466); // 27, 
	ma.AddCloudPoint(0.609, 0.582); // 28, 
	ma.AddCloudPoint(0.593, 0.709); // 29, 
	ma.AddCloudPoint(0.638, 0.375); // 30, 
	ma.AddCloudPoint(0.659, 0.225); // 31, 
	ma.AddCloudPoint(0.752, 0.415); // 32, 
	ma.AddCloudPoint(0.721, 0.520); // 33, 
	ma.AddCloudPoint(0.741, 0.624); // 34, 
	ma.AddCloudPoint(0.760, 0.275); // 35, 
	ma.AddCloudPoint(0.705, 0.134); // 36, 
	ma.AddCloudPoint(0.803, 0.140); // 37, 
	ma.AddCloudPoint(0.863, 0.262); // 38, 
	ma.AddCloudPoint(0.881, 0.430); // 39, 
	ma.AddCloudPoint(0.867, 0.580); // 40, 

	cout << ma.Generate() << " elements generated" << endl;
	ma.WriteVTU("talphashape01.vtu");
	cout << "File <talphashape01.vtu> generated" << endl;

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
