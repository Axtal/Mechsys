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
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/util/array.h>
#include <mechsys/util/colors.h>
#include <mechsys/util/stopwatch.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    int seed = 123;
    if (argc>1) seed = atoi(argv[1]);
    srand(seed);
    Array<String> clrs;
    Colors::GetRandom (10, clrs);
    for (size_t i=0; i<clrs.Size(); ++i) cout << clrs[i] << endl;
    return 0;
}
MECHSYS_CATCH
