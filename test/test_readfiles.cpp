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
#include <mechsys/inpfile.h>
#include <mechsys/matfile.h>
#include <mechsys/datfiles.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    InpFile  inp;
    MatFile  mat;
    DatFiles dat;

    inp.Read ("main.inp");
    mat.Read ("materials.mat");
    dat.Read (Array<String>("FO1.CTP.01.dat", /*justone*/true));

    cout << "\n############################# main.inp #############################\n";
    cout << inp << endl;

    cout << "\n########################## materials.mat ###########################\n";
    cout << mat << endl;

    cout << "\n######################### FO1.CTP.01.dat ###########################\n";
    cout << "sx = \n" << dat.sx[0] << endl;

    return 0;
}
MECHSYS_CATCH
