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
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    cout << "\n/////////////////////////////////////////////////////////////////////////////////\n" << endl;

    int a,b,c,d;
    a = 10;
    b = 2;
    c = 7;
    d = 1;
    cout << "Sorting: " << a << " " << b << " " << c << " " << d << endl;
    Util::Sort(a,b,c,d);
    cout << "Result:  " << a << " " << b << " " << c << " " << d << endl;

    cout << "\n/////////////////////////////////////////////////////////////////////////////////\n" << endl;

    String keys("ux uy uz  fx fy fz");
    cout << "Keys: " << keys << endl;
    cout << "Is 'fx' in Keys ? => " << (Util::HasKey(keys,"fx") ? " YES " : " NO ") << endl;
    cout << "Is 'FX' in Keys ? => " << (Util::HasKey(keys,"FX") ? " YES " : " NO ") << endl;
    cout << "Is 'Fx' in Keys ? => " << (Util::HasKey(keys,"Fx") ? " YES " : " NO ") << endl;
    cout << "Is 'fX' in Keys ? => " << (Util::HasKey(keys,"fX") ? " YES " : " NO ") << endl;
    cout << "Is 'uy' in Keys ? => " << (Util::HasKey(keys,"uy") ? " YES " : " NO ") << endl;

    return 0;
}
MECHSYS_CATCH
