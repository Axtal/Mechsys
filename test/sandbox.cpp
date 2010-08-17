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
#include <mechsys/util/stopwatch.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    Util::Stopwatch stopwatch;

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
    Array<String> array;
    Util::Keys2Array (keys,array);
    cout << "Keys: " << keys << endl;
    cout << "Is 'fx' in Keys ? => " << (Util::HasKey(keys,"fx") ? " YES " : " NO ") << endl;
    cout << "Is 'FX' in Keys ? => " << (Util::HasKey(keys,"FX") ? " YES " : " NO ") << endl;
    cout << "Is 'Fx' in Keys ? => " << (Util::HasKey(keys,"Fx") ? " YES " : " NO ") << endl;
    cout << "Is 'fX' in Keys ? => " << (Util::HasKey(keys,"fX") ? " YES " : " NO ") << endl;
    cout << "Is 'uy' in Keys ? => " << (Util::HasKey(keys,"uy") ? " YES " : " NO ") << endl;
    cout << "Array of keys = " << array << endl;

    cout << "\n/////////////////////////////////////////////////////////////////////////////////\n" << endl;

    Array<int> A(7); A = 1, 4, 2, 6, 0, 3, -3;
    long pos = A.Find(8);
    cout << "A         = " << A << endl;
    cout << "Is 8 in A = " << (pos<0 ? "false" : "true: pos = ");
    if (pos>=0) cout << pos << endl;
    else cout << endl;
    pos = A.Find(6);
    cout << "Is 6 in A = " << (pos<0 ? "false" : "true: pos = ");
    if (pos>=0) cout << pos << endl;
    else cout << endl;
    cout << "min(A)    = " << A.TheMin() << endl;
    cout << "First(A)  = " << (*A.GetPtr()) << endl;
    cout << "Last(A)   = " << A.Last() << endl;
    cout << endl;

    return 0;
}
MECHSYS_CATCH
