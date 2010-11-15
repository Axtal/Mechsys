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
#include <mechsys/util/maps.h>
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

    cout << "\n/////////////////////////////////////////////////////////////////////////////////\n" << endl;

    Table tab;
    tab.SetZero ("A B C D", /*nrows*/4);
    tab("A") = 0.0, 0.1, 0.2, 0.3;
    tab("B") = 1.0, 1.1, 1.2, 1.3;
    tab("C") = 2.0, 2.1, 2.2, 2.3;
    tab("D") = 3.0, 3.1, 3.2, 3.3;
    tab.Write ("table_sandbox.dat", /*NF*/"%4.1f");
    cout << "File <table_sandbox.dat> written" << endl;
    Table tab2;
    tab2.Read ("table_sandbox.dat");
    cout << "Reading back Table from file:\n" << tab2 << endl;
    double error = 0.0;
    for (size_t i=0; i<4; ++i)
    {
        error += fabs(tab2("A",i) - tab("A",i));
        error += fabs(tab2("B",i) - tab("B",i));
        error += fabs(tab2("C",i) - tab("C",i));
        error += fabs(tab2("D",i) - tab("D",i));
    }
    if (error>1.0e-15) throw new Fatal("Table read/write failed");

    cout << "\n/////////////////////////////////////////////////////////////////////////////////\n" << endl;

    SDPair p1, p2;
    p1.Set ("k l m", 10., 20., 30.);
    p1.Set ("K L M", 100., 200., 300.);
    //p2.Set ("A B M", 11., 22., 33.); // error => values cannot be summed up (by design)
    p2.Set ("A B", 11., 22.); // error => values cannot be summed up (by design)
    p1 += p2;
    Dict D1, D2, D3;
    D1.Set (-1, "a b c", 1., 2., 3.);
    D1.Set (-1, "M N P", 7., 8., 9.);
    D1.Set (-2, "d e f", 4., 5., 6.);
    D2.Set (-3, "aa bb", 11., 22.);
    //D3.Set (-1, "MM N", 33., 44.); // error
    D3.Set (-1, "MM NN", 33., 44.);
    D1 += D2;
    D1 += D3;
    cout << "p1 : " << p1 << endl;
    cout << "p2 : " << p2 << endl;
    cout << "D1 =\n" << D1 << endl;
    cout << "D2 =\n" << D2 << endl;
    cout << "D3 =\n" << D3 << endl;

    cout << endl;
    return 0;
}
MECHSYS_CATCH
