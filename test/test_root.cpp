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

// Std Lib
#include <iostream>
#include <fstream>

// MechSys
#include <mechsys/numerical/root.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using Util::_8s;

class Prob1
{
public:
    double F    (double x, void * UserData) { return pow(x,3.0) - 0.165*pow(x,2.0) + 3.993e-4; }
    double dFdx (double x, void * UserData) { return 3.0*pow(x,2.0) - 0.33*x; }
    void   Cout () { cout << TERM_CLR1 << "Prob1: x^3 - 0.165*x^2 + 3.993e-4" << TERM_RST << endl; }
};

class Prob2
{
public:
    double F    (double x, void * UserData) { return pow(x-1.0, 3.0) + 0.512; }
    double dFdx (double x, void * UserData) { return 3.0*pow(x-1.0,2.0); }
    void   Cout () { cout << TERM_CLR1 << "Prob1: (x-1)^3 + 0.512" << TERM_RST << endl; }
};

class Prob3
{
public:
    double F    (double x, void * UserData) { return pow(x,3.0) - 0.03*pow(x,2.0) + 2.4e-6; }
    double dFdx (double x, void * UserData) { return 3.0*pow(x,2.0) - 0.06*x; }
    void   Cout () { cout << TERM_CLR1 << "Prob1: x^3 - 0.03*x^2 + 2.4e-6" << TERM_RST << endl; }
};

int main(int argc, char **argv) try
{
    cout << endl << TERM_BLACK_WHITE << "####################################### Prob1 ##########################################" << TERM_RST << endl;

    // solve
    Prob1 prob1;
    Numerical::Root<Prob1> root(&prob1, &Prob1::F, &Prob1::dFdx);
    double xsol = root.Solve (0.0, 0.11);
    cout << "\n============================= Brent ===============================\n";
    prob1.Cout ();
    cout << "xsol = " << xsol    << endl;
    cout << "It   = " << root.It << endl;

    cout << "\n============================ Newton ===============================\n";
    root.Scheme  = "Newton";
    root.Verbose = true;
    double x_guess = 0.5;
    cout << "\n---------------- with guess: x=" << x_guess << " ----------------\n";
    prob1.Cout ();
    xsol = root.Solve (0.0, 0.11, &x_guess);
    cout << "xsol = " << xsol    << endl;
    cout << "It   = " << root.It << endl;

    x_guess = 0.05;
    cout << "\n---------------- with guess: x=" << x_guess << " ----------------\n";
    prob1.Cout ();
    xsol = root.Solve (0.0, 0.11, &x_guess);
    cout << "xsol = " << xsol    << endl;
    cout << "It   = " << root.It << endl;

    cout << "\n----------------- no guess -----------------------\n";
    prob1.Cout ();
    xsol = root.Solve (0.0, 0.11);
    cout << "xsol = " << xsol    << endl;
    cout << "It   = " << root.It << endl;

    cout << endl << TERM_BLACK_WHITE << "####################################### Prob2 ##########################################" << TERM_RST << endl;

    Prob2 prob2;
    Numerical::Root<Prob2> root2(&prob2, &Prob2::F, &Prob2::dFdx);
    cout << "\n============================= Brent ===============================\n";
    prob2.Cout ();
    xsol = root2.Solve (-2.0, 6.0);
    cout << "xsol = " << xsol     << endl;
    cout << "It   = " << root2.It << endl;

    cout << "\n============================ Newton ===============================\n";
    root2.Scheme  = "Newton";
    root2.Verbose = true;
    x_guess = 0.5;
    cout << "\n---------------- with guess: x=" << x_guess << " ----------------\n";
    prob2.Cout ();
    xsol = root2.Solve (-2.0, 6.0, &x_guess);
    cout << "xsol = " << xsol    << endl;
    cout << "It   = " << root2.It << endl;

    x_guess = 5.0;
    cout << "\n---------------- with guess: x=" << x_guess << " ----------------\n";
    prob2.Cout ();
    xsol = root2.Solve (-2.0, 6.0, &x_guess);
    cout << "xsol = " << xsol    << endl;
    cout << "It   = " << root2.It << endl;

    cout << endl << TERM_BLACK_WHITE << "####################################### Prob3 ##########################################" << TERM_RST << endl;

    Prob3 prob3;
    Numerical::Root<Prob3> root3(&prob3, &Prob3::F, &Prob3::dFdx);
    cout << "\n============================= Brent ===============================\n";
    prob3.Cout ();
    xsol = root3.Solve (0.0, 0.025);
    cout << "xsol = " << xsol     << endl;
    cout << "It   = " << root3.It << endl;

    cout << "\n============================ Newton ===============================\n";
    root3.Scheme  = "Newton";
    root3.Verbose = true;
    x_guess = 0.01999;
    cout << "\n---------------- with guess: x=" << x_guess << " ----------------\n";
    prob3.Cout ();
    xsol = root3.Solve (0.0, 0.025, &x_guess);
    cout << "xsol = " << xsol    << endl;
    cout << "It   = " << root3.It << endl;

    x_guess = 0.02;
    cout << "\n---------------- with guess: x=" << x_guess << " ----------------\n";
    prob3.Cout ();
    xsol = root3.Solve (0.0, 0.025, &x_guess);
    cout << "xsol = " << xsol    << endl;
    cout << "It   = " << root3.It << endl;

    // end
    return 0;
}
MECHSYS_CATCH
