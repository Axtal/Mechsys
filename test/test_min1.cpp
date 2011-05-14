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
#include <mechsys/numerical/min.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using Util::_8s;

// Paraboloid centered at (xc,yc) with scale factors (a,b) and minimum zmin
class Prob1
{
public:
    Prob1 () : xc(1.0), yc(2.0), a(10.0), b(20.0), zmin(30.0) {}
    double F  (double const X[]) { return a * pow(X[0] - xc, 2.0) + b * pow(X[1] - yc, 2.0) + zmin; }
    void   dF (double const X[], double dFdX[])
    {
        dFdX[0] = 2.0 * a * (X[0] - xc);
        dFdX[1] = 2.0 * b * (X[1] - yc);
    }
    double xc,yc,a,b,zmin;
};

int main(int argc, char **argv) try
{
    // solve
    Prob1 p1;
    Numerical::Min<Prob1> sol(&p1, &Prob1::F, &Prob1::dF, /*ndim*/2);
    sol.Debug  = true;

    // with grads (ConjFR)
    sol.SetScheme ("ConjFR");
    cout << "\nWith grads (ConjFR)\n";
    double x[]  = {5.0, 7.0};
    double fmin = sol.Find (x);
    bool   err1 = (fabs(x[0]-1.0)>1.e-15 || fabs(x[1]-2.0)>1.e-16 || fabs(fmin-30.0)>1.e-16);
    cout << "xmin = (" << x[0] << ", " << x[1] << ")  =>  fmin = " << fmin << endl;
    cout << "err_x = " << fabs(x[0]-1.0) << ", err_y = " << fabs(x[1]-2.0) << ", err_fmin = " << fabs(fmin-30.0) << endl;

    // no grads (NMSimplex)
    sol.SetScheme ("NMSimplex");
    x[0] = 5.0;   x[1] = 7.0;
    cout << "\nNo grads (NMSimplex)\n";
    fmin = sol.Find (x, /*stepsize*/1.0, /*tol*/1.e-2);
    bool err2 = (fabs(x[0]-1.0)>1.e-2 || fabs(x[1]-2.0)>1.e-2 || fabs(fmin-30.0)>1.e-3);
    cout << "xmin = (" << x[0] << ", " << x[1] << ")  =>  fmin = " << fmin << endl;
    cout << "err_x = " << fabs(x[0]-1.0) << ", err_y = " << fabs(x[1]-2.0) << ", err_fmin = " << fabs(fmin-30.0) << endl;

    // with grads (BFGS)
    sol.SetScheme ("BFGS");
    x[0] = 5.0;   x[1] = 7.0;
    cout << "\nWith grads (BFGS)\n";
    fmin = sol.Find (x);
    bool err3 = (fabs(x[0]-1.0)>1.e-13 || fabs(x[1]-2.0)>1.0e-12 || fabs(fmin-30.0)>1.0e-16);
    cout << "xmin = (" << x[0] << ", " << x[1] << ")  =>  fmin = " << fmin << endl;
    cout << "err_x = " << fabs(x[0]-1.0) << ", err_y = " << fabs(x[1]-2.0) << ", err_fmin = " << fabs(fmin-30.0) << endl;

    // end
    if (err1) { cout << TERM_CLR_RED_H << "Error with ConjFR"    << endl; return 1; }
    if (err2) { cout << TERM_CLR_RED_H << "Error with NMSimplex" << endl; return 1; }
    if (err3) { cout << TERM_CLR_RED_H << "Error with BFGS"      << endl; return 1; }
    return 0;
}
MECHSYS_CATCH
