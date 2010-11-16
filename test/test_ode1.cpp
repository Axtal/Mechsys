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
#include <mechsys/numerical/odesolver.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;
using Util::_8s;

/* Solving:
 *             x''(t) + \mu x'(t) (x(t)^2 - 1) + x(t) = 0
 * Using:
 *             x' = y
 *             y' = -x + \mu y (1-x^2)
 */

class ODE
{
public:
    ODE (double Mu) : mu(Mu) {}
    int Fun (double t, double const Y[], double dYdt[])
    {
        dYdt[0] =  Y[1];
        dYdt[1] = -Y[0] - mu*Y[1]*(Y[0]*Y[0] - 1);
        return GSL_SUCCESS;
    }
    double mu;
};

int main(int argc, char **argv) try
{
    // data
    double mu    = 2.0;
    double y0ini = 0.0;
    double y1ini = 3.0;
    double tf    = 10.209022;
    double dtOut = 0.1;

    // output file
    std::ofstream of("test_ode1.dat", std::ios::out);
    of << _8s<<"t" << _8s<<"y0"  << _8s<<"y1"  << endl;
    of << _8s<<0.0 << _8s<<y0ini << _8s<<y1ini << endl;

    // solve
    ODE ode(mu);
    Numerical::ODESolver<ODE> sol(&ode, &ODE::Fun, /*neq*/2);
    sol.t    = 0.0;
    sol.Y[0] = y0ini;
    sol.Y[1] = y1ini;
    double tout = 0.0 + dtOut;
    while (sol.t<tf)
    {
        sol.Evolve (tout);
        of << _8s<<sol.t << _8s<<sol.Y[0] << _8s<<sol.Y[1] << endl;
        tout += dtOut;
    }

    // close file
    of.close();
    cout << "\nFile <" << TERM_CLR_BLUE_H << "test_ode1.dat" << TERM_RST << "> written" << endl;
    cout << endl;

    // end
    return 0;
}
MECHSYS_CATCH
