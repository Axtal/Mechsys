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
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// MechSys
#include <mechsys/numerical/numdiff.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>

using std::cout;
using std::endl;
using Util::_8s;
using Util::PI;

class Problem
{
public:
    double yFun (double x) const
    {
        return pow(x,1.5);
    }
};

int main(int argc, char **argv) try
{
    Problem                  prob;
    Numerical::Diff<Problem> nd(&prob);

    double result = nd.DyDx (&Problem::yFun, /*At*/2.0);
    printf ("f(x) = x^(3/2)\n");
    printf ("x = 2.0\n");
    printf ("f'(x) = %.10f +/- %.10f\n", result, nd.LastAbsErr);
    printf ("exact = %.10f\n\n", 1.5 * sqrt(2.0));

    nd.Method = 1;
    result = nd.DyDx (&Problem::yFun, /*At*/0.0);
    printf ("x = 0.0\n");
    printf ("f'(x) = %.10f +/- %.10f\n", result, nd.LastAbsErr);
    printf ("exact = %.10f\n\n", 0.);

    // end
    return 0;
}
MECHSYS_CATCH
