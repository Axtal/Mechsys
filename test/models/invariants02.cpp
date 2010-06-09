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
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;
using Util::SQ2;
using Util::SQ6;
using Util::PI;

int main(int argc, char **argv) try
{
    int tst = 2;
    if (argc>1) tst = atoi(argv[1]);

    Vec_t eps(4);
    switch (tst)
    {
        case 1: { eps = -0.1, -0.1, -0.1, 0.0*SQ2;  break; }
        case 2: { eps = -0.1, -0.1, -0.2, 0.0*SQ2;  break; }
        case 3: { eps = -0.2, -0.1, -0.1, 0.0*SQ2;  break; }
        default: throw new Fatal("main: Test = %d is not available",tst);
    }

    double ev  = Calc_ev    (eps);
    double ed  = Calc_ed    (eps);
    double evo = Calc_evoct (eps);
    double edo = Calc_edoct (eps);

    Vec3_t L;
    pqth2L (-evo, edo, -150.0*PI/180.0, L, "oct");
    double ex = L(2);
    double ey = L(1);
    double ez = L(0);

    String buf;
    buf.Printf("eps                = [%g, %g, %g, %g]  \n",eps(0),eps(1),eps(2),eps(3)); cout<<buf;
    buf.Printf("ev,    ed          = %g, %g            \n",ev,  ed);                     cout<<buf;
    buf.Printf("evoct, edoct       = %g, %g            \n",evo, edo);                    cout<<buf;
    buf.Printf("ex,  ey,  ez       = %g, %g, %g        \n",ex, ey, ez);                  cout<<buf;


    return 0;
}
MECHSYS_CATCH
