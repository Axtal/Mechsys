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
#include <mechsys/fem/solver.h>

using std::cout;
using std::endl;
using Util::_8s;

int main(int argc, char **argv) try
{
    int    a   = 7;
    int    b   = 12;
    double L   = 100.0;
    double m   = 2.0;
    if (argc>1) a = atoi(argv[1]);
    if (argc>2) b = atoi(argv[2]);
    if (argc>3) L = atof(argv[3]);
    if (argc>4) m = atof(argv[4]);

    std::ofstream of("test_timesteps.res", std::ios::out);
    of << _8s<< "i" << _8s<< "dt_sch0" << _8s<<"dt_sch1" << endl;
    for (int i=0; i<b; ++i)
        of << _8s<< i << _8s<< FEM::Solver::Timestep(i,a,L,/*sch*/0,m)
                      << _8s<< FEM::Solver::Timestep(i,a,L,/*sch*/1,m) << endl;
    of.close();

    cout << "\nFile <" << TERM_CLR_BLUE_H << "test_timesteps.res" << TERM_RST << "> written" << endl;
    cout << endl;
    return 0;
}
MECHSYS_CATCH
