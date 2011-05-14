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
#include <fstream>

// MechSys
#include<mechsys/linalg/matvec.h>
#include<mechsys/models/unsatflow.h>
#include<mechsys/numerical/odesolver.h>

using std::cout;
using std::endl;
using Util::_6_3;
using Util::_8s;
using Util::_10_6;
using Util::SQ2;
using Util::SQ3;
using Util::SQ6;
using Util::PI;
const double TRUE  = 1.0;
const double FALSE = 0.0;

int main(int argc, char **argv) try
{
    // input
    if (argc<2)
    {
        cout << "Usage:\n";
        cout << "         " << argv[0] << " WRC {pcf}\n";
        cout << "in which\n";
        cout << "           WRC: 0 = BC\n";
        cout << "                1 = HZ\n";
        cout << "                2 = ZI\n";
        return 0;
    }
    int    wrc  = atoi(argv[1]);
    double pc   = 0.0;
    double Sw   = 1.0;
    double pcf  = 100.0;
    int    ndiv = 100;
    if (argc>2) pcf = atof(argv[2]);

    SDPair prms, inis;
    prms.Set ("gamW kwsat akw WRC", 10.0, 1.0, 2.2, (double)wrc);
    inis.Set ("pw Sw n", -pc, Sw, 0.3);
    UnsatFlow      mdl(/*ndim*/3, prms, /*anothermdl*/NULL);
    UnsatFlowState sta(/*ndim*/3);
    mdl.InitIvs (inis, &sta);

    std::ofstream of("test_wrcs.res",std::ios::out);
    of << _8s<<"pc"    << _8s<<"Sw"    << _8s<<"kw"                            << endl;
    of << _8s<< sta.pc << _8s<< sta.Sw << _8s<< mdl.rw(sta.Sw)*mdl.kwsatb*10.0 << endl;

    double dpc = (pcf-pc)/ndiv;
    while (sta.pc<pcf)
    {
        mdl.Update (-dpc, /*dev*/0.0, &sta);
        of << _8s<< sta.pc << _8s<< sta.Sw << _8s<< mdl.rw(sta.Sw)*mdl.kwsatb*10.0 << endl;
    }

    of.close();
    cout << "file <" << TERM_CLR_BLUE_H << "test_wrcs.res" << TERM_RST << "> written" << endl;

    // end
    return 0;
}
MECHSYS_CATCH
