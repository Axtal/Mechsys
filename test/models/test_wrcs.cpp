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
        cout << "         " << argv[0] << " WRC {pwfin}\n";
        cout << "in which\n";
        cout << "           WRC: BC\n";
        cout << "                HZ\n";
        cout << "                ZI\n";
        cout << "                PW\n";
        return 0;
    }
    String wrc   = argv[1];
    double pw0   =   0.0;
    double pwfin = -10.0;
    int    ndiv  = 200;
    if (argc>2) pwfin = atof(argv[2]);
    double porosity = 0.3;

    // model
    SDPair prms, inis;
    String keys("gamW kwsat akw por ");
    keys.append(wrc);
    prms.Set (keys.CStr(), 10.0, 1.0, 2.2, porosity, TRUE);
    inis.Set ("pw", pw0);
    UnsatFlow      mdl(/*ndim*/3, prms, /*anothermdl*/NULL);
    UnsatFlowState sta(/*ndim*/3);
    mdl.InitIvs (inis, &sta);

    // generate data
    cout << "\nModel: " << mdl.Name << endl;
    //Array<double> pcs(-pw0, 10.0, -pw0);
    Array<double> pcs(-pw0, 2.6, 1.0, 10.0, -pw0);
    mdl.GenCurve (pcs, "test_wrcs", ndiv);

    // end
    return 0;
}
MECHSYS_CATCH
