/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;
using std::ofstream;
using DEM::Domain;

int main(int argc, char **argv) try
{
    // set the simulation domain ////////////////////////////////////////////////////////////////////////////
    
    Domain d;
    d.CamPos = Vec3_t(0, 0.0, 2.0*16.0); // position of camera

    // particle
    d.AddVoroPack (-1, 0.1, 16,16,1, 16,16,1, 1.0, false, 1000, 1.0, "cylinder");
    d.FreeParticles = d.Particles;
    d.WritePOV("cylinder");
    d.WriteBPY("cylinder");

    Array<double> X,Y,D;
    d.GetGSD(X,Y,D);

    ofstream fg("granulometry.txt");
    for (size_t i = 0;i < X.Size() ; i++ )
    {
        fg << Util::_10_6 << X[i] << Util::_8s << Y[i] << std::endl;
    }
    fg.close();

    return 0;
}
MECHSYS_CATCH
