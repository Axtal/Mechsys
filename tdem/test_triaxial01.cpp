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
#include "dem/domain.h"
#include "util/fatal.h"

using std::cout;
using std::endl;
using DEM::Domain;

int main(int argc, char **argv) try
{

    // Setting the stress path with the invariants p,q,alpha for the p cte paths
    double p = -0.1;
    double q = 0.3;
    double alpha = 15.0*M_PI/180.0;
    
    // the stress rates are functions of these invariants
    double srx = q*sqrt(2.0/3.0)*cos(alpha)/40;
    double sry = q*sqrt(1.0/2.0)*sin(alpha)/40-q*sqrt(1.0/6.0)*cos(alpha)/40;
    double srz = -q*sqrt(1.0/2.0)*sin(alpha)/40-q*sqrt(1.0/6.0)*cos(alpha)/40;

    Domain d;
    // Creating the Voronoi packing of particles
    d.AddVoroPack(/*Tag*/-1,/*R*/0.1,/*Lx*/4,/*Ly*/4,/*Lz*/4,/*nx*/4,/*ny*/4,/*nz*/4,/*Periodic?*/true,/*rho*/1.0);
    //d.AddRice(-1,Vec3_t(0.0,0.0,0.0),2.0,0.1,1.0);
    d.GenBox(/*InitialTag*/-2,/*Lx*/6,/*Ly*/6,/*Lz*/6,/*R*/0.1);

    // First stage compression
    d.SetTxTest(Vec3_t(/*ex*/0.0,/*ey*/0.0,/*ez*/0.0),Vec3_t(/*sx*/-0.1,/*sy*/-0.1,/*sz*/-0.1));

    d.WriteBPY("test_triaxial01");

    d.Solve (/*tf*/10, 0.001, /*dtOut*/0.1, "test_triaxial01a", /*CamPos*/Vec3_t(0,35,0));

    //Second stage Stress Path

    d.SetTxTest(Vec3_t(/*ex*/0.0,/*ey*/0.0,/*ez*/0.0),Vec3_t(/*sx*/-0.1,/*sy*/-0.1,/*sz*/-0.1),Vec3_t(/*srx*/srx,/*sry*/sry,/*srz*/srz));

    d.Solve (/*tf*/40, 0.001, /*dtOut*/0.1, "test_triaxial01b", /*CamPos*/Vec3_t(0,35,0));
}
MECHSYS_CATCH
