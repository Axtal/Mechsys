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


int main(int argc, char **argv) try
{
    Domain d;
    // Creating the Voronoi packing of particles
    d.AddVoronoiPacking(/*Tag*/-1,/*R*/0.1,/*Lx*/6,/*Ly*/6,/*Lz*/6,/*nx*/6,/*ny*/6,/*nz*/6,/*Periodic?*/true,/*rho*/1.0);
    //d.AddRice(-1,Vec3_t(0.0,0.0,0.0),2.0,0.1,1.0);
    d.GenBox(/*InitialTag*/-2,/*Lx*/8,/*Ly*/8,/*Lz*/8,/*R*/0.1);

    // First stage compression
    d.SetTriaxialTest(Vec3_t(/*sx*/0.1,/*sy*/0.1,/*sz*/0.1),Vec3_t(/*ex*/0.0,/*ey*/0.0,/*ez*/0.0));

    d.WriteBPY("test_triaxial01");

    d.Solve (/*tf*/10, 0.001, /*dtOut*/0.1, "test_triaxial01a", /*CamPos*/Vec3_t(0,35,0));

    //Second stage monotonic load

    d.SetTriaxialTest(Vec3_t(/*sx*/0.1,/*sy*/0.1,/*sz*/0.0),Vec3_t(/*ex*/0.0,/*ey*/0.0,/*ez*/0.005));

    d.Solve (/*tf*/40, 0.001, /*dtOut*/0.1, "test_triaxial01b", /*CamPos*/Vec3_t(0,35,0));
}
MECHSYS_CATCH
