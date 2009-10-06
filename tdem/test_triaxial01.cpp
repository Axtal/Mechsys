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
    d.AddVoronoiPacking(-1,0.1,4,4,5,4,4,5,1.0);
    d.GenBox(-2,8.0,8.0,10.0,0.1);

    // First stage compression
    Dict A;
    A.Set(-2,"fx fy fz",-10.0,0.0,0.0)
     .Set(-3,"fx fy fz",10.0,0.0,0.0)
     .Set(-4,"fx fy fz",0.0,-10.0,0.0)
     .Set(-5,"fx fy fz",0.0,10.0,0.0)
     .Set(-6,"fx fy fz",0.0,0.0,-10.0)
     .Set(-7,"fx fy fz",0.0,0.0,10.0);

    d.SetProps(A);

    d.WriteBPY("test_triaxial01");

    d.Solve (/*tf*/30, 0.001, /*dtOut*/0.1, "test_triaxial01a", /*CamPos*/Vec3_t(0,25,0));

    //Second stage monotonic load
    Dict B;
    B.Set(-2,"fx fy fz",-10.0,0.0,0.0)
     .Set(-3,"fx fy fz",10.0,0.0,0.0)
     .Set(-4,"fx fy fz",0.0,-10.0,0.0)
     .Set(-5,"fx fy fz",0.0,10.0,0.0)
     .Set(-6,"vx vy vz",0.0,0.0,-0.05)
     .Set(-7,"vx vy vz",0.0,0.0,0.05);

    d.SetProps(B);

    d.Solve (/*tf*/30, 0.001, /*dtOut*/0.1, "test_triaxial01b", /*CamPos*/Vec3_t(0,25,0));
}
MECHSYS_CATCH
