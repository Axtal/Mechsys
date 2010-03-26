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

// Std lib
#include <math.h>

// MechSys
#include <mechsys/sph/domain.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    SPHDomain dom;
    dom.AddBox(Vec3_t(5.5,0.0,5.0),1,10,10,0.5,1.0,true);
    dom.AddBox(Vec3_t(-5.5,0.0,5.0),1,10,10,0.5,1.0,true);
    //dom.AddBox(Vec3_t(0.0,5.5,5.0),10,1,10,0.5,1.0,true);
    //dom.AddBox(Vec3_t(0.0,-5.5,5.0),10,1,10,0.5,1.0,true);
    //dom.AddBox(Vec3_t(0.0,0.0,-0.5),10,10,1,0.5,1.0,true);
    //dom.AddBox(Vec3_t(-3.5,0.0,5.0),3,1,10,0.5,1.0,false);
    dom.WriteBPY("test01");

    double dt = 0.001;

}
MECHSYS_CATCH
