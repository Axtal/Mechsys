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
    Vec3_t Xmax( 20.5, 20.5,0.0);
    Vec3_t Xmin(-20.5,-20.5,0.0);
    iVec3_t Div(200,200,1);
    SPHDomain dom(Div,Xmin,Xmax);
    dom.CamPos = 0.0,60.0,0.0;
    dom.Gravity = 0.0,-0.05,0.0;
    
    dom.AddBox(Vec3_t (-20.5, 0.0, 0.0), 1,40,1,0.5,0.5,20.0,true);
    dom.AddBox(Vec3_t ( 20.5, 0.0, 0.0), 1,40,1,0.5,0.5,20.0,true);
    dom.AddBox(Vec3_t ( 0.0, -20.5,0.0), 42,1,1,0.5,0.5,20.0,true);
    dom.AddRandomBox(Vec3_t(-10.0,0.0,0.0),20,40,0,20,40,1,10.0,0.5);

    dom.Solve(/*tf*/400.0,/*dt*/0.001,/*dtOut*/4,"test01");
    return 0;
}
MECHSYS_CATCH
