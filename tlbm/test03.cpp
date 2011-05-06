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
// Magnus effect

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
};

void Report (Domain & Dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //std::cout << Dom.Particles[0]->X(0) << " " << Dom.Particles[0]->X(1) << " "
              //<< Dom.Particles[1]->X(0) << " " << Dom.Particles[1]->X(1) << std::endl;
}


int main(int argc, char **argv) try
{
    size_t nx = 100;
    size_t ny = 100;
    double nu = 0.01;
    Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), /*dx*/1.0, /*dt*/1.0);
    Dom.AddDisk(0,Vec3_t(15.0,20.0,0.0),Vec3_t( 0.01,0.0,0.0),Vec3_t(0.0,0.0, 0.01),3.0,10.0,1.0);
    Dom.AddDisk(0,Vec3_t(85.0,20.0,0.0),Vec3_t(-0.01,0.01,0.0),Vec3_t(0.0,0.0,-0.01),3.0,10.0,1.0);
    UserData dat;
    Dom.UserData = &dat;
    double rho0 = 0.001;
    Vec3_t v0(0.0,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat.Cells.Size();i++)
    {
        Dom.Lat.Cells[i]->Initialize(rho0, v0);
    }

     //Set solid boundaries
    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat.GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat.GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }
    for (size_t i=0;i<ny;i++)
    {
        Dom.Lat.GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        Dom.Lat.GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    }

    //Solving
    Dom.Solve(10000.0,20.0,NULL,Report,"test03");
}
MECHSYS_CATCH

