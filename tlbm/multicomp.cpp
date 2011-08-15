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
#include <stdlib.h>

// MechSys
#include <mechsys/lbm/Domain.h>

using std::cout;
using std::endl;
struct UserData
{
    Vec3_t             g;
};

void Setup(Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t j=0;j<dom.Lat.Size();j++)
    for (size_t i=0;i<dom.Lat[j].Cells.Size();i++)
    {
        Cell * c = dom.Lat[j].Cells[i];
        c->BForcef = c->Density()*dat.g;
    }

}


int main(int argc, char **argv) try
{
    Array<double> nu(2);
    nu[0] = 1.0/6.0;
    nu[1] = 1.0/6.0;

    size_t nx = 50, ny = 50;

    // Setting top and bottom wall as solid
    Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;
    dat.g           = 0.0,-0.001,0.0;
    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[1].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }
    //for (size_t i=0;i<ny;i++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(0   ,i,0))->IsSolid = true;
        //Dom.Lat[1].GetCell(iVec3_t(nx-1,i,0))->IsSolid = true;
    //}

    // Set inner drop
    int obsX = nx/2, obsY = ny/2;
    int radius =  5;

	for (size_t i=0; i<nx; ++i)
	for (size_t j=0; j<ny; ++j)
    {
		Vec3_t V;  V = 0.0, 0.0, 0.0;
		if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
		{
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(1300.0,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(100.0,V);
		}
		else
		{
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(100.0,V);
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(1300.0,V);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
		}
    }

    // Set parameters
    Dom.Lat[0].G = -200.0;
    Dom.Lat[0].Gs= -400.0;
    //Dom.Lat[0].G = -0.0;
    //Dom.Lat[0].Gs= -0.0;
    Dom.Lat[1].G =  0.0;
    Dom.Lat[1].Gs=  0.0;
    Dom.Gmix     =  0.001;

    Dom.Solve(5000,50.0,Setup,NULL,"multicomp");


    return 0;
}
MECHSYS_CATCH
