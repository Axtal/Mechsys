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

// FLow in a pipe with obstacle

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    Array<Cell *> Left;
    Array<Cell *> Right;
    double        vmax;
    double        rho;
};

void Setup (Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //for (size_t i=0;i<dat.Left.Size();i++)
    //{
        //dat.Left [i]->Initialize(dat.Left [i]->RhoBC,dat.Left [i]->VelBC);
        //dat.Right[i]->Initialize(dat.Right[i]->RhoBC,dat.Right[i]->VelBC);
    //}

	// Cells with prescribed velocity
	for (size_t i=0; i<dat.Left.Size(); ++i)
	{
		Cell * c = dat.Left[i];
		if (c->IsSolid) continue;
		double rho = (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[3]+c->F[6]+c->F[7]))/(1.0-c->VelBC(0));
		c->F[1] = c->F[3] + (2.0/3.0)*rho*c->VelBC(0);
		c->F[5] = c->F[7] + (1.0/6.0)*rho*c->VelBC(0) + 0.5*rho*c->VelBC[1] - 0.5*(c->F[2]-c->F[4]);
		c->F[8] = c->F[6] + (1.0/6.0)*rho*c->VelBC(0) - 0.5*rho*c->VelBC[1] + 0.5*(c->F[2]-c->F[4]);
	}

	// Cells with prescribed density
	for (size_t i=0; i<dat.Right.Size(); ++i)
	{
		Cell * c = dat.Right[i];
		if (c->IsSolid) continue;
		double vx = -1.0 + (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[1]+c->F[5]+c->F[8]))/c->RhoBC;
		c->F[3] = c->F[1] - (2.0/3.0)*c->RhoBC*vx; 
		c->F[7] = c->F[5] - (1.0/6.0)*c->RhoBC*vx + 0.5*(c->F[2]-c->F[4]);
		c->F[6] = c->F[8] - (1.0/6.0)*c->RhoBC*vx - 0.5*(c->F[2]-c->F[4]);
	}
}

int main(int argc, char **argv) try
{
    double u_max  = 0.1;                 // Poiseuille's maximum velocity
    double Re     = 100;                 // Reynold's number
    size_t nx = 400;
    size_t ny = 100;
    int radius = ny/10 + 1;           // radius of inner circle (obstacle)
    double nu     = u_max*(2*radius)/Re; // viscocity
    Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;

    dat.vmax = u_max;
    //Assigning the left and right cells
    for (size_t i=0;i<ny;i++)
    {
        dat.Left .Push(Dom.Lat[0].GetCell(iVec3_t(0   ,i,0)));
        dat.Right.Push(Dom.Lat[0].GetCell(iVec3_t(nx-1,i,0)));
        
        // set parabolic profile
        double L  = ny - 2;                       // channel width in cell units
        double yp = i - 1.5;                      // ordinate of cell
        double vx = dat.vmax*4/(L*L)*(L*yp - yp*yp); // horizontal velocity
        double vy = 0.0;                          // vertical velocity
		Vec3_t v(vx, vy, 0.0);                    // velocity vector
        dat.Left [i]->VelBC = v;
        dat.Left [i]->RhoBC = 1.0;
        //dat.Right[i]->VelBC = 0.08,0.0,0.0;
        dat.Right[i]->VelBC = v;
        dat.Right[i]->RhoBC = 1.0;

    }
    dat.rho  = 1.0;

	// set inner obstacle
	int obsX   = ny/2;   // x position
	int obsY   = ny/2+3; // y position
    Dom.AddDisk(0,Vec3_t(obsX,obsY,0),Vec3_t(0.0,0.0,0.0),Vec3_t(0.0,0.0,0.0),   3.0,radius,1.0);
    Dom.AddDisk(0,Vec3_t(nx/2.0,obsY,0),Vec3_t(0.0,0.0,0.0),Vec3_t(0.0,0.0,0.0),30.0,2*radius,1.0);



    //Assigning solid boundaries at top and bottom
    for (size_t i=0;i<nx;i++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,0))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,0))->IsSolid = true;
    }

    double rho0 = 1.0;
    Vec3_t v0(0.08,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Cells.Size();i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
    }

    //Solving
    Dom.Time = 0.0;
    Dom.Solve(10000.0,20.0,Setup,NULL,"test01");
 
}
MECHSYS_CATCH

