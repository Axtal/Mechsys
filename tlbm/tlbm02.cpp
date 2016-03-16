/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
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
// Vorticity test

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

void Setup (LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dat.Left.Size();i++)
    {
        dat.Left [i]->Initialize(dat.Left [i]->RhoBC,dat.Left [i]->VelBC);
        dat.Right[i]->Initialize(dat.Right[i]->RhoBC,dat.Right[i]->VelBC);
    }

	// Cells with prescribed velocity
	//for (size_t i=0; i<dat.Left.Size(); ++i)
	//{
		//Cell * c = dat.Left[i];
		//if (c->IsSolid) continue;
		//double rho = (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[3]+c->F[6]+c->F[7]))/(1.0-c->VelBC(0));
		//c->F[1] = c->F[3] + (2.0/3.0)*rho*c->VelBC(0);
		//c->F[5] = c->F[7] + (1.0/6.0)*rho*c->VelBC(0) + 0.5*rho*c->VelBC[1] - 0.5*(c->F[2]-c->F[4]);
		//c->F[8] = c->F[6] + (1.0/6.0)*rho*c->VelBC(0) - 0.5*rho*c->VelBC[1] + 0.5*(c->F[2]-c->F[4]);
	//}

	// Cells with prescribed density
	//for (size_t i=0; i<dat.Right.Size(); ++i)
	//{
		//Cell * c = dat.Right[i];
		//if (c->IsSolid) continue;
		//double vx = -1.0 + (c->F[0]+c->F[2]+c->F[4] + 2.0*(c->F[1]+c->F[5]+c->F[8]))/c->RhoBC;
		//c->F[3] = c->F[1] - (2.0/3.0)*c->RhoBC*vx; 
		//c->F[7] = c->F[5] - (1.0/6.0)*c->RhoBC*vx + 0.5*(c->F[2]-c->F[4]);
		//c->F[6] = c->F[8] - (1.0/6.0)*c->RhoBC*vx - 0.5*(c->F[2]-c->F[4]);
	//}
}

int main(int argc, char **argv) try
{
    double u_max  = 0.1;                 // Poiseuille's maximum velocity
    double Re     = 100;                 // Reynold's number
    size_t nx = 200;
    size_t ny = 50;
    size_t nz = 50;
    int radius = ny/10 + 1;           // radius of inner circle (obstacle)
    double nu     = u_max*(2*radius)/Re; // viscocity
    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), 1.0, 1.0);
    UserData dat;
    Dom.UserData = &dat;

    dat.vmax = u_max;
    //Assigning the left and right cells
    size_t n = 0;
    for (size_t i=0;i<ny;i++)
    for (size_t j=0;j<nz;j++)
    {
        dat.Left .Push(Dom.Lat[0].GetCell(iVec3_t(0   ,i,j)));
        dat.Right.Push(Dom.Lat[0].GetCell(iVec3_t(nx-1,i,j)));
        
        // set parabolic profile
        double L  = ny - 2;                       // channel width in cell units
        double yp = i - 1.5;                      // coordinate of cell
        double zp = j - 1.5;                      // coordinate of cell
        double vx = dat.vmax*4/(L*L*L*L)*(L*yp - yp*yp)*(L*zp - zp*zp); // horizontal velocity
        double vy = 0.0;                          // vertical velocity
		Vec3_t v(vx, vy, 0.0);                    // velocity vector
        dat.Left [n]->VelBC = v;
        dat.Left [n]->RhoBC = 1.0;
        dat.Right[n]->VelBC = 0.08,0.0,0.0;
        //dat.Right[i]->VelBC = v;
        dat.Right[n]->RhoBC = 1.0;
        n++;
    }
    dat.rho  = 1.0;

	// set inner obstacle
	int obsX   = ny/2;   // x position
	int obsY   = ny/2; // y position
	int obsZ   = ny/2; // z position
    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        if (pow(i-obsX,2.0)+pow(j-obsY,2.0)+pow(k-obsZ,2.0)<radius*radius) Dom.Lat[0].GetCell(iVec3_t(i,j,k))->IsSolid = true;
    }


    //Assigning solid boundaries at top and bottom
    for (size_t i=0;i<nx;i++)
    for (size_t j=0;j<ny;j++)
    {
        Dom.Lat[0].GetCell(iVec3_t(i,0   ,j))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,ny-1,j))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,j,0   ))->IsSolid = true;
        Dom.Lat[0].GetCell(iVec3_t(i,j,ny-1))->IsSolid = true;
    }

    double rho0 = 1.0;
    Vec3_t v0(0.08,0.0,0.0);

    //Initializing values
    for (size_t i=0;i<Dom.Lat[0].Ncells;i++)
    {
        Dom.Lat[0].Cells[i]->Initialize(rho0, v0);
    }

    //Solving
    Dom.Solve(10000.0,20.0,Setup,NULL,"test02");
 
}
MECHSYS_CATCH

