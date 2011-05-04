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

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    Vec3_t g;
};

void Setup(Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dom.Lat.Cells.Size();i++)
    {
        Cell * c = dom.Lat.Cells[i];
        c->BForcef = c->Density()*dat.g;
    }
    for (size_t i=0;i<dom.Particles.Size();i++)
    {
        dom.Particles[i]->Ff = dom.Particles[i]->M*dat.g;
    }
}


int main(int argc, char **argv) try
{
    size_t nx = 100;
    size_t ny = 100;
    double nu = 1.0/6.0;
    Domain Dom(D2Q9, nu, iVec3_t(nx,ny,1), /*dx*/1.0, /*dt*/1.0);
    UserData dat;
    Dom.UserData = &dat;
    Dom.Lat.G    = -200.0;
    Dom.Lat.Gs   = -200.0;
    dat.g        = 0.0,-0.001,0.0;

    //Initializing values
    //for (size_t i=0;i<Dom.Lat.Cells.Size();i++)
    //{
        //double rho0 = 200.0 + (1.0*rand())/RAND_MAX;
        //Vec3_t v0(0.0,0.0,0.0);
        //Dom.Lat.Cells[i]->Initialize(rho0, v0);
    //}
    
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

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        if (j<ny/2.0) Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize(1300.0,v0);
        else          Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize(  50.0,v0);

        //if (pow(i-39,2)+pow(j-50,2)<100)
        //{
            //Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize(1300.0,v0);
        //}
        //else
        //{
            //Dom.Lat.GetCell(iVec3_t(i,j,0))->Initialize(  50.0,v0);
        //}
    }

    Dom.AddDisk(0,Vec3_t(50.0,70.0,0.0),Vec3_t(0.0,-0.01,0.0),OrthoSys::O,4000.0,10.0,1.0);
    //Dom.Particles[Dom.Particles.Size()-1]->FixVelocity();

    //Solving
    Dom.Solve(10000.0,20.0,Setup,NULL,"test05");
}
MECHSYS_CATCH

