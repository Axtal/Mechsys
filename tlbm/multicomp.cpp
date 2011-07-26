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
//#include <mechsys/lbm/mixture.h>

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
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(1300.0,V);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(100.0,V);
		}
		else
		{
            //Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
            //Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(100.0,V);
            Dom.Lat[0].GetCell(iVec3_t(i,j,0))->Initialize(1300.0,V);
            Dom.Lat[1].GetCell(iVec3_t(i,j,0))->Initialize(0.1,V);
		}
    }

    // Set parameters
    Dom.Lat[0].G = -200.0;
    Dom.Lat[0].Gs= -400.0;
    Dom.Lat[1].G =  0.0;
    Dom.Lat[1].Gs=  0.0;
    Dom.Gmix     =  0.001;

    Dom.Solve(5000,50.0,Setup,NULL,"multicomp");


    return 0;
	//Now you have to enter the true viscosity as an array for each component
	//double nu[2];
	//nu[0]=1.0;
	//nu[1]=3.0;

	// Allocate lattice
	//LBM::Mixture m( /*FileKey*/ "multi", /*Is3D*/false, /*NComp*/2,nu/*true viscosity*/,/*Nx*/50, /*Ny*/50,1./*Nz*/,1./*space step*/,1./*time step*/);

	// Set walls (top and bottom)
	// Lattice 0
	//for (size_t i=0; i<m.GetLattice(0)->Top()   .Size(); ++i) m.GetLattice(0)->Top()[i]->SetSolid();
	//for (size_t i=0; i<m.GetLattice(0)->Bottom().Size(); ++i) m.GetLattice(0)->Bottom()[i]->SetSolid();
	// Lattice 1
	//for (size_t i=0; i<m.GetLattice(1)->Top()   .Size(); ++i) m.GetLattice(1)->Top()[i]->SetSolid();
	//for (size_t i=0; i<m.GetLattice(1)->Bottom().Size(); ++i) m.GetLattice(1)->Bottom()[i]->SetSolid();

	//// Initialize
	//for (size_t i=0; i<m.Nx(); ++i)
	//for (size_t j=0; j<m.Ny(); ++j)
	//{
	//	//double rho0 = 1.4 +(.02*rand())/RAND_MAX;
		//Vec3_t V; V = 0.0, 0.0, 0.0;
		//m.GetLattice(0)->GetCell(i,j)->Initialize (0.1, V);
		//m.GetLattice(1)->GetCell(i,j)->Initialize (1.0, V);
	//}

	// lattice[0] : Water
	// lattice[1] : Oil

	// Properties

	// Set inner drop
	//int obsX = m.Nx()/2;   // x position
	//int obsY = m.Ny()/2;   // y position
	//int radius = m.Nx()/5; // Inital drop radius
//
	//for (size_t i=0; i<m.Nx(); ++i)
	//for (size_t j=0; j<m.Ny(); ++j)
	//{
		//Vec3_t V;  V = 0.0, 0.0, 0.0;
		//if (pow((int)(i)-obsX,2.0) + pow((int)(j)-obsY,2.0) <= pow(radius,2.0)) // circle equation
		//{
			//m.GetLattice(0)->GetCell(i,j)->Initialize (0.01, V,m.GetLattice(0)->Cs());
			//m.GetLattice(1)->GetCell(i,j)->Initialize (1.0, V,m.GetLattice(1)->Cs());
		//}
		//else
		//{
			//m.GetLattice(0)->GetCell(i,j)->Initialize (1.0, V,m.GetLattice(0)->Cs());
			//m.GetLattice(1)->GetCell(i,j)->Initialize (0.01, V,m.GetLattice(1)->Cs());
		//}
	//}

	//// Properties
	//m.GetLattice(0)->SetG(-3.0)->SetGSolid(-0.0);
	//m.GetLattice(1)->SetG(-1.0)->SetGSolid(-0.0);
	
	//m.SetGravity(0.0, -0.0005, 0.0);
	//m.SetMixG(0.5);

	// Solve
	//m.Solve(/*tIni*/0.0, /*tFin*/5000.0, /*dtOut*/25.0);
}
MECHSYS_CATCH
