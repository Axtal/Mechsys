/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2014 Sergio Galindo                                    *
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

/////////////////////// Test 02 The dielectric medium

// MechSys
#include <mechsys/emlbm2/Domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
    size_t nproc = 1; 
    if (argc==2) nproc=atoi(argv[1]);
    int nx = 10;
    int ny = 10;
    int nz = 200;
    EMLBM::Domain Dom(iVec3_t(nx,ny,nz), 1.0, 1.0);
    double E0 = 0.001;
    double B0 = sqrt(2.0)*E0;
    double alpha = 0.01;
    double z0 = 40;

    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Vec3_t E(-E0*exp(-alpha*(k-z0)*(k-z0)),0.0,0.0);
        Vec3_t B(0.0,B0*exp(-alpha*(k-z0)*(k-z0)),0.0);
        Dom.Lat.GetCell(iVec3_t(i,j,k))->Initialize(0.0,OrthoSys::O,E,B);
        Dom.Lat.GetCell(iVec3_t(i,j,k))->Eps = 0.75*tanh(k-nz/2)+1.75;
    }
    //Dom.WriteXDMF("test");
    Dom.Solve(200.0,2.0,NULL,NULL,"temlbm02",true,nproc);

    return 0;
}
MECHSYS_CATCH

