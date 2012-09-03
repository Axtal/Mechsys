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
// Stokes law

//STD
#include<iostream>

// MechSys
#include <mechsys/lbm/Domain.h>

struct UserData
{
    std::ofstream      oss_ss;       ///< file for particle data
    Vec3_t                acc;
    Array<Cell *>        xmin;
    Array<Cell *>        xmax;
};

void Setup(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    for (size_t i=0;i<dom.Lat[0].Cells.Size();i++)
    {
        Cell * c   = dom.Lat[0].Cells[i];
        c->BForcef = c->Rho*dat.acc;
    }
    //for (size_t i=0;i<dat.xmin.Size();i++)
    //{
        //Cell * c = dat.xmin[i];
        //if(c->IsSolid) continue;
        //c->F[1] = 1.0/3.0 *(-2*c->F[0] - 4*c->F[10] - 4*c->F[12]-4*c->F[14]-c->F[2]-2*c->F[3]-2*c->F[4]-2*c->F[5]-2*c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->F[7] = 1.0/24.0*(-2*c->F[0] - 4*c->F[10] - 4*c->F[12]-4*c->F[14]-4*c->F[2]  +c->F[3]-5*c->F[4]  +c->F[5]-5*c->F[6]+20*c->F[8]+2*c->RhoBC);
        //c->F[9] = 1.0/24.0*(-2*c->F[0] + 20*c->F[10] - 4*c->F[12]-4*c->F[14]-4*c->F[2]+c->F[3]-5*c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->F[11]= 1.0/24.0*(-2*c->F[0] - 4*c->F[10] + 20*c->F[12]-4*c->F[14]-4*c->F[2]-5*c->F[3]+c->F[4]  +c->F[5]-5*c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->F[13]= 1.0/24.0*(-2*c->F[0] - 4*c->F[10] - 4 *c->F[12]+20*c->F[14]-4*c->F[2]-5*c->F[3]+  c->F[4]-5*c->F[5]+c->F[6]-4*c->F[8]+2*c->RhoBC);
        //c->Rho = c->VelDen(c->Vel);
    //}
    //for (size_t i=0;i<dat.xmax.Size();i++)
    //{
        //Cell * c = dat.xmax[i];
        //if(c->IsSolid) continue;
        //c->F[2] = 1/3.0* (-2*c->F[0]-c->F[1]-2*(2*c->F[11]+2*c->F[13]+c->F[3]+c->F[4]+c->F[5]+c->F[6]+2*c->F[7]+2*c->F[9]-c->RhoBC));
        //c->F[8] = 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] - 5*c->F[5] + c->F[6] +20*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        //c->F[10]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] - 4*c->F[13] - 5*c->F[3] + c->F[4] + c->F[5] - 5*c->F[6] - 4*c->F[7] + 20*c->F[9] + 2*c->RhoBC) ;
        //c->F[12]= 1/24.0*(-2*c->F[0] - 4*c->F[1] + 20*c->F[11] - 4*c->F[13] + c->F[3] - 5*c->F[4] - 5*c->F[5] + c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        //c->F[14]= 1/24.0*(-2*c->F[0] - 4*c->F[1] - 4*c->F[11] + 20*c->F[13] + c->F[3] - 5*c->F[4] + c->F[5] - 5*c->F[6] -  4*c->F[7] - 4*c->F[9] + 2*c->RhoBC);
        //c->Rho = c->VelDen(c->Vel);
    //}
}


void Report(LBM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("%s_force.res",dom.FileKey.CStr());
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "Fx" << Util::_8s << "V" << Util::_8s << "Rho \n";
    }
    if (!dom.Finished) 
    {
        Vec3_t Flux = OrthoSys::O;
        double M    = 0.0;
        size_t nc   = 0;
        for (size_t i=0;i<dom.Lat[0].Cells.Size();i++)
        {
            Cell * c = dom.Lat[0].Cells[i];
            if (c->IsSolid||c->Gamma>1.0e-8) continue;
            Flux += c->Rho*c->Vel;
            M += c->Rho;
            nc++;
        }
        Flux/=M;
        dat.oss_ss << Util::_10_6 << dom.Time << Util::_8s << dom.Particles[0]->F(0) << Util::_8s << Flux(0) << Util::_8s << M/nc << std::endl;
    }
    else
    {
        dat.oss_ss.close();
    }
}

int main(int argc, char **argv) try
{
    size_t Nproc = 1; 
    size_t nx = 100;
    size_t ny = 50;
    size_t nz = 50;
    double nu = 0.01;
    double dx = 1.0;
    double dt = 1.0;
    double Dp = 0.1;
    double R  = 10.0;
    double Tf = 40000.0;
    if (argc>=2)
    {
        Nproc = atoi(argv[1]);
        nx    = atoi(argv[2]);
        //ny = nz = nx/2;
        ny = nz = nx;
        nu    = atof(argv[3]);
        Dp    = atof(argv[4]);
        R     = atof(argv[5]);
        Tf    = atof(argv[6]);
    }

    

    LBM::Domain Dom(D3Q15, nu, iVec3_t(nx,ny,nz), dx, dt);
    UserData dat;
    Dom.UserData = &dat;
    dat.acc      = Vec3_t(Dp,0.0,0.0);

    //for (size_t i=0;i<nx;i++)
    //for (size_t j=0;j<ny;j++)
    //{
        //Dom.Lat[0].GetCell(iVec3_t(i,0   ,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,ny-1,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,j,0   ))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(i,j,ny-1))->IsSolid = true;
    //}

    //for (size_t i=0;i<ny;i++)
    //for (size_t j=0;j<nz;j++)
    //{
        //dat.xmin.Push(Dom.Lat[0].GetCell(iVec3_t(0   ,i,j)));
        //dat.xmax.Push(Dom.Lat[0].GetCell(iVec3_t(nx-1,i,j)));
        //Dom.Lat[0].GetCell(iVec3_t(0   ,i,j))->IsSolid = true;
        //Dom.Lat[0].GetCell(iVec3_t(nx-1,i,j))->IsSolid = true;
    //}


    for (int i=0;i<nx;i++)
    for (int j=0;j<ny;j++)
    for (int k=0;k<nz;k++)
    {
        Vec3_t v0(0.0,0.0,0.0);
        Dom.Lat[0].GetCell(iVec3_t(i,j,k))->Initialize(1.0,v0);
    }

    //Setting boundary conditions
    //for (size_t i=0;i<dat.xmin.Size();i++)
    //{
        //dat.xmin[i]->RhoBC = 1.0 + Dp;
        //dat.xmax[i]->RhoBC = 1.0;
        //dat.xmin[i]->Initialize(1.0+Dp,OrthoSys::O);
        //dat.xmax[i]->Initialize(1.0,    OrthoSys::O);
    //}

    Dom.AddSphere(-1,Vec3_t(0.5*nx*dx,0.5*ny*dx,0.5*nz*dx),R,3.0);
    Dom.Particles[0]->FixVeloc();
    //Dom.Particles[0]->v = Vec3_t(vel,0.0,0.0);

    //Solving
    Dom.Solve(Tf,0.01*Tf,Setup,Report,"test08",true,Nproc);
}
MECHSYS_CATCH

