/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2016 Sergio Galindo                                    *
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

/////////////////////////////LBM OpenCL implementation////////////////////

typedef struct lbm_aux
{
    size_t     Nl;          ///< Numver of Lattices
    size_t     Nneigh;      ///< Number of Neighbors
    size_t     NCPairs;     ///< Number of cell pairs
    size_t     Nx;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Ny;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Nz;          ///< Integer vector with the dimensions of the LBM domain
    size_t     Op[27];      ///< Array with the opposite directions for bounce back calculation
    double3    C[27];       ///< Collection of discrete velocity vectors
    double     EEk[27];     ///< Dyadic product of discrete velocities for LES calculation
    double     W[27];       ///< Collection of discrete weights
    double     Tau[2];      ///< Collection of characteristic collision times
    double     G[2];        ///< Collection of cohesive constants for multiphase simulation
    double     Rhoref[2];   ///< Collection of cohesive constants for multiphase simulation
    double     Psi[2];      ///< Collection of cohesive constants for multiphase simulation
    double     Gmix;        ///< Repulsion constant for multicomponent simulation
    double     Cs;          ///< Lattice speed
    
} d_lbm_aux;

//ulong Pt2idx(uint3 iv,uint3 Dim)
//{
    //return iv.x + iv.y*Dim.x + iv.z*Dim.x*Dim.y;
//}
//
//ulong3 idx2Pt(size_t n,uint3 Dim)
//{
    //ulong3 tmp;
    //tmp.x = n%Dim.x;
    //tmp.y = (n/Dim.x)%Dim.y;
    //tmp.z = n/(Dim.x*Dim.y);
    //return tmp;
//}

void kernel CheckUpLoad (global struct lbm_aux * lbmaux)
{
    //printf("Nl          %d \n",  lbmaux[0].Nl     );
    //printf("Nneigh      %d \n",  lbmaux[0].Nneigh );
    //printf("NCP         %d \n",  lbmaux[0].NCPairs);
    //printf("Dim      %d %d \n",0,lbmaux[0].Nx );
    //printf("Dim      %d %d \n",1,lbmaux[0].Ny );
    //printf("Dim      %d %d \n",2,lbmaux[0].Nz );
    //for (size_t i=0;i < lbmaux[0].Nneigh;i++)
    //{
        //printf("C      %d %f %f %f \n",i,lbmaux[0].C[i].x,lbmaux[0].C[i].y,lbmaux[0].C[i].z);
    //}
    //for (size_t i=0;i < lbmaux[0].Nneigh;i++)
    //{
        //printf("EEk    %d %f       \n",i,lbmaux[0].EEk[i]);
    //}
}

void kernel ApplyForcesSC(global const bool * IsSolid, global double3* BForce, global const double * Rho, global const struct lbm_aux * lbmaux)
{
    //uint3 ic = idx2Pt(get_global_id(0),lbmaux[0].Ndim);
    //printf("%d %d %d %d \n",get_global_id(0),ic.s1,ic.s2,ic.s3);
    size_t ic  = get_global_id(0);
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = ic/(lbmaux[0].Nx*lbmaux[0].Ny);

    //BForce[ic] = (double3)(0.0,0.0,0.0);
    for (size_t k=1;k<lbmaux[0].Nneigh;k++)
    {
        size_t inx = (size_t)((long)icx + (long)lbmaux[0].C[k].x + (long)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((long)icy + (long)lbmaux[0].C[k].y + (long)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((long)icz + (long)lbmaux[0].C[k].z + (long)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;
        
        //if (ic==0) printf("ic: %lu %lu %lu %lu k: %lu in: %lu %lu %lu %lu \n",ic,icx,icy,icz,k,in,inx,iny,inz);
        
        double psic = 0.0;
        double psin = 0.0;
        if (!IsSolid[ic]) psic = lbmaux[0].Psi[0]*exp(-lbmaux[0].Rhoref[0]/Rho[ic]);
        if (!IsSolid[in]) psin = lbmaux[0].Psi[0]*exp(-lbmaux[0].Rhoref[0]/Rho[in]);
        
        BForce[ic] += -lbmaux[0].G[0]*lbmaux[0].W[k]*psic*psin*lbmaux[0].C[k];

        //if (ic==0) printf("k: %lu %f %f %f \n",k,BForce[ic].x,BForce[ic].y,BForce[ic].z);
    }

    //printf("%lu %lu %lu %lu \n",ic,icx,icy,icz);
    //printf("%d Dim %d %d %d \n",ic,lbmaux[0].Nx,lbmaux[0].Ny,lbmaux[0].Nz);
    //printf("Dim      %d %d \n",ic,icx );
    //printf("Dim      %d %d \n",ic,icy );
    //printf("Dim      %d %d \n",ic,icz );
}
                                                                                                              
void kernel CollideSC    (global const bool * IsSolid, global double * F, global double * Ftemp, global const double3* BForce, global const double3* Vel, global const double * Rho, global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    if (!IsSolid[ic])
    {
        double3 vel  = Vel[ic]+BForce[ic]/Rho[ic];
        double rho   = Rho[ic];
        double VdotV = dot(vel,vel);
        double Cs    = lbmaux[0].Cs;
        double tau   = lbmaux[0].Tau[0];

        bool valid = true;
        double alpha = 1.0;
        while (valid)
        {
            valid = false;
            for (size_t k=0;k<lbmaux[0].Nneigh;k++)
            {
                double VdotC = dot(vel,lbmaux[0].C[k]);
                double Feq   = lbmaux[0].W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k] - alpha*(F[ic*lbmaux[0].Nneigh + k] - Feq)/tau;
                if (Ftemp[ic*lbmaux[0].Nneigh + k]<-1.0e-12)
                {
                    double temp = tau*F[ic*lbmaux[0].Nneigh + k]/(F[ic*lbmaux[0].Nneigh + k] - Feq);
                    if (temp<alpha) alpha = temp;
                    valid = true;
                }
            }
        }
    }
    else
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            Ftemp[ic*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + lbmaux[0].Op[k]]; 
        }
    }
    for (size_t k=0;k<lbmaux[0].Nneigh;k++)
    {
        F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k]; 
    }
}

void kernel Stream1     (global double * F, global double * Ftemp, global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    size_t icx = ic%lbmaux[0].Nx;
    size_t icy = (ic/lbmaux[0].Nx)%lbmaux[0].Ny;
    size_t icz = ic/(lbmaux[0].Nx*lbmaux[0].Ny);
    for (size_t k=1;k<lbmaux[0].Nneigh;k++)
    {
        size_t inx = (size_t)((long)icx + (long)lbmaux[0].C[k].x + (long)lbmaux[0].Nx)%lbmaux[0].Nx;
        size_t iny = (size_t)((long)icy + (long)lbmaux[0].C[k].y + (long)lbmaux[0].Ny)%lbmaux[0].Ny;
        size_t inz = (size_t)((long)icz + (long)lbmaux[0].C[k].z + (long)lbmaux[0].Nz)%lbmaux[0].Nz;
        size_t in  = inx + iny*lbmaux[0].Nx + inz*lbmaux[0].Nx*lbmaux[0].Ny;
        Ftemp[in*lbmaux[0].Nneigh + k] = F[ic*lbmaux[0].Nneigh + k];
    }
}

void kernel Stream2     (global const bool * IsSolid, global double * F, global double * Ftemp, global double3* BForce, global double3* Vel, global double * Rho, global const struct lbm_aux * lbmaux)
{
    size_t ic  = get_global_id(0);
    Rho   [ic] = 0.0;
    Vel   [ic] = (double3)(0.0,0.0,0.0);
    BForce[ic] = (double3)(0.0,0.0,0.0);
    if (!IsSolid[ic])
    {
        for (size_t k=0;k<lbmaux[0].Nneigh;k++)
        {
            F[ic*lbmaux[0].Nneigh + k] = Ftemp[ic*lbmaux[0].Nneigh + k];
            Rho[ic] += F[ic*lbmaux[0].Nneigh + k];
            Vel[ic] += F[ic*lbmaux[0].Nneigh + k]*lbmaux[0].C[k];
            //if (ic==0) printf("k: %lu %f %f %f %f \n",k,Rho[ic],Vel[ic].x,Vel[ic].y,Vel[ic].z);
        }
        Vel[ic] /= Rho[ic];
    }
}
