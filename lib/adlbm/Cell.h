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

#ifndef MECHSYS_ADLBM_CELL_H
#define MECHSYS_ADLBM_CELL_H

// Std lib
#ifdef USE_THREAD
    #include <pthread.h>
#endif

//OpenMP
#ifdef USE_OMP
    #include <omp.h>
#endif

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>



class Cell
{
public:
	static const Vec3_t  C   [9]; ///< Local velocities (D2Q9)
	static const double  W   [9]; ///< Local weights (D2Q9)
    static       double  Cs;      ///< Speed parameter
    static const size_t  Op  [9]; ///< Index of opposing cells for bounce back
    static const size_t  Nneigh;  ///< Number of neighbors
   
    //Constructor
    Cell (size_t ID, iVec3_t Indexes, iVec3_t Ndim, double Cs, double Dt); ///< Constructor, it receives the grid type, ht einteger position and the total integer dimension of the domain and the spatial and time steps
    
    // Methods
    double       Feq(size_t k);                         ///< Calculate the equilibrium distribution function F
    double       Geq(size_t k);                         ///< Calculate the equilibrium distribution function G
    void         CalcProp();                                       ///< Calculate the vectorial properties with the new distributions functions
    void         Initialize(double TheRho, double TheTemp, Vec3_t & TheV);       ///< Initialize cell with a given velocity and density

    // Data
    size_t       ID;       ///< Tag for the particle
    double       Dt;       ///< Time step
    
    bool         IsSolid;  ///< It is a solid node
    iVec3_t      Index;    ///< Vector of indexes

    double        *      F; ///< Distribution functions for vector potentials
    double        *  Ftemp; ///< Temporary distribution functions
    double        *      G; ///< Distribution functions for vector potentials
    double        *  Gtemp; ///< Temporary distribution functions
    size_t        * Neighs; ///< Array of neighbors indexes
    double             Dif; ///< Diffusion coefficient
    double            Tauc; ///< Local Relaxation Time
    double             Rho; ///< Fluid density
    double            Temp; ///< Temperature
    Vec3_t             Vel; ///< Fluid velocity
    Vec3_t            Flux; ///< Heat Flux

};

inline Cell::Cell(size_t TheID, iVec3_t TheIndexes, iVec3_t TheNdim, double TheCs, double TheDt)
{
    ID      = TheID;
    Index   = TheIndexes;
    Dt      = TheDt;
    F      = new double [Nneigh];
    Ftemp  = new double [Nneigh];
    G      = new double [Nneigh];
    Gtemp  = new double [Nneigh];
    Initialize(1.0,1.0,OrthoSys::O);

    Neighs  = new size_t [Nneigh];

    //Set neighbors
    for (size_t k=0;k<Nneigh;k++)
    {
        //iVec3_t nindex = Index + C[k];
        blitz::TinyVector<int,3>   nindex;
        nindex[0] = static_cast<int>(Index[0])+ C[k][0];
        nindex[1] = static_cast<int>(Index[1])+ C[k][1];
        nindex[2] = static_cast<int>(Index[2])+ C[k][2];
        if (nindex[0]==                          -1) nindex[0] = TheNdim[0]-1;
        if (nindex[0]==static_cast<int>(TheNdim[0])) nindex[0] = 0;
        if (nindex[1]==                          -1) nindex[1] = TheNdim[1]-1;
        if (nindex[1]==static_cast<int>(TheNdim[1])) nindex[1] = 0;
        if (nindex[2]==                          -1) nindex[2] = TheNdim[2]-1;
        if (nindex[2]==static_cast<int>(TheNdim[2])) nindex[2] = 0;

        Neighs[k] =  nindex[0] + nindex[1]*TheNdim[0] + nindex[2]*TheNdim[0]*TheNdim[1];
    }

}


inline void Cell::CalcProp()
{
    Rho  = 0.0;
    Temp = 0.0;
    Vel  = OrthoSys::O;
    Flux = OrthoSys::O;
    if (!IsSolid)
    {
        for (size_t i=0;i<Nneigh;i++)
        {
            Rho  += F[i];
            Vel  += Cs*F[i]*C[i];
        }
        Vel /= Rho;
    }
    for (size_t i=0;i<Nneigh;i++)
    {
        Temp += G[i];
        Flux += (F[i]-Feq(i))*C[i];
    }
    Flux = 6.0*Flux/(Dt*Cs*Cs*(1.0+2.0*Tauc))+Temp*Vel;
}

inline double Cell::Feq(size_t k)
{   
    double VdotC = dot(Vel,C[k]);
    double VdotV = dot(Vel,Vel);
    return W[k]*Rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

inline double Cell::Geq(size_t k)
{
    double VdotC = dot(Vel,C[k]);
    double VdotV = dot(Vel,Vel);
    return W[k]*Temp*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

inline void Cell::Initialize(double TheRho, double TheTemp, Vec3_t & TheVel)
{
    Rho = TheRho;
    Temp= TheTemp;
    Vel = TheVel;

    for (size_t i=0;i<Nneigh;i++)
    {
        F[i] = Feq(i);
        G[i] = Geq(i);
    }
    if (IsSolid)
    {
        Rho = 0.0;
        Vel = OrthoSys::O;
    }
}



const Vec3_t Cell::C   [9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
const double Cell::W   [9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const size_t Cell::Op  [9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 }; 
const size_t Cell::Nneigh  = 9;
      double Cell::Cs      = 1.0;



#endif // MECHSYS_ADLBM_CELL_H
