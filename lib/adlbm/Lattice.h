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

#ifndef MECHSYS_ADLBM_LATTICE_H
#define MECHSYS_ADLBM_LATTICE_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// MechSys
#include <mechsys/adlbm/Cell.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>

class Lattice
{
public:
    //typedefs
    typedef void (*ptFun_t) (Lattice & Lat, void * UserData);

    //Constructors
    Lattice () {};            //Default
    Lattice (double Thenu, double TheDif, iVec3_t TheNdim, double Thedx, double Thedt);

    //Methods
    void Stream1    (size_t Np);                                      ///< Stream the velocity distributions
    void Stream2    (size_t Np);                                      ///< Stream the velocity distributions
    void CalcProps  (size_t Np);                                      ///< Calculate the fluid properties
    Cell * GetCell(iVec3_t const & v);                                ///< Get pointer to cell at v


    //Data
    iVec3_t                                   Ndim;             // Dimensions of the lattice
    size_t                                    idx_out;          // The discrete time step
    size_t                                    Ncells;           // Number of cells per lattice
    double                                    Time;             // The current time
    double                                    dx;               // grid space
    double                                    dt;               // time step
    double                                    Nu;               // Real viscosity
    double                                    Dif;              // Diffusion coefficient
    double                                    Tau;              // Relaxation time
    double                                    Tauc;             // Relaxation time for diffusion 
    Cell                                   ** Cells;            // Array of pointer cells
};

inline Lattice::Lattice(double Thenu, double TheDif, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Ndim = TheNdim;
    dx   = Thedx;
    dt   = Thedt;
    Nu   = Thenu;
    Dif  = TheDif;
    Tau  = 3.0*Nu *dt/(dx*dx) + 0.5;
    Tauc = 3.0*Dif*dt/(dx*dx) + 0.5;



    //Cells.Resize(Ndim[0]*Ndim[1]*Ndim[2]);
    Cells = new Cell * [Ndim[0]*Ndim[1]*Ndim[2]];
    Ncells = Ndim[0]*Ndim[1]*Ndim[2];
    size_t n = 0;
    for (size_t k=0;k<Ndim[2];k++)
    for (size_t j=0;j<Ndim[1];j++)
    for (size_t i=0;i<Ndim[0];i++)
    {
        //Cells[n] =  new Cell(n,TheMethod,iVec3_t(i,j,k),Ndim,dx/dt,Tau);
        Cells[n] = new Cell(n,iVec3_t(i,j,k),Ndim,dx/dt,dt);
        n++;
    } 
    Cells[0]->Cs = dx/dt;
}

inline void Lattice::Stream1(size_t Np)
{
    // Assign temporal distributions
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Ncells;i++)
    for (size_t j=0;j<Cell::Nneigh;j++)
    {
        Cells[Cells[i]->Neighs[j]]->Ftemp[j] = Cells[i]->F[j];
        Cells[Cells[i]->Neighs[j]]->Gtemp[j] = Cells[i]->G[j];
    }
}

inline void Lattice::Stream2(size_t Np)
{
    //Swap the distribution values
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Ncells;i++)
    {
        double * Ftemp  = Cells[i]->F;
        Cells[i]->F     = Cells[i]->Ftemp;
        Cells[i]->Ftemp = Ftemp;
        double * Gtemp  = Cells[i]->G;
        Cells[i]->G     = Cells[i]->Gtemp;
        Cells[i]->Gtemp = Gtemp;
    }
}

inline void Lattice::CalcProps(size_t Np)
{
    //Calculate fields and densities
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Ncells;i++)
    {
        Cells[i]->CalcProp();
    }
}



inline Cell * Lattice::GetCell(iVec3_t const & v)
{
    return Cells[v[0] + v[1]*Ndim[0] + v[2]*Ndim[0]*Ndim[1]];
}

#endif
