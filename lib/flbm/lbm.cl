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
    uint3      Ndim;        ///< Integer vector with the dimensions of the LBM domain
    uint3      CellPairs;   ///< Vector with the cell pair system for for calculation
    double3    C[27];       ///< Collection of discrete velocity vectors
    double     Ekk[27];     ///< Dyadic product of discrete velocities for LES calculation
    double     W[27];       ///< Collection of discrete weights
    double     Tau[2];      ///< Collection of characteristic collision times
    double     G[2];        ///< Collection of cohesive constants for multiphase simulation
    double     Rhoref[2];   ///< Collection of cohesive constants for multiphase simulation
    double     Psi[2];      ///< Collection of cohesive constants for multiphase simulation
    double     Gmix;        ///< Repulsion constant for multicomponent simulation
    double     Cs;          ///< Lattice speed
    
} d_lbm_aux;


void kernel ApplyForces (global bool * IsSolid, global double * F, global double * Ftemp, global double3* BForce, global double3* Vel, global double * Rho, global struct lbm_aux * lbmaux)
{                                                                                                             
                                                                                                              
}                                                                                                             
                                                                                                              
void kernel Collide     (global bool * IsSolid, global double * F, global double * Ftemp, global double3* BForce, global double3* Vel, global double * Rho, global struct lbm_aux * lbmaux)
{                                                                                                             
                                                                                                              
}                                                                                                             
void kernel Stream      (global bool * IsSolid, global double * F, global double * Ftemp, global double3* BForce, global double3* Vel, global double * Rho, global struct lbm_aux * lbmaux)
{
    
}
