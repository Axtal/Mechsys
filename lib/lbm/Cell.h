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

#ifndef MECHSYS_LBM_CELL_H
#define MECHSYS_LBM_CELL_H

// MechSys
#include <mechsys/util/array.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>

enum LBMethod
{
    D2Q9,     ///< 2D 9 velocities
    D3Q15,    ///< 3D 15 velocities
    D3Q27     ///< 3D 27 velocities
};

class Cell
{
public:
	static const double   WEIGHTSD2Q9   [ 9]; ///< Weights for the equilibrium distribution functions (D2Q9)
	static const double   WEIGHTSD3Q15  [15]; ///< Weights for the equilibrium distribution functions (D3Q15)
	//static const double   WEIGHTSD3Q27  [27]; ///< Weights for the equilibrium distribution functions (D3Q27)
	static const Vec3_t   LVELOCD2Q9    [ 9]; ///< Local velocities (D2Q9) 
	static const Vec3_t   LVELOCD3Q15   [15]; ///< Local velocities (D3Q15)
	//static const Vec3_t   LVELOCD3Q27   [27]; ///< Local velocities (D3Q27)
	static const size_t   OPPOSITED2Q9  [ 9]; ///< Opposite directions (D2Q9) 
	static const size_t   OPPOSITED3Q15 [15]; ///< Opposite directions (D3Q15)
	//static const size_t   OPPOSITED3Q27 [27]; ///< Opposite directions (D3Q27)
   
    //Constructor
    Cell (int Tag, LBMethod Method, iVec3_t Indexes, iVec3_t Ndim, double Cs, double Tau); ///< Constructor, it receives the grid type, ht einteger position and the total integer dimension of the domain and the spatial and time steps
    
    // Methods
    double       Density();                                        ///< Calculate density
    void         Velocity(Vec3_t & V);                             ///< Calculate velocity
    double       VelDen(Vec3_t & V);                               ///< Calculate both (faster)
    double       Feq(size_t k, Vec3_t const & V, double Rho);      ///< Calculate the equilibrium distribution function
    void         Initialize(double Rho, Vec3_t const & V);         ///< Initialize cell with a given velocity and density

    // Data
    LBMethod     Method;   ///< Is 2D, 3D and how many velocities it has
    bool         IsSolid;  ///< It is a solid node
    double       Tau;      ///< Relaxation Time
    double       Gamma;    ///< Solid/Fluid ratio
    double       Gs;       ///< Interaction constant between solid and fluid
    size_t       Nneigh;   ///< Number of neighbors
    int          Tag;      ///< Tag for the particle
    double       Cs;       ///< Velocity of the grid
    iVec3_t      Index;    ///< Vector of indexes
    Vec3_t       VelP;     ///< Velocity of the contact particle
    Vec3_t       VelBC;    ///< Velocity at boundary
    Vec3_t       BForce;   ///< Applied body force
    Vec3_t       BForcef;  ///< Fixed Applied body force
    double       RhoBC;    ///< Density at boundary

    size_t const  * Op;    ///< Pointer to opposite velocities
    double const  * W;     ///< Pointer to array of weights
    Vec3_t const  * C;     ///< Pointer to velocity constants
    Array<double>   F;     ///< Distribution functions
    Array<double>   Ftemp; ///< Temporary distribution functions
    Array<double>   Omeis; ///< Array of collision operators
    Array<size_t>   Neighs;///< Array of neighbors indexes
};

inline Cell::Cell(int TheTag, LBMethod TheMethod, iVec3_t TheIndexes, iVec3_t TheNdim, double TheCs, double TheTau)
{
    Tag    = TheTag;
    Method = TheMethod;
    Index  = TheIndexes;
    Cs     = TheCs;
    Gamma  = 0.0;
    Tau    = TheTau;
    if (Method==D2Q9)
    {
        Nneigh = 9;
        W      = WEIGHTSD2Q9;
        C      = LVELOCD2Q9;
        Op     = OPPOSITED2Q9;
    }
    if (Method==D3Q15)
    {
        Nneigh = 15;
        W      = WEIGHTSD3Q15;
        C      = LVELOCD3Q15;
        Op     = OPPOSITED3Q15;
    }
    F.     Resize(Nneigh);
    Ftemp. Resize(Nneigh);
    Neighs.Resize(Nneigh);
    Omeis .Resize(Nneigh);
    BForcef = 0.0,0.0,0.0;

    //Set neighbors
    for (size_t k=0;k<Nneigh;k++)
    {
        //iVec3_t nindex = Index + C[k];
        blitz::TinyVector<int,3>  nindex = Index + C[k];
        if (nindex[0]==        -1) nindex[0] = TheNdim[0]-1;
        if (nindex[0]==TheNdim[0]) nindex[0] = 0;
        if (nindex[1]==        -1) nindex[1] = TheNdim[1]-1;
        if (nindex[1]==TheNdim[1]) nindex[1] = 0;
        if (nindex[2]==        -1) nindex[2] = TheNdim[2]-1;
        if (nindex[2]==TheNdim[2]) nindex[2] = 0;

        Neighs[k] =  nindex[0] + nindex[1]*TheNdim[0] + nindex[2]*TheNdim[0]*TheNdim[1];
    }
}

inline double Cell::Density()
{
    if (IsSolid) return 0.0;
    double rho = 0.0;
    for (size_t k=0;k<Nneigh;k++) rho += F[k];
    return rho;
}

inline void Cell::Velocity(Vec3_t & V)
{
    V = 0.0, 0.0, 0.0;
    if (IsSolid) return;
    double rho = Density();
    if (rho<1.0e-12) return;
    for (size_t k=0;k<Nneigh;k++) V += F[k]*C[k];
    V *= Cs/rho;
}

inline double Cell::VelDen(Vec3_t & V)
{
    V = 0.0, 0.0, 0.0;
    if (IsSolid) return 0.0;
    double rho = Density();
    if (rho<1.0e-12) return 0.0;
    for (size_t k=0;k<Nneigh;k++) V += F[k]*C[k];
    V *= Cs/rho;
    return rho;
}

inline double Cell::Feq(size_t k, Vec3_t const & V, double Rho)
{
    double VdotC = dot(V,C[k]);
    double VdotV = dot(V,V);
    return W[k]*Rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}

inline void Cell::Initialize(double Rho, Vec3_t const & V)
{
    for (size_t k=0;k<Nneigh;k++) F[k] = Feq(k,V,Rho);
}


const double Cell::WEIGHTSD2Q9 [ 9] = { 4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36. };
const double Cell::WEIGHTSD3Q15[15] = { 2./9., 1./9., 1./9., 1./9., 1./9.,  1./9.,  1./9., 1./72., 1./72. , 1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
const Vec3_t Cell::LVELOCD2Q9  [ 9] = { {0,0,0}, {1,0,0}, {0,1,0}, {-1,0,0}, {0,-1,0}, {1,1,0}, {-1,1,0}, {-1,-1,0}, {1,-1,0} };
const Vec3_t Cell::LVELOCD3Q15 [15] =
{
	{ 0, 0, 0}, { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, 
	{ 0, 0, 1}, { 0, 0,-1}, { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, 
	{-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1} 
};
const size_t Cell::OPPOSITED2Q9  [ 9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };                       ///< Opposite directions (D2Q9) 
const size_t Cell::OPPOSITED3Q15 [15] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13}; ///< Opposite directions (D3Q15)
#endif // MECHSYS_LBM_CELL_H
