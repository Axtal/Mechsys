/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2016 Sergio Galindo and Zi Li                          *
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

//=========================================================================================================
// This is a reduced/improved version of sergio's original Domain.h file
// GPU part is cut/reserved for the moment   
// Contact force + EFS                                                                                     // --- 2017.03.08 
//=========================================================================================================










//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& MECHSYS_FLBM &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

//======================================== Head files =====================================================
#ifndef MECHSYS_FLBM_DOMAIN_H
#define MECHSYS_FLBM_DOMAIN_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

//************************************---p1---*********************************

// Std lib
#ifdef USE_OMP
#include <omp.h>
#endif

// STD
#include<iostream>

// Mechsys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>

enum LBMethod
{
    D2Q5,                     // 2D 5  velocities
    D2Q9,                     // 2D 9  velocities
    D3Q15,                    // 3D 15 velocities
    D3Q19,                    // 3D 19 velocities
};                            // please note that there is <;> 
//======================================== Head files =====================================================










//######################################## namespace FLBM #################################################
namespace FLBM
{

//====================================== inline functions =================================================

//************************************---p2---*********************************
inline size_t Pt2idx(iVec3_t iv, iVec3_t & Dim)  { return iv(0) + iv(1)*Dim(0) + iv(2)*Dim(0)*Dim(1); }    // Calculates the index of the cell at coordinates iv for a cubic lattice of dimensions Dim

inline void   idx2Pt(size_t n, iVec3_t & iv, iVec3_t & Dim)                                                // Calculates the coordinates from the index --- interesting! put return in the statement
{
    iv(0) = n%Dim(0);
    iv(1) =(n/Dim(0))%(Dim(1));
    iv(2) = n/(Dim(0)*Dim(1));
}
//======================================= inline functions ================================================










//==================================== class_Domain statements ============================================
class Domain
{
typedef void (*ptDFun_t) (Domain & Dom, void * UserData);                                                  // --- what does this line means?  --- 20170302

public:
    static const double   WEIGHTSD2Q5    [ 5];              // Weights for the equilibrium distribution functions (D2Q5)
    static const double   WEIGHTSD2Q9    [ 9];              // Weights for the equilibrium distribution functions (D2Q9)
    static const double   WEIGHTSD3Q15   [15];              // Weights for the equilibrium distribution functions (D3Q15)
    static const double   WEIGHTSD3Q19   [19];              // Weights for the equilibrium distribution functions (D3Q19)
 
    static const size_t   OPPOSITED2Q5   [ 5];              // Opposite directions (D2Q5) 
    static const size_t   OPPOSITED2Q9   [ 9];              // Opposite directions (D2Q9) 
    static const size_t   OPPOSITED3Q15  [15];              // Opposite directions (D3Q15)
    static const size_t   OPPOSITED3Q19  [19];              // Opposite directions (D3Q19)

    static const Vec3_t   LVELOCD2Q5     [ 5];              // Local velocities (D2Q5) 
    static const Vec3_t   LVELOCD2Q9     [ 9];              // Local velocities (D2Q9) 
    static const Vec3_t   LVELOCD3Q15    [15];              // Local velocities (D3Q15)
    static const Vec3_t   LVELOCD3Q19    [19];              // Local velocities (D3Q19)

    #ifdef USE_OMP
    omp_lock_t   lck;                                       // to protect variables in multithreading
    #endif

    //Data variables
    double       dt;                                        // Time Step
    double       dx;                                        // Grid size
    double       Cs;                                        // Lattice Velocity              //it is not sound speed, and be very careful           ---20161021 ---20170301

    iVec3_t      Ndim;                                      // Lattice Dimensions            // nx,ny,nz, only 3 numbers
    size_t       Nl;                                        // Number of lattices (fluids)   //single/two component(s)                              ---20161021
    size_t       Ncells;                                    // Number of cells
    size_t       Nneigh;                                    // Number of Neighbors, depends on the scheme

    size_t       NCellPairs;                                // Number of cell pairs --- Array for pair calculation
    iVec3_t   *  CellPairs;                                 // Pairs of cells for molecular force calculation -- fluid-fluid interaction 
    size_t       NCellPaira;                                // Number of cell pairs --- Array for pair calculation                                  ---20161101
    iVec3_t   *  CellPaira;                                 // Pairs of cells for molecular force calculation -- fluid-fluid interaction            ---20161101

    double ***** F;                                         // The array containing the pdf functions with the lattice fluid, the t,x,y,z coordinates and the direction of velocity space 
    double ***** Ftemp;                                     // A similar array to hold provisional data
    double ***** F2;                                        // The array                                                                            ---20170301 
    double ***** F2temp;                                    // The array                                                                            ---20170301
    double ***** F3;                                        // The array                                                                            ---20170301
    double ***** F3temp;                                    // The array                                                                            ---20170301

    Vec3_t const * C;                                       // The array of lattice velocities
    double *     EEk;                                       // Dyadic product of the velocity vectors     //it should be a dot product--20161021 
    size_t const * Op;                                      // An array containing the indexes of the opposite direction for bounce back conditions
    double const *  W;                                      // An array with the direction weights

    double *     G;                                         // The attractive constant for multiphase simulations
    double *     Psi;                                       // Parameters for the Shan Chen pseudo potential
    double *     Rhoref;                                    // Parameters for the Shan Chen pseudo potential
    double *     Gs;                                        // The attractive constant for solid phase

    double       Gmix;                                      // The mixing constant for multicomponent simulations
    double       Gf;                                        // The mixing constant for exponentially-decreasing fluid-fluid interaction             --- 20161125
    double *     Gw;                                        // Parameters for the exponentially-decreasing fluid-solid interaction                  --- 20161101
    double       Elta;                                      // Length scale of the exponentially-decreasing fluid-solid interaction                 --- 20161101

    Vec3_t ****  BForce;                                    // Body Force for each cell   --3D array
    double *     Tau;                                       // The characteristic time of the lattice
    Vec3_t ****  Vel;                                       // The fluid velocities       --3D array
    double ****  Rho;                                       // The fluid densities

    bool   ****  IsSolid;                                   // An identifier to see if the cell is a solid cell
    bool   ****  IsSolidEdge;                               // An identifier to see if the cell is a solid cell at the surface                      --- 20161208
    bool         IsFirstTime;                               // Bool variable checking if it is the first time function Setup is called

    size_t       Nproc;                                     // Number of processors for openmp
    size_t       idx_out;                                   // The discrete time step for output
    String       FileKey;                                   // File Key for output files
    void *       UserData;                                  // User Data

    double       Sc;                                        // Smagorinsky constant                                    //what is it?                ---20161021
    size_t       Step;                                      // Length of averaging cube to save data
    double       Time;                                      // Simulation time variable

    //Special constructor with only one component and having the same parameters as above
    Domain (LBMethod      Method,                           // Type of array, for example D2Q9
    double                nu,                               // Viscosity for each fluid
    iVec3_t               Ndim,                             // Cell divisions per side
    double                dx,                               // Space spacing
    double                dt);                              // Time step

    //Methods--function--process--20161021
    void   WriteXDMF(char const * FileKey);                                    // 1< Write the domain data in xdmf file
    void   Initialize(size_t k, iVec3_t idx, double Rho, Vec3_t & Vel);        // 2< Initialize each cell with a given density and velocity
    void   SolidPhase(Vec3_t const & XC, double RC);                           // 3< Same function as SolidDisk() in /tlbm (old Mechsys)                                        -- 20161222

    double Feq(size_t k, double Rho, Vec3_t & Vel);                            // 4< The equilibrium function
    void   Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL,               // 5< Solve the Domain dynamics    
                 ptDFun_t ptReport=NULL,  char const * FileKey=NULL, 
                 bool RenderVideo=true,   size_t Nproc=1);                     

    void   ApplyForcesSC();                                                    // 6<  Apply the molecular forces for the single component case
    void   CollideSC();                                                        // 7<  The collide step of LBM for single component simulations
    void   StreamSC();                                                         // 8<  The stream step of LBM SC

    void   ApplyForces_contact();                                              // 9<  Contact force
    void   CollideEFS_SRT();                                                   // 10< Explicit forcing scheme (EFS)
    void   StreamEFS_SRT();                                                    // 11< The stream step of LBM SC


//************************************---p3---*********************************

//************************************---p4---*********************************
};                                                                             // class Domain and please note that there is <;> 
//======================================= class_Domain statements =========================================










//=========================================== Constant weights ============================================
const double Domain::WEIGHTSD2Q5   [ 5] = { 2./6. , 1./6. , 1./6. , 1./6. , 1./6. };
const double Domain::WEIGHTSD2Q9   [ 9] = { 4./9. , 1./9. , 1./9. , 1./9. , 1./9. , 1./36., 1./36., 1./36., 1./36.};
const double Domain::WEIGHTSD3Q15  [15] = { 2./9. , 1./9. , 1./9. , 1./9. , 1./9. , 1./9. , 1./9. , 1./72., 1./72., 
                                            1./72., 1./72., 1./72., 1./72., 1./72., 1./72.};
const double Domain::WEIGHTSD3Q19  [19] = { 1./3. , 1./18., 1./18., 1./18., 1./18., 1./18., 1./18., 1./36., 1./36., 1./36., 
                                            1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36., 1./36.};

const size_t Domain::OPPOSITED2Q5  [ 5] = { 0, 3, 4, 1, 2};                                                                // Opposite directions (D2Q5) 
const size_t Domain::OPPOSITED2Q9  [ 9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6};                                                    // Opposite directions (D2Q9) 
const size_t Domain::OPPOSITED3Q15 [15] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13};                             // Opposite directions (D3Q15)
const size_t Domain::OPPOSITED3Q19 [19] = { 0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};             // Opposite directions (D3Q19)

const Vec3_t Domain::LVELOCD2Q5    [ 5] = { { 0, 0, 0}, { 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}};                   // method of defining and citing 2D matrix   
const Vec3_t Domain::LVELOCD2Q9    [ 9] = { { 0, 0, 0}, { 1, 0, 0}, { 0, 1, 0}, {-1, 0, 0}, { 0,-1, 0}, 
                                            { 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}};
const Vec3_t Domain::LVELOCD3Q15   [15] = { { 0, 0, 0}, { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, 
					    { 0, 0, 1}, { 0, 0,-1}, { 1, 1, 1}, {-1,-1,-1}, { 1, 1,-1}, 
					    {-1,-1, 1}, { 1,-1, 1}, {-1, 1,-1}, { 1,-1,-1}, {-1, 1, 1}};
const Vec3_t Domain::LVELOCD3Q19   [19] = { { 0, 0, 0}, 
					    { 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1}, 
					    { 1, 1, 0}, {-1,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 1, 0, 1}, {-1, 0,-1},
					    { 1, 0,-1}, {-1, 0, 1}, { 0, 1, 1}, { 0,-1,-1}, { 0, 1,-1}, { 0,-1, 1}};
//=========================================== Constant weights ============================================










//=========================================== Single_component_Domain =====================================
inline Domain::Domain(LBMethod TheMethod, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)       
{
    Array<double> nu(1);
    nu[0] = Thenu;
    
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (nu.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one fluid please");
    if (TheNdim(2) >1&&(TheMethod==D2Q9 ||TheMethod==D2Q5 ))  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (TheNdim(2)==1&&(TheMethod==D3Q15||TheMethod==D3Q19))  throw new Fatal("LBM::Domain: Ndim(2) is 1. Either change the method to D2Q9 or increase the z-dimension");
   
    if (TheMethod==D2Q5)
    {
        Nneigh = 5;
        W      = WEIGHTSD2Q5;
        C      = LVELOCD2Q5;
        Op     = OPPOSITED2Q5;
    }
    if (TheMethod==D2Q9)
    {
        Nneigh = 9;
        W      = WEIGHTSD2Q9;
        C      = LVELOCD2Q9;
        Op     = OPPOSITED2Q9;

    }
    if (TheMethod==D3Q15)
    {
        Nneigh = 15;
        W      = WEIGHTSD3Q15;
        C      = LVELOCD3Q15;
        Op     = OPPOSITED3Q15;
    }
    if (TheMethod==D3Q19)
    {
        Nneigh = 19;
        W      = WEIGHTSD3Q19;
        C      = LVELOCD3Q19;
        Op     = OPPOSITED3Q19;
    }
    
    Time        = 0.0;
    dt          = Thedt;
    dx          = Thedx;
    Cs          = dx/dt;
    Step        = 1;
    Sc          = 0.17;
    Nl          = 1;
    Ndim        = TheNdim;
    Ncells      = Ndim(0)*Ndim(1)*Ndim(2);
    IsFirstTime = true;

    Tau    = new double [Nl];
    G      = new double [Nl];
    Gs     = new double [Nl];
    Rhoref = new double [Nl];
    Psi    = new double [Nl];
    Gw     = new double [Nl];                                                                              // --- 20161101

    Elta   = 0.0;                                                                                          // --- 20161129
    Gmix   = 0.0;
    Gf     = 0.0;

    F           = new double **** [Nl];
    Ftemp       = new double **** [Nl];
    F2          = new double **** [Nl];
    F2temp      = new double **** [Nl];
    F3          = new double **** [Nl];
    F3temp      = new double **** [Nl];

    Vel         = new Vec3_t ***  [Nl];
    BForce      = new Vec3_t ***  [Nl];
    Rho         = new double ***  [Nl];
    IsSolid     = new bool   ***  [Nl];
    IsSolidEdge = new bool   ***  [Nl];                                                                    // --- 20161208

    for (size_t i=0;i<Nl;i++)
    {
        Tau     [i]    = 3.0*nu[i]*dt/(dx*dx)+0.5;
        G       [i]    = 0.0;
        Gs      [i]    = 0.0;
        Rhoref  [i]    = 200.0;
        Psi     [i]    = 4.0;
        Gw      [i]    = 0.0;                                                                              // --- 20161101 

        F           [i]    = new double *** [Ndim(0)];
        Ftemp       [i]    = new double *** [Ndim(0)];
        F2          [i]    = new double *** [Ndim(0)];
        F2temp      [i]    = new double *** [Ndim(0)];
        F3          [i]    = new double *** [Ndim(0)];
        F3temp      [i]    = new double *** [Ndim(0)];

        Vel         [i]    = new Vec3_t **  [Ndim(0)];
        BForce      [i]    = new Vec3_t **  [Ndim(0)];
        Rho         [i]    = new double **  [Ndim(0)];
        IsSolid     [i]    = new bool   **  [Ndim(0)];
        IsSolidEdge [i]    = new bool   **  [Ndim(0)];                                                     // --- 20161208

        for (size_t nx=0;nx<Ndim(0);nx++)
        {
            F           [i][nx]    = new double ** [Ndim(1)];
            Ftemp       [i][nx]    = new double ** [Ndim(1)];
            F2          [i][nx]    = new double ** [Ndim(1)];
            F2temp      [i][nx]    = new double ** [Ndim(1)];
            F3          [i][nx]    = new double ** [Ndim(1)];
            F3temp      [i][nx]    = new double ** [Ndim(1)];

            Vel         [i][nx]    = new Vec3_t *  [Ndim(1)];
            BForce      [i][nx]    = new Vec3_t *  [Ndim(1)];
            Rho         [i][nx]    = new double *  [Ndim(1)];
            IsSolid     [i][nx]    = new bool   *  [Ndim(1)];
            IsSolidEdge [i][nx]    = new bool   *  [Ndim(1)];                                              // --- 20161208

            for (size_t ny=0;ny<Ndim(1);ny++)
            {
                F           [i][nx][ny]    = new double * [Ndim(2)];
                Ftemp       [i][nx][ny]    = new double * [Ndim(2)];
                F2          [i][nx][ny]    = new double * [Ndim(2)];
                F2temp      [i][nx][ny]    = new double * [Ndim(2)];
                F3          [i][nx][ny]    = new double * [Ndim(2)];
                F3temp      [i][nx][ny]    = new double * [Ndim(2)];

                Vel         [i][nx][ny]    = new Vec3_t   [Ndim(2)];
                BForce      [i][nx][ny]    = new Vec3_t   [Ndim(2)];
                Rho         [i][nx][ny]    = new double   [Ndim(2)];
                IsSolid     [i][nx][ny]    = new bool     [Ndim(2)];
                IsSolidEdge [i][nx][ny]    = new bool     [Ndim(2)];                                       // --- 20161208

                for (size_t nz=0;nz<Ndim(2);nz++)
                {
                    F        [i][nx][ny][nz]    = new double [Nneigh];
                    Ftemp    [i][nx][ny][nz]    = new double [Nneigh];
                    F2       [i][nx][ny][nz]    = new double [Nneigh];
                    F2temp   [i][nx][ny][nz]    = new double [Nneigh];
                    F3       [i][nx][ny][nz]    = new double [Nneigh];
                    F3temp   [i][nx][ny][nz]    = new double [Nneigh];

                    IsSolid    [i][nx][ny][nz]  = false;
                    IsSolidEdge[i][nx][ny][nz]  = false;                                                   // --- 20161208

                    for (size_t nn=0;nn<Nneigh;nn++)
                    {
                        F     [i][nx][ny][nz][nn] = 0.0;
                        Ftemp [i][nx][ny][nz][nn] = 0.0;
                        F2    [i][nx][ny][nz][nn] = 0.0;
                        F2temp[i][nx][ny][nz][nn] = 0.0;
                        F3    [i][nx][ny][nz][nn] = 0.0;
                        F3temp[i][nx][ny][nz][nn] = 0.0;
                    }
                }
            }
        }
    }

    EEk = new double [Nneigh];                                                                             // Useful line. 
    for (size_t k=0;k<Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(C[k][n]*C[k][m]);                                                               // This line tells me how to cite 2D matrix. --- 20170228
        }
    }

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Nl*Ncells,TERM_RST);
    #ifdef USE_OMP
    omp_init_lock(&lck);
    #endif
}
//=========================================== Single_component_Domain =====================================










//===================================== function 1 -- WriteXDMF ===========================================
inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Ndim[0]/Step;
    size_t  Ny = Ndim[1]/Step;
    size_t  Nz = Ndim[2]/Step;

    for (size_t j=0;j<Nl;j++)
    {
        // Creating data sets
        double * Density   = new double[  Nx*Ny*Nz];
        double * Gamma     = new double[  Nx*Ny*Nz];
        double * Vvec      = new double[3*Nx*Ny*Nz];

        size_t i=0;
        #ifdef USE_OMP
        for (size_t m=0;m<Ndim(2);m+=Step)
        for (size_t l=0;l<Ndim(1);l+=Step)
        for (size_t n=0;n<Ndim(0);n+=Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            Vec3_t vel    = OrthoSys::O;

            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho    += Rho    [j][n+ni][l+li][m+mi];
                gamma  += IsSolid[j][n+ni][l+li][m+mi] ? 1.0: 0.0;
                vel    += Vel    [j][n+ni][l+li][m+mi];
            }
            rho  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            vel  /= Step*Step*Step;
            Density [i]  = (double) rho;
            Gamma   [i]  = (double) gamma;
            Vvec[3*i  ]  = (double) vel(0);
            Vvec[3*i+1]  = (double) vel(1);
            Vvec[3*i+2]  = (double) vel(2);
            i++;
        }        
        #endif

        //Writing data to h5 file
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Density );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Gamma   );
        }
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity_%d",j);
        H5LTmake_dataset_double(file_id,dsname.CStr(),1,dims,Vvec    );
        dims[0] = 1;
        int N[1];
        N[0] = Nx;
        dsname.Printf("Nx");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Ny;
        dsname.Printf("Ny");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Nz;
        dsname.Printf("Nz");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);

        delete [] Density ;
        delete [] Gamma   ;
        delete [] Vvec    ;
    }

    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    //Writing xmf file
    std::ostringstream oss;

    if (Ndim(2)==1)  // 2-D
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Ndim(1) << " " << Ndim(0) << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";

        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }

        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Ndim(0) << " " << Ndim(1) << " " << Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    
    else            // 3-D
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"LBM_Mesh\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*dx << " " << Step*dx  << " " << Step*dx  << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";

        for (size_t j=0;j<Nl;j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }

        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}
//===================================== function 1 -- WriteXDMF ===========================================










//===================================== function 2 -- SolidPhase ==========================================
inline void Domain::SolidPhase(Vec3_t const & XC, double RC)                                               // impose to be solid nodes   --- 20170302 
{
    #ifdef USE_OMP
    for (size_t m=0; m<Ndim(0); m++)
    for (size_t n=0; n<Ndim(1); n++)
    for (size_t l=0; l<Ndim(2); l++)  
    {
        Vec3_t XX(m,n,l);
        if (norm(XX-XC) < RC)    IsSolid [0][m][n][l] = true ;                     
    }
    #endif
}
//===================================== function 2 -- SolidPhase ==========================================










//======================================== function 3 -- Feq ==============================================
inline double Domain::Feq(size_t k, double Rho, Vec3_t & V)
{
    double VdotC = dot(V,C[k]);
    double VdotV = dot(V,V);
    return W[k]*Rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
}
//======================================== function 3 -- Feq ==============================================










//===================================== function 4 -- Initialize ==========================================
inline void Domain::Initialize(size_t il, iVec3_t idx, double TheRho, Vec3_t & TheVel)                     // --- pointwise specification                                    --- 20170302 
{
    size_t nx = Ndim(0);                                                                                   // --- newly added                                                --- 20161215
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    size_t ix = idx(0);
    size_t iy = idx(1);
    size_t iz = idx(2);

    //BForce[il][ix][iy][iz] = OrthoSys::O;                                                                // BForce is specified externally and should not be added here     --- 20170301 

    if (!IsSolid[il][ix][iy][iz])
    {
        Vel[il][ix][iy][iz] = TheVel;
        Rho[il][ix][iy][iz] = TheRho;
    }

    else
    {
        Vel[il][ix][iy][iz] = OrthoSys::O;
        Rho[il][ix][iy][iz] = 0.0;
    }

    #ifdef USE_OMP
    for (size_t k=0;k<Nneigh;k++)
    {
        F [il][ix][iy][iz][k] = Feq(k,TheRho,TheVel);
        F2[il][ix][iy][iz][k] = 3.0/TheRho * Feq(k,TheRho,TheVel) * dot(BForce[il][ix][iy][iz],(C[k]-TheVel));
        F3[il][ix][iy][iz][k] = F[il][ix][iy][iz][k] - 0.5*dt*F2[il][ix][iy][iz][k];                       // dt has been specified in domian()                               --- 20170301 
    }
    #endif
                                                                        
    if (IsSolid[il][ix][iy][iz] == true)                                                                   // --- newly added IsSolidEdge identification                      --- 20161215
    {
        #ifdef USE_OMP
        for (size_t k=1;k<Nneigh;k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) ) ;
            size_t niy = (size_t)((int)iy + (int)C[k](1) ) ;
            size_t niz = (size_t)((int)iz + (int)C[k](2) ) ;

            if (nix>=0 && nix<nx &&  niy>=0 && niy<ny && niz>=0 && niz<nz && IsSolid[il][nix][niy][niz] == false)   
            //if (IsSolid[i][nix][niy][niz] == false)  
            {
                IsSolidEdge[il][ix][iy][iz]  = true;  
            }
        }
        #endif
    }
}
//===================================== function 4 -- Initialize ==========================================










//========================================= function 5 -- Solve ===========================================
inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{

//************************************---p6---*********************************

    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Nproc = TheNproc;

    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1    , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                    , TERM_RST);
    for (size_t i=0;i<Nl;i++)
    {
    printf("%s  Tau of Lattice %zd               =  %g%s\n"       ,TERM_CLR2, i, Tau[i]                             , TERM_RST);
    }

//************************************---p7---*********************************   

    {
        size_t nx      = Ndim(0);
        size_t ny      = Ndim(1);
        size_t nz      = Ndim(2);
        //size_t Naround = nx*ny*nz-1;                                                                     // --- each particle has (nx*ny*nz-1) partners in Gw-Elta force setup --- 20161101

//=============================================================================================== 
//      Array< Array<iVec3_t> > CPPT(Nproc) ;   // CPP  is defined here and shown first time             // --- the particle pair in G/Gmix/Gs-neighbour force setup           --- 20161125 

//      #ifdef USE_OMP
//      #pragma omp parallel for schedule(static) num_threads(Nproc)
//      #endif
//	for (size_t ix=0;ix<nx;ix++)
//      for (size_t iy=0;iy<ny;iy++)
//	for (size_t iz=0;iz<nz;iz++)
//        {
//            size_t nc = Pt2idx(iVec3_t(ix,iy,iz),Ndim);                                                  // --- Gs/G/Gmix-neighbour and Gw/Gf-Elta force setps center around central particle
//            for (size_t k=1;k<Nneigh;k++)                                                                // ---                                                                --- 20161125 
//            {
//                size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0)) % Ndim(0);
//                size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1)) % Ndim(1);
//                size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2)) % Ndim(2);

//                size_t nb  = Pt2idx(iVec3_t(nix,niy,niz),Ndim);                                          // --- nc-nb is number of centeral-neighbour particle                 --- 20161125
//                if (nb>nc)                                                                               // --- take a half of the whole as the fluid pair                     --- 20161125
//                {
//                    CPPT[omp_get_thread_num()].Push(iVec3_t(nc,nb,k));
//                }
//            }
//        }

//	Array<iVec3_t> CPP(0);
//	for (size_t nt=0;nt<Nproc;nt++)
//	for (size_t it=0;it<CPPT[nt].Size();it++)                                                          //--- use of size             --- 20170309 
//	{
//		CPP.Push(CPPT[nt][it]);
//	}
   
//        NCellPairs = CPP.Size();                                                                         // --- the particle pair in G/Gmix/Gs-neighbour force setup           --- 20161125 
//        CellPairs  = new iVec3_t [NCellPairs];

//        for (size_t n=0;n<NCellPairs;n++)
//        {
//            CellPairs[n] = CPP[n];
//        }
//=============================================================================================== 

        Array< Array<iVec3_t> > CPPaT(Nproc);   // CPPa is defined here and shown first time               // --- newly defined for the particle pair in Gw-Elta force setup     --- 20161101

        #ifdef USE_OMP
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        #endif
        for (size_t ix=0;ix<nx;ix++)  //nx
        for (size_t iy=0;iy<ny;iy++)  //ny
        for (size_t iz=0;iz<nz;iz++)
        {
            //if (IsSolid [0][ix][iy][iz])   
             
            if (IsSolidEdge [0][ix][iy][iz])                                                               // --- 4 limitations to save memory                                    --- 20170307
            {   size_t nc = Pt2idx(iVec3_t(ix,iy,iz),Ndim);                                                // --- Gs/G/Gmix-neighbour and Gw/Gf-Elta force setps center around central particle  
                double Ns      = 15;                                                                       // --- newly added module --- 20170214
                double k       =  0;

                double ixa_min =  0;
                double ixa_max =  0; 
                double iya_min =  0; 
                double iya_max =  0; 
                double iza_min =  0; 
                double iza_max =  0; 

                (double)(ix-Ns)> 0 ? ixa_min=(ix-Ns) : ixa_min= 0;
                (double)(ix+Ns)<nx ? ixa_max=(ix+Ns) : ixa_max=nx;
                (double)(iy-Ns)> 0 ? iya_min=(iy-Ns) : iya_min= 0;
                (double)(iy+Ns)<ny ? iya_max=(iy+Ns) : iya_max=ny;
                (double)(iz-Ns)> 0 ? iza_min=(iz-Ns) : iza_min= 0;
                (double)(iz+Ns)<nz ? iza_max=(iz+Ns) : iza_max=nz;

                for (size_t ixa=ixa_min;ixa<ixa_max;ixa++)
                for (size_t iya=iya_min;iya<iya_max;iya++)
                for (size_t iza=iza_min;iza<iza_max;iza++)
                {                                                                                          // --- newly added module --- 20161101
                    size_t  na = Pt2idx(iVec3_t(ixa,iya,iza),Ndim); 
                    //if (na < nc)
                    if (!IsSolid [0][ixa][iya][iza])
                    {
                     double Distance = pow(((double)ixa-(double)ix)*((double)ixa-(double)ix)+((double)iya-(double)iy)*((double)iya-(double)iy)+((double)iza-(double)iz)*((double)iza-(double)iz),0.5);
                     double Weight   = exp(-Distance/Elta);

                       if (Weight > 0.0002)                                                                // --- interaction distance is  8lu for elta=1.0.                     --- 20161216
                     //if (Weight > 1.0e-7)                                                                // --- interaction distance is 16lu for elta=1.0.                     --- 20161216
                       {
                          CPPaT[omp_get_thread_num()].Push(iVec3_t(nc,na,k));                              // --- I do not understand k in iVec3_t(nc,na,k) now                  --- 20161125 
                       }  
                    }  
                    k++;
                }
            }
        }

	Array<iVec3_t> CPPa(0);    
	for (size_t nt=0;nt<Nproc;nt++)
        //for (size_t nt=0;nt<1 ;nt++)
	for (size_t it=0;it<CPPaT[nt].Size();it++)
	{
	    CPPa.Push(CPPaT[nt][it]);
	}

        NCellPaira = CPPa.Size();                                                                          // --- the particle pair in Gw-Elta force setup                       --- 20161125 
        CellPaira  = new iVec3_t [NCellPaira];

        for (size_t n=0;n<NCellPaira;n++)
        {
            CellPaira[n] = CPPa[n];
        }
    }

//************************************---p8---*********************************
     
    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);                                                 //externally specified BC is imported                                     --- 20170302 
        if (Time >= tout)       //The LBM dynamics
        {

//************************************---p9---*********************************

            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);                                             //idx_out: the number of output frames                                      --- 20170301
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    WriteXDMF(fn.CStr());
                    //#else
                    //WriteVTK (fn.CStr());
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }   
        
//************************************---p10---*********************************

        #ifdef USE_OMP
        if (Nl==1)                                                                                         //Nl=1 or 2 to identify single-component or two-component model               --- 20161125 
        {
            //if (fabs(G[0])>1.0e-12)         
            //ApplyForcesSC();

            //CollideSC();
            //StreamSC();

            ApplyForces_contact();
            CollideEFS_SRT();
            StreamEFS_SRT();
        }
        #endif

        Time += dt;
        //std::cout << Time << std::endl;
    }
}
//==================================== function 5 -- Solve ================================================










//==================================== function 6 -- ApplyForcesSC ========================================
inline void Domain::ApplyForcesSC()
{
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t n=0;n<NCellPairs;n++)
    {
        iVec3_t idxc,idxn;
        size_t k = CellPairs[n](2);                                                                        // CellPairs[n](2) --- k is extracted                                        --- 20161125 
        idx2Pt(CellPairs[n](0),idxc,Ndim);
        idx2Pt(CellPairs[n](1),idxn,Ndim);                                                               

        size_t ixc = idxc(0);
        size_t iyc = idxc(1);
        size_t izc = idxc(2);
        size_t ixn = idxn(0);
        size_t iyn = idxn(1);
        size_t izn = idxn(2);

        double psic = 0.0;
        double psin = 0.0;
        double Gt   = G[0]; 
 
        IsSolid[0][ixc][iyc][izc] ? psic = 0.0 : psic = Psi[0]*exp(-Rhoref[0]/Rho[0][ixc][iyc][izc]);
        IsSolid[0][ixn][iyn][izn] ? psin = 0.0 : psin = Psi[0]*exp(-Rhoref[0]/Rho[0][ixn][iyn][izn]);

        Vec3_t bforce = -Gt*W[k]*C[k]*psic*psin;                                                           // CellPairs[n](2) --- k is used here                                        --- 20161125 
        BForce[0][ixc][iyc][izc] += bforce;
        BForce[0][ixn][iyn][izn] -= bforce;
    }
}
//==================================== function 6 -- ApplyForcesSC ========================================










//====================================== function 7 -- CollideSC ==========================================
inline void Domain::CollideSC()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[0][ix][iy][iz])
        {
            double NonEq[Nneigh];
            double Q = 0.0;
            double tau = Tau[0];
            double rho = Rho[0][ix][iy][iz];

        //iVec3_t bforce = BForce[0][ix][iy][iz];
        //std::cout  << "---"<< BForce[0][ix][iy][iz] << "---"<<std::endl;

            Vec3_t vel = Vel[0][ix][iy][iz]+dt*tau*BForce[0][ix][iy][iz]/rho;
            double VdotV = dot(vel,vel);
            for (size_t k=0;k<Nneigh;k++)
            {
                double VdotC = dot(vel,C[k]);
                double Feq   = W[k]*rho*(1.0 + 3.0*VdotC/Cs + 4.5*VdotC*VdotC/(Cs*Cs) - 1.5*VdotV/(Cs*Cs));
                NonEq[k] = F[0][ix][iy][iz][k] - Feq;
                Q +=  NonEq[k]*NonEq[k]*EEk[k];
            }
            Q = sqrt(2.0*Q);
            tau = 0.5*(tau+sqrt(tau*tau + 6.0*Q*Sc/rho));

            bool valid = true;
            double alpha = 1.0;
            while (valid)
            {
                valid = false;
                for (size_t k=0;k<Nneigh;k++)
                {
                    Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][k] - alpha*NonEq[k]/tau;
                    if (Ftemp[0][ix][iy][iz][k]<-1.0e-12)
                    {
                        //std::cout << Ftemp[0][ix][iy][iz][k] << std::endl;
                        double temp =  tau*F[0][ix][iy][iz][k]/NonEq[k];
                        if (temp<alpha) alpha = temp;
                        valid = true;
                    }
                    if (std::isnan(Ftemp[0][ix][iy][iz][k]))
                    {
                        std::cout << "CollideSC: Nan found, resetting" << std::endl;
                        std::cout << " " << alpha << " " << iVec3_t(ix,iy,iz) << " " << k << " " << std::endl;
                        throw new Fatal("Domain::CollideSC: Distribution funcitons gave nan value, check parameters");
                    }
                }
            }
        }

        else
        
        //if    (IsSolid[0][ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Ftemp[0][ix][iy][iz][k] = F[0][ix][iy][iz][Op[k]];
            }
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;
}
//====================================== function 7 -- CollideSC ==========================================










//===================================== function 8 -- StreamSC ============================================
inline void Domain::StreamSC()
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        for (size_t k=0;k<Nneigh;k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);                          // periodic BC scheme 
            Ftemp[0][nix][niy][niz][k] = F[0][ix][iy][iz][k];                                              // 
        }
    }

    double ***** tmp = F;
    F = Ftemp;
    Ftemp = tmp;

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        //BForce[0][ix][iy][iz] = OrthoSys::O;
        Vel   [0][ix][iy][iz] = OrthoSys::O;
        Rho   [0][ix][iy][iz] = 0.0;
        if (!IsSolid[0][ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Rho[0][ix][iy][iz] +=  F[0][ix][iy][iz][k];                                                // calculating density
                Vel[0][ix][iy][iz] +=  F[0][ix][iy][iz][k]*C[k];
            }
            Vel[0][ix][iy][iz] /= Rho[0][ix][iy][iz];                                                      // calculating velocity
        }
    }
}
//===================================== function 8 -- StreamSC ============================================










//==================================== function 9 -- ApplyForces_contact =====================================
inline void Domain::ApplyForces_contact() 
{
    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t n=0;n<NCellPaira;n++)                                                                                  
    {
        iVec3_t idxc,idxn;
      //size_t k = CellPaira[n](2);                                                                                    // CellPairs[n](2) --- k (direction) related weight is not used here--- 20161125 
        idx2Pt(CellPaira[n](0),idxc,Ndim);
        idx2Pt(CellPaira[n](1),idxn,Ndim);

        size_t ixc = idxc(0);
        size_t iyc = idxc(1);
        size_t izc = idxc(2);
        size_t ixn = idxn(0);
        size_t iyn = idxn(1);
        size_t izn = idxn(2);

        double psic = 0.0;
        double psin = 0.0;
        double Gt   = G[0]; 

        IsSolidEdge[0][ixc][iyc][izc] ? psic = 1.0, Gt = Gw[0] : psic = Rho[0][ixc][iyc][izc];                         // Ideal fluid --- only solid surface exert a force on fluid      --- 20161216
        IsSolidEdge[0][ixn][iyn][izn] ? psin = 1.0, Gt = Gw[0] : psin = Rho[0][ixn][iyn][izn];

        double  Dist = pow(((double)ixn-(double)ixc)*((double)ixn-(double)ixc)+((double)iyn-(double)iyc)*((double)iyn-(double)iyc)+((double)izn-(double)izc)*((double)izn-(double)izc),0.5);
        double  Wk   = exp(-Dist/Elta);
        Vec3_t  Ck   = {((double)ixn-(double)ixc)/Dist,((double)iyn-(double)iyc)/Dist,((double)izn-(double)izc)/Dist}; // --- 20161101 --- grid independence is awful                    --- 20161220

        Vec3_t bforce = -Gt*Wk*Ck*psic*psin;
        BForce[0][ixc][iyc][izc] += bforce;
        BForce[0][ixn][iyn][izn] -= bforce;
    }
}
//===================================== function 9 -- ApplyForces_contact ====================================










//===================================== function 10 -- CollideEFS_SRT ======================================
inline void Domain::CollideEFS_SRT()                                                                         // explicit forcing scheme --- 20170301
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        if (!IsSolid[0][ix][iy][iz])
        {
            double tau   = Tau[0];
            double rho   = Rho[0][ix][iy][iz];
            Vec3_t vel   = Vel[0][ix][iy][iz];

            for (size_t k=0;k<Nneigh;k++)
            {
                F3temp[0][ix][iy][iz][k] = 1.0/tau*Feq(k,rho,vel) + (1.0-0.5/tau)*dt*F2[0][ix][iy][iz][k] + (1.0-1.0/tau)*F3[0][ix][iy][iz][k];   
            }
        }

        else
        //if    (IsSolid[0][ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                F3temp[0][ix][iy][iz][k] = F3[0][ix][iy][iz][Op[k]];                                        // bounce-back BC scheme --- 20170228
            }
        }
    }

    double ***** tmp = F3;
    F3 = F3temp;
    F3temp = tmp;
}
//===================================== function 10 -- CollideEFS_SRT ======================================










//===================================== function 11 -- StreamEFS_SRT =======================================
inline void Domain::StreamEFS_SRT()      
{
    size_t nx = Ndim(0);
    size_t ny = Ndim(1);
    size_t nz = Ndim(2);

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        for (size_t k=0;k<Nneigh;k++)
        {
            size_t nix = (size_t)((int)ix + (int)C[k](0) + (int)Ndim(0))%Ndim(0);
            size_t niy = (size_t)((int)iy + (int)C[k](1) + (int)Ndim(1))%Ndim(1);
            size_t niz = (size_t)((int)iz + (int)C[k](2) + (int)Ndim(2))%Ndim(2);                          // periodic BC scheme 
            F3temp[0][nix][niy][niz][k] = F3[0][ix][iy][iz][k];                                            //
        }
    }

    double ***** tmp = F3;
    F3  = F3temp;
    F3temp = tmp;

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vel   [0][ix][iy][iz] = OrthoSys::O;                                                               // in fact the solid nodes are specified zero values. Good.   --- 20170302 
        Rho   [0][ix][iy][iz] = 0.0;
        if (!IsSolid[0][ix][iy][iz])
        {
            for (size_t k=0;k<Nneigh;k++)
            {
                Rho[0][ix][iy][iz] +=  F3[0][ix][iy][iz][k];                                               // calculating density
                Vel[0][ix][iy][iz] +=  F3[0][ix][iy][iz][k]*C[k];
            }
            Vel[0][ix][iy][iz] += 0.5*dt*BForce[0][ix][iy][iz];  // calculating velocity
            Vel[0][ix][iy][iz] /= Rho[0][ix][iy][iz];
        }
    }

    #ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    #endif
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        double rho   = Rho[0][ix][iy][iz];
        Vec3_t vel   = Vel[0][ix][iy][iz];

        for (size_t k=0;k<Nneigh;k++)
        {
            F2 [0][ix][iy][iz][k] = 3.0/rho * Feq(k,rho,vel) * dot(BForce[0][ix][iy][iz],(C[k]-vel));
        }
    }
}
//===================================== function 11 -- StreamEFS_SRT =======================================










}                                                                                                          // namespace FLBM
//###################################### namespace FLBM ###################################################

#endif                                                                                                     // #ifndef MECHSYS_FLBM_DOMAIN_H 
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& MECHSYS_FLBM &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
