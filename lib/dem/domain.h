/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
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

/** @file dem/domain.h .*/

#ifndef MECHSYS_DEM_DOMAIN_H
#define MECHSYS_DEM_DOMAIN_H

// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>
#include <list>
#include <utility> // for std::pair

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// Voro++
#include "src/voro++.cc" // old
/*
#include "voro++.hh" // new
using namespace voro;
using namespace std;
*/

// VTK
#ifdef USE_VTK
#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkLineSource.h>
#include <vtkPolyData.h>
#include <vtkVersion.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#endif // USE_VTK

// MechSys
#include <mechsys/dem/interacton.h>
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/tree.h>

namespace DEM
{

class Domain
{
public:
    // typedefs
    typedef void (*ptFun_t) (Domain & Dom, void * UserData);

    // Constructor
    Domain(void * UserData=NULL);

    // Destructor
    ~Domain();

    // Particle generation
    void GenSpheres      (int Tag, double L, size_t N, double rho, char const * Type,
                          size_t Randomseed, double fraction, double RminFraction = 1.0);                                        ///< General spheres
    void GenSpheresBox (int Tag, Vec3_t const & X0, Vec3_t const & X1,                                                           ///< Generate spheres within a rectangular box defined by the vectors X0 and X1
                        double R, double rho, char const * Type, size_t Randomseed, double fraction, double RminFraction);
    void GenRice         (int Tag, double L, size_t N, double R, double rho, size_t Randomseed, double fraction);                ///< General rices
    void GenBox          (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, bool Cohesion=false);            ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenOpenBox      (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf);                                 ///< Generate five walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBoundingBox  (int InitialTag, double R, double Cf,bool Cohesion=false);                                              ///< Generate o bounding box enclosing the previous included particles.
    void GenBoundingPlane(int InitialTag, double R, double Cf,bool Cohesion=false);                                              ///< Same as GenBounding but only generates one pair of planes.
    void GenFromMesh     (Mesh::Generic & M, double R, double rho, bool cohesion=false, bool MC=true, double thickness = 0.0);   ///< Generate particles from a FEM mesh generator
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bool Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                        ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni
    // Single particle addition
    void AddSphere   (int Tag, Vec3_t const & X, double R, double rho);                                                          ///< Add sphere
    void AddCube     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddRecBox   (int Tag, Vec3_t const & X, Vec3_t const & L, double R, double rho, double Angle=0, Vec3_t * Axis=NULL);    ///< Add a rectangular box with dimensions given by the vector L
    void AddTetra    (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a tetrahedron at position X with spheroradius R, side of length L and density rho
    void AddDrill    (int Tag, Vec3_t const & X, double R, double Lt, double Ll, double rho);                                    ///< A drill made as a combination of a cube and a pyramid.
    void AddRice     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a rice at position X with spheroradius R, side of length L and density rho
    void AddPlane    (int Tag, Vec3_t const & X, double R, double Lx,double Ly, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddVoroCell (int Tag, voronoicell_neighbor & VC, double R, double rho, bool Erode, Vec3_t nv = iVec3_t(1.0,1.0,1.0));   ///< Add a single voronoi cell, it should be built before tough
    void AddTorus    (int Tag, Vec3_t const & X, Vec3_t const & N, double Rmax, double R, double rho);                           ///< Add a single torus at position X with a normal N, circunference Rmax and spheroradius R
    void AddCylinder (int Tag, Vec3_t const & X0, double R0, Vec3_t const & X1, double R1, double R, double rho);                ///< Add a cylinder formed by the connection of two circles at positions X0 and X1 and radii R0 and R1


    // 
    //void AddParticle (DEM::Particle * Pa);                                                                                       ///< Add a particle as an exact copy of particle Pa

    // Methods
    void SetProps          (Dict & D);                                                                          ///< Set the properties of individual grains by dictionaries
    void Initialize        (double dt=0.0);                                                                     ///< Set the particles to a initial state and asign the possible insteractions
    void Solve             (double tf, double dt, double dtOut, ptFun_t ptSetup=NULL, ptFun_t ptReport=NULL,
                            char const * FileKey=NULL, size_t VOut=3, size_t Nproc=1,double minEkin=0.0);       ///< Run simulation the simulation up to time tf, with dt and dtOut the time and report steps. The funstion Setup and Report are used to control the workflow form outside, filekey is used to name the report files. VOut has the options 0 no visualization, 1 povray, 2 xmdf and 3 both. minEkin is a minimun of kinetic energy before the simulation stops
    void WritePOV          (char const * FileKey);                                                              ///< Write POV file
    void WriteBPY          (char const * FileKey);                                                              ///< Write BPY (Blender) file
#ifdef USE_HDF5    
    void WriteBF           (char const * FileKey);                                                              ///< Save a h5 with branch and force information
    void WriteXDMF         (char const * FileKey);                                                              ///< Save a xdmf file for visualization
    void Save              (char const * FileKey);                                                              ///< Save the current domain
    void Load              (char const * FileKey);                                                              ///< Load the domain form a file
#endif

#ifdef USE_VTK
    void WriteVTKContacts  (char const * FileKey);                                                              ///< Save a vtk - vtp file for conatcs visualization
#endif // USE_VTK

#ifdef USE_THREAD
    void UpdateLinkedCells ();                                                                                  ///< Update the linked cells
#endif
    void BoundingBox       (Vec3_t & minX, Vec3_t & maxX);                                                      ///< Defines the rectangular box that encloses the particles.
    void Center            (Vec3_t C = Vec3_t(0.0,0.0,0.0));                                                    ///< Centers the domain around C
    void ClearInteractons  ();                                                                                  ///< Reset the interactons
    void ResetInteractons  ();                                                                                  ///< Reset the interactons
    void ResetDisplacements();                                                                                  ///< Reset the displacements
    double MaxDisplacement ();                                                                                  ///< Calculate maximun displacement
    void ResetContacts     ();                                                                                  ///< Reset the displacements
    void ResetBoundaries   ();                                                                                  ///< Reset the Boundary particles
    void EnergyOutput      (size_t IdxOut, std::ostream & OutFile);                                             ///< Output of the energy variables
    void GetGSD            (Array<double> & X, Array<double> & Y, Array<double> & D, size_t NDiv=10) const;     ///< Get the Grain Size Distribution
    void Clusters          ();                                                                                  ///< Check the bounded particles in the domain and how many connected clusters are still present
    void DelParticles      (Array<int> const & Tags);                                                           ///< Delete particle

    // Access methods
    Particle       * GetParticle  (int Tag, bool Check=true);       ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    Particle const & GetParticle  (int Tag, bool Check=true) const; ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    void             GetParticles (int Tag, Array<Particle*> & P);  ///< Find all particles with Tag

    // Auxiliar methods
    void   LinearMomentum  (Vec3_t & L);                    ///< Return total momentum of the system
    void   AngularMomentum (Vec3_t & L);                    ///< Return total angular momentum of the system
    double CalcEnergy      (double & Ekin, double & Epot);  ///< Return total energy of the system

#ifdef USE_THREAD
    pthread_mutex_t lck;                                                           ///< to protect variables in multithreading
    Array<pair<size_t, size_t> >                      ListPosPairs;                ///< List of all possible particles pairs
    iVec3_t                                           LCellDim;                    ///< Dimensions of the linked cell array
    Array<Array <size_t> >                            LinkedCell;                  ///< Linked Cell array for optimization.
    Vec3_t                                            LCxmin;                      ///< Bounding box low   limit for the linked cell array
    Vec3_t                                            LCxmax;                      ///< Bounding box upper limit for the linked cell array
#endif

    // Data
    bool                                              Initialized;                 ///< System (particles and interactons) initialized ?
    bool                                              Finished;                    ///< Has the simulation finished
    bool                                              Dilate;                      ///< True if eroded particles should be dilated for visualization
    Array<size_t>                                     FreePar;                     ///< Particles that are free
    Array<size_t>                                     NoFreePar;                   ///< Particles that are not free
    Array<Particle*>                                  Particles;                   ///< All particles in domain
    Array<Particle*>                                  ParXmax;                     ///< Particles that are on the Xmax boundary for periodic boudary conditions
    Array<Interacton*>                                Interactons;                 ///< All interactons
    Array<CInteracton*>                               CInteractons;                ///< Contact interactons
    Array<BInteracton*>                               BInteractons;                ///< Cohesion interactons
    Array<Interacton*>                                PInteractons;                ///< Interactons for periodic conditions
    Array<CInteracton*>                               CPInteractons;               ///< Contact interacton for periodic conditions
    Vec3_t                                            CamPos;                      ///< Camera position for POV
    double                                            Time;                        ///< Current time
    double                                            Dt;                          ///< Time step
    double                                            Evis;                        ///< Energy dissipated by the viscosity of the grains
    double                                            Efric;                       ///< Energy dissipated by friction
    double                                            Wext;                        ///< Work done by external forces
    double                                            Vs;                          ///< Volume occupied by the grains
    double                                            Ms;                          ///< Total mass of the particles
    double                                            Alpha;                       ///< Verlet distance
    double                                            Beta;                        ///< Binmultiplier
    double                                            Xmax;                        ///< Maximun distance along the X axis (Periodic Boundary)
    double                                            Xmin;                        ///< Minimun distance along the X axis (Periodic Boundary)
    double                                            MaxDmax;                     ///< Maximun value for the radious of the spheres surronding each particle
    void *                                            UserData;                    ///< Some user data
    String                                            FileKey;                     ///< File Key for output files
    size_t                                            idx_out;                     ///< Index of output
    set<pair<Particle *, Particle *> >                Listofpairs;                 ///< List of pair of particles associated per interacton for memory optimization
    set<pair<Particle *, Particle *> >                PListofpairs;                ///< List of pair of particles associated per interacton for memory optimization under periodic boundary conditions
    Array<Array <int> >                               Listofclusters;              ///< List of particles belonging to bounded clusters (applies only for cohesion simulations)
    
    
#ifdef USE_BOOST_PYTHON
    void PyAddSphere (int Tag, BPy::tuple const & X, double R, double rho)                                                         { AddSphere (Tag,Tup2Vec3(X),R,rho); }
    void PyAddCube   (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddCube  (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddTetra  (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddTetra (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddRice   (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddRice  (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddPlane  (int Tag, BPy::tuple const & X, double R, double Lx,double Ly, double rho, double Ang, BPy::tuple const & Ax) { Vec3_t a(Tup2Vec3(Ax)); AddPlane (Tag,Tup2Vec3(X),R,Lx,Ly,rho,Ang,&a); }
    void PySetCamPos (BPy::tuple const & PyCamPos)                                                                                 { CamPos = Tup2Vec3(PyCamPos); }
    void PyGetParticles(BPy::list & P)
    {
        for (size_t i=0; i<Particles.Size(); ++i)
        {
            BPy::list p,V,E,F;
            double radius = Particles[i]->PyGetFeatures (V, E, F);
            p.append (radius);
            p.append (V);
            p.append (E);
            p.append (F);
            P.append (p);
        }
    }
    void PyGetGSD (BPy::list & X, BPy::list & Y, BPy::list & D, int NDiv=10)
    {
        Array<double> x, y, d;
        GetGSD (x, y, d, NDiv);
        for (size_t i=0; i<x.Size(); ++i)
        {
            X.append (x[i]);
            Y.append (y[i]);
        }
    }
#endif

};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

#ifdef USE_THREAD

struct MtData   /// A structure for the multi-tread data
{
    size_t                  ProcRank; ///< Rank of the thread
    size_t                    N_Proc; ///< Total number of threads
    DEM::Domain *                Dom; ///< Pointer to the lbm domain
    double                       Dmx; ///< Maximun displacement
    Array<pair<size_t,size_t> >   LC; ///< A temporal list of new contacts
    Array<size_t>                LCI; ///< A temporal array of posible Cinteractions
    Array<size_t>                LCB; ///< A temporal array of posible Binteractions
    Array<pair<size_t,size_t> >  LPC; ///< A temporal list of new contacts for periodic boundary conditions
    Array<size_t>               LPCI; ///< A temporal array of posible Cinteractions for periodic boundary conditions
    Array<size_t>                LBP; ///< A temporal array of possible boundary particles
    Array<pair<iVec3_t,size_t> > LLC; ///< A temporal array of possible linked cells locations
};

void * GlobalIni(void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
	for (size_t i=In;i<Fn;i++)
	{
        // set the force and torque to the fixed values
        (*P)[i]->F = (*P)[i]->Ff;
        (*P)[i]->T = (*P)[i]->Tf;
        for (size_t n=0;n<3;n++)
        {
            for (size_t m=0;m<3;m++)  
            {
                (*P)[i]->M(n,m)=0.0;
                (*P)[i]->B(n,m)=0.0;
            }
        }

        // initialize the coordination (number of contacts per particle) number and the Bdry flag
        (*P)[i]->Comp = 0.0;
        (*P)[i]->Cn   = 0.0;
        (*P)[i]->Bdry = false;
    }
    return NULL;
}

void * GlobalForce(void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Interacton * > * I = &dat.Dom->Interactons;
	size_t Ni = I->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = I->Size() : Fn = (dat.ProcRank+1)*Ni;
	for (size_t i=In;i<Fn;i++)
	{
		if ((*I)[i]->CalcForce(dat.Dom->Dt))
        {
            dat.Dom->Save     ("error");
            dat.Dom->WriteXDMF("error");
            std::cout << "Maximun overlap detected between particles at time " << dat.Dom->Time << std::endl;
            sleep(1);
            throw new Fatal("Maximun overlap detected between particles");
        }
	}
    return NULL;
}

void * GlobalMove(void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.Dmx = 0.0;
	for (size_t i=In;i<Fn;i++)
	{
        //std::cout << "1" << std::endl;
		(*P)[i]->Translate(dat.Dom->Dt);
        //std::cout << "2" << std::endl;
		(*P)[i]->Rotate(dat.Dom->Dt);
        //std::cout << "3" << std::endl;
        if ((*P)[i]->MaxDisplacement()>dat.Dmx) dat.Dmx = (*P)[i]->MaxDisplacement();
	}
    return NULL;
}

void * GlobalResetDisplacement(void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.Dmx = 0.0;
    dat.LLC.Resize(0);
	for (size_t i=In;i<Fn;i++)
    {
        (*P)[i]->ResetDisplacements();
        if ((*P)[i]->IsFree())
        {
            iVec3_t idx = ((*P)[i]->x - dat.Dom->LCxmin)/(2.0*dat.Dom->Beta*dat.Dom->MaxDmax);
            dat.LLC.Push(make_pair(idx,i));
            //if (i==345)
            //{
                //std::cout << idx << " " << (*P)[i]->x << " " << dat.Dom->Time << " " << dat.ProcRank << std::endl;
                //std::cout << dat.LLC.Last().first << " " << dat.LLC.Last().second << std::endl;
            //}
        }
    }
    return NULL;
}

void * GlobalResetContacts1 (void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
	size_t Ni = dat.Dom->ListPosPairs.Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = dat.Dom->ListPosPairs.Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LC.Resize(0);

    //std::cout << dat.Dom->ListPosPairs.Size() << " " << Ni << std::endl;

    for (size_t n=In;n<Fn;n++)
    {
        size_t i = dat.Dom->ListPosPairs[n].first;
        size_t j = dat.Dom->ListPosPairs[n].second;
        bool pi_has_vf = !dat.Dom->Particles[i]->IsFree();
        bool pj_has_vf = !dat.Dom->Particles[j]->IsFree();

        bool close = (Distance(dat.Dom->Particles[i]->x,dat.Dom->Particles[j]->x)<=dat.Dom->Particles[i]->Dmax+dat.Dom->Particles[j]->Dmax+2*dat.Dom->Alpha);
        if ((pi_has_vf && pj_has_vf) || !close) continue;
        //if (!close) continue;
        
        // checking if the interacton exist for that pair of particles
        set<pair<Particle *, Particle *> >::iterator it = dat.Dom->Listofpairs.find(make_pair(dat.Dom->Particles[i],dat.Dom->Particles[j]));
        if (it != dat.Dom->Listofpairs.end())
        {
            continue;
        }
        dat.LC.Push(make_pair(i,j));
    }
    return NULL;
}

void * GlobalResetContacts2 (void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
	size_t Ni = dat.Dom->CInteractons.Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = dat.Dom->CInteractons.Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LCI.Resize(0);
    for (size_t n=In;n<Fn;n++)
    {
        if(dat.Dom->CInteractons[n]->UpdateContacts(dat.Dom->Alpha)) dat.LCI.Push(n);
    }
	Ni = dat.Dom->BInteractons.Size()/dat.N_Proc;
    In = dat.ProcRank*Ni;
    dat.ProcRank == dat.N_Proc-1 ? Fn = dat.Dom->BInteractons.Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LCB.Resize(0);
    for (size_t n=In;n<Fn;n++)
    {
        if(dat.Dom->BInteractons[n]->UpdateContacts(dat.Dom->Alpha)) dat.LCB.Push(n);
    }
    return NULL;
}

void * GlobalResetBoundaries1 (void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LBP.Resize(0);
	for (size_t i=In;i<Fn;i++)
    {
        Particle * Pa = (*P)[i];
        if ((Pa->MinX()>dat.Dom->Xmax)&&Pa->IsFree())
        {
            Vec3_t v(dat.Dom->Xmin-dat.Dom->Xmax,0.0,0.0);
            Pa->Translate(v);
        }
        if ((Pa->MaxX()<dat.Dom->Xmin)&&Pa->IsFree())
        {
            Vec3_t v(dat.Dom->Xmax-dat.Dom->Xmin,0.0,0.0);
            Pa->Translate(v);
        }
        if ((Pa->MaxX()>dat.Dom->Xmax-2.0*dat.Dom->Alpha-2.0*dat.Dom->MaxDmax)&&Pa->IsFree())
        {
            dat.LBP.Push(i);
        }
    }
    return NULL;
}

void * GlobalResetBoundaries2 (void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->ParXmax;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LPC.Resize(0);
    Vec3_t v(dat.Dom->Xmin-dat.Dom->Xmax,0.0,0.0);
	for (size_t i=In;i<Fn;i++)
    {
        Particle * P1 = (*P)[i];
        P1->Translate(v);
        for (size_t j=0; j<dat.Dom->Particles.Size(); j++)
        {
            Particle * P2 = dat.Dom->Particles[j];
            if (P1==P2||!P2->IsFree()||P->Has(P2)) continue;
            bool close = (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*dat.Dom->Alpha);
            if (!close) continue;
            set<pair<Particle *, Particle *> >::iterator it = dat.Dom->PListofpairs.find(make_pair(P1,P2));
            if (it != dat.Dom->PListofpairs.end())
            {
                continue;
            }
            dat.LPC.Push(make_pair(i,j));
        }
    }
    return NULL;
}

void * GlobalResetBoundaries3 (void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
	size_t Ni = dat.Dom->CPInteractons.Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = dat.Dom->CPInteractons.Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LPCI.Resize(0);
    for (size_t n=In;n<Fn;n++)
    {
        if(dat.Dom->CPInteractons[n]->UpdateContacts(dat.Dom->Alpha)) dat.LPCI.Push(n);
    }
    return NULL;
}

void * GlobalPerTranslate(void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->ParXmax;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    Vec3_t v(dat.Dom->Xmin-dat.Dom->Xmax,0.0,0.0);
	for (size_t i=In;i<Fn;i++)
    {
        (*P)[i]->Translate(v);
    }
    return NULL;
}

void * GlobalPerTranslateBack(void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->ParXmax;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    Vec3_t v(dat.Dom->Xmax-dat.Dom->Xmin,0.0,0.0);
	for (size_t i=In;i<Fn;i++)
    {
        (*P)[i]->Translate(v);
    }
    return NULL;
}

void * GlobalPerForce(void * Data)
{
    DEM::MtData & dat = (*static_cast<DEM::MtData *>(Data));
    Array<Interacton * > * I = &dat.Dom->PInteractons;
	size_t Ni = I->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = I->Size() : Fn = (dat.ProcRank+1)*Ni;
	for (size_t i=In;i<Fn;i++)
	{
		if ((*I)[i]->CalcForce(dat.Dom->Dt))
        {
            dat.Dom->Save     ("error");
            dat.Dom->WriteXDMF("error");
            std::cout << "Maximun overlap detected between particles at Time" << dat.Dom->Time << std::endl;
            sleep(1);
            throw new Fatal("Maximun overlap detected between particles");
        }
	}
    return NULL;
}

#endif

// Constructor & Destructor

inline Domain::Domain (void * UD)
    :  Initialized(false), Dilate(false), Time(0.0), Alpha(0.05), Beta(1.0), UserData(UD)
{
    Xmax = Xmin = 0.0;
    CamPos = 1.0, 2.0, 3.0;
#ifdef USE_THREAD
    pthread_mutex_init(&lck,NULL);
#endif
}

inline Domain::~Domain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    for (size_t i=0; i<Interactons.Size(); ++i) if (Interactons[i]!=NULL) delete Interactons[i];
}

// Particle generation

inline void Domain::GenSpheres (int Tag, double L, size_t N, double rho,char const * Type, size_t Randomseed, double fraction, double RminFraction)
{
    // find radius from the edge's length
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating packing of spheres -----------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    srand(Randomseed);
    double R = L/(2.0*N);
    if (strcmp(Type,"Normal")==0)
    {
        for (size_t n=0; n<N*N*N; ++n) 
        {
            Vec3_t pos(-L/2.0+R, -L/2.0+R, -L/2.0+R);
            size_t i = (n%N);
            size_t j = (n/N)%N;
            size_t k = (n/(N*N));
            pos += Vec3_t(2.0*i*R, 2.0*j*R, 2.0*k*R);
            if (rand()<fraction*RAND_MAX) AddSphere (Tag,pos,R*RminFraction+(1.0*rand())/RAND_MAX*(R-R*RminFraction),rho);
        }
    }
    else if (strcmp(Type,"HCP")==0)
    {
        size_t nx = N;
        size_t ny = int(L/(sqrt(3.0)*R));
        size_t nz = int(L/(sqrt(8.0/3.0)*R));
        for (size_t k = 0; k < nz; k++)
        {
            for (size_t j = 0; j < ny; j++)
            {
                Vec3_t X;
                if (k%2==0) X = Vec3_t(-2*R-L/2.0,R-L/2.0,2*R-L/2.0+k*sqrt(8.0/3.0)*R);
                else X = Vec3_t(-R-L/2.0,R+sqrt(1.0/3.0)*R-L/2.0,2*R-L/2.0+k*sqrt(8.0/3.0)*R);
                if (j%2==0) X += Vec3_t(R,j*sqrt(3.0)*R,0.0);
                else X += Vec3_t(0.0,j*sqrt(3.0)*R,0.0);
                for (size_t i = 0; i < nx; i++)
                {
                    X += Vec3_t(2*R,0.0,0.0);
                    if (rand()<fraction*RAND_MAX) AddSphere(Tag,X,R*RminFraction+(1.0*rand())/RAND_MAX*(R-R*RminFraction),rho);
                }
            }
        }
    }
    else throw new Fatal ("Right now there are only two possible packings available the Normal and the HCP, packing %s is not implemented yet",Type);
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::GenSpheresBox (int Tag, Vec3_t const & X0, Vec3_t const & X1, double R, double rho, char const * Type, size_t Randomseed, double fraction, double RminFraction)
{
    // find radius from the edge's length
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating packing of spheres -----------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    srand(Randomseed);
    if (strcmp(Type,"Normal")==0)
    {
        size_t nx = 0.5*(X1(0)-X0(0))/R;
        size_t ny = 0.5*(X1(1)-X0(1))/R;
        size_t nz = 0.5*(X1(2)-X0(2))/R;
        for (size_t i = 0; i < nx; i++)
        for (size_t j = 0; j < ny; j++)
        for (size_t k = 0; k < nz; k++)
        {
            Vec3_t pos(-(X1(0)-X0(0))/2.0+R, -(X1(1)-X0(1))/2.0+R, -(X1(2)-X0(2))/2.0+R);
            pos += Vec3_t(2.0*i*R, 2.0*j*R, 2.0*k*R);
            if (rand()<fraction*RAND_MAX) AddSphere (Tag,pos,R*RminFraction+(1.0*rand())/RAND_MAX*(R-R*RminFraction),rho);
        }
    }

    else if (strcmp(Type,"HCP")==0)
    {
        size_t nx = 0.5*(X1(0)-X0(0))/R-1;
        size_t ny = int((X1(1)-X0(1))/(sqrt(3.0)*R));
        size_t nz = int((X1(2)-X0(2))/(sqrt(8.0/3.0)*R));
        for (size_t k = 0; k < nz; k++)
        {
            for (size_t j = 0; j < ny; j++)
            {
                Vec3_t X;
                if (k%2==0) X = Vec3_t(-R,R,2*R+k*sqrt(8.0/3.0)*R) + X0;
                else X = Vec3_t(0.0,R+sqrt(1.0/3.0)*R,2*R+k*sqrt(8.0/3.0)*R) + X0;
                if (j%2==0) X += Vec3_t(R,j*sqrt(3.0)*R,0.0);
                else X += Vec3_t(0.0,j*sqrt(3.0)*R,0.0);
                for (size_t i = 0; i < nx; i++)
                {
                    X += Vec3_t(2*R,0.0,0.0);
                    //std::cout << X << X0 << X1 << std::endl;
                    if ((X(0)<X0(0))||(X(0)>X1(0))) continue;
                    if ((X(1)<X0(1))||(X(1)>X1(1))) continue;
                    if ((X(2)<X0(2))||(X(2)>X1(2))) continue;
                    if (rand()<fraction*RAND_MAX) AddSphere (Tag,X,R*RminFraction+(1.0*rand())/RAND_MAX*(R-R*RminFraction),rho);
                }
            }
        }
    }
    
    else throw new Fatal ("Right now there are only two possible packings available the Normal and the HCP, packing %s is not implemented yet",Type);
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::GenRice (int Tag, double L, size_t N, double R, double rho, size_t Randomseed, double fraction)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating packing of 'rices' -----------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    srand(Randomseed);
    double dL = L/N;
    for (size_t n=0; n<N*N*N; ++n) 
    {
        Vec3_t pos(-L/2.0+dL, -L/2.0+dL, -L/2.0+dL);
        size_t i = (n%N);
        size_t j = (n/N)%N;
        size_t k = (n/(N*N));
        pos += Vec3_t(2.0*i*dL, 2.0*j*dL, 2.0*k*dL);
        if (rand()<fraction*RAND_MAX) AddRice (Tag, pos, R, dL-2*R, rho);
    }
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::GenBox (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, bool Cohesion)
{
    /*                         +----------------+
     *                       ,'|              ,'|
     *                     ,'  |  ___       ,'  |
     *     z             ,'    |,'4,'  [1],'    |
     *     |           ,'      |~~~     ,'      |
     *    ,+--y      +'===============+'  ,'|   |
     *  x'           |   ,'|   |      |   |2|   |
     *               |   |3|   |      |   |,'   |
     *               |   |,'   +- - - | +- - - -+
     *               |       ,'       |       ,'
     *               |     ,' [0]  ___|     ,'
     *               |   ,'      ,'5,'|   ,'
     *               | ,'        ~~~  | ,'
     *               +----------------+'
     */

    
    // add faces of box
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    size_t IIndex = Particles.Size();  // First face index
    AddPlane (InitialTag,   Vec3_t(Lx/2.0,0.0,0.0),  R, Cf*Lz, Cf*Ly, 1.0, M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-1, Vec3_t(-Lx/2.0,0.0,0.0), R, Cf*Lz, Cf*Ly, 1.0, 3.0*M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-2, Vec3_t(0.0,Ly/2.0,0.0),  R, Cf*Lx, Cf*Lz, 1.0, 3.0*M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-3, Vec3_t(0.0,-Ly/2.0,0.0), R, Cf*Lx, Cf*Lz, 1.0, M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-4, Vec3_t(0.0,0.0,Lz/2.0),  R, Cf*Lx, Cf*Ly, 1.0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-5, Vec3_t(0.0,0.0,-Lz/2.0), R, Cf*Lx, Cf*Ly, 1.0, M_PI, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);

    // define some tolerance for comparissions
    if (Cohesion)
    {
        double tol1 = 1.0e-8;
        double tol2 = 1.0e-3;
        for (size_t i=0;i<IIndex;i++)
        {
            Particle * P1 = Particles[i];
            for (size_t j=IIndex;j<Particles.Size();j++)
            {
                Particle * P2 = Particles[j];
                for (size_t k=0;k<P1->Faces.Size();k++)
                {
                    Face * F1 = P1->Faces[k];
                    Vec3_t n1,c1;
                    F1->Normal  (n1);
                    F1->Centroid(c1);
                    Face * F2 = P2->Faces[0];
                    Vec3_t n2,c2;
                    F2->Normal  (n2);
                    F2->Centroid(c2);
                    Vec3_t n = 0.5*(n1-n2);
                    n/=norm(n);
                    if ((fabs(dot(n1,n2)+1.0)<tol1)
                       &&(fabs(dot(c2-c1,n)-2*R)<tol2))
                    {
                        BInteractons.Push(new BInteracton(P1,P2,k,1));
                        break;
                    }
                }
            }        
        }
    }
}

inline void Domain::GenOpenBox (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf)
{
    /*                         +----------------+
     *                       ,'|              ,'|
     *                     ,'  |  ___       ,'  |
     *     z             ,'    |,'N,'  [1],'    |
     *     |           ,'      |~~~     ,'      |
     *    ,+--y      +'===============+'  ,'|   |
     *  x'           |   ,'|   |      |   |2|   |
     *               |   |3|   |      |   |,'   |
     *               |   |,'   +- - - | +- - - -+
     *               |       ,'       |       ,'
     *               |     ,' [0]  ___|     ,'
     *               |   ,'      ,'4,'|   ,'
     *               | ,'        ~~~  | ,'
     *               +----------------+'
     */


    // Creates an open box without the top lid, acts as a container
    
    // add faces of box
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    AddPlane (InitialTag,   Vec3_t(Lx/2.0,0.0,0.0),  R, Cf*Lz, Cf*Ly, 1.0, M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-1, Vec3_t(-Lx/2.0,0.0,0.0), R, Cf*Lz, Cf*Ly, 1.0, 3.0*M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-2, Vec3_t(0.0,Ly/2.0,0.0),  R, Cf*Lx, Cf*Lz, 1.0, 3.0*M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-3, Vec3_t(0.0,-Ly/2.0,0.0), R, Cf*Lx, Cf*Lz, 1.0, M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-4, Vec3_t(0.0,0.0,-Lz/2.0), R, Cf*Lx, Cf*Ly, 1.0, M_PI, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
}

inline void Domain::GenBoundingBox (int InitialTag, double R, double Cf,bool Cohesion)
{
    Center();
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    GenBox(InitialTag, maxX(0)-minX(0)+2*R, maxX(1)-minX(1)+2*R, maxX(2)-minX(2)+2*R, R, Cf,Cohesion);
}

inline void Domain::GenBoundingPlane (int InitialTag, double R, double Cf,bool Cohesion)
{
    Center();
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    size_t IIndex = Particles.Size();  // First face index
    double Lx = maxX(0)-minX(0)+2*R;
    double Ly = maxX(1)-minX(1)+2*R;
    double Lz = maxX(2)-minX(2)+2*R;
    AddPlane (InitialTag  , Vec3_t(0.0,Ly/2.0,0.0),  R, Cf*Lx, Cf*Lz, 1.0, 3.0*M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-1, Vec3_t(0.0,-Ly/2.0,0.0), R, Cf*Lx, Cf*Lz, 1.0, M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);

    // define some tolerance for comparissions
    if (Cohesion)
    {
        double tol1 = 1.0e-8;
        double tol2 = 1.0e-3;
        for (size_t i=0;i<IIndex;i++)
        {
            Particle * P1 = Particles[i];
            for (size_t j=IIndex;j<Particles.Size();j++)
            {
                Particle * P2 = Particles[j];
                for (size_t k=0;k<P1->Faces.Size();k++)
                {
                    Face * F1 = P1->Faces[k];
                    Vec3_t n1,c1;
                    F1->Normal  (n1);
                    F1->Centroid(c1);
                    Face * F2 = P2->Faces[0];
                    Vec3_t n2,c2;
                    F2->Normal  (n2);
                    F2->Centroid(c2);
                    Vec3_t n = 0.5*(n1-n2);
                    n/=norm(n);
                    if ((fabs(dot(n1,n2)+1.0)<tol1)
                       &&(fabs(dot(c2-c1,n)-2*R)<tol2))
                    {
                        BInteractons.Push(new BInteracton(P1,P2,k,1));
                        break;
                    }
                }
            }        
        }
    }
}

inline void Domain::GenFromMesh (Mesh::Generic & M, double R, double rho, bool Cohesion, bool MC, double thickness)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating particles from mesh ----------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    size_t IIndex = Particles.Size();

    for (size_t i=0; i<M.Cells.Size(); ++i)
    {

        Array<Vec3_t> V;             // Array of vertices
        Array<Array <int> > E;       // Array of edges
        Array<Array <int> > F;       // array of faces
        if (M.NDim==3)
        {
            if (thickness > 0.0) throw new Fatal("Domain::GenFromMesh: Thickness should not be used in a 3D mesh");
            Array<Mesh::Vertex*> const & verts = M.Cells[i]->V;
            size_t nverts = verts.Size();

            // verts
            V.Resize(nverts);
            for (size_t j=0; j<nverts; ++j)
            {
                V[j] = verts[j]->C;
            }

            // edges
            size_t nedges = Mesh::NVertsToNEdges3D[nverts];
            E.Resize(nedges);
            for (size_t j=0; j<nedges; ++j)
            {
                E[j].Push (Mesh::NVertsToEdge3D[nverts][j][0]);
                E[j].Push (Mesh::NVertsToEdge3D[nverts][j][1]);
            }

            size_t nfaces = Mesh::NVertsToNFaces3D[nverts];
            size_t nvperf = Mesh::NVertsToNVertsPerFace3D[nverts];
            F.Resize(nfaces);
            for (size_t j=0; j<nfaces; ++j)
            {
                for (size_t k=0; k<nvperf; ++k)
                {
                    // TODO: check if face is planar or not
                    F[j].Push(Mesh::NVertsToFace3D[nverts][j][k]);
                }
            }
        }
        else if (M.NDim==2)
        {
            if (thickness <= 0.0) throw new Fatal("Domain::GenFromMesh: Thickness should be positive in a 2D mesh");
            Array<Mesh::Vertex*> const & verts = M.Cells[i]->V;
            size_t nverts = verts.Size();
            V.Resize(2*nverts);
            for (size_t j=0; j<nverts; ++j)
            {
                V[j] = verts[j]->C;
                V[j+nverts] = verts[j]->C + Vec3_t(0.0,0.0,thickness);
            }
            size_t nedges = 3*nverts;
            E.Resize(nedges);
            for (size_t j=0; j<nverts; ++j)
            {
                E[j].Push (Mesh::NVertsToEdge2D[nverts][j][0]);
                E[j].Push (Mesh::NVertsToEdge2D[nverts][j][1]);
                E[j+nverts].Push (Mesh::NVertsToEdge2D[nverts][j][0]+nverts);
                E[j+nverts].Push (Mesh::NVertsToEdge2D[nverts][j][1]+nverts);
                E[j+2*nverts].Push(j);
                E[j+2*nverts].Push(j+nverts);
            }
            size_t nfaces = nverts+2;
            F.Resize(nfaces);
            for (size_t j=0; j<nverts; ++j)
            {
                F[j].Push (Mesh::NVertsToEdge2D[nverts][j][0]);
                F[j].Push (Mesh::NVertsToEdge2D[nverts][j][1]);
                F[j].Push (Mesh::NVertsToEdge2D[nverts][j][1]+nverts);
                F[j].Push (Mesh::NVertsToEdge2D[nverts][j][0]+nverts);
                F[nverts].Push(nverts-1-j);
                F[nverts+1].Push(j+nverts);
            }
        }

        double vol; // volume of the polyhedron
        Vec3_t CM;  // Center of mass of the polyhedron
        Mat3_t It;  // Inertia tensor of the polyhedron
        PolyhedraMP(V,F,vol,CM,It);
        Erosion(V,E,F,R);

        // add particle
        Particles.Push (new Particle(M.Cells[i]->Tag, V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
        Particles[Particles.Size()-1]->Eroded = true;
        Particles[Particles.Size()-1]->Index = Particles.Size()-1;
        if (!MC)
        {
            Particles[Particles.Size()-1]->x       = CM;
            Particles[Particles.Size()-1]->Props.V = vol;
            Particles[Particles.Size()-1]->Props.m = vol*rho;
            Vec3_t I;
            Quaternion_t Q;
            Vec3_t xp,yp,zp;
            Eig(It,I,xp,yp,zp);
            CheckDestroGiro(xp,yp,zp);
            I *= rho;
            Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
            Q(1) = (yp(2)-zp(1))/(4*Q(0));
            Q(2) = (zp(0)-xp(2))/(4*Q(0));
            Q(3) = (xp(1)-yp(0))/(4*Q(0));
            Q = Q/norm(Q);
            Particles[Particles.Size()-1]->I     = I;
            Particles[Particles.Size()-1]->Q     = Q;
            double Dmax = Distance(CM,V[0])+R;
            for (size_t i=1; i<V.Size(); ++i)
            {
                if (Distance(CM,V[i])+R > Dmax) Dmax = Distance(CM,V[i])+R;
            }
            Particles[Particles.Size()-1]->Ekin = 0.0;
            Particles[Particles.Size()-1]->Erot = 0.0;
            Particles[Particles.Size()-1]->Dmax  = Dmax;
            Particles[Particles.Size()-1]->PropsReady = true;
        }
    }

    Array<Array <int> > Neigh(Particles.Size()-IIndex);
    Array<Array <int> > FNeigh(Particles.Size()-IIndex);
    if(Cohesion)
    {
        M.FindNeigh();
        //std::cout << M;
        for (size_t i=0; i<M.Cells.Size(); ++i) 
        {
            for (Mesh::Neighs_t::const_iterator p=M.Cells[i]->Neighs.begin(); p!=M.Cells[i]->Neighs.end(); ++p)
            {
                Neigh[i].Push(p->second.second->ID);
                FNeigh[i].Push(p->second.first);
            }           
        }
        for (size_t i=0; i<Neigh.Size(); ++i)
        {
            for (size_t j=0; j<Neigh[i].Size(); ++j)
            {
                size_t index = Neigh[Neigh[i][j]].Find(i);
                if ((size_t)Neigh[i][j]>i) BInteractons.Push(new BInteracton(Particles[i+IIndex],Particles[Neigh[i][j]+IIndex],FNeigh[i][j],FNeigh[Neigh[i][j]][index]));
            }
        }
        
    }

    // info
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::AddVoroPack (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho
                                 , bool Cohesion, bool Periodic,size_t Randomseed, double fraction, Vec3_t qin)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Adding Voronoi particles packing --------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    srand(Randomseed);
    //const double x_min=-(nx*Lx/2.0), x_max=nx*Lx/2.0;
    //const double y_min=-(ny*Ly/2.0), y_max=ny*Ly/2.0;
    //const double z_min=-(nz*Lz/2.0), z_max=nz*Lz/2.0;
    const double x_min=-(nx/2.0), x_max=nx/2.0;
    const double y_min=-(ny/2.0), y_max=ny/2.0;
    const double z_min=-(nz/2.0), z_max=nz/2.0;
    container con(x_min,x_max,y_min,y_max,z_min,z_max,nx,ny,nz, Periodic,Periodic,Periodic,8);
    int n = 0;
    for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = x_min+(i+0.5*qin(0)+(1-qin(0))*double(rand())/RAND_MAX)*(x_max-x_min)/nx;
                double y = y_min+(j+0.5*qin(1)+(1-qin(1))*double(rand())/RAND_MAX)*(y_max-y_min)/ny;
                double z = z_min+(k+0.5*qin(2)+(1-qin(2))*double(rand())/RAND_MAX)*(z_max-z_min)/nz;
                con.put (n,x,y,z);
                n++;
            }
        }
    }

    double x,y,z,px,py,pz;
    container *cp = & con;
    voropp_loop l1(cp);
    //c_loop_all l1(con); // new: trying to get new voro working
    int q,s;
    voronoicell_neighbor c;
    s=l1.init(con.ax,con.bx,con.ay,con.by,con.az,con.bz,px,py,pz);

    Array<Array <size_t> > ListBpairs(n);
    size_t IIndex = Particles.Size();
    do 
    {
        for(q=0;q<con.co[s];q++) 
        {
            x=con.p[s][con.sz*q]+px;y=con.p[s][con.sz*q+1]+py;z=con.p[s][con.sz*q+2]+pz;
            if(x>con.ax&&x<con.bx&&y>con.ay&&y<con.by&&z>con.az&&z<con.bz) 
            {
                if(con.compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z)) 
                {

                    if (rand()<fraction*RAND_MAX)
                    {
                        AddVoroCell(Tag,c,R,rho,true,Vec3_t(Lx/nx,Ly/ny,Lz/nz));
                        Vec3_t trans(Lx*x/nx,Ly*y/ny,Lz*z/nz);
                        Particle * P = Particles[Particles.Size()-1];
                        P->Translate(trans);
                    }
                }
            }
        }
    } while((s=l1.inc(px,py,pz))!=-1);

    // info
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);


    if (Cohesion)
    {
        //if (fraction<1.0) throw new Fatal("Domain::AddVoroPack: With the Cohesion all particles should be considered, plese change the fraction to 1.0");

        // define some tolerance for comparissions
        double tol1 = 1.0e-8;
        double tol2 = 1.0e-3;
        for (size_t i=IIndex;i<Particles.Size()-1;i++)
        {
            Particle * P1 = Particles[i];
            for (size_t j=i+1;j<Particles.Size();j++)
            {
                Particle * P2 = Particles[j];
                if (Distance(P1->x,P2->x)<P1->Dmax+P2->Dmax)
                {
                    for (size_t k=0;k<P1->Faces.Size();k++)
                    {
                        Face * F1 = P1->Faces[k];
                        Vec3_t n1,c1;
                        F1->Normal  (n1);
                        F1->Centroid(c1);
                        bool found = false;
                        for (size_t l=0;l<P2->Faces.Size();l++)
                        {
                            Face * F2 = P2->Faces[l];
                            Vec3_t n2,c2;
                            F2->Normal  (n2);
                            F2->Centroid(c2);
                            Vec3_t n = 0.5*(n1-n2);
                            n/=norm(n);
                            if ((fabs(dot(n1,n2)+1.0)<tol1)
                               &&(fabs(Distance(c1,*F2)-2*R)<tol2)
                               &&(fabs(Distance(c2,*F1)-2*R)<tol2))
                            {
                                BInteractons.Push(new BInteracton(P1,P2,k,l));
                                found = true;
                                break;
                            }
                        }
                        if (found) break;
                    }
                }
            }
        }
    }
}

// Sihgle particle addition

inline void Domain::AddSphere (int Tag,Vec3_t const & X, double R, double rho)
{
    // vertices
    Array<Vec3_t> V(1);
    V[0] = X;

    // edges
    Array<Array <int> > E(0); // no edges

    // faces
    Array<Array <int> > F(0); // no faces

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    Quaternion_t q;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = (4.0/3.0)*M_PI*R*R*R;
    Particles[Particles.Size()-1]->Props.m    = rho*(4.0/3.0)*M_PI*R*R*R;
    Particles[Particles.Size()-1]->I          = (2.0/5.0)*Particles[Particles.Size()-1]->Props.m*R*R;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;

}

inline void Domain::AddCube (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(8);
    double l = L/2.0;
    V[0] = -l, -l, -l;
    V[1] =  l, -l, -l;
    V[2] =  l,  l, -l;
    V[3] = -l,  l, -l;
    V[4] = -l, -l,  l;
    V[5] =  l, -l,  l;
    V[6] =  l,  l,  l;
    V[7] = -l,  l,  l;

    // edges
    Array<Array <int> > E(12);
    for (size_t i=0; i<12; ++i) E[i].Resize(2);
    E[ 0] = 0, 1;
    E[ 1] = 1, 2;
    E[ 2] = 2, 3;
    E[ 3] = 3, 0;
    E[ 4] = 4, 5;
    E[ 5] = 5, 6;
    E[ 6] = 6, 7;
    E[ 7] = 7, 4;
    E[ 8] = 0, 4;
    E[ 9] = 1, 5;
    E[10] = 2, 6;
    E[11] = 3, 7;

    // faces
    Array<Array <int> > F(6);
    for (size_t i=0; i<6; i++) F[i].Resize(4);
    F[0] = 4, 7, 3, 0;
    F[1] = 1, 2, 6, 5;
    F[2] = 0, 1, 5, 4;
    F[3] = 2, 3, 7, 6;
    F[4] = 0, 3, 2, 1;
    F[5] = 4, 5, 6, 7;

    // calculate the rotation
    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    // clean up
    if (!ThereisanAxis) delete Axis;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = L*L*L;
    Particles[Particles.Size()-1]->Props.m    = rho*L*L*L;
    Particles[Particles.Size()-1]->I          = L*L, L*L, L*L;
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m/6.0;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(3.0*L*L/4.0)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    

}

inline void Domain::AddRecBox (int Tag, Vec3_t const & X, Vec3_t const & L, double R, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(8);
    double lx = L(0)/2.0;
    double ly = L(1)/2.0;
    double lz = L(2)/2.0;
    V[0] = -lx, -ly, -lz;
    V[1] =  lx, -ly, -lz;
    V[2] =  lx,  ly, -lz;
    V[3] = -lx,  ly, -lz;
    V[4] = -lx, -ly,  lz;
    V[5] =  lx, -ly,  lz;
    V[6] =  lx,  ly,  lz;
    V[7] = -lx,  ly,  lz;

    // edges
    Array<Array <int> > E(12);
    for (size_t i=0; i<12; ++i) E[i].Resize(2);
    E[ 0] = 0, 1;
    E[ 1] = 1, 2;
    E[ 2] = 2, 3;
    E[ 3] = 3, 0;
    E[ 4] = 4, 5;
    E[ 5] = 5, 6;
    E[ 6] = 6, 7;
    E[ 7] = 7, 4;
    E[ 8] = 0, 4;
    E[ 9] = 1, 5;
    E[10] = 2, 6;
    E[11] = 3, 7;

    // faces
    Array<Array <int> > F(6);
    for (size_t i=0; i<6; i++) F[i].Resize(4);
    F[0] = 4, 7, 3, 0;
    F[1] = 1, 2, 6, 5;
    F[2] = 0, 1, 5, 4;
    F[3] = 2, 3, 7, 6;
    F[4] = 0, 3, 2, 1;
    F[5] = 4, 5, 6, 7;

    // calculate the rotation
    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    // clean up
    if (!ThereisanAxis) delete Axis;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = L(0)*L(1)*L(2);
    Particles[Particles.Size()-1]->Props.m    = rho*L(0)*L(1)*L(2);
    Particles[Particles.Size()-1]->I          = L(1)*L(1) + L(2)*L(2), L(0)*L(0) + L(2)*L(2), L(0)*L(0) + L(1)*L(1);
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m/12.0;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = 0.5*norm(L)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    

}

inline void Domain::AddTetra (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    double sq8 = sqrt(8.0);
    Array<Vec3_t> V(4);
    V[0] =  L/sq8,  L/sq8, L/sq8;
    V[1] = -L/sq8, -L/sq8, L/sq8;
    V[2] = -L/sq8,  L/sq8,-L/sq8;
    V[3] =  L/sq8, -L/sq8,-L/sq8;

    // edges
    Array<Array <int> > E(6);
    for (size_t i=0; i<6; ++i) E[i].Resize(2);
    E[0] = 0, 1;
    E[1] = 1, 2;
    E[2] = 2, 0;
    E[3] = 0, 3;
    E[4] = 1, 3;
    E[5] = 2, 3;

    // face
    Array<Array <int> > F;
    F.Resize(4);
    for (size_t i=0; i<4; ++i) F[i].Resize(3);
    F[0] = 0, 3, 2;
    F[1] = 0, 1, 3;
    F[2] = 0, 2, 1;
    F[3] = 1, 2, 3;

    // calculate the rotation
    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    // clean up
    if (!ThereisanAxis) delete Axis;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = sqrt(2.0)*L*L*L/12.0;
    Particles[Particles.Size()-1]->Props.m    = rho*sqrt(2.0)*L*L*L/12.0;
    Particles[Particles.Size()-1]->I          = L*L, L*L, L*L;
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m/20.0;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(3.0*L*L/8.0)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
}

inline void Domain::AddDrill (int Tag, const Vec3_t & X, double R, double Lt, double Ll, double rho)
{
    Array<Vec3_t> V(9);
    V[0] =  Lt/2.0,  Lt/2.0, Ll/2.0;
    V[1] = -Lt/2.0,  Lt/2.0, Ll/2.0;
    V[2] = -Lt/2.0, -Lt/2.0, Ll/2.0;
    V[3] =  Lt/2.0, -Lt/2.0, Ll/2.0;
    V[4] =  Lt/2.0,  Lt/2.0,    0.0;
    V[5] = -Lt/2.0,  Lt/2.0,    0.0;
    V[6] = -Lt/2.0, -Lt/2.0,    0.0;
    V[7] =  Lt/2.0, -Lt/2.0,    0.0;
    V[8] =     0.0,     0.0,-Ll/2.0;
    
    Array<Array <int> > E(12);
    for (size_t i=0; i<12; ++i) E[i].Resize(2);
    for (size_t i=0;i<4;i++)
    {
        E[i]   = i  , (i+1)%4    ;
        E[i+4] = i+4, (i+5)%4 + 4;
        E[i+8] = i+4, 8;
    }

    Array<Array <int> > F;
    F.Resize(9);
    F[0].Resize(4);
    F[0] = 0, 1, 2, 3;
    F[1].Resize(4);
    F[1] = 0, 4, 5, 1;
    F[2].Resize(4);
    F[2] = 1, 5, 6, 2;
    F[3].Resize(4);
    F[3] = 2, 6, 7, 3;
    F[4].Resize(4);
    F[4] = 3, 7, 4, 0;
    F[5].Resize(3);
    F[5] = 4, 8, 5;
    F[6].Resize(3);
    F[6] = 5, 8, 6;
    F[7].Resize(3);
    F[7] = 6, 8, 7;
    F[8].Resize(3);
    F[8] = 7, 8, 4;


    double vol; // volume of the polyhedron
    Vec3_t CM;  // Center of mass of the polyhedron
    Mat3_t It;  // Inertia tensor of the polyhedron
    PolyhedraMP(V,F,vol,CM,It);
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    Particles[Particles.Size()-1]->x       = CM;
    Particles[Particles.Size()-1]->Props.V = vol;
    Particles[Particles.Size()-1]->Props.m = vol*rho;
    Vec3_t I;
    Quaternion_t Q;
    Vec3_t xp,yp,zp;
    Eig(It,I,xp,yp,zp);
    CheckDestroGiro(xp,yp,zp);
    I *= rho;
    Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    Q(1) = (yp(2)-zp(1))/(4*Q(0));
    Q(2) = (zp(0)-xp(2))/(4*Q(0));
    Q(3) = (xp(1)-yp(0))/(4*Q(0));
    Q = Q/norm(Q);
    Particles[Particles.Size()-1]->I     = I;
    Particles[Particles.Size()-1]->Q     = Q;
    double Dmax = Distance(CM,V[0])+R;
    for (size_t i=1; i<V.Size(); ++i)
    {
        if (Distance(CM,V[i])+R > Dmax) Dmax = Distance(CM,V[i])+R;
    }
    Particles[Particles.Size()-1]->Ekin = 0.0;
    Particles[Particles.Size()-1]->Erot = 0.0;
    Particles[Particles.Size()-1]->Dmax  = Dmax;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index = Particles.Size()-1;

    Vec3_t Y = X;
    Particles[Particles.Size()-1]->Translate(Y);

}

inline void Domain::AddRice (int Tag, const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(2);
    V[0] = 0.0, 0.0,  L/2;
    V[1] = 0.0, 0.0, -L/2;

    // edges
    Array<Array <int> > E(1);
    E[0].Resize(2);
    E[0] = 0, 1;

    // faces
    Array<Array <int> > F(0); // no faces

    // calculate the rotation
    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = (4.0/3.0)*M_PI*R*R*R + M_PI*L*R*R;
    Particles[Particles.Size()-1]->Props.m    = rho*Particles[Particles.Size()-1]->Props.V;
    Particles[Particles.Size()-1]->I          = (1.0/3.0)*M_PI*rho*R*R*R*L*L + (1.0/12.0)*M_PI*rho*R*R*L*L*L + (3.0/4.0)*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R,
                                                (1.0/3.0)*M_PI*rho*R*R*R*L*L + (1.0/12.0)*M_PI*rho*R*R*L*L*L + (3.0/4.0)*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R,
                                                0.5*M_PI*rho*R*R*R*R*L + (8.0/15.0)*M_PI*rho*R*R*R*R*R;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = 0.5*L + R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;

    // clean up
    if (!ThereisanAxis) delete Axis;
}

inline void Domain::AddPlane (int Tag, const Vec3_t & X, double R, double Lx, double Ly, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(4);
    double lx = Lx/2.0, ly = Ly/2.0;
    V[0] = -lx, -ly, 0.0;
    V[1] =  lx, -ly, 0.0;
    V[2] =  lx,  ly, 0.0;
    V[3] = -lx,  ly, 0.0;

    // edges
    Array<Array <int> > E(4);
    for (size_t i=0; i<4; ++i) E[i].Resize(2);
    E[ 0] = 0, 1;
    E[ 1] = 1, 2;
    E[ 2] = 2, 3;
    E[ 3] = 3, 0;

    // faces
    Array<Array <int> > F(1);
    F[0].Resize(4);
    F[0] = 0, 3, 2, 1;

    bool ThereisanAxis = true;
    if (Axis==NULL)
    {
        Angle   = 0.;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
        ThereisanAxis = false;
    }
    Quaternion_t q;
    NormalizeRotation (Angle,(*Axis),q);
    for (size_t i=0; i<V.Size(); i++)
    {
        Vec3_t t;
        Rotation (V[i],q,t);
        V[i] = t+X;
    }

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = Lx*Ly*2*R;
    Particles[Particles.Size()-1]->Props.m    = rho*Lx*Ly*2*R;
    Particles[Particles.Size()-1]->I          = (1.0/12.0)*(Ly*Ly+4*R*R),(1.0/12.0)*(Lx*Lx+4*R*R),(1.0/12.0)*(Lx*Lx+Ly*Ly);
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(Lx*Lx+Ly*Ly)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    // clean up
    if (!ThereisanAxis) delete Axis;
}

inline void Domain::AddVoroCell (int Tag, voronoicell_neighbor & VC, double R, double rho, bool Erode, Vec3_t nv)
{
    Array<Vec3_t> V(VC.p);
    Array<Array <int> > E;
    Array<int> Eaux(2);
    for(int i=0;i<VC.p;i++) 
    {
        V[i] = Vec3_t(0.5*VC.pts[3*i]*nv(0),0.5*VC.pts[3*i+1]*nv(1),0.5*VC.pts[3*i+2]*nv(2));
        for(int j=0;j<VC.nu[i];j++) 
        {
            int k=VC.ed[i][j];
            if (VC.ed[i][j]<i) 
            {
                Eaux[0] = i;
                Eaux[1] = k;
                E.Push(Eaux);
            }
        }
    }
    Array<Array <int> > F;
    Array<int> Faux;
    for(int i=0;i<VC.p;i++) 
    {
        for(int j=0;j<VC.nu[i];j++) 
        {
            int k=VC.ed[i][j];
            if (k>=0) 
            {
                Faux.Push(i);
                VC.ed[i][j]=-1-k;
                int l=VC.cycle_up(VC.ed[i][VC.nu[i]+j],k);
                do 
                {
                    Faux.Push(k);
                    int m=VC.ed[k][l];
                    VC.ed[k][l]=-1-m;
                    l=VC.cycle_up(VC.ed[k][VC.nu[k]+l],m);
                    k=m;
                } while (k!=i);
                Array<int> Faux2(Faux.Size());
                for (size_t l = 0; l < Faux.Size();l++)
                {
                    Faux2[l] = Faux[Faux.Size()-1-l];
                }

                F.Push(Faux2);
                Faux.Clear();
                Faux2.Clear();
            }
        }
    }
    VC.reset_edges();
    double vol; // volume of the polyhedron
    Vec3_t CM;  // Center of mass of the polyhedron
    Mat3_t It;  // Inertia tensor of the polyhedron
    PolyhedraMP(V,F,vol,CM,It);
    if (Erode) Erosion(V,E,F,R);
    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    if (Erode) Particles[Particles.Size()-1]->Eroded = true;
    Particles[Particles.Size()-1]->x       = CM;
    Particles[Particles.Size()-1]->Props.V = vol;
    Particles[Particles.Size()-1]->Props.m = vol*rho;
    Vec3_t I;
    Quaternion_t Q;
    Vec3_t xp,yp,zp;
    Eig(It,I,xp,yp,zp);
    CheckDestroGiro(xp,yp,zp);
    I *= rho;
    Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    Q(1) = (yp(2)-zp(1))/(4*Q(0));
    Q(2) = (zp(0)-xp(2))/(4*Q(0));
    Q(3) = (xp(1)-yp(0))/(4*Q(0));
    Q = Q/norm(Q);
    Particles[Particles.Size()-1]->I     = I;
    Particles[Particles.Size()-1]->Q     = Q;
    double Dmax = Distance(CM,V[0])+R;
    for (size_t i=1; i<V.Size(); ++i)
    {
        if (Distance(CM,V[i])+R > Dmax) Dmax = Distance(CM,V[i])+R;
    }
    Particles[Particles.Size()-1]->Ekin = 0.0;
    Particles[Particles.Size()-1]->Erot = 0.0;
    Particles[Particles.Size()-1]->Dmax  = Dmax;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index = Particles.Size()-1;
}

inline void Domain::AddTorus (int Tag, Vec3_t const & X, Vec3_t const & N, double Rmax, double R, double rho)
{
    // Normalize normal vector
    Vec3_t n = N/norm(N);

    // Create the 2 vertices that define the torus
    Vec3_t P1 = OrthoSys::e0 - dot(OrthoSys::e0,n)*n;
    if (norm(P1)<1.0e-12) P1 = OrthoSys::e1 - dot(OrthoSys::e1,n)*n;
    P1       /= norm(P1);
    Vec3_t P2 = cross(n,P1);
    P2       /= norm(P2);
    Array<Vec3_t > V(2);
    V[0] = X + Rmax*P1;
    V[1] = X + Rmax*P2;

    // the torus has no edges or faces
    Array<Array <int> > E(0);
    Array<Array <int> > F(0);

    //Add the particle just with two vertices
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    //Input all the mass properties
    Vec3_t xp,yp,zp;
    xp = P1;
    zp = P2;
    yp = cross(zp,xp);
    CheckDestroGiro(xp,yp,zp);
    Quaternion_t q;
    q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    q(1) = (yp(2)-zp(1))/(4*q(0));
    q(2) = (zp(0)-xp(2))/(4*q(0));
    q(3) = (xp(1)-yp(0))/(4*q(0));
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = 2*M_PI*M_PI*R*R*Rmax;
    Particles[Particles.Size()-1]->Props.m    = rho*2*M_PI*M_PI*R*R*Rmax;
    Particles[Particles.Size()-1]->I          = 0.125*(4*Rmax*Rmax + 5*R*R), (Rmax*Rmax+0.75*R*R), 0.125*(4*Rmax*Rmax + 5*R*R);
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = Rmax + R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    
    Particles[Particles.Size()-1]->Tori.Push(new Torus(&Particles[Particles.Size()-1]->x,Particles[Particles.Size()-1]->Verts[0],Particles[Particles.Size()-1]->Verts[1]));

}

inline void Domain::AddCylinder (int Tag, Vec3_t const & X0, double R0, Vec3_t const & X1, double R1, double R, double rho)
{
    // Define all key variables of the cylinder
    Vec3_t n = X1 - X0;
    n /= norm(n);
    Vec3_t P1 = OrthoSys::e0 - dot(OrthoSys::e0,n)*n;
    if (norm(P1)<1.0e-12) P1 = OrthoSys::e1 - dot(OrthoSys::e1,n)*n;
    P1       /= norm(P1);
    Vec3_t P2 = cross(n,P1);
    P2       /= norm(P2);

    //The cylinder is defined by 6 vertices
    Array<Vec3_t > V(6);
    V[0] = X0 + R0*P1;
    V[1] = X0 + R0*P2;
    V[2] = X0 - R0*P1;
    V[3] = X1 + R1*P1;
    V[4] = X1 + R1*P2;
    V[5] = X1 - R1*P1;

    // It has no edges or faces
    Array<Array <int> > E(0);
    Array<Array <int> > F(0);

    //Add the particle just with the vertices
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    //Input all the mass properties
    Vec3_t xp,yp,zp;
    xp = P1;
    zp = P2;
    yp = cross(zp,xp);
    CheckDestroGiro(xp,yp,zp);
    Quaternion_t q;
    q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
    q(1) = (yp(2)-zp(1))/(4*q(0));
    q(2) = (zp(0)-xp(2))/(4*q(0));
    q(3) = (xp(1)-yp(0))/(4*q(0));
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = 4.0*M_PI*R0*R*norm(X1-X0);
    Particles[Particles.Size()-1]->Props.m    = rho*4.0*M_PI*R0*R*norm(X1-X0);
    Particles[Particles.Size()-1]->I          = 1.0, 1.0, 1.0;
    Particles[Particles.Size()-1]->x          = 0.5*(X0 + X1);
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = sqrt(0.25*dot(X1-X0,X1-X0)+max(R0,R1)*max(R0,R1)) + R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;
    Particles[Particles.Size()-1]->Tori.Push     (new Torus(&X0,Particles[Particles.Size()-1]->Verts[0],Particles[Particles.Size()-1]->Verts[1]));
    Particles[Particles.Size()-1]->Tori.Push     (new Torus(&X1,Particles[Particles.Size()-1]->Verts[3],Particles[Particles.Size()-1]->Verts[4]));
    Particles[Particles.Size()-1]->Cylinders.Push(new Cylinder(Particles[Particles.Size()-1]->Tori[0],Particles[Particles.Size()-1]->Tori[1],Particles[Particles.Size()-1]->Verts[2],Particles[Particles.Size()-1]->Verts[5]));
}

// Methods

inline void Domain::SetProps (Dict & D)
{
    for (size_t i =0 ; i<Particles.Size(); i++)
    {
        for (size_t j=0; j<D.Keys.Size(); ++j)
        {
            int tag = D.Keys[j];
            if (tag==Particles[i]->Tag)
            {
                SDPair const & p = D(tag);
                if (p.HasKey("Gn"))
                {
                    Particles[i]->Props.Gn = p("Gn");
                }
                if (p.HasKey("Gt"))
                {
                    Particles[i]->Props.Gt = p("Gt");
                }
                if (p.HasKey("Kn"))
                {
                    Particles[i]->Props.Kn = p("Kn");
                }
                if (p.HasKey("Kt"))
                {
                    Particles[i]->Props.Kt = p("Kt");
                }
                if (p.HasKey("Bn"))
                {
                    Particles[i]->Props.Bn = p("Bn");
                }
                if (p.HasKey("Bt"))
                {
                    Particles[i]->Props.Bt = p("Bt");
                }
                if (p.HasKey("Bm"))
                {
                    Particles[i]->Props.Bm = p("Bm");
                }
                if (p.HasKey("Mu"))
                {
                    Particles[i]->Props.Mu = p("Mu");
                }
                if (p.HasKey("Eps"))
                {
                    Particles[i]->Props.eps = p("Eps");
                }
                if (p.HasKey("Beta"))
                {
                    Particles[i]->Props.Beta = p("Beta");
                }
                if (p.HasKey("Eta"))
                {
                    Particles[i]->Props.Eta = p("Eta");
                }
            }
        }
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        BInteractons[i]->UpdateParameters();
    }
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        CInteractons[i]->UpdateParameters();
    }
}

inline void Domain::Initialize (double dt)
{
    if (!Initialized)
    {
        // initialize all particles
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->Initialize(i);
            Particles[i]->InitializeVelocity(dt);
        }
        //Initializing the energies
        Evis = 0.0;
        Efric = 0.0;
        Wext = 0.0;

        // initialize
        //ResetInteractons();
        // info
        Util::Stopwatch stopwatch;
        printf("\n%s--- Initializing particles ------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
        // set flag
        Initialized = true;

        // info
        double Ekin, Epot, Etot;
        Etot = CalcEnergy (Ekin, Epot);
        printf("%s  Kinematic energy   = %g%s\n",TERM_CLR4, Ekin, TERM_RST);
        printf("%s  Potential energy   = %g%s\n",TERM_CLR4, Epot, TERM_RST);
        printf("%s  Total energy       = %g%s\n",TERM_CLR2, Etot, TERM_RST);
    }
    else
    {
        for (size_t i=0; i<Particles.Size(); i++)
        {
            if (Particles[i]->vxf) Particles[i]->xb(0) = Particles[i]->x(0) - Particles[i]->v(0)*dt;
            if (Particles[i]->vyf) Particles[i]->xb(1) = Particles[i]->x(1) - Particles[i]->v(1)*dt;
            if (Particles[i]->vzf) Particles[i]->xb(2) = Particles[i]->x(2) - Particles[i]->v(2)*dt;
        }
    }

}


inline void Domain::Solve (double tf, double dt, double dtOut, ptFun_t ptSetup, ptFun_t ptReport, char const * TheFileKey, size_t VOut, size_t Nproc, double minEkin)
{
    if (VOut > 3) throw new Fatal("Domain::Solve The visualization argument can only have 4 values: 0 None, 1 povray visualization, 2 xdmf visualization and 3 both options");
    // Assigning some domain particles especifically to the output
    FileKey.Printf("%s",TheFileKey);
    idx_out = 0;
    
    //Assigning the vlaue for the time step to the domain variable
    Dt = dt;

    // initialize particles
    Initialize (Dt);


    // calc the total volume of particles (solids)
    FreePar.Resize(0);
    NoFreePar.Resize(0);
    Vs = 0.0;
    Ms = 0.0;
    MaxDmax = 0.0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        if (Particles[i]->IsFree())
        {
            Vs += Particles[i]->Props.V;
            Ms += Particles[i]->Props.m;
            if (Particles[i]->Dmax>MaxDmax) MaxDmax = Particles[i]->Dmax;
            FreePar.Push(i);
        }
        else NoFreePar.Push(i);
    }

    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    printf("%s  Total mass   of free particles   = %g%s\n",TERM_CLR4, Ms              , TERM_RST);
    printf("%s  Total volume of free particles   = %g%s\n",TERM_CLR4, Vs              , TERM_RST);
    printf("%s  Total number of particles        = %zd%s\n",TERM_CLR2, Particles.Size(), TERM_RST);

    // solve
    double t0   = Time;     // initial time
    double tout = t0; // time position for output

    Finished = false;

    // string to output energy data, if user gives the FileKey
    std::ostringstream oss_energy; 
    EnergyOutput (idx_out, oss_energy);

#ifdef USE_THREAD
    DEM::MtData MTD[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = this;
        MTD[i].Dmx      = 0.0;
    }
    pthread_t thrs[Nproc];
    
    //for (size_t i=0; i<Particles.Size()-1; i++)
    //for (size_t j=i+1; j<Particles.Size(); j++)
    //{
        //ListPosPairs.Push(make_pair(i,j));
    //}
    LinkedCell.Resize(0);
    BoundingBox(LCxmin,LCxmax);
    LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax) + iVec3_t(1,1,1);
    LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));
    //std::cout << LCellDim << std::endl;
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_create(&thrs[i], NULL, GlobalResetDisplacement, &MTD[i]);
    }
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_join(thrs[i], NULL);
        for (size_t j=0;j<MTD[i].LLC.Size();j++)
        {
            size_t idx = Pt2idx(MTD[i].LLC[j].first,LCellDim);
            LinkedCell[idx].Push(MTD[i].LLC[j].second);
        }
    }

    UpdateLinkedCells();
    //std::cout << LinkedCell.Size() << std::endl;
    //std::cout << ListPosPairs.Size() << endl;

    //std::cout << "2 " << CInteractons.Size() << std::endl;
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_create(&thrs[i], NULL, GlobalResetContacts1, &MTD[i]);
    }
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_join(thrs[i], NULL);
        //std::cout << MTD[i].LC.Size() << std::endl;
        for (size_t j=0;j<MTD[i].LC.Size();j++)
        {
        //std::cout << MTD[i].LC.Size() << std::endl;
            size_t n = MTD[i].LC[j].first;
            size_t m = MTD[i].LC[j].second;
            Listofpairs.insert(make_pair(Particles[n],Particles[m]));
            if (Particles[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
            {
                CInteractons.Push (new CInteractonSphere(Particles[n],Particles[m]));
            }
            else
            {
                CInteractons.Push (new CInteracton(Particles[n],Particles[m]));
            }
        }
    }
    //std::cout << ListPosPairs.Size() << endl;
    //std::cout << "2 " << CInteractons.Size() << std::endl;
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_create(&thrs[i], NULL, GlobalResetContacts2, &MTD[i]);
    }
    Interactons.Resize(0);
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_join(thrs[i], NULL);
        for (size_t j=0;j<MTD[i].LCI.Size();j++)
        {
            Interactons.Push(CInteractons[MTD[i].LCI[j]]);
        }
        for (size_t j=0;j<MTD[i].LCB.Size();j++)
        {
            Interactons.Push(BInteractons[MTD[i].LCB[j]]);
        }
    }

    if (Xmax-Xmin>Alpha)
    {
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalResetBoundaries1, &MTD[i]);
        }
        ParXmax.Resize(0);
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
            for (size_t j=0;j<MTD[i].LBP.Size();j++)
            {
                ParXmax.Push(Particles[MTD[i].LBP[j]]);
            }
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalResetBoundaries2, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
            for (size_t j=0;j<MTD[i].LPC.Size();j++)
            {
                size_t n = MTD[i].LPC[j].first;
                size_t m = MTD[i].LPC[j].second;
                PListofpairs.insert(make_pair(ParXmax[n],Particles[m]));
                if (ParXmax[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
                {
                    CPInteractons.Push (new CInteractonSphere(ParXmax[n],Particles[m]));
                }
                else
                {
                    CPInteractons.Push (new CInteracton(ParXmax[n],Particles[m]));
                }
            }
        }
        //std::cout << "2 " << CPInteractons.Size() << std::endl;
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalResetBoundaries3, &MTD[i]);
        }
        PInteractons.Resize(0);
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
            for (size_t j=0;j<MTD[i].LPCI.Size();j++)
            {
                PInteractons.Push(CPInteractons[MTD[i].LPCI[j]]);
            }
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalPerTranslateBack, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << ParXmax.Size() << std::endl;
    }

#else

    // set the displacement of the particles to zero (for the Halo)
    ResetDisplacements();

    // build the map of possible contacts (for the Halo)
    ResetContacts();
    if (fabs(Xmax-Xmin)>Alpha) ResetBoundaries();

#endif

    // run
    while (Time<tf)
    {

        // output
        if (Time>=tout)
        {
            double Ekin,Epot;
            CalcEnergy(Ekin,Epot);
            if (Ekin<minEkin&&Time>0.1*tf)
            {
                printf("\n%s--- Minimun energy reached ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
                break;
            }
            if (BInteractons.Size()>0) Clusters();
            if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
#ifdef USE_HDF5
                if (VOut==2||VOut==3)    WriteXDMF    (fn.CStr());
#endif
                if (VOut==1||VOut==3)    WritePOV     (fn.CStr());
                //EnergyOutput (idx_out, oss_energy);
            }
            idx_out++;
            tout += dtOut;
        }
#ifdef USE_THREAD
        //std::cout << "1" << std::endl;
        //Initialize particles
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalIni, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }

        //Calculate forces
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalForce, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }

        if (Xmax-Xmin>Alpha)
        {
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalPerTranslate, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalPerForce, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalPerTranslateBack, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
            }
        }

        //std::cout << "2" << std::endl;
        // tell the user function to update its data
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);

        // Move Particles
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalMove, &MTD[i]);
        }
        double maxdis = 0.0;
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
            if (maxdis<MTD[i].Dmx) maxdis = MTD[i].Dmx;
        }
        //std::cout << "3" << std::endl;

        if (maxdis>Alpha)
        {
            LinkedCell.Resize(0);
            BoundingBox(LCxmin,LCxmax);
            LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax) + iVec3_t(1,1,1);
            LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));
            //std::cout << "1" << std::endl;
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalResetDisplacement, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
                for (size_t j=0;j<MTD[i].LLC.Size();j++)
                {
                    size_t idx = Pt2idx(MTD[i].LLC[j].first,LCellDim);
                    LinkedCell[idx].Push(MTD[i].LLC[j].second);
                }
            }
            UpdateLinkedCells();
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalResetContacts1, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
                for (size_t j=0;j<MTD[i].LC.Size();j++)
                {
                    size_t n = MTD[i].LC[j].first;
                    size_t m = MTD[i].LC[j].second;
                    Listofpairs.insert(make_pair(Particles[n],Particles[m]));
                    if (Particles[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
                    {
                        CInteractons.Push (new CInteractonSphere(Particles[n],Particles[m]));
                    }
                    else
                    {
                        CInteractons.Push (new CInteracton(Particles[n],Particles[m]));
                    }
                }
            }
            //std::cout << "2 " << CInteractons.Size() << std::endl;
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalResetContacts2, &MTD[i]);
            }
            Interactons.Resize(0);
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
                for (size_t j=0;j<MTD[i].LCI.Size();j++)
                {
                    Interactons.Push(CInteractons[MTD[i].LCI[j]]);
                }
                for (size_t j=0;j<MTD[i].LCB.Size();j++)
                {
                    Interactons.Push(BInteractons[MTD[i].LCB[j]]);
                }
            }


            ///////////////// Periodic Boundaries ////////////////////////
            //std::cout << "1" << std::endl;
            if (Xmax-Xmin>Alpha)
            {
                for (size_t i=0;i<Nproc;i++)
                {
                    pthread_create(&thrs[i], NULL, GlobalResetBoundaries1, &MTD[i]);
                }
                ParXmax.Resize(0);
                for (size_t i=0;i<Nproc;i++)
                {
                    pthread_join(thrs[i], NULL);
                    for (size_t j=0;j<MTD[i].LBP.Size();j++)
                    {
                        ParXmax.Push(Particles[MTD[i].LBP[j]]);
                    }
                }
                for (size_t i=0;i<Nproc;i++)
                {
                    pthread_create(&thrs[i], NULL, GlobalResetBoundaries2, &MTD[i]);
                }
                for (size_t i=0;i<Nproc;i++)
                {
                    pthread_join(thrs[i], NULL);
                    for (size_t j=0;j<MTD[i].LPC.Size();j++)
                    {
                        size_t n = MTD[i].LPC[j].first;
                        size_t m = MTD[i].LPC[j].second;
                        PListofpairs.insert(make_pair(ParXmax[n],Particles[m]));
                        if (ParXmax[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
                        {
                            CPInteractons.Push (new CInteractonSphere(ParXmax[n],Particles[m]));
                        }
                        else
                        {
                            CPInteractons.Push (new CInteracton(ParXmax[n],Particles[m]));
                        }
                    }
                }
                //std::cout << "2 " << CPInteractons.Size() << std::endl;
                for (size_t i=0;i<Nproc;i++)
                {
                    pthread_create(&thrs[i], NULL, GlobalResetBoundaries3, &MTD[i]);
                }
                PInteractons.Resize(0);
                for (size_t i=0;i<Nproc;i++)
                {
                    pthread_join(thrs[i], NULL);
                    for (size_t j=0;j<MTD[i].LPCI.Size();j++)
                    {
                        PInteractons.Push(CPInteractons[MTD[i].LPCI[j]]);
                    }
                }
                for (size_t i=0;i<Nproc;i++)
                {
                    pthread_create(&thrs[i], NULL, GlobalPerTranslateBack, &MTD[i]);
                }
                for (size_t i=0;i<Nproc;i++)
                {
                    pthread_join(thrs[i], NULL);
                }
                //std::cout << ParXmax.Size() << std::endl;
            }
            //std::cout << "3 " << PInteractons.Size() << std::endl;
        }


        //std::cout << "4 " << Time << std::endl;
#else 


        // initialize forces and torques
        for (size_t i=0; i<Particles.Size(); i++)
        {
            // set the force and torque to the fixed values
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Tf;
            for (size_t n=0;n<3;n++)
            {
                for (size_t m=0;m<3;m++)  
                {
                    Particles[i]->M(n,m)=0.0;
                    Particles[i]->B(n,m)=0.0;
                }
            }

            // initialize the coordination (number of contacts per particle) number and the flag of the particle begin in contact with the container
            Particles[i]->Cn   = 0.0;
            Particles[i]->Bdry = false;

            // external work added to the system by the fixed forces Ff
            Wext += dot(Particles[i]->Ff,Particles[i]->v)*Dt;
        }

        // calc contact forces: collision and bonding (cohesion)
        for (size_t i=0; i<Interactons.Size(); i++)
        {
            Interactons[i]->CalcForce (Dt);
        }

        if (fabs(Xmax-Xmin)>Alpha)
        {
            Vec3_t vmax(Xmin-Xmax,0.0,0.0);
            Vec3_t vmin(Xmax-Xmin,0.0,0.0);
            //for (size_t i=0;i<ParXmin.Size();i++) ParXmin[i]->Translate(vmin);
            for (size_t i=0;i<ParXmax.Size();i++) ParXmax[i]->Translate(vmax);

            for (size_t i=0; i<PInteractons.Size(); i++)
            {   
                Interacton * PI = PInteractons[i];
                PI->CalcForce(Dt);
            }

            //for (size_t i=0;i<ParXmin.Size();i++) ParXmin[i]->Translate(vmax);
            for (size_t i=0;i<ParXmax.Size();i++) ParXmax[i]->Translate(vmin);
        }

        // calculate the collision energy
        for (size_t i=0; i<CInteractons.Size(); i++)
        {
            Evis  += CInteractons[i]-> dEvis;
            Efric += CInteractons[i]-> dEfric;
        }

        // tell the user function to update its data
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);

        // move particles
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->Translate (Dt);
            Particles[i]->Rotate    (Dt);
        }

        double maxdis = MaxDisplacement();

        // update the Halos
        if (maxdis>Alpha)
        {
            ResetDisplacements();
            ResetContacts();
            if (fabs(Xmax-Xmin)>Alpha) ResetBoundaries();
        }
#endif
        

        // next time position
        Time += Dt;
    }

    // last output
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    // save energy data
    //if (TheFileKey!=NULL)
    //{
        //String fn;
        //fn.Printf("%s_energy.res",TheFileKey);
        //std::ofstream fe(fn.CStr());
        //fe << oss_energy.str();
        //fe.close();
    //}

    // info
    double Ekin, Epot, Etot;
    Etot = CalcEnergy (Ekin, Epot);
    printf("%s  Kinematic energy   = %g%s\n",TERM_CLR4, Ekin, TERM_RST);
    printf("%s  Potential energy   = %g%s\n",TERM_CLR4, Epot, TERM_RST);
    printf("%s  Total energy       = %g%s\n",TERM_CLR2, Etot, TERM_RST);
}

inline void Domain::WritePOV (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".pov");
    std::ofstream of(fn.CStr(), std::ios::out);
    POVHeader (of);
    POVSetCam (of, CamPos, OrthoSys::O);
    Array <String> Colors(10);
    Colors = "Gray","Blue","Yellow","Gold","Green","Blue","Orange","Salmon","Copper","Aquamarine";
    for (size_t i=0; i<Particles.Size(); i++)
    {
        if (!Particles[i]->IsFree()) Particles[i]->Draw(of,"Col_Glass_Bluish");
        else
        {
            bool found = false;
            for (size_t j=0;j<Listofclusters.Size();j++)
            {
                if (Listofclusters[j].Has(i)&&BInteractons.Size()>0)
                {
                    Particles[i]->Draw(of,Colors[j%10].CStr());
                    found = true;
                    break;
                }
            }
            if (!found) Particles[i]->Draw(of,"Red");
        }
    }

}

inline void Domain::WriteBPY (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".bpy");
    std::ofstream of(fn.CStr(), std::ios::out);
    BPYHeader(of);
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Draw (of,"",true);
    of.close();
}

#ifdef USE_HDF5

inline void Domain::WriteBF (char const * FileKey)
{

    size_t n_fn = 0;

    for (size_t i=0;i<CInteractons.Size();i++)
    {
        if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree())) n_fn++;
        //if (norm(CInteractons[i]->Fnet)>1.0e-12) n_fn++;
    }

    if (n_fn==0) return;
    
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    float  *  Fnnet = new float[3*n_fn];
    float  *  Ftnet = new float[3*n_fn];
    float  * Branch = new float[3*n_fn];
    float  *   Orig = new float[3*n_fn];
    float  *    Fn  = new float[  n_fn];
    int    *    ID1 = new   int[  n_fn];
    int    *    ID2 = new   int[  n_fn];

    size_t idx = 0;

    for (size_t i=0;i<CInteractons.Size();i++)
    {
        if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree()))
        //if (norm(CInteractons[i]->Fnet)>1.0e-12)
        {
            Fnnet [3*idx  ] = float(CInteractons[i]->Fnet  (0));
            Fnnet [3*idx+1] = float(CInteractons[i]->Fnet  (1));
            Fnnet [3*idx+2] = float(CInteractons[i]->Fnet  (2));
            Ftnet [3*idx  ] = float(CInteractons[i]->Ftnet (0));
            Ftnet [3*idx+1] = float(CInteractons[i]->Ftnet (1));
            Ftnet [3*idx+2] = float(CInteractons[i]->Ftnet (2));
            Branch[3*idx  ] = float(CInteractons[i]->P1->x(0)-CInteractons[i]->P2->x(0));
            Branch[3*idx+1] = float(CInteractons[i]->P1->x(1)-CInteractons[i]->P2->x(1)); 
            Branch[3*idx+2] = float(CInteractons[i]->P1->x(2)-CInteractons[i]->P2->x(2)); 
            //Orig  [3*idx  ] = 0.5*float(CInteractons[i]->P1->x(0)+CInteractons[i]->P2->x(0));
            //Orig  [3*idx+1] = 0.5*float(CInteractons[i]->P1->x(1)+CInteractons[i]->P2->x(1)); 
            //Orig  [3*idx+2] = 0.5*float(CInteractons[i]->P1->x(2)+CInteractons[i]->P2->x(2)); 
            Orig  [3*idx  ] = float(CInteractons[i]->P2->x(0));
            Orig  [3*idx+1] = float(CInteractons[i]->P2->x(1)); 
            Orig  [3*idx+2] = float(CInteractons[i]->P2->x(2)); 
            Fn    [idx]     = float(norm(CInteractons[i]->Fnet));
            ID1   [idx]     = int  (CInteractons[i]->P1->Index);
            ID2   [idx]     = int  (CInteractons[i]->P2->Index);
            idx++;
        }
    }

    hsize_t dims[1];
    dims[0] = 3*n_fn;
    String dsname;
    dsname.Printf("Normal");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Fnnet );
    dsname.Printf("Tangential");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ftnet );
    dsname.Printf("Branch");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Branch);
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Orig);
    dims[0] = n_fn;
    dsname.Printf("Fn");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Fn);
    dsname.Printf("ID1");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,ID1   );
    dsname.Printf("ID2");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,ID2   );


    delete [] Fnnet;
    delete [] Ftnet;
    delete [] Branch;
    delete [] Orig;
    delete [] Fn;
    delete [] ID1;
    delete [] ID2;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"BranchForce\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << n_fn << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << n_fn << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Normal\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Normal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tangential\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tangential \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Branch\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Branch \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Fn\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Fn \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ID1\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ID1 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"ID2\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << n_fn << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/ID2 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::WriteXDMF (char const * FileKey)
{
    size_t N_Faces = 0;
    size_t N_Verts = 0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        for (size_t j=0;j<Particles[i]->Faces.Size();j++)
        {
            N_Faces += Particles[i]->Faces[j]->Edges.Size();
        }
        N_Verts += Particles[i]->Verts.Size() + Particles[i]->Faces.Size();
    }

    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (N_Faces>0)
    {

        //Geometric information
        float  * Verts   = new float [3*N_Verts];
        int    * FaceCon = new int   [3*N_Faces];
        
        //Atributes
        int    * Tags    = new int   [  N_Faces];
        int    * Clus    = new int   [  N_Faces];
        float  * Vel     = new float [  N_Faces];
        float  * Ome     = new float [  N_Faces];
        float  * FComp   = new float [  N_Faces];
        //float  * Stress  = new float [9*N_Faces];

        size_t n_verts = 0;
        size_t n_faces = 0;
        size_t n_attrs = 0;
        //size_t n_attrv = 0;
        //size_t n_attrt = 0;
        for (size_t i=0;i<Particles.Size();i++)
        {
            Particle * Pa = Particles[i];
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(Pa->Verts.Size());
            Array<Vec3_t> Vres (Pa->Verts.Size());
            for (size_t j=0;j<Pa->Verts.Size();j++)
            {
                Vtemp[j] = *Pa->Verts[j];
                Vres [j] = *Pa->Verts[j];
            }
            double multiplier = 0.0;
            if (Dilate&&Pa->Eroded&&Pa->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,Pa->EdgeCon,Pa->FaceCon,Vres,Pa->Props.R);
                multiplier = 1.0;
            }
            for (size_t j=0;j<Pa->Verts.Size();j++)
            {
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(0);
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(1);
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(2);
                Verts[n_verts++] = float(Vres[j](0));
                Verts[n_verts++] = float(Vres[j](1));
                Verts[n_verts++] = float(Vres[j](2));
            }
            size_t n_reff = n_verts/3;
            for (size_t j=0;j<Pa->FaceCon.Size();j++)
            {
                Vec3_t C,N;
                Pa->Faces[j]->Centroid(C);
                Pa->Faces[j]->Normal(N);
                Verts[n_verts++] = float(C(0) + multiplier*Pa->Props.R*N(0));
                Verts[n_verts++] = float(C(1) + multiplier*Pa->Props.R*N(1));
                Verts[n_verts++] = float(C(2) + multiplier*Pa->Props.R*N(2));
                //Verts[n_verts++] = (float) C(0);
                //Verts[n_verts++] = (float) C(1);
                //Verts[n_verts++] = (float) C(2);
                for (size_t k=0;k<Pa->FaceCon[j].Size();k++)
                {
                    size_t nin = Pa->FaceCon[j][k];
                    size_t nen = Pa->FaceCon[j][(k+1)%Pa->FaceCon[j].Size()];
                    FaceCon[n_faces++] = int(n_reff + j);  
                    FaceCon[n_faces++] = int(n_refv + nin);
                    FaceCon[n_faces++] = int(n_refv + nen);

                    //Writing the attributes
                    Tags  [n_attrs] = int(Pa->Tag);
                    Clus  [n_attrs] = size_t(Pa->Cluster);
                    Vel   [n_attrs] = float(norm(Pa->v));
                    Ome   [n_attrs] = float(norm(Pa->w));
                    FComp [n_attrs] = float(Pa->Comp);
                    n_attrs++;

                    //Vel [n_attrv  ] = (float) Pa->v(0);
                    //Vel [n_attrv+1] = (float) Pa->v(1);
                    //Vel [n_attrv+2] = (float) Pa->v(2);
                    //n_attrv += 3;

                    //Stress[n_attrt  ] = (float) Pa->M(0,0);
                    //Stress[n_attrt+1] = (float) Pa->M(1,0);
                    //Stress[n_attrt+2] = (float) Pa->M(2,0);
                    //Stress[n_attrt+3] = (float) Pa->M(0,1);
                    //Stress[n_attrt+4] = (float) Pa->M(1,1);
                    //Stress[n_attrt+5] = (float) Pa->M(2,1);
                    //Stress[n_attrt+6] = (float) Pa->M(0,2);
                    //Stress[n_attrt+7] = (float) Pa->M(1,2);
                    //Stress[n_attrt+8] = (float) Pa->M(2,2);
                    //n_attrt += 9;
                }
            }
        }

        //Write the data
        hsize_t dims[1];
        String dsname;
        dims[0] = 3*N_Verts;
        dsname.Printf("Verts");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Verts);
        dims[0] = 3*N_Faces;
        dsname.Printf("FaceCon");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,FaceCon);
        dims[0] = N_Faces;
        dsname.Printf("Tag");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tags   );
        dims[0] = N_Faces;
        dsname.Printf("Cluster");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Clus   );
        dims[0] = N_Faces;
        dsname.Printf("Velocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vel);
        dims[0] = N_Faces;
        dsname.Printf("AngVelocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ome);
        dims[0] = N_Faces;
        dsname.Printf("FaceComp");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,FComp);
        //dims[0] = 9*N_Faces;
        //dsname.Printf("Stress");
        //H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Stress);
        
        //Erasing the data
        delete [] Verts;
        delete [] FaceCon;
        delete [] Tags;
        delete [] Clus;
        delete [] Vel;
        delete [] Ome;
        delete [] FComp;
        //delete [] Stress;
    }
    // Storing center of mass data
    
    float * Radius = new float[  Particles.Size()];
    float * Posvec = new float[3*Particles.Size()];
    float * Velvec = new float[3*Particles.Size()];
    float * Omevec = new float[3*Particles.Size()];
    int   * Tag    = new int  [  Particles.Size()];
    float * Comp   = new float[  Particles.Size()];
    for (size_t i=0;i<Particles.Size();i++)
    {
        Vec3_t Ome;
        Rotation(Particles[i]->w,Particles[i]->Q,Ome);
        Radius[i]     = float(Particles[i]->Dmax);
        Posvec[3*i  ] = float(Particles[i]->x(0));
        Posvec[3*i+1] = float(Particles[i]->x(1));
        Posvec[3*i+2] = float(Particles[i]->x(2));
        Velvec[3*i  ] = float(Particles[i]->v(0));
        Velvec[3*i+1] = float(Particles[i]->v(1));
        Velvec[3*i+2] = float(Particles[i]->v(2));
        Omevec[3*i  ] = float(Ome(0));
        Omevec[3*i+1] = float(Ome(1)); 
        Omevec[3*i+2] = float(Ome(2)); 
        Tag   [i]     = int  (Particles[i]->Tag);  
        //Comp  [i]     = float(Particles[i]->M(0,0) + Particles[i]->M(1,1) + Particles[i]->M(2,2));
        Comp[i]       = float(Particles[i]->Comp);
    }

    hsize_t dims[1];
    dims[0] = 3*Particles.Size();
    String dsname;
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("PVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dsname.Printf("PAngVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Omevec);
    dims[0] = Particles.Size();
    dsname.Printf("Radius");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
    dsname.Printf("PTag");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tag   );
    dsname.Printf("PComp");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Comp  );


    delete [] Radius;
    delete [] Posvec;
    delete [] Velvec;
    delete [] Omevec;
    delete [] Tag;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
    

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    if(N_Faces>0)
    {
    oss << "   <Grid Name=\"DEM_Faces\">\n";
    oss << "     <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << N_Faces << "\">\n";
    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << N_Faces << " 3\">\n";
    oss << "        " << fn.CStr() <<":/FaceCon \n";
    oss << "       </DataItem>\n";
    oss << "     </Topology>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << N_Verts << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Verts \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Cluster\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Cluster \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/AngVelocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"FaceComp\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/FaceComp \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    //oss << "     </Attribute>\n";
    //oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Cell\">\n";
    //oss << "       <DataItem Dimensions=\"" << N_Faces << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    //oss << "        " << fn.CStr() <<":/AngVelocity \n";
    //oss << "       </DataItem>\n";
    //oss << "     </Attribute>\n";
    //oss << "     <Attribute Name=\"Stress\" AttributeType=\"Tensor\" Center=\"Cell\">\n";
    //oss << "       <DataItem Dimensions=\"" << N_Faces << " 3 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    //oss << "        " << fn.CStr() <<":/Stress \n";
    //oss << "       </DataItem>\n";
    //oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << "   <Grid Name=\"CMCenter\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Compression\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PComp \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVel\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";


    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::Save (char const * FileKey)
{

    // Opening the file for writing
    String fn(FileKey);
    fn.append(".hdf5");
    if (Util::FileExists(fn))
    {
        String command;
        command.Printf("rm %s",fn.CStr());
        system(command.CStr());
    }
    hid_t file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Storing the number of particles in the domain
    int data[1];
    data[0]=Particles.Size();
    hsize_t dims[1];
    dims[0]=1;
    H5LTmake_dataset_int(file_id,"/NP",1,dims,data);

    for (size_t i=0; i<Particles.Size(); i++)
    {
        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gcreate(file_id, par.CStr(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


        // Storing some scalar variables
        double dat[1];
        dat[0] = Particles[i]->Props.R;
        H5LTmake_dataset_double(group_id,"SR",1,dims,dat);
        dat[0] = Particles[i]->Props.rho;
        H5LTmake_dataset_double(group_id,"Rho",1,dims,dat);
        dat[0] = Particles[i]->Props.m;
        H5LTmake_dataset_double(group_id,"m",1,dims,dat);
        dat[0] = Particles[i]->Props.V;
        H5LTmake_dataset_double(group_id,"V",1,dims,dat);
        dat[0] = Particles[i]->Diam;
        H5LTmake_dataset_double(group_id,"Diam",1,dims,dat);
        dat[0] = Particles[i]->Dmax;
        H5LTmake_dataset_double(group_id,"Dmax",1,dims,dat);
        int datint[1];
        datint[0] = Particles[i]->Index;
        H5LTmake_dataset_int(group_id,"Index",1,dims,datint);


        int tag[1];
        tag[0] = Particles[i]->Tag;
        H5LTmake_dataset_int(group_id,"Tag",1,dims,tag);

        // Storing vectorial variables
        double cd[3];
        hsize_t dd[1];
        dd[0] = 3;

        cd[0]=Particles[i]->x(0);
        cd[1]=Particles[i]->x(1);
        cd[2]=Particles[i]->x(2);
        H5LTmake_dataset_double(group_id,"x",1,dd,cd);

        cd[0]=Particles[i]->xb(0);
        cd[1]=Particles[i]->xb(1);
        cd[2]=Particles[i]->xb(2);
        H5LTmake_dataset_double(group_id,"xb",1,dd,cd);

        cd[0]=Particles[i]->v(0);
        cd[1]=Particles[i]->v(1);
        cd[2]=Particles[i]->v(2);
        H5LTmake_dataset_double(group_id,"v",1,dd,cd);

        cd[0]=Particles[i]->w(0);
        cd[1]=Particles[i]->w(1);
        cd[2]=Particles[i]->w(2);
        H5LTmake_dataset_double(group_id,"w",1,dd,cd);

        cd[0]=Particles[i]->wb(0);
        cd[1]=Particles[i]->wb(1);
        cd[2]=Particles[i]->wb(2);
        H5LTmake_dataset_double(group_id,"wb",1,dd,cd);

        cd[0]=Particles[i]->I(0);
        cd[1]=Particles[i]->I(1);
        cd[2]=Particles[i]->I(2);
        H5LTmake_dataset_double(group_id,"I",1,dd,cd);

        double cq[4];
        dd[0] = 4;
        cq[0]=Particles[i]->Q(0);
        cq[1]=Particles[i]->Q(1);
        cq[2]=Particles[i]->Q(2);
        cq[3]=Particles[i]->Q(3);
        H5LTmake_dataset_double(group_id,"Q",1,dd,cq);




        // Storing the number of vertices of each particle
        data[0] = Particles[i]->Verts.Size();
        H5LTmake_dataset_int(group_id,"n_vertices",1,dims,data);
        hid_t gv_id;
        gv_id = H5Gcreate(group_id,"Verts", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Storing each vertex 
        for (size_t j=0;j<Particles[i]->Verts.Size();j++)
        {
            String parv;
            parv.Printf("Verts_%08d",j);
            double cod[3];
            cod[0]=(*Particles[i]->Verts[j])(0);
            cod[1]=(*Particles[i]->Verts[j])(1);
            cod[2]=(*Particles[i]->Verts[j])(2);
            hsize_t dim[1];
            dim[0]=3;
            H5LTmake_dataset_double(gv_id,parv.CStr(),1,dim,cod);
        }

        // Number of edges of the particle
        data[0] = Particles[i]->Edges.Size();
        H5LTmake_dataset_int(group_id,"n_edges",1,dims,data);
        gv_id = H5Gcreate(group_id,"Edges", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        // Edges
        for (size_t j=0;j<Particles[i]->Edges.Size();j++)
        {
            String parv;
            parv.Printf("Edges_%08d",j);
            int co[2];
            co[0] = Particles[i]->EdgeCon[j][0];
            co[1] = Particles[i]->EdgeCon[j][1];
            hsize_t dim[1];
            dim[0] =2;
            H5LTmake_dataset_int(gv_id,parv.CStr(),1,dim,co);
        }
        
        // Number of faces of the particle
        data[0] = Particles[i]->Faces.Size();
        H5LTmake_dataset_int(group_id,"n_faces",1,dims,data);
        gv_id = H5Gcreate(group_id,"Faces", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        // Faces
        for (size_t j=0;j<Particles[i]->Faces.Size();j++)
        {
            String parv;
            parv.Printf("Faces_%08d",j);
            int co[Particles[i]->FaceCon[j].Size()];
            hsize_t dim[1];
            dim[0]= Particles[i]->FaceCon[j].Size();
            for (size_t k=0;k<Particles[i]->FaceCon[j].Size();k++)
            {
                co[k]=Particles[i]->FaceCon[j][k];
            }
            H5LTmake_dataset_int(gv_id,parv.CStr(),1,dim,co);
        }
        
    }

    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
    //sleep(5);
}

inline void Domain::Load (char const * FileKey)
{

    // Opening the file for reading
    String fn(FileKey);
    fn.append(".hdf5");
    if (!Util::FileExists(fn)) throw new Fatal("File <%s> not found",fn.CStr());
    printf("\n%s--- Loading file %s --------------------------------------------%s\n",TERM_CLR1,fn.CStr(),TERM_RST);
    hid_t file_id;
    file_id = H5Fopen(fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Number of particles in the domain
    int data[1];
    H5LTread_dataset_int(file_id,"/NP",data);
    size_t NP = data[0];

    // Loading the particles
    for (size_t i=0; i<NP; i++)
    {

        // Creating the string and the group for each particle
        hid_t group_id;
        String par;
        par.Printf("/Particle_%08d",i);
        group_id = H5Gopen(file_id, par.CStr(),H5P_DEFAULT);

        // Finding the particle's position for the domain decomposition
        double X[3];
        H5LTread_dataset_double(group_id,"x",X);


        // Loading the Vertices
        H5LTread_dataset_int(group_id,"n_vertices",data);
        size_t nv = data[0];
        hid_t gv_id;
        gv_id = H5Gopen(group_id,"Verts", H5P_DEFAULT);
        Array<Vec3_t> V;

        for (size_t j=0;j<nv;j++)
        {
            String parv;
            parv.Printf("Verts_%08d",j);
            double cod[3];
            H5LTread_dataset_double(gv_id,parv.CStr(),cod);
            V.Push(Vec3_t(cod[0],cod[1],cod[2]));
        }
        
        // Loading the edges
        H5LTread_dataset_int(group_id,"n_edges",data);
        size_t ne = data[0];
        gv_id = H5Gopen(group_id,"Edges", H5P_DEFAULT);
        Array<Array <int> > E;

        for (size_t j=0;j<ne;j++)
        {
            String parv;
            parv.Printf("Edges_%08d",j);
            int cod[2];
            H5LTread_dataset_int(gv_id,parv.CStr(),cod);
            Array<int> Ep(2);
            Ep[0]=cod[0];
            Ep[1]=cod[1];
            E.Push(Ep);
        }

        // Loading the faces

        // Number of faces of the particle
        H5LTread_dataset_int(group_id,"n_faces",data);
        size_t nf = data[0];
        gv_id = H5Gopen(group_id,"Faces", H5P_DEFAULT);
        Array<Array <int> > F;
        
        // Faces
        for (size_t j=0;j<nf;j++)
        {
            String parv;
            parv.Printf("Faces_%08d",j);
            hsize_t dim[1];
            H5LTget_dataset_info(gv_id,parv.CStr(),dim,NULL,NULL);
            size_t ns = (size_t)dim[0];
            int co[ns];
            Array<int> Fp(ns);

            H5LTread_dataset_int(gv_id,parv.CStr(),co);
            
            for (size_t k=0;k<ns;k++)
            {
                Fp[k] = co[k];
            }

            F.Push(Fp);

        }

        Particles.Push (new Particle(-1,V,E,F,OrthoSys::O,OrthoSys::O,0.1,1.0));

        // Loading vectorial variables
        Particles[Particles.Size()-1]->x = Vec3_t(X[0],X[1],X[2]);
        double cd[3];
        H5LTread_dataset_double(group_id,"xb",cd);
        Particles[Particles.Size()-1]->xb = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"v",cd);
        Particles[Particles.Size()-1]->v = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"w",cd);
        Particles[Particles.Size()-1]->w = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"wb",cd);
        Particles[Particles.Size()-1]->wb = Vec3_t(cd[0],cd[1],cd[2]);
        H5LTread_dataset_double(group_id,"I",cd);
        Particles[Particles.Size()-1]->I = Vec3_t(cd[0],cd[1],cd[2]);

        double cq[4];
        H5LTread_dataset_double(group_id,"Q",cq);
        Particles[Particles.Size()-1]->Q = Quaternion_t(cq[0],cq[1],cq[2],cq[3]);

        // Loading the scalar quantities of the particle
        double dat[1];
        H5LTread_dataset_double(group_id,"SR",dat);
        Particles[Particles.Size()-1]->Props.R = dat[0];
        H5LTread_dataset_double(group_id,"Rho",dat);
        Particles[Particles.Size()-1]->Props.rho = dat[0];
        H5LTread_dataset_double(group_id,"m",dat);
        Particles[Particles.Size()-1]->Props.m = dat[0];
        H5LTread_dataset_double(group_id,"V",dat);
        Particles[Particles.Size()-1]->Props.V = dat[0];
        H5LTread_dataset_double(group_id,"Diam",dat);
        Particles[Particles.Size()-1]->Diam = dat[0];
        H5LTread_dataset_double(group_id,"Dmax",dat);
        Particles[Particles.Size()-1]->Dmax = dat[0];
        int datint[1];
        H5LTread_dataset_int(group_id,"Index",datint);
        Particles[Particles.Size()-1]->Index = datint[0];
        int tag[1];
        H5LTread_dataset_int(group_id,"Tag",tag);
        Particles[Particles.Size()-1]->Tag = tag[0];
        Particles[Particles.Size()-1]->PropsReady = true;

    }


    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);
}

#endif


#ifdef USE_VTK

void Domain::WriteVTKContacts  (char const * FileKey)
{
    size_t ncontacts = 0;

    for (size_t i=0;i<CInteractons.Size();i++)
    {
        if ((norm(CInteractons[i]->Fnet)>0.0)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree())) ncontacts++;
    }

    if (ncontacts==0) return;
    
    // Create a vtkPoints object and store the points in it
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetDataTypeToDouble(); // set the precision to double
    // add the points     // WARNING: not considering periodic conditions or walls
    for (size_t i=0;i<CInteractons.Size();i++)
    {
        if ((norm(CInteractons[i]->Fnet)>0.0)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree()))
        {
            Vec3_t R1 = CInteractons[i]->P1->x; 
            Vec3_t R2 = CInteractons[i]->P2->x; 
            points->InsertNextPoint(R1(0), R1(1), R1(2)); 
            points->InsertNextPoint(R2(0), R2(1), R2(2)); 
        }
    } 
    const auto npoints = points->GetNumberOfPoints();
    
    // Create a cell array to store the lines in and add the lines to it
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    // create each line according to contact
    for(int ip = 0; ip < npoints; ip += 2) 
    {
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0, ip);
        line->GetPointIds()->SetId(1, ip+1);
        lines->InsertNextCell(line);
    }
    
    // Create a polydata to store everything in
    vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
    // Add the points to the dataset
    polyData->SetPoints(points);
    // Add the lines to the dataset
    polyData->SetLines(lines);
    
    // ADD CELL DATA
    // SEE: http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/MiscCellData
    // SCALARS
    /*// uid1
    vtkSmartPointer<vtkIntArray> uid1 = 
      vtkSmartPointer<vtkIntArray>::New();
    uid1->SetNumberOfComponents(1); 
    uid1->SetNumberOfTuples(ncontacts); 
    uid1->SetName("uid1");
    for (int ic = 0; ic < ncontacts; ++ic) {
      uid1->SetValue(ic, contacts[ic].uid1_); 
    }
    polyData->GetCellData()->AddArray(uid1);*/
    // VECTORS
    // Normal
    vtkSmartPointer<vtkDoubleArray> NormalForce = 
    vtkSmartPointer<vtkDoubleArray>::New();
    NormalForce->SetNumberOfComponents(1); 
    NormalForce->SetNumberOfTuples(ncontacts); 
    NormalForce->SetName("Normal Force");
    size_t ic = 0;
    for (size_t i=0;i<CInteractons.Size();i++)
    {
        if ((norm(CInteractons[i]->Fnet)>0.0)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree()))
        {
            Vec3_t P = CInteractons[i]->Fnet;
            double data[1] = {norm(P)}; 
            NormalForce->SetTupleValue(ic, data); 
            ic++;
        }
    }
    polyData->GetCellData()->AddArray(NormalForce);
    
    // Write the file
    String fn(FileKey);
    fn.append(".vtp");
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(fn.CStr());
#if VTK_MAJOR_VERSION <= 5
    writer->SetInput(polyData);
#else
    writer->SetInputData(polydata);
#endif
    //// Optional - set the mode.
    writer->SetDataModeToBinary(); // default
    //writer->SetDataModeToAscii();
    writer->SetCompressorTypeToZLib(); // or ToNone
    // write
    writer->Write();    
}

#endif // USE_VTK


#ifdef USE_THREAD

inline void Domain::UpdateLinkedCells()
{
    //std::cout << "a ) Updating linked cells "  << Time << std::endl;
    //
    
    size_t n = 0;
    for (size_t i=0;i<LinkedCell.Size();i++)
    {
        //for (size_t j=0;j<LinkedCell[i].Size();j++)
        //{
            //std::cout << LinkedCell[i][j] << " ";
        //}
        //std::cout << std::endl;
        n += LinkedCell[i].Size();
    }

    //std::cout << "Number of linked cells   " << LinkedCell.Size() << " " << LCellDim << std::endl;
    //std::cout << "Limits                   " << LCxmin            << " " << LCxmax   << " " << MaxDmax << std::endl;
    //std::cout << "Particles in linked cells" << n                 << std::endl;
    //std::cout << "Free Particles           " << FreePar.Size()    << std::endl;
    if (n!=FreePar.Size())
    {
        //ofstream lin("linked.txt");
        //ofstream pos("position.txt");
        //for (size_t i=0;i<LinkedCell.Size();i++)
        //{
            //lin << i << std::endl;
            //for (size_t j=0;j<LinkedCell[i].Size();j++)
            //{
                //lin << LinkedCell[i][j] << " ";
            //}
            //lin << std::endl;
        //}
        //for (size_t i=0;i<Particles.Size();i++)
        //{
            //iVec3_t cell = (Particles[i]->x - LCxmin)/(2.0*Beta*MaxDmax);
            //size_t  idx  = Pt2idx(cell,LCellDim);
            //pos << i << " " << Particles[i]->x << " " << cell << " " << idx << std::endl;
        //}
        //lin.close();
        //pos.close();
        
        //for (size_t i=0;i<FreePar.Size();i++)
        //{
            //bool found = false;
            //for (size_t j=0;j<LinkedCell.Size();j++)
            //{
                //if (LinkedCell[j].Has(FreePar[i])) found = true;
            //}
            //if (!found) 
            //{
                //iVec3_t cell = (Particles[FreePar[i]]->x - LCxmin)/(2.0*Beta*MaxDmax);
                //size_t  idx  = Pt2idx(cell,LCellDim);
                //if (FreePar[i]==345) std::cout << FreePar[i] << " " << Particles[FreePar[i]]->x << " " << cell << " " << idx << " " << LinkedCell[idx].Size() <<  std::endl;
            //}
        //}
        throw new Fatal("Domain::UpdateLinkedCells: Linked cells dont match");
    }
    
    ListPosPairs.Resize(0);
    for (size_t i=0;i<FreePar.Size();i++)
    for (size_t j=0;j<NoFreePar.Size();j++)
    {
        ListPosPairs.Push(make_pair(FreePar[i],NoFreePar[j]));
    }

    for (size_t k=0;k<LCellDim(2);k++)
    for (size_t j=0;j<LCellDim(1);j++)
    for (size_t i=0;i<LCellDim(0);i++)
    {
        iVec3_t Pt(i,j,k);
        size_t idx = Pt2idx(Pt,LCellDim);
        //std::cout << Pt << " " << idx << " " << LinkedCell[idx].Size() << std::endl;
        if (LinkedCell[idx].Size()==0) continue;

        for (size_t n=0  ;n<LinkedCell[idx].Size()-1;n++)
        for (size_t m=n+1;m<LinkedCell[idx].Size()  ;m++)
        {
            size_t i1 = LinkedCell[idx][n];
            size_t i2 = LinkedCell[idx][m];
            if (i1==i2) continue;
            //if (i1==i2) std::cout << i1 << " " << i2 << " " << n << " " << m << " " << idx << std::endl;
            ListPosPairs.Push(make_pair(i1,i2));
        }
        //std::cout << k << " " << std::min(LCellDim(0),k+1) << std::endl;;
        for (size_t knb=std::max(0,int(k)-1);knb<=std::min(LCellDim(2)-1,k+1);knb++)
        for (size_t jnb=std::max(0,int(j)-1);jnb<=std::min(LCellDim(1)-1,j+1);jnb++)
        for (size_t inb=std::max(0,int(i)-1);inb<=std::min(LCellDim(0)-1,i+1);inb++)
        {
            //std ::cout << inb << " " << jnb << " " << knb << std::endl;
            iVec3_t Ptnb(inb,jnb,knb);
            size_t idxnb = Pt2idx(Ptnb,LCellDim);
            //std::cout << "a" << Pt << " " << idx << " " << LinkedCell[idx].Size() << std::endl;
            //std::cout << "b" << Ptnb << " " << idxnb << " " << LinkedCell[idxnb].Size() << std::endl;
            if (idxnb>idx)
            {
                for (size_t n=0;n<LinkedCell[idx].Size()  ;n++)
                {
                    for (size_t m=0;m<LinkedCell[idxnb].Size()  ;m++)
                    {
                        //size_t i1 = LinkedCell[idx  ][n];
                        //size_t i2 = LinkedCell[idxnb][m];
                        size_t i1 = std::min(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        size_t i2 = std::max(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        if (i1==i2) continue;
                        //if (i1==i2) std::cout << i1 << " " << i2 << " " << n << " " << m << " " << idx << " " << idxnb << std::endl;
                        ListPosPairs.Push(make_pair(i1,i2));
                    }
                }
            }
        }
    }
    //std::cout << "b ) Updating linked cells "  << Time << " " << ListPosPairs.Size() << std::endl;
    //for (size_t i=0;i<ListPosPairs.Size();i++)
    //{
        //std::cout << ListPosPairs[i].first << " " << ListPosPairs[i].second << std::endl;
    //}
}

#endif

inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    if (Particles.Size()==0) throw new Fatal("DEM::Domain::BoundingBox: There are no particles to build the bounding box");
    minX = Vec3_t(Particles[0]->MinX(), Particles[0]->MinY(), Particles[0]->MinZ());
    maxX = Vec3_t(Particles[0]->MaxX(), Particles[0]->MaxY(), Particles[0]->MaxZ());
    for (size_t i=1; i<Particles.Size(); i++)
    {
        if (minX(0)>Particles[i]->MinX()&&Particles[i]->IsFree()) minX(0) = Particles[i]->MinX();
        if (minX(1)>Particles[i]->MinY()&&Particles[i]->IsFree()) minX(1) = Particles[i]->MinY();
        if (minX(2)>Particles[i]->MinZ()&&Particles[i]->IsFree()) minX(2) = Particles[i]->MinZ();
        if (maxX(0)<Particles[i]->MaxX()&&Particles[i]->IsFree()) maxX(0) = Particles[i]->MaxX();
        if (maxX(1)<Particles[i]->MaxY()&&Particles[i]->IsFree()) maxX(1) = Particles[i]->MaxY();
        if (maxX(2)<Particles[i]->MaxZ()&&Particles[i]->IsFree()) maxX(2) = Particles[i]->MaxZ();
    }
}

inline void Domain::Center(Vec3_t C)
{
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    Vec3_t Transport(-0.5*(maxX+minX));
    Transport += C;
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Translate(Transport);
}

inline void Domain::ClearInteractons()
{
    // delete old interactors
    for (size_t i=0; i<CInteractons.Size(); ++i)
    {
        if (CInteractons[i]!=NULL) delete CInteractons[i];
    }
    CInteractons.Resize(0);
    for (size_t i=0; i<BInteractons.Size(); ++i)
    {
        if (BInteractons[i]!=NULL) delete BInteractons[i];
    }
    BInteractons.Resize(0);
    Interactons.Resize(0);
    Listofpairs.clear();

    for (size_t i=0; i<CPInteractons.Size(); ++i)
    {
        if (CPInteractons[i]!=NULL) delete CPInteractons[i];
    }
    CPInteractons.Resize(0);
    PInteractons.Resize(0);
    PListofpairs.clear();
}

inline void Domain::ResetInteractons()
{
    // delete old interactors
    for (size_t i=0; i<CInteractons.Size(); ++i)
    {
        if (CInteractons[i]!=NULL) delete CInteractons[i];
    }

    // new interactors
    CInteractons.Resize(0);
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            bool pj_has_vf = !Particles[j]->IsFree();


            // if both particles have any component specified or they are far away, don't create any intereactor
            bool close = (Distance(Particles[i]->x,Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close ) continue;
            Listofpairs.insert(make_pair(Particles[i],Particles[j]));

            // if both particles are spheres (just one vertex)
            if (Particles[i]->Verts.Size()==1 && Particles[j]->Verts.Size()==1)
            {
                CInteractons.Push (new CInteractonSphere(Particles[i],Particles[j]));
            }

            // normal particles
            else
            {
                CInteractons.Push (new CInteracton(Particles[i],Particles[j]));
            }
        }
    }
}

inline void Domain::ResetDisplacements()
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->ResetDisplacements();
    }
}

inline double Domain::MaxDisplacement()
{
    double md = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        double mpd = Particles[i]->MaxDisplacement();
        if (mpd > md) md = mpd;
    }
    return md;
}

inline void Domain::ResetContacts()
{   
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            bool pj_has_vf = !Particles[j]->IsFree();

            bool close = (Distance(Particles[i]->x,Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close) continue;
            
            // checking if the interacton exist for that pair of particles
            set<pair<Particle *, Particle *> >::iterator it = Listofpairs.find(make_pair(Particles[i],Particles[j]));
            if (it != Listofpairs.end())
            {
                continue;
            }
            Listofpairs.insert(make_pair(Particles[i],Particles[j]));
            
            // if both particles are spheres (just one vertex)
            if (Particles[i]->Verts.Size()==1 && Particles[j]->Verts.Size()==1)
            {
                CInteractons.Push (new CInteractonSphere(Particles[i],Particles[j]));
            }

            // normal particles
            else
            {
                CInteractons.Push (new CInteracton(Particles[i],Particles[j]));
            }
        }
    }
    Interactons.Resize(0);
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        if(CInteractons[i]->UpdateContacts(Alpha)) Interactons.Push(CInteractons[i]);
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        if(BInteractons[i]->UpdateContacts(Alpha)) Interactons.Push(BInteractons[i]);
    }
}

inline void Domain::ResetBoundaries()
{
    ParXmax.Resize(0);
    for (size_t i=0;i<Particles.Size();i++)
    {
        Particle * Pa = Particles[i];
        if ((Pa->MinX()>Xmax)&&Pa->IsFree())
        {
            Vec3_t v(Xmin-Xmax,0.0,0.0);
            Pa->Translate(v);
        }
        if ((Pa->MaxX()<Xmin)&&Pa->IsFree())
        {
            Vec3_t v(Xmax-Xmin,0.0,0.0);
            Pa->Translate(v);
        }
        if ((Pa->MaxX()>Xmax-2.0*Alpha-2.0*MaxDmax)&&Pa->IsFree())
        {
            ParXmax.Push(Pa);
        }
    }
    
    for (size_t i=0;i<ParXmax.Size();i++) 
    {
        Particle * P1 = ParXmax[i];
        Vec3_t v(Xmin-Xmax,0.0,0.0);
        P1->Translate(v);
        for (size_t j=0; j<Particles.Size(); j++)
        {
            Particle * P2 = Particles[j];
            if (P1==P2||!P2->IsFree()||ParXmax.Has(P2)) continue;
            bool close = (Distance(P1->x,P2->x)<=P1->Dmax+P2->Dmax+2*Alpha);
            if (!close) continue;
            set<pair<Particle *, Particle *> >::iterator it = PListofpairs.find(make_pair(P1,P2));
            if (it != PListofpairs.end())
            {
                continue;
            }
            PListofpairs.insert(make_pair(P1,P2));
            // if both particles are spheres (just one vertex)
            if (P1->Verts.Size()==1 && P2->Verts.Size()==1)
            {
                CPInteractons.Push (new CInteractonSphere(P1,P2));
            }

            // normal particles
            else
            {
                CPInteractons.Push (new CInteracton(P1,P2));
            }
        }
    }
    PInteractons.Resize(0);
    for (size_t i=0; i<CPInteractons.Size(); i++)
    {
        if(CPInteractons[i]->UpdateContacts(Alpha)) PInteractons.Push(CPInteractons[i]);
    }
    Vec3_t vmin(Xmax-Xmin,0.0,0.0);
    for (size_t i=0;i<ParXmax.Size();i++) ParXmax[i]->Translate(vmin);
}
// Auxiliar methods

inline void Domain::LinearMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        L += Particles[i]->Props.m*Particles[i]->v;
    }
}

inline void Domain::AngularMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Vec3_t t1,t2;
        t1 = Particles[i]->I(0)*Particles[i]->w(0),Particles[i]->I(1)*Particles[i]->w(1),Particles[i]->I(2)*Particles[i]->w(2);
        Rotation (t1,Particles[i]->Q,t2);
        L += Particles[i]->Props.m*cross(Particles[i]->x,Particles[i]->v)+t2;
    }
}

inline double Domain::CalcEnergy (double & Ekin, double & Epot)
{
    // kinematic energy
    Ekin = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Ekin += 0.5*Particles[i]->Props.m*dot(Particles[i]->v,Particles[i]->v)
                + 0.5*(Particles[i]->I(0)*Particles[i]->w(0)*Particles[i]->w(0)
                      +Particles[i]->I(1)*Particles[i]->w(1)*Particles[i]->w(1)
                      +Particles[i]->I(2)*Particles[i]->w(2)*Particles[i]->w(2));
    }

    // potential energy
    Epot = 0.0;
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        Epot += CInteractons[i]->Epot;
    }

    // total energy
    return Ekin + Epot;
}

inline void Domain::EnergyOutput (size_t IdxOut, std::ostream & OF)
{
    // header
    if (IdxOut==0)
    {
        OF << Util::_10_6 << "Time" << Util::_8s << "Ekin" << Util::_8s << "Epot" << Util::_8s << "Evis" << Util::_8s << "Efric" << Util::_8s << "Wext" << std::endl;
    }
    double Ekin,Epot;
    CalcEnergy(Ekin,Epot);
    OF << Util::_10_6 << Time << Util::_8s << Ekin << Util::_8s << Epot << Util::_8s << Evis << Util::_8s << Efric << Util::_8s << Wext << std::endl;
}

inline void Domain::GetGSD (Array<double> & X, Array<double> & Y, Array<double> & D, size_t NDiv) const
{
    // calc GSD information
    Array<double> Vg;
    double Vs = 0.0;

    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particle * P = Particles[i];
        double Diam = sqrt((P->MaxX()-P->MinX())*(P->MaxX()-P->MinX())+(P->MaxY()-P->MinY())*(P->MaxY()-P->MinY())+(P->MaxZ()-P->MinZ())*(P->MaxX()-P->MinX()));
        Vs += Particles[i]->Props.V;
        Vg.Push(Particles[i]->Props.V);
        D.Push(Diam);
    }
    double Dmin  = D[D.TheMin()]; // minimum diameter
    double Dmax  = D[D.TheMax()]; // maximum diameter
    double Dspan = (Dmax-Dmin)/NDiv;
    for (size_t i=0; i<=NDiv; i++)
    {
        X.Push (i*Dspan+Dmin);
        double cumsum = 0;
        for (size_t j=0; j<D.Size(); j++)
        {
            if (D[j]<=i*Dspan+Dmin) cumsum++;
        }
        Y.Push (cumsum/Particles.Size());
    }
}

inline void Domain::Clusters ()
{
    Array<int> connections;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (BInteractons[i]->valid)
        {
            connections.Push(BInteractons[i]->P1->Index);
            connections.Push(BInteractons[i]->P2->Index);
        }
    }

    Util::Tree tree(connections);
    tree.GetClusters(Listofclusters);
    for (size_t i=0;i<Listofclusters.Size();i++)
    {
        for (size_t j=0;j<Listofclusters[i].Size();j++)
        {
            Particles[Listofclusters[i][j]]->Cluster = i;
        }
    }
}

inline void Domain::DelParticles (Array<int> const & Tags)
{
    Array<int> idxs; // indices to be deleted
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        for (size_t j=0; j<Tags.Size(); ++j)
        {
            if (Particles[i]->Tag==Tags[j]) idxs.Push(i);
        }
    }
    if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particle to be deleted");
    Particles.DelItems (idxs);
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        Particles[i]->Index = i;
    }
}

inline Particle * Domain::GetParticle (int Tag, bool Check)
{
    size_t idx   = 0;
    size_t count = 0;
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag)
        {
            if (!Check) return Particles[i];
            idx = i;
            count++;
        }
    }
    if      (count==0) throw new Fatal("Domain::GetParticle: Could not find Particle with Tag==%d",Tag);
    else if (count>1)  throw new Fatal("Domain::GetParticle: There are more than one particle with Tag==%d",Tag);
    return Particles[idx];
}

inline Particle const & Domain::GetParticle (int Tag, bool Check) const
{
    size_t idx   = 0;
    size_t count = 0;
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag)
        {
            if (!Check) return (*Particles[i]);
            idx = i;
            count++;
        }
    }
    if      (count==0) throw new Fatal("Domain::GetParticle: Could not find Particle with Tag==%d",Tag);
    else if (count>1)  throw new Fatal("Domain::GetParticle: There are more than one particle with Tag==%d",Tag);
    return (*Particles[idx]);
}

inline void Domain::GetParticles (int Tag, Array<Particle*> & P)
{
    P.Resize(0);
    for (size_t i=0; i<Particles.Size(); ++i)
    {
        if (Particles[i]->Tag==Tag) P.Push(Particles[i]);
    }
    if (P.Size()==0) throw new Fatal("Domain::GetParticles: Could not find any Particle with Tag==%d",Tag);
}

}; // namespace DEM

#endif // MECHSYS_DEM_DOMAIN_H
