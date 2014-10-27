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


#ifndef MECHSYS_LBM_DOMAIN_H
#define MECHSYS_LBM_DOMAIN_H

// STD
#include <map>
#include <vector>
#include <utility>
#include <set>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/lbm/Lattice.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/util/util.h>
#include <mechsys/util/maps.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/util/tree.h>

using std::set;
using std::map;
using std::pair;
using std::make_pair;

namespace LBM
{

struct ParticleCellPair
{
    size_t ICell;         ///< Index of the cell
    size_t IPar;          ///< Index of the particle
    Array<size_t> IGeo;   ///< Array of index of the geometric feature
};

struct MtData;

class Domain
{
public:
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructors
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    Array<double>         nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Special constructor with only one component, the parameters are the same as above
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    double                nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step
    
    //Utility methods
    void BoundingBox       (Vec3_t & minX, Vec3_t & maxX);                                                      ///< Defines the rectangular box that encloses the particles.
    void Center            (Vec3_t C = Vec3_t(0.0,0.0,0.0));                                                    ///< Centers the domain around C
    void SetProps          (Dict & D);                                                                          ///< Set the properties of individual grains by dictionaries
    void Clusters          ();                                                                                  ///< Check the bounded particles in the domain and how many connected clusters are still present
    
    //Methods for adding particles
    void AddSphere       (int Tag, Vec3_t const & X, double R, double rho);                                                                              ///< Add sphere
    void GenSpheresBox   (int Tag, Vec3_t const & X0, Vec3_t const & X1, double R, double rho, size_t Randomseed, double fraction, double RminFraction); ///< Create an array of spheres
    void AddCube         (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);                                ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddOcta         (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);                                ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddTetra        (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);                                ///< Add a octahedron at position X with spheroradius R, side of length L and density rho
    void AddPlane        (int Tag, Vec3_t const & X, double R, double Lx,double Ly , double rho, double Angle=0, Vec3_t * Axis=NULL);                    ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void GenBox          (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, double rho, bool Cohesion=false);                        ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBoundingBox  (int InitialTag, double R, double Cf, double rho, bool Cohesion=false);                                                         ///< Generate o bounding box enclosing the previous included particles.
    void GenFromMesh     (Mesh::Generic & M, double R, double rho, bool cohesion=false, bool MC=true, double thickness = 0.0);                           ///< Generate particles from a FEM mesh generator
    void AddVoroCell     (int Tag, voro::voronoicell & VC, double R, double rho, bool Erode, Vec3_t nv = iVec3_t(1.0,1.0,1.0));                          ///< Add a single voronoi cell, it should be built before tough
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bool Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                                                ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bVec3_t Periodic,size_t Randomseed, double fraction, Vec3_t q = OrthoSys::O);                                             ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni, Periodic conditions are chosen for each particle
    // Access methods
    DEM::Particle       * GetParticle  (int Tag, bool Check=true);       ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    DEM::Particle const & GetParticle  (int Tag, bool Check=true) const; ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    void                 DelParticles  (Array<int> const & Tags);        ///< Delete particle
    
    //Methods
#ifdef USE_HDF5
    void WriteBF           (char const * FileKey);  ///< Save a h5 with branch and force information
    void WriteFrac         (char const * FileKey);  ///< Save a xdmf file for fracture visualization
    void WriteXDMF         (char const * FileKey);  ///< Write the domain data in xdmf file
    void Load              (char const * FileKey);  ///< Load particle data from Mechsys DEM
#endif
    void UpdateLinkedCells ();                                                                                  ///< Update the linked cells

    void Initialize     (double dt=0.0);                                                                                              ///< Set the particles to a initial state and asign the possible insteractions
    void ApplyForce     (size_t n = 0, size_t Np = 1, bool MC=false);                                                                 ///< Apply the interaction forces and the collision operator
    void Collide        (size_t n = 0, size_t Np = 1);                                                                                ///< Apply the interaction forces and the collision operator
    void ImprintLattice (size_t n = 0, size_t Np = 1);                                                                                ///< Imprint the DEM particles into the lattices
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                                ///< Solve the Domain dynamics
    void ResetContacts();                                                                                                             ///< Reset contacts for verlet method DEM
    void ResetDisplacements();                                                                                                        ///< Reset the displacements for the verlet method DEM
    double  MaxDisplacement();                                                                                                        ///< Give the maximun displacement of DEM particles

#ifdef USE_THREAD
    Array<pair<size_t, size_t> >                ListPosPairs;         ///< List of all possible particles pairs
#elif USE_OMP
    Array<pair<size_t, size_t> >                ListPosPairs;         ///< List of all possible particles pairs
    iVec3_t                                     LCellDim;             ///< Dimensions of the linked cell array
    Array<Array <size_t> >                      LinkedCell;           ///< Linked Cell array for optimization.
    Vec3_t                                      LCxmin;               ///< Bounding box low   limit for the linked cell array
    Vec3_t                                      LCxmax;               ///< Bounding box upper limit for the linked cell array
#endif
    //Data
    bool                                         Initialized;         ///< System (particles and interactons) initialized ?
    bool                                              PrtVec;         ///< Print Vector data into the xdmf-h5 files
    bool                                            Finished;         ///< Has the simulation finished
    bool                                              Dilate;         ///< True if eroded particles should be dilated for visualization
    Array<size_t>                                    FreePar;         ///< Particles that are free
    Array<size_t>                                  NoFreePar;         ///< Particles that are not free
    String                                           FileKey;         ///< File Key for output files
    Array <Lattice>                                      Lat;         ///< Fluid Lattices
    Array <DEM::Particle *>                        Particles;         ///< Array of Disks
    Array <DEM::Interacton *>                    Interactons;         ///< Array of insteractons
    Array <DEM::CInteracton *>                  CInteractons;         ///< Array of valid  collision interactons
    Array <DEM::BInteracton *>                  BInteractons;         ///< Cohesion interactons
    Array <iVec3_t>                                CellPairs;         ///< Pairs of cells
    Array <ParticleCellPair>                    ParCellPairs;         ///< Pairs of cells and particles
    set<pair<DEM::Particle *, DEM::Particle *> > Listofpairs;         ///< List of pair of particles associated per interacton for memory optimization
    double                                              Time;         ///< Time of the simulation
    double                                                dt;         ///< Timestep
    double                                             Alpha;         ///< Verlet distance
    double                                              Beta;         ///< Binmultiplier
    double                                              Gmix;         ///< Interaction constant for the mixture
    double                                            Voltot;         ///< toala particle volume
    double                                           MaxDmax;         ///< Maximun value for the radious of the spheres surronding each particle
    double                                                Ms;         ///< Total mass of the particles
    double                                                Vs;         ///< Volume occupied by the grains
    double                                                Sc;         ///< Smagorinsky constant
    double                                             Fconv;         ///< Force conversion factor
    Array <double>                                       EEk;         ///< Diadic velocity tensor trace
    void *                                          UserData;         ///< User Data
    size_t                                             Nproc;         ///< Number of cores for multithreading
    size_t                                           idx_out;         ///< The discrete time step
    size_t                                              Step;         ///< The space step to reduce the size of the h5 file for visualization
    Array<Array <int> >                       Listofclusters;         ///< List of particles belonging to bounded clusters (applies only for cohesion simulations)
    MtData *                                             MTD;         ///< Multithread data
};

struct MtData
{
    size_t                        ProcRank; ///< Rank of the thread
    size_t                          N_Proc; ///< Total number of threads
    LBM::Domain *                      Dom; ///< Pointer to the lbm domain
    double                             Dmx; ///< Maximun displacement
    double                              dt; ///< Time step
    Array<pair<size_t,size_t> >         LC; ///< A temporal list of new contacts
    Array<size_t>                      LCI; ///< A temporal array of posible Cinteractions
    Array<size_t>                      LCB; ///< A temporal array of posible Binteractions
    Array<ParticleCellPair>            LPC; ///< A temporal array of possible particle cell contacts
    Array<std::pair<iVec3_t,size_t> >  LLC; ///< A temporal array of possible linked cells locations
    Array<std::pair<size_t,size_t> >   LPP; ///< A temporal array of possible partcle types
};

#ifdef USE_THREAD

void * GlobalIni(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].SetZeroGamma(dat.ProcRank, dat.N_Proc);
    }
    Array<DEM::Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    for (size_t i=In;i<Fn;i++)
    {
        (*P)[i]->F = (*P)[i]->Ff;
        (*P)[i]->T = (*P)[i]->Tf;
    }
    return NULL;
}

void * GlobalImprint(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    dat.Dom->ImprintLattice(dat.ProcRank, dat.N_Proc);
    return NULL;
}

void * GlobalForce(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    Array<DEM::Interacton * > * I = &dat.Dom->Interactons;
	size_t Ni = I->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = I->Size() : Fn = (dat.ProcRank+1)*Ni;
	for (size_t i=In;i<Fn;i++)
	{
		if ((*I)[i]->CalcForce(dat.dt))
        {
            dat.Dom->WriteXDMF("error");
            std::cout << "Maximun overlap detected between particles at time " << dat.Dom->Time << std::endl;
            sleep(1);
            throw new Fatal("Maximun overlap detected between particles");
        }
        if (isnan(norm((*I)[i]->F1)))
        {
            std::cout << "Calculated force is not a number at time " << dat.Dom->Time    << std::endl;
            std::cout << "Particle 1 " << (*I)[i]->P1->Index << " "  << (*I)[i]->P1->Tag << std::endl;
            std::cout << *(*I)[i]->P1 << std::endl;
            std::cout << "Particle 2 " << (*I)[i]->P2->Index << " "  << (*I)[i]->P2->Tag << std::endl;
            std::cout << *(*I)[i]->P2 << std::endl;
            sleep(1);
            throw new Fatal("Force is not a number detected between particles");
        }


#ifdef USE_THREAD
        pthread_mutex_lock  (&(*I)[i]->P1->lck);
#endif
        (*I)[i]->P1->F += (*I)[i]->F1;
        (*I)[i]->P1->T += (*I)[i]->T1;
#ifdef USE_THREAD
        pthread_mutex_unlock(&(*I)[i]->P1->lck);
        pthread_mutex_lock  (&(*I)[i]->P2->lck);
#endif
        (*I)[i]->P2->F += (*I)[i]->F2;
        (*I)[i]->P2->T += (*I)[i]->T2;
#ifdef USE_THREAD
        pthread_mutex_unlock(&(*I)[i]->P2->lck);
#endif
	}
    return NULL;
}

void * GlobalMove(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    Array<DEM::Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.Dmx = 0.0;
	for (size_t i=In;i<Fn;i++)
	{
		(*P)[i]->Translate(dat.dt);
		(*P)[i]->Rotate   (dat.dt);
        if ((*P)[i]->MaxDisplacement()>dat.Dmx) dat.Dmx = (*P)[i]->MaxDisplacement();
	}
    return NULL;
}

void * GlobalApplyForce (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    bool MC = false;
    if (dat.Dom->Lat.Size()==2)
    {
        if (fabs(dat.Dom->Lat[0].G)<1.0e-9&&fabs(dat.Dom->Lat[1].G)<1.0e-9) MC = true;
    }
    dat.Dom->ApplyForce(dat.ProcRank, dat.N_Proc, MC);
    return NULL;
}

void * GlobalCollide (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    dat.Dom->Collide(dat.ProcRank, dat.N_Proc);
    return NULL;
}

void * GlobalBounceBack (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].BounceBack(dat.ProcRank, dat.N_Proc);
    }
    return NULL;
}

void * GlobalStream (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].Stream(dat.ProcRank, dat.N_Proc);
    }
    return NULL;
}

void * GlobalStream1 (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].Stream1(dat.ProcRank, dat.N_Proc);
    }
    return NULL;
}

void * GlobalStream2 (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].Stream2(dat.ProcRank, dat.N_Proc);
    }
    return NULL;
}

void * GlobalResetDisplacement(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    Array<DEM::Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.Dmx = 0.0;
	for (size_t i=In;i<Fn;i++)
    {
        (*P)[i]->ResetDisplacements();
    }
    return NULL;
}

void * GlobalResetContacts1 (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
	size_t Ni = dat.Dom->ListPosPairs.Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = dat.Dom->ListPosPairs.Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LC.Resize(0);

    for (size_t n=In;n<Fn;n++)
    {
        size_t i = dat.Dom->ListPosPairs[n].first;
        size_t j = dat.Dom->ListPosPairs[n].second;
        bool pi_has_vf = !dat.Dom->Particles[i]->IsFree();
        bool pj_has_vf = !dat.Dom->Particles[j]->IsFree();

        bool close = (norm(dat.Dom->Particles[i]->x-dat.Dom->Particles[j]->x)<=dat.Dom->Particles[i]->Dmax+dat.Dom->Particles[j]->Dmax+2*dat.Dom->Alpha);
        if ((pi_has_vf && pj_has_vf) || !close) continue;
        
        // checking if the interacton exist for that pair of particles
        set<pair<DEM::Particle *, DEM::Particle *> >::iterator it = dat.Dom->Listofpairs.find(std::make_pair(dat.Dom->Particles[i],dat.Dom->Particles[j]));
        if (it != dat.Dom->Listofpairs.end())
        {
            continue;
        }
        dat.LC.Push(make_pair(i,j));
    }

    Array<DEM::Particle * > * P = &dat.Dom->Particles;
	Ni = P->Size()/dat.N_Proc;
    In = dat.ProcRank*Ni;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
//
	//Ni = dat.Dom->Lat[0].Cells.Size()/dat.N_Proc;
    //In = dat.ProcRank*Ni;
    //dat.ProcRank == dat.N_Proc-1 ? Fn = dat.Dom->Lat[0].Cells.Size() : Fn = (dat.ProcRank+1)*Ni;
    
    dat.LPC.Resize(0);

    if (dat.Dom->Lat[0].Ndim(2)==1)
    {
        // TODO 2D case
    }
    else
    {
	    for (size_t i=In;i<Fn;i++)
        {
            //Cell  * cell = dat.Dom->Lat[0].Cells[i];
            DEM::Particle * Pa = (*P)[i];
            //std::cout << std::max(0.0,double(Pa->x(0)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " 
                      //<< std::max(0.0,double(Pa->x(1)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " 
                      //<< std::max(0.0,double(Pa->x(2)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " << std::endl;
            //std::cout << std::min(double(dat.Dom->Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " "
                      //<< std::min(double(dat.Dom->Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " "
                      //<< std::min(double(dat.Dom->Lat[0].Ndim(2)-1),double(Pa->x(2)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " << std::endl;

            for (size_t n=std::max(0.0,double(Pa->x(0)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx);n<=std::min(double(dat.Dom->Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx);n++)
            for (size_t m=std::max(0.0,double(Pa->x(1)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx);m<=std::min(double(dat.Dom->Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx);m++)
            for (size_t l=std::max(0.0,double(Pa->x(2)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx);l<=std::min(double(dat.Dom->Lat[0].Ndim(2)-1),double(Pa->x(2)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx);l++)
            //for (size_t n=0;n<dat.Dom->Particles.Size();n++)
            {
                //DEM::Particle * Pa = dat.Dom->Particles[n];
                Cell  * cell = dat.Dom->Lat[0].GetCell(iVec3_t(n,m,l));
                double x     = dat.Dom->Lat[0].dx*(cell->Index(0));
                double y     = dat.Dom->Lat[0].dx*(cell->Index(1));
                double z     = dat.Dom->Lat[0].dx*(cell->Index(2));
                Vec3_t  C(x,y,z);
                if ((norm(C-Pa->x)>2.0*dat.Dom->Alpha+2.0*dat.Dom->Lat[0].dx+Pa->Dmax)||(Pa->Bdry==true)) continue;
                ParticleCellPair NewPCP;
                NewPCP.IPar = i;
                //NewPCP.IPar = n;
                NewPCP.ICell= cell->ID;
                bool valid = false;
                if (Pa->Faces.Size()>0)
                {
                    for (size_t j=0;j<Pa->Faces.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Faces[j])<2.0*dat.Dom->Alpha+Pa->Props.R)
                        {
                            if (Pa->Faces[j]->Area()<2.0*M_PI*Pa->Props.R*Pa->Props.R)
                            {
                                //std::cout << "it did" << std::endl;
                                continue;
                            }
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                    }
                }
                else if (Pa->Edges.Size()>0)
                {
                    for (size_t j=0;j<Pa->Edges.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Edges[j])<2.0*dat.Dom->Alpha+Pa->Props.R) 
                        {
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                           
                    }
                }
                else if (Pa->Verts.Size()>0)
                {
                    for (size_t j=0;j<Pa->Verts.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Verts[j])<2.0*dat.Dom->Alpha+Pa->Props.R)
                        {
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                    }
                }
                if (Pa->IsInsideFaceOnly(C)) valid = true;

                if (valid) dat.LPC.Push(NewPCP);
            }
        }
    }
    return NULL;
}

void * GlobalResetContacts2 (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
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

#endif

inline Domain::Domain(LBMethod Method, Array<double> nu, iVec3_t Ndim, double dx, double Thedt)
{
    Initialized = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (nu.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one fluid please");
    if (Ndim(2) >1&&Method==D2Q9)  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (Ndim(2)==1&&Method==D3Q15) throw new Fatal("LBM::Domain: Ndim(2) is 1. Either change the method to D2Q9 or increase the z-dimension");
    for (size_t i=0;i<nu.Size();i++)
    {
        Lat.Push(Lattice(Method,nu[i],Ndim,dx,Thedt));
    }
    Time   = 0.0;
    dt     = Thedt;
    Alpha  = 10.0;
    Beta   = 1.0;
    Step   = 1;
    Sc     = 0.17;
    PrtVec = true;
    Fconv  = 1.0;


    EEk.Resize(Lat[0].Cells[0]->Nneigh);
    for (size_t k=0;k<Lat[0].Cells[0]->Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(Lat[0].Cells[0]->C[k][n]*Lat[0].Cells[0]->C[k][m]);
        }
    }

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Size()*Lat[0].Ncells,TERM_RST);
}

inline Domain::Domain(LBMethod Method, double nu, iVec3_t Ndim, double dx, double Thedt)
{
    Initialized = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    Lat.Push(Lattice(Method,nu,Ndim,dx,Thedt));
    if (Ndim(2) >1&&Method==D2Q9)  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (Ndim(2)==1&&Method==D3Q15) throw new Fatal("LBM::Domain: Ndim(2) is greater than 1. Either change the method to D2Q9 or increse the z-dimension");
    Time   = 0.0;
    dt     = Thedt;
    Alpha  = 10.0;
    Beta   = 1.0;
    Step   = 1;
    Sc     = 0.17;
    PrtVec = true;
    Fconv  = 1.0;

    EEk.Resize(Lat[0].Cells[0]->Nneigh);
    for (size_t k=0;k<Lat[0].Cells[0]->Nneigh;k++)
    {
        EEk[k]    = 0.0;
        for (size_t n=0;n<3;n++)
        for (size_t m=0;m<3;m++)
        {
            EEk[k] += fabs(Lat[0].Cells[0]->C[k][n]*Lat[0].Cells[0]->C[k][m]);
        }
    }

    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Size()*Lat[0].Ncells,TERM_RST);
}

#ifdef USE_HDF5

inline void Domain::WriteBF (char const * FileKey)
{

    size_t n_fn = 0;

    for (size_t i=0;i<CInteractons.Size();i++)
    {
        //if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree())) n_fn++;
        if (norm(CInteractons[i]->Fnet)>1.0e-12) n_fn++;
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

    // Saving Collision forces
    for (size_t i=0;i<CInteractons.Size();i++)
    {
        //if ((norm(CInteractons[i]->Fnet)>1.0e-12)&&(CInteractons[i]->P1->IsFree()&&CInteractons[i]->P2->IsFree()))
        if (norm(CInteractons[i]->Fnet)>1.0e-12)
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

    //Saving Cohesive forces
    if (BInteractons.Size()>0)
    {
    float  *  Bnnet = new float[3*BInteractons.Size()];
    float  *  Btnet = new float[3*BInteractons.Size()];
    float  *  BOrig = new float[3*BInteractons.Size()];
    int    *  BVal  = new   int[  BInteractons.Size()];
    int    *  BID1  = new   int[  BInteractons.Size()];
    int    *  BID2  = new   int[  BInteractons.Size()];

    idx = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        Bnnet [3*idx  ] = float(BInteractons[i]->Fnet  (0));
        Bnnet [3*idx+1] = float(BInteractons[i]->Fnet  (1));
        Bnnet [3*idx+2] = float(BInteractons[i]->Fnet  (2));
        Btnet [3*idx  ] = float(BInteractons[i]->Ftnet (0));
        Btnet [3*idx+1] = float(BInteractons[i]->Ftnet (1));
        Btnet [3*idx+2] = float(BInteractons[i]->Ftnet (2));
        BOrig [3*idx  ] = float(BInteractons[i]->xnet(0));
        BOrig [3*idx+1] = float(BInteractons[i]->xnet(1)); 
        BOrig [3*idx+2] = float(BInteractons[i]->xnet(2)); 
        BVal  [idx]     = int  (BInteractons[i]->valid);
        BID1  [idx]     = int  (BInteractons[i]->P1->Index);
        BID2  [idx]     = int  (BInteractons[i]->P2->Index);
        idx++;

    }
    hsize_t dims[1];
    dims[0] = 3*BInteractons.Size();
    String dsname;
    dsname.Printf("BNormal");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Bnnet );
    dsname.Printf("BTangential");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Btnet );
    dsname.Printf("BPosition");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,BOrig );
    dims[0] = BInteractons.Size();
    dsname.Printf("BVal");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BVal  );
    dsname.Printf("BID1");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BID1  );
    dsname.Printf("BID2");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,BID2  );

    delete [] Bnnet;
    delete [] Btnet;
    delete [] BOrig;
    delete [] BVal ;
    delete [] BID1 ;
    delete [] BID2 ;


    }


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
    if (BInteractons.Size()>0)
    {
    oss << "   <Grid Name=\"CohesiveForce\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << BInteractons.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << BInteractons.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/BPosition \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"BNormal\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BNormal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BTangential\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BTangential \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BVal\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BVal \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BID1\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BID1 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"BID2\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << BInteractons.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/BID2 \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Domain::WriteFrac (char const * FileKey)
{

    // Counting the number of non valid cohesive interactons
    size_t N_Faces = 0;
    size_t N_Verts = 0;
    size_t nvbi = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (!BInteractons[i]->valid) 
        {
            DEM::Particle * P1 = BInteractons[i]->P1;
            DEM::Particle * P2 = BInteractons[i]->P2;
            DEM::Face     * F1 = P1->Faces[BInteractons[i]->IF1];
            DEM::Face     * F2 = P2->Faces[BInteractons[i]->IF2];
            nvbi++;
            N_Faces += F1->Edges.Size();
            N_Faces += F2->Edges.Size();
            N_Verts += F1->Edges.Size() + 1;
            N_Verts += F2->Edges.Size() + 1;
        }
    }

    //std::cout << "1 " << nvbi << std::endl;

    if (nvbi==0) return;

    //Geometric information
    float  * Verts   = new float [3*N_Verts];
    int    * FaceCon = new int   [3*N_Faces];

    //Atributes
    int    * Tags    = new int   [  N_Faces];
       
    size_t n_verts = 0;
    size_t n_faces = 0;
    size_t n_attrs = 0;
    for (size_t i=0;i<BInteractons.Size();i++)
    {
        if (BInteractons[i]->valid) continue;
        DEM::Particle * P1 = BInteractons[i]->P1;
        DEM::Particle * P2 = BInteractons[i]->P2;
        size_t    IF1 = BInteractons[i]->IF1;
        size_t    IF2 = BInteractons[i]->IF2;
        //std::cout << P1 << " " << P2 <<std::endl;

        //For P1
        {
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(P1->Verts.Size());
            Array<Vec3_t> Vres (P1->Verts.Size());
            for (size_t j=0;j<P1->Verts.Size();j++)
            {
                Vtemp[j] = *P1->Verts[j];
                Vres [j] = *P1->Verts[j];
            }
            double multiplier = 0.0;
            if (P1->Eroded&&P1->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,P1->EdgeCon,P1->FaceCon,Vres,P1->Props.R);
                multiplier = 1.0;
            }
            for (size_t j=0;j<P1->FaceCon[IF1].Size();j++)
            {
                size_t k = P1->FaceCon[IF1][j];
                Verts[n_verts++] = float(Vres[k](0));
                Verts[n_verts++] = float(Vres[k](1));
                Verts[n_verts++] = float(Vres[k](2));
            }
            size_t n_reff = n_verts/3;
            Vec3_t C,N;
            P1->Faces[IF1]->Centroid(C);
            P1->Faces[IF1]->Normal(N);
            Verts[n_verts++] = float(C(0) + multiplier*P1->Props.R*N(0));
            Verts[n_verts++] = float(C(1) + multiplier*P1->Props.R*N(1));
            Verts[n_verts++] = float(C(2) + multiplier*P1->Props.R*N(2));
            for (size_t j=0;j<P1->FaceCon[IF1].Size();j++)
            {
                FaceCon[n_faces++] = int(n_reff);  
                FaceCon[n_faces++] = int(n_refv + j);
                FaceCon[n_faces++] = int(n_refv + (j+1)%(P1->FaceCon[IF1].Size()));
                Tags   [n_attrs]   = int(P1->Tag);
                n_attrs++;
            }
        }
        //std::cout << "2" << std::endl;
        //For P2
        {
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(P2->Verts.Size());
            Array<Vec3_t> Vres (P2->Verts.Size());
            for (size_t j=0;j<P2->Verts.Size();j++)
            {
                Vtemp[j] = *P2->Verts[j];
                Vres [j] = *P2->Verts[j];
            }
            //std::cout << "3" << std::endl;
            double multiplier = 0.0;
            if (P2->Eroded&&P2->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,P2->EdgeCon,P2->FaceCon,Vres,P2->Props.R);
                multiplier = 1.0;
            }
            //std::cout << "4" << std::endl;
            for (size_t j=0;j<P2->FaceCon[IF2].Size();j++)
            {
                size_t k = P2->FaceCon[IF2][j];
                Verts[n_verts++] = float(Vres[k](0));
                Verts[n_verts++] = float(Vres[k](1));
                Verts[n_verts++] = float(Vres[k](2));
            }
            //std::cout << "5" << std::endl;
            size_t n_reff = n_verts/3;
            Vec3_t C,N;
            P2->Faces[IF2]->Centroid(C);
            P2->Faces[IF2]->Normal(N);
            Verts[n_verts++] = float(C(0) + multiplier*P2->Props.R*N(0));
            Verts[n_verts++] = float(C(1) + multiplier*P2->Props.R*N(1));
            Verts[n_verts++] = float(C(2) + multiplier*P2->Props.R*N(2));
            //std::cout << "6" << std::endl;
            for (size_t j=0;j<P2->FaceCon[IF2].Size();j++)
            {
                FaceCon[n_faces++] = int(n_reff);  
                FaceCon[n_faces++] = int(n_refv + j);
                FaceCon[n_faces++] = int(n_refv + (j+1)%(P2->FaceCon[IF2].Size()));
                Tags   [n_attrs]   = int(P2->Tag);
                n_attrs++;
            }
            //std::cout << "7" << std::endl;
        }
    }
    //std::cout << n_faces << " " << N_Faces << std::endl;
    //Write the data
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
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

    //Erasing the data
    delete [] Verts;
    delete [] FaceCon;
    delete [] Tags;

    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
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
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";

    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Lat[0].Ndim[0]/Step;
    size_t  Ny = Lat[0].Ndim[1]/Step;
    size_t  Nz = Lat[0].Ndim[2]/Step;
    for (size_t j=0;j<Lat.Size();j++)
    {
        // Creating data sets
        float * Density   = new float[  Nx*Ny*Nz];
        float * Gamma     = new float[  Nx*Ny*Nz];
        float * Vvec      = new float[3*Nx*Ny*Nz];

        size_t i=0;
        for (size_t m=0;m<Lat[0].Ndim(2);m+=Step)
        for (size_t l=0;l<Lat[0].Ndim(1);l+=Step)
        for (size_t n=0;n<Lat[0].Ndim(0);n+=Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            Vec3_t vel    = OrthoSys::O;

            for (size_t ni=0;ni<Step;ni++)
            for (size_t li=0;li<Step;li++)
            for (size_t mi=0;mi<Step;mi++)
            {
                rho  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Rho;
                gamma+= Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Gamma;
                vel  += Lat[j].GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel;
            }
            rho  /= Step*Step*Step;
            gamma/= Step*Step*Step;
            vel  /= Step*Step*Step;
            Density [i]  = (float) rho;
            Gamma   [i]  = (float) Lat[j].Cells[i]->IsSolid&&Step==1? 1.0: gamma;
            Vvec[3*i  ]  = (float) vel(0);
            Vvec[3*i+1]  = (float) vel(1);
            Vvec[3*i+2]  = (float) vel(2);
            i++;
        }


        //for (size_t i=0;i<Lat[j].Cells.Size();i++)
        //{
            //double rho   = Lat[j].Cells[i]->Rho;
            //Vec3_t vel   = Lat[j].Cells[i]->Vel;
            //Density [i]  = (float) rho;
            //Gamma   [i]  = (float) Lat[j].Cells[i]->IsSolid? 1.0:Lat[j].Cells[i]->Gamma;
            //Vvec[3*i  ]  = (float) vel(0);
            //Vvec[3*i+1]  = (float) vel(1);
            //Vvec[3*i+2]  = (float) vel(2);
        //}

        //Write the data
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Gamma   );
        }
        if (PrtVec)
        {
            dims[0] = 3*Nx*Ny*Nz;
            dsname.Printf("Velocity_%d",j);
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vvec    );
        }
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


    size_t N_Faces = 0;
    size_t N_Verts = 0;
    //Writing particle data
    if (Particles.Size()>0)
    {
        for (size_t i=0; i<Particles.Size(); i++) 
        { 
            for (size_t j=0;j<Particles[i]->Faces.Size();j++)
            {
                N_Faces += Particles[i]->Faces[j]->Edges.Size();
            }
            N_Verts += Particles[i]->Verts.Size() + Particles[i]->Faces.Size();
        }
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
            //float  * Stress  = new float [9*N_Faces];

            size_t n_verts = 0;
            size_t n_faces = 0;
            size_t n_attrs = 0;
            //size_t n_attrv = 0;
            //size_t n_attrt = 0;
            for (size_t i=0;i<Particles.Size();i++)
            {
                DEM::Particle * Pa = Particles[i];
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
                    Verts[n_verts++] = (float) Vres[j](0);
                    Verts[n_verts++] = (float) Vres[j](1);
                    Verts[n_verts++] = (float) Vres[j](2);
                }
                size_t n_reff = n_verts/3;
                for (size_t j=0;j<Pa->FaceCon.Size();j++)
                {
                    Vec3_t C,N;
                    Pa->Faces[j]->Centroid(C);
                    Pa->Faces[j]->Normal(N);
                    Verts[n_verts++] = (float) C(0) + multiplier*Pa->Props.R*N(0);
                    Verts[n_verts++] = (float) C(1) + multiplier*Pa->Props.R*N(1);
                    Verts[n_verts++] = (float) C(2) + multiplier*Pa->Props.R*N(2);
                    //Verts[n_verts++] = (float) C(0);
                    //Verts[n_verts++] = (float) C(1);
                    //Verts[n_verts++] = (float) C(2);
                    for (size_t k=0;k<Pa->FaceCon[j].Size();k++)
                    {
                        size_t nin = Pa->FaceCon[j][k];
                        size_t nen = Pa->FaceCon[j][(k+1)%Pa->FaceCon[j].Size()];
                        FaceCon[n_faces++] = (int) n_reff + j;  
                        FaceCon[n_faces++] = (int) n_refv + nin;
                        FaceCon[n_faces++] = (int) n_refv + nen;

                        //Writing the attributes
                        Tags  [n_attrs] = (int)   Pa->Tag;
                        Clus  [n_attrs] = size_t(Pa->Cluster);
                        Vel   [n_attrs] = (float) norm(Pa->v);
                        Ome   [n_attrs] = (float) norm(Pa->w);
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
            //delete [] Stress;
        }
        //Creating data sets
        float * Radius = new float[  Particles.Size()];
        float * Posvec = new float[3*Particles.Size()];
        float * Velvec = new float[3*Particles.Size()];
        float * Omevec = new float[3*Particles.Size()];
        int   * Tags   = new int  [  Particles.Size()];
        for (size_t i=0;i<Particles.Size();i++)
        {
            Vec3_t Ome;
            Rotation(Particles[i]->w,Particles[i]->Q,Ome);
            Radius[i]     = (float) Particles[i]->Dmax;
            Posvec[3*i  ] = (float) Particles[i]->x(0);
            Posvec[3*i+1] = (float) Particles[i]->x(1);
            Posvec[3*i+2] = (float) Particles[i]->x(2);
            Velvec[3*i  ] = (float) Particles[i]->v(0);
            Velvec[3*i+1] = (float) Particles[i]->v(1);
            Velvec[3*i+2] = (float) Particles[i]->v(2);
            Omevec[3*i  ] = (float)  Ome(0);
            Omevec[3*i+1] = (float)  Ome(1);
            Omevec[3*i+2] = (float)  Ome(2);
            Tags  [i]     = (int)   Particles[i]->Tag;
        }

        hsize_t dims[1];
        dims[0] = 3*Particles.Size();
        String dsname;
        dsname.Printf("Position");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
        dsname.Printf("PVelocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
        dsname.Printf("PAngVel");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Omevec);
        dims[0] = Particles.Size();
        dsname.Printf("Radius");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
        dsname.Printf("PTag");
        H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tags  );


        delete [] Radius;
        delete [] Posvec;
        delete [] Velvec;
        delete [] Omevec;
        delete [] Tags  ;
    }
    
    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    //std::cout << "2" << std::endl;

    if (Lat[0].Ndim(2)==1)
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"2DCoRectMesh\" Dimensions=\"" << Lat[0].Ndim(1) << " " << Lat[0].Ndim(0) << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDY\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"2\"> 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Lat.Size();j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        if (PrtVec)
        {
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    else
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
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*Lat[0].dx << " " << Step*Lat[0].dx  << " " << Step*Lat[0].dx  << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Lat.Size();j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        if (PrtVec)
        {
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        }
        oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        if(Particles.Size()>0)
        {
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
        oss << "   <Grid Name=\"DEM_Center\" GridType=\"Uniform\">\n";
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
        oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/PVelocity\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"AngVel\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/PAngVel\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";
        }
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
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

        Particles.Push (new DEM::Particle(-1,V,E,F,OrthoSys::O,OrthoSys::O,0.1,1.0));

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
        //Particles[Particles.Size()-1]->Index = datint[0];
        Particles[Particles.Size()-1]->Index = Particles.Size()-1;
        int tag[1];
        H5LTread_dataset_int(group_id,"Tag",tag);
        Particles[Particles.Size()-1]->Tag = tag[0];
        Particles[Particles.Size()-1]->PropsReady = true;

    }


    H5Fclose(file_id);
    printf("\n%s--- Done --------------------------------------------%s\n",TERM_CLR2,TERM_RST);
}

#endif

inline void Domain::ApplyForce(size_t n, size_t Np, bool MC)
{
    size_t Ni = CellPairs.Size()/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = CellPairs.Size() : Fn = (n+1)*Ni;

    if (MC)
    {
#ifdef USE_OMP
        In = 0;
        Fn = CellPairs.Size();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
        for (size_t i=In;i<Fn;i++)
        {
            size_t ind1 = CellPairs[i](0);
            size_t ind2 = CellPairs[i](1);
            size_t vec  = CellPairs[i](2);
            Cell * c = Lat[0].Cells[ind1];
            Cell *nb = Lat[1].Cells[ind2];
            double nb_psi,psi,G=Gmix/(dt*dt);
            if (c->IsSolid)
            {
                psi    = 1.0;
                G      = Lat[1].Gs*c->Gs;
            }
            else psi   = c ->Rho;
            if (nb->IsSolid)
            {
                nb_psi = 1.0;
                G      = Lat[0].Gs*nb->Gs;
            }
            else nb_psi = nb->Rho;
            Vec3_t  BF = -G*psi*nb_psi*c->W[vec]*c->C[vec];

#ifdef USE_THREAD
            pthread_mutex_lock  (&c ->lck);
#elif USE_OMP
            omp_set_lock        (&c->lck);
#endif
            c ->BForce += BF;
#ifdef USE_THREAD
            pthread_mutex_unlock(&c ->lck);
            pthread_mutex_lock  (&nb->lck);
#elif USE_OMP
            omp_unset_lock      (&c ->lck);
            omp_set_lock        (&nb->lck);
#endif
            nb->BForce -= BF;
#ifdef USE_THREAD
            pthread_mutex_unlock(&nb->lck);
#elif USE_OMP
            omp_unset_lock      (&nb->lck);
#endif

            c  = Lat[1].Cells[ind1];
            nb = Lat[0].Cells[ind2];
            G  = Gmix/(dt*dt);
            if (c->IsSolid)
            {
                psi    = 1.0;
                G      = Lat[0].Gs;
            }
            else psi   = c ->Rho;
            if (nb->IsSolid)
            {
                nb_psi = 1.0;
                G      = Lat[1].Gs;
            }
            else nb_psi = nb->Rho;
            BF = -G*psi*nb_psi*c->W[vec]*c->C[vec];

#ifdef USE_THREAD
            pthread_mutex_lock  (&c ->lck);
#elif USE_OMP
            omp_set_lock        (&c->lck);
#endif
            c ->BForce += BF;
#ifdef USE_THREAD
            pthread_mutex_unlock(&c ->lck);
            pthread_mutex_lock  (&nb->lck);
#elif USE_OMP
            omp_unset_lock      (&c ->lck);
            omp_set_lock        (&nb->lck);
#endif
            nb->BForce -= BF;
#ifdef USE_THREAD
            pthread_mutex_unlock(&nb->lck);
#elif USE_OMP
            omp_unset_lock      (&nb->lck);
#endif
        }
    }
    else
    {
#ifdef USE_OMP
        In = 0;
        Fn = CellPairs.Size();
        #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
        for (size_t i=In;i<Fn;i++)
        {
            size_t ind1 = CellPairs[i](0);
            size_t ind2 = CellPairs[i](1);
            size_t vec  = CellPairs[i](2);
            for (size_t j=0;j<Lat.Size();j++)
            {
                bool solid = false;
                Cell * c = Lat[j].Cells[ind1];
                double psi;
                if (c->IsSolid||c->Gamma>0.0)
                {
                    psi   = 1.0;
                    solid = true;
                }
                else fabs(Lat[j].G)>1.0e-12 ? psi = Lat[j].Psi(c->Rho) : psi = c->Rho;
                for (size_t k=0;k<Lat.Size();k++)
                {
                    Cell * nb = Lat[k].Cells[ind2];
                    double nb_psi;
                    if (nb->IsSolid||nb->Gamma>0.0)
                    {
                        nb_psi = 1.0;
                        solid  = true;
                    }
                    else fabs(Lat[j].G)>1.0e-12 ? nb_psi = Lat[k].Psi(nb->Rho) : nb_psi = nb->Rho;
                    double G;
                    solid ? G = Lat[j].Gs*2.0*ReducedValue(nb->Gs,c->Gs) : G = Lat[j].G; 
                    Vec3_t BF(OrthoSys::O);
                    if (j==k)
                    {
                        BF += -G*psi*nb_psi*c->W[vec]*c->C[vec];
                    }
                    else if(!solid)
                    {
                        BF += -Gmix*c->Rho*nb->Rho*c->W[vec]*c->C[vec];
                    }
#ifdef USE_THREAD
                    pthread_mutex_lock  (&c ->lck);
#elif USE_OMP
                    omp_set_lock        (&c ->lck);
#endif
                    c ->BForce += BF;
#ifdef USE_THREAD
                    pthread_mutex_unlock(&c ->lck);
                    pthread_mutex_lock  (&nb->lck);
#elif USE_OMP
                    omp_unset_lock      (&c ->lck);
                    omp_set_lock        (&nb->lck);
#endif
                    nb->BForce -= BF;
#ifdef USE_THREAD
                    pthread_mutex_unlock(&nb->lck);
#elif USE_OMP
                    omp_unset_lock      (&nb->lck);
#endif
                }
            }
        }
    }
}

void Domain::Collide (size_t n, size_t Np)
{
	size_t Ni = Lat[0].Ncells/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Lat[0].Ncells : Fn = (n+1)*Ni;
#ifdef USE_OMP
    In = 0;
    Fn = Lat[0].Ncells;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
    for (size_t i=In;i<Fn;i++)
    {
        Vec3_t num(0.0,0.0,0.0);
        double den = 0.0;
        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double tau = Lat[j].Tau;
            num += c->Vel*c->Rho/tau;
            den += c->Rho/tau;
        }
        Vec3_t Vmix = num/den;


        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double rho = c->Rho;
            //if (rho<1.0e-12) continue;
            //if (c->IsSolid||rho<1.0e-12) continue;
            if (fabs(c->Gamma-1.0)<1.0e-12&&fabs(Lat[j].G)>1.0e-12) continue;
            //if (fabs(c->Gamma-1.0)<1.0e-12) continue;
            if (!c->IsSolid)
            {
                //Calculate Smagorinsky LES model
                Array<double> NonEq(c->Nneigh);
                Vec3_t DV  = Vmix + c->BForce*dt/rho;
                double Tau = Lat[j].Tau;
                double Q = 0.0;
                for (size_t k=0;k<c->Nneigh;k++)
                {
                    double FDeqn = c->Feq(k,DV,rho);
                    NonEq[k] = c->F[k] - FDeqn;
                    Q += NonEq[k]*NonEq[k]*EEk[k];
                }
                Q = sqrt(2.0*Q);
                Tau = 0.5*(Tau + sqrt(Tau*Tau + 6.0*Q*Sc/rho));
                //std::cout << Tau << std::endl;

                double Bn;
                rho<10e-12 ? Bn =0.0 : Bn = (c->Gamma*(Lat[j].Tau-0.5))/((1.0-c->Gamma)+(Lat[j].Tau-0.5));
                bool valid  = true;
                double alphal = 1.0;
                double alphat = 1.0;
                size_t num  = 0;
                //double newrho = 0.0;
                while (valid)
                {
                    //newrho = 0.0;
                    valid = false;
                    alphal  = alphat;
                    for (size_t k=0;k<c->Nneigh;k++)
                    {
                        //double FDeqn = c->Feq(k,DV,rho);
                        //c->Ftemp[k] = c->F[k] - alphal*((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]);
                        c->Ftemp[k] = c->F[k] - alphal*((1 - Bn)*(NonEq[k])/Tau - Bn*c->Omeis[k]);
                        //newrho += c->Ftemp[k];
                        if (c->Ftemp[k]<0.0&&num<1)
                        {
                            //double temp = fabs(c->F[k]/((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]));
                            double temp = fabs(c->F[k]/((1 - Bn)*(NonEq[k])/Tau - Bn*c->Omeis[k]));
                            if (temp<alphat) alphat = temp;
                            valid = true;
                        }
                    }
                    num++;
                    if (num>2) 
                    {
                        throw new Fatal("Domain::Collide: Redefine your time step, the current value ensures unstability");
                    }
                }
                for (size_t k=0;k<c->Nneigh;k++)
                {
                    if (std::isnan(c->Ftemp[k]))
                    {
                        c->Gamma = 2.0;
                        #ifdef USE_HDF5
                        WriteXDMF("error");
                        #endif
                        std::cout << "Nan found, resetting" << std::endl;
                        std::cout << c->Density() << " " << c->BForce << " " << num << " " << alphat << " " << c->Index << " " << c->IsSolid << " " << j << " " << k << " " << std::endl;
                        c->Ftemp[k] = 1.0e-15;
                        //c->Ftemp[k] = c->F[k];
                        //throw new Fatal("Domain::Collide: Body force gives nan value, check parameters");
                    }
                    c->F[k] = fabs(c->Ftemp[k]);
                    //c->F[k] = fabs(c->Ftemp[k])*c->Rho/newrho;
                    //std::cout << newrho << std::endl;
                }
            }
            else
            {
                for (size_t j = 1;j<c->Nneigh;j++)
                {
                    c->Ftemp[j] = c->F[j];
                }
                for (size_t j = 1;j<c->Nneigh;j++)
                {
                    c->F[j]     = c->Ftemp[c->Op[j]];
                }
            }
        }
    }   
}

void Domain::ImprintLattice (size_t n,size_t Np)
{
    
	size_t Ni = ParCellPairs.Size()/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = ParCellPairs.Size() : Fn = (n+1)*Ni;
    //std::cout << "Im proccess = " << n << std::endl;
    // 2D imprint
    //if (Lat[0].Ndim(2)==1)
    //{
        //for (size_t i = In;i<Fn;i++)
        //{
            //DEM::Particle * Pa = Particles[i];
            //for (size_t n=std::max(0.0,double(Pa->x(0)-Pa->Props.R-Lat[0].dx)/Lat[0].dx);n<=std::min(double(Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Props.R+Lat[0].dx)/Lat[0].dx);n++)
            //for (size_t m=std::max(0.0,double(Pa->x(1)-Pa->Props.R-Lat[0].dx)/Lat[0].dx);m<=std::min(double(Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Props.R+Lat[0].dx)/Lat[0].dx);m++)
            //{
                //Cell  * cell = Lat[0].GetCell(iVec3_t(n,m,0));
                //double x     = Lat[0].dx*cell->Index(0);
                //double y     = Lat[0].dx*cell->Index(1);
                //double z     = Lat[0].dx*cell->Index(2);
                //Vec3_t  C(x,y,z);
                //Array<Vec3_t> P(4);
//
                //P[0] = C - 0.5*Lat[0].dx*OrthoSys::e0 - 0.5*Lat[0].dx*OrthoSys::e1;
                //P[1] = C + 0.5*Lat[0].dx*OrthoSys::e0 - 0.5*Lat[0].dx*OrthoSys::e1;
                //P[2] = C + 0.5*Lat[0].dx*OrthoSys::e0 + 0.5*Lat[0].dx*OrthoSys::e1;
                //P[3] = C - 0.5*Lat[0].dx*OrthoSys::e0 + 0.5*Lat[0].dx*OrthoSys::e1;
                //double dmin = 2*Pa->Props.R;
                //double dmax = 0.0;
                //for (size_t j=0;j<P.Size();j++)
                //{
                    //double dist = norm(P[j] - Pa->x);
                    //if (dmin>dist) dmin = dist;
                    //if (dmax<dist) dmax = dist;
                //}
                //if (dmin > Pa->Props.R + Lat[0].dx) continue;
//
                //double len = 0.0;
//
                //if (dmax < Pa->Props.R)
                //{
                    //len = 4.0;
                //}
                //else
                //{
                    //for (size_t j=0;j<4;j++)
                    //{
                        //Vec3_t D = P[(j+1)%4] - P[j];
                        //double a = dot(D,D);
                        //double b = 2*dot(P[j]-Pa->x,D);
                        //double c = dot(P[j]-Pa->x,P[j]-Pa->x) - Pa->Props.R*Pa->Props.R;
                        //if (b*b-4*a*c>0.0)
                        //{
                            //double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
                            //double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
                            //if (ta>1.0&&tb>1.0) continue;
                            //if (ta<0.0&&tb<0.0) continue;
                            //if (ta<0.0) ta = 0.0;
                            //if (tb>1.0) tb = 1.0;
                            //len += norm((tb-ta)*D);
                        //}
                    //}
                //}
//
                //for (size_t j=0;j<Lat.Size();j++)
                //{
                    //cell = Lat[j].GetCell(iVec3_t(n,m,0));
                    //cell->Gamma   = std::max(len/(4.0*Lat[0].dx),cell->Gamma);
                    //if (fabs(cell->Gamma-1.0)<1.0e-12&&fabs(Lat[0].G)>1.0e-12) continue;
                    //Vec3_t B      = C - Pa->x;
                    //Vec3_t VelP   = Pa->x + cross(Pa->x,B);
                    //double rho = cell->Rho;
                    //double Bn  = (cell->Gamma*(cell->Tau-0.5))/((1.0-cell->Gamma)+(cell->Tau-0.5));
                    //for (size_t k=0;k<cell->Nneigh;k++)
                    //{
                        //double Fvpp    = cell->Feq(cell->Op[k],VelP,rho);
                        //double Fvp     = cell->Feq(k          ,VelP,rho);
                        //cell->Omeis[k] = cell->F[cell->Op[k]] - Fvpp - (cell->F[k] - Fvp);
                        //Vec3_t Flbm    = -Bn*cell->Omeis[k]*cell->C[k];
                        //Pa->F          += Flbm;
                        //Pa->T          += cross(B,Flbm);
                    //}
                //}
            //}
        //}
    //}

    //3D imprint
#ifdef USE_OMP
    In = 0;
    Fn = ParCellPairs.Size();
    #pragma omp parallel for schedule(static) num_threads(Nproc)
#endif
    for (size_t i = In;i<Fn;i++)
    {
        DEM::Particle  * Pa   = Particles[ParCellPairs[i].IPar];
        Cell           * cell = Lat[0].Cells[ParCellPairs[i].ICell];
        double x              = Lat[0].dx*(cell->Index(0));
        double y              = Lat[0].dx*(cell->Index(1));
        double z              = Lat[0].dx*(cell->Index(2));
        Vec3_t  C(x,y,z);
        Vec3_t  Xtemp,Xs,Xstemp;
        double len,minl = Pa->Dmax;

        //std::cout << "1" << std::endl;
        if (Pa->IsInsideFaceOnly(C)) len = 12.0*Lat[0].dx;
        else if (ParCellPairs[i].IGeo.Size()==0) continue;
        else
        {
            if (Pa->Faces.Size()>0)
            {
                for (size_t j=0;j<ParCellPairs[i].IGeo.Size();j++)
                {
                    DEM::Distance(C,*Pa->Faces[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                    if (norm(Xtemp-Xstemp) < minl)
                    {
                        minl = norm(Xtemp-Xstemp);
                        Xs   = Xstemp;
                    }
                }
            }
            else if (Pa->Edges.Size()>0)
            {
                for (size_t j=0;j<ParCellPairs[i].IGeo.Size();j++)
                {
                    DEM::Distance(C,*Pa->Edges[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                    if (norm(Xtemp-Xstemp) < minl)
                    {
                        minl = norm(Xtemp-Xstemp);
                        Xs   = Xstemp;
                    }
                }
            }
            else if (Pa->Verts.Size()>0)
            {
                for (size_t j=0;j<ParCellPairs[i].IGeo.Size();j++)
                {
                    DEM::Distance(C,*Pa->Verts[ParCellPairs[i].IGeo[j]],Xtemp,Xstemp);
                    if (norm(Xtemp-Xstemp) < minl)
                    {
                        minl = norm(Xtemp-Xstemp);
                        Xs   = Xstemp;
                    }
                }
            }
            len = DEM::SphereCube(Xs,C,Pa->Props.R,Lat[0].dx);
        }
        
        //std::cout << "2" << std::endl;
        if (fabs(len)<1.0e-12) continue;

        for (size_t j=0;j<Lat.Size();j++)
        {
            double Tau = Lat[j].Tau;
            cell = Lat[j].Cells[ParCellPairs[i].ICell];
            cell->Gamma   = std::max(len/(12.0*Lat[0].dx),cell->Gamma);
            //if (fabs(cell->Gamma-1.0)<1.0e-12)
            if (fabs(cell->Gamma-1.0)<1.0e-12&&fabs(Lat[0].G)>1.0e-12) 
            {
                continue;
            }
            Vec3_t B      = C - Pa->x;
            Vec3_t tmp;
            Rotation(Pa->w,Pa->Q,tmp);
            Vec3_t VelP   = Pa->v + cross(tmp,B);
            double rho = cell->Rho;
            double Bn  = (cell->Gamma*(Tau-0.5))/((1.0-cell->Gamma)+(Tau-0.5));
            for (size_t k=0;k<cell->Nneigh;k++)
            {
                double Fvpp    = cell->Feq(cell->Op[k],VelP,rho);
                double Fvp     = cell->Feq(k          ,VelP,rho);
                cell->Omeis[k] = cell->F[cell->Op[k]] - Fvpp - (cell->F[k] - Fvp);
                Vec3_t Flbm    = -Fconv*Bn*cell->Omeis[k]*cell->C[k]*cell->Cs*cell->Cs*Lat[0].dx*Lat[0].dx;
                Vec3_t T,Tt;
                Tt =           cross(B,Flbm);
                Quaternion_t q;
                Conjugate    (Pa->Q,q);
                Rotation     (Tt,q,T);
                //std::cout << "1" << std::endl;
#ifdef USE_THREAD
                pthread_mutex_lock(&Pa->lck);
#elif USE_OMP
                omp_set_lock      (&Pa->lck);
#endif
                Pa->F          += Flbm;
                Pa->T          += T;
#ifdef USE_THREAD
                pthread_mutex_unlock(&Pa->lck);
#elif USE_OMP
                omp_unset_lock    (&Pa->lck);
#endif
                //std::cout << "2" << std::endl;
//#ifdef USE_THREAD
                //pthread_mutex_lock  (&cell->lck);
//#endif
                //cell->Omeis[k] = Omeis;
                //cell->Gamma    = Gamma;
//#ifdef USE_THREAD
                //pthread_mutex_unlock(&cell->lck);
//#endif
                //std::cout << "3" << std::endl;
            }
        }
        //std::cout << "3" << std::endl;
    }
}

inline void Domain::ResetDisplacements()
{
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].Dmx = 0.0;
        MTD[i].LLC.Resize(0);
    }
    //std::cout << "1" << std::endl;
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Particles.Size();i++)
    {
        Particles[i]->ResetDisplacements();
        if(Particles[i]->IsFree())
        {
            iVec3_t idx = (Particles[i]->x - LCxmin)/(2.0*Beta*MaxDmax);
            MTD[omp_get_thread_num()].LLC.Push(std::make_pair(idx,i));
        }
    }
    //std::cout << "2" << std::endl;
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LLC.Size();j++)
        {
            //std::cout << i << " " << j << std::endl;
            size_t idx = DEM::Pt2idx(MTD[i].LLC[j].first,LCellDim);
            LinkedCell[idx].Push(MTD[i].LLC[j].second);
        }
    }
    //std::cout << "3" << std::endl;
#else
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->ResetDisplacements();
    }
#endif
}

inline double Domain::MaxDisplacement()
{
    double md = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        double mdp = Particles[i]->MaxDisplacement();
        if (mdp>md) md = mdp;
    }
    return md;
}

#ifdef USE_OMP
inline void Domain::UpdateLinkedCells()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].LPP.Resize(0);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<FreePar.Size()  ;i++)
    for (size_t j=0;j<NoFreePar.Size();j++)
    {
        size_t i1 = std::min(FreePar[i],NoFreePar[j]);
        size_t i2 = std::max(FreePar[i],NoFreePar[j]);
        MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t idx=0;idx<LinkedCell.Size();idx++)
    {
        if (LinkedCell[idx].Size()==0) continue;
        iVec3_t index;
        DEM::idx2Pt(idx,index,LCellDim);
        for (size_t n=0  ;n<LinkedCell[idx].Size()-1;n++)
        for (size_t m=n+1;m<LinkedCell[idx].Size()  ;m++)
        {
            size_t i1 = LinkedCell[idx][n];
            size_t i2 = LinkedCell[idx][m];
            if (i1==i2) continue;
            MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
        }
        size_t i = index(0);
        size_t j = index(1);
        size_t k = index(2);
        for (size_t knb=std::max(0,int(k)-1);knb<=std::min(LCellDim(2)-1,k+1);knb++)
        for (size_t jnb=std::max(0,int(j)-1);jnb<=std::min(LCellDim(1)-1,j+1);jnb++)
        for (size_t inb=std::max(0,int(i)-1);inb<=std::min(LCellDim(0)-1,i+1);inb++)
        {
            iVec3_t Ptnb(inb,jnb,knb);
            size_t idxnb = DEM::Pt2idx(Ptnb,LCellDim);
            if (idxnb>idx)
            {
                for (size_t n=0;n<LinkedCell[idx].Size()  ;n++)
                {
                    for (size_t m=0;m<LinkedCell[idxnb].Size()  ;m++)
                    {
                        size_t i1 = std::min(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        size_t i2 = std::max(LinkedCell[idx  ][n],LinkedCell[idxnb][m]);
                        if (i1==i2) continue;
                        MTD[omp_get_thread_num()].LPP.Push(std::make_pair(i1,i2));
                    }
                }
            }
        }
    }
    size_t Npp = 0;
    for (size_t i=0;i<Nproc;i++)
    {
        Npp += MTD[i].LPP.Size();
    }
    ListPosPairs.Resize(Npp);
    size_t idx = 0;
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LPP.Size();j++)
        {
            ListPosPairs[idx] = MTD[i].LPP[j];
            idx++;
        }
    }
}
#endif

inline void Domain::ResetContacts()
{
#ifdef USE_OMP
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].LC.Resize(0);
        MTD[i].LCI.Resize(0);
        MTD[i].LCB.Resize(0);
        MTD[i].LPC.Resize(0);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<ListPosPairs.Size();n++)
    {
        size_t i = ListPosPairs[n].first;
        size_t j = ListPosPairs[n].second;
        bool pi_has_vf = !Particles[i]->IsFree();
        bool pj_has_vf = !Particles[j]->IsFree();
        bool close = (DEM::Distance(Particles[i]->x,Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
        if ((pi_has_vf && pj_has_vf) || !close) continue;
        std::set<std::pair<DEM::Particle *, DEM::Particle *> >::iterator it = Listofpairs.find(std::make_pair(Particles[i],Particles[j]));
        if (it != Listofpairs.end())
        {
            continue;
        }
        MTD[omp_get_thread_num()].LC.Push(std::make_pair(i,j));
    }
    for (size_t i=0;i<Nproc;i++)
    {
        //std::cout << MTD[i].LC.Size() << std::endl;
        for (size_t j=0;j<MTD[i].LC.Size();j++)
        {
        //std::cout << MTD[i].LC.Size() << std::endl;
            size_t n = MTD[i].LC[j].first;
            size_t m = MTD[i].LC[j].second;
            Listofpairs.insert(std::make_pair(Particles[n],Particles[m]));
            if (Particles[n]->Verts.Size()==1 && Particles[m]->Verts.Size()==1)
            {
                CInteractons.Push (new DEM::CInteractonSphere(Particles[n],Particles[m]));
            }
            else
            {
                CInteractons.Push (new DEM::CInteracton(Particles[n],Particles[m]));
            }
        }
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<CInteractons.Size();n++)
    {
        if(CInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LCI.Push(n);
    }
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t n=0;n<BInteractons.Size();n++)
    {
        if(BInteractons[n]->UpdateContacts(Alpha)) MTD[omp_get_thread_num()].LCB.Push(n);
    }
    Interactons.Resize(0);
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LCI.Size();j++)
        {
            Interactons.Push(CInteractons[MTD[i].LCI[j]]);
        }
        for (size_t j=0;j<MTD[i].LCB.Size();j++)
        {
            Interactons.Push(BInteractons[MTD[i].LCB[j]]);
        }
    }
    if (Lat[0].Ndim(2)==1)
    {
        // TODO 2D case
    }
    else
    {
    #pragma omp parallel for schedule(static) num_threads(Nproc)
	    for (size_t i=0;i<Particles.Size();i++)
        {
            //Cell  * cell = dat.Dom->Lat[0].Cells[i];
            DEM::Particle * Pa = Particles[i];
            //std::cout << std::max(0.0,double(Pa->x(0)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " 
                      //<< std::max(0.0,double(Pa->x(1)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " 
                      //<< std::max(0.0,double(Pa->x(2)-Pa->Dmax-2.0*dat.Dom->Alpha-dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " << std::endl;
            //std::cout << std::min(double(dat.Dom->Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " "
                      //<< std::min(double(dat.Dom->Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " "
                      //<< std::min(double(dat.Dom->Lat[0].Ndim(2)-1),double(Pa->x(2)+Pa->Dmax+2.0*dat.Dom->Alpha+dat.Dom->Lat[0].dx)/dat.Dom->Lat[0].dx) << " " << std::endl;

            for (size_t n=std::max(0.0,double(Pa->x(0)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);n<=std::min(double(Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);n++)
            for (size_t m=std::max(0.0,double(Pa->x(1)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);m<=std::min(double(Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);m++)
            for (size_t l=std::max(0.0,double(Pa->x(2)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);l<=std::min(double(Lat[0].Ndim(2)-1),double(Pa->x(2)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);l++)
            //for (size_t n=0;n<dat.Dom->Particles.Size();n++)
            {
                //DEM::Particle * Pa = dat.Dom->Particles[n];
                Cell  * cell = Lat[0].GetCell(iVec3_t(n,m,l));
                double x     = Lat[0].dx*(cell->Index(0));
                double y     = Lat[0].dx*(cell->Index(1));
                double z     = Lat[0].dx*(cell->Index(2));
                Vec3_t  C(x,y,z);
                if ((norm(C-Pa->x)>2.0*Alpha+2.0*Lat[0].dx+Pa->Dmax)||(Pa->Bdry==true)) continue;
                ParticleCellPair NewPCP;
                NewPCP.IPar = i;
                //NewPCP.IPar = n;
                NewPCP.ICell= cell->ID;
                bool valid = false;
                if (Pa->Faces.Size()>0)
                {
                    for (size_t j=0;j<Pa->Faces.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Faces[j])<2.0*Alpha+Pa->Props.R)
                        {
                            if (Pa->Faces[j]->Area()<2.0*M_PI*Pa->Props.R*Pa->Props.R)
                            {
                                //std::cout << "it did" << std::endl;
                                continue;
                            }
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                    }
                }
                else if (Pa->Edges.Size()>0)
                {
                    for (size_t j=0;j<Pa->Edges.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Edges[j])<2.0*Alpha+Pa->Props.R) 
                        {
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                           
                    }
                }
                else if (Pa->Verts.Size()>0)
                {
                    for (size_t j=0;j<Pa->Verts.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Verts[j])<2.0*Alpha+Pa->Props.R)
                        {
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                    }
                }
                if (Pa->IsInsideFaceOnly(C)) valid = true;

                if (valid) MTD[omp_get_thread_num()].LPC.Push(NewPCP);
            }
        }
    }
    ParCellPairs.Resize(0);
    for (size_t i=0;i<Nproc;i++)
    {
        for (size_t j=0;j<MTD[i].LPC.Size();j++)
        {
            ParCellPairs.Push(MTD[i].LPC[j]);
        }
    }
#else
	if (Particles.Size()==0) return;
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            bool pj_has_vf = !Particles[j]->IsFree();

            bool close = (norm(Particles[i]->x-Particles[j]->x)<=Particles[i]->Dmax+Particles[j]->Dmax+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close) continue;
            
            // checking if the interacton exist for that pair of particles
            set<std::pair<DEM::Particle *, DEM::Particle *> >::iterator it = Listofpairs.find(std::make_pair(Particles[i],Particles[j]));
            if (it != Listofpairs.end())
            {
                continue;
            }
            Listofpairs.insert(std::make_pair(Particles[i],Particles[j]));
            if (Particles[i]->Verts.Size()==1 && Particles[j]->Verts.Size()==1)
            {
                CInteractons.Push (new DEM::CInteractonSphere(Particles[i],Particles[j]));
            }
            else
            {
                CInteractons.Push (new DEM::CInteracton(Particles[i],Particles[j]));
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

    ParCellPairs.Resize(0);

    if (Lat[0].Ndim(2)==1)
    {
        // TODO 2D case
    }
    else
    {
        for (size_t i=0; i<Particles.Size(); i++)
        {
            DEM::Particle * Pa = Particles[i];
            for (size_t n=std::max(0.0,double(Pa->x(0)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);n<=std::min(double(Lat[0].Ndim(0)-1),double(Pa->x(0)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);n++)
            for (size_t m=std::max(0.0,double(Pa->x(1)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);m<=std::min(double(Lat[0].Ndim(1)-1),double(Pa->x(1)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);m++)
            for (size_t l=std::max(0.0,double(Pa->x(2)-Pa->Dmax-2.0*Alpha-Lat[0].dx)/Lat[0].dx);l<=std::min(double(Lat[0].Ndim(2)-1),double(Pa->x(2)+Pa->Dmax+2.0*Alpha+Lat[0].dx)/Lat[0].dx);l++)
            {
                Cell  * cell = Lat[0].GetCell(iVec3_t(n,m,l));
                double x     = Lat[0].dx*(cell->Index(0));
                double y     = Lat[0].dx*(cell->Index(1));
                double z     = Lat[0].dx*(cell->Index(2));
                Vec3_t  C(x,y,z);
                ParticleCellPair NewPCP;
                NewPCP.IPar = i;
                NewPCP.ICell= cell->ID;
                bool valid = false;
                if (Pa->Faces.Size()>0)
                {
                    for (size_t j=0;j<Pa->Faces.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Faces[j])<2.0*Alpha+Pa->Props.R)                        
                        {
                            if (Pa->Faces[j]->Area()<2.0*M_PI*Pa->Props.R*Pa->Props.R)
                            {
                                //std::cout << "it did" << std::endl;
                                continue;
                            }
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                    }
                }
                else if (Pa->Edges.Size()>0)
                {
                    for (size_t j=0;j<Pa->Edges.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Edges[j])<2.0*Alpha+Pa->Props.R)
                        {
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                    }
                }
                else if (Pa->Verts.Size()>0)
                {
                    for (size_t j=0;j<Pa->Verts.Size();j++)
                    {
                        if (DEM::Distance(C,*Pa->Verts[j])<2.0*Alpha+Pa->Props.R)
                        {
                            NewPCP.IGeo.Push(j);
                            valid = true;
                        }
                    }
                }
                if (Pa->IsInsideFaceOnly(C)) valid = true;
                if (valid) ParCellPairs.Push(NewPCP);
            }
        }
    }
#endif
}

//Utility methods
inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    if (Particles.Size()==0)
    {
        minX = OrthoSys::O;
        maxX = OrthoSys::O;
    }
    else
    {
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
}

inline void Domain::Center(Vec3_t C)
{
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    Vec3_t Transport(-0.5*(maxX+minX));
    Transport += C;
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Translate(Transport);
}

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

//DEM particle methods
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
    Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

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

inline void Domain::GenSpheresBox (int Tag, Vec3_t const & X0, Vec3_t const & X1, double R, double rho, size_t Randomseed, double fraction, double RminFraction)
{
    // find radius from the edge's length
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating packing of spheres -----------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    srand(Randomseed);

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
                if (rand()<fraction*RAND_MAX) AddSphere (Tag,X,R*RminFraction+(1.0*rand())/RAND_MAX*(R-R*RminFraction),rho);
            }
        }
    }
    
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
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
    Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

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
    Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

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

inline void Domain::AddOcta (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
{
    // vertices
    Array<Vec3_t> V(6);
    double l = L/sqrt(2.0);
    V[0] =  l,0.0,0.0;
    V[1] = -l,0.0,0.0;
    V[2] =0.0,  l,0.0;
    V[3] =0.0, -l,0.0;
    V[4] =0.0,0.0,  l;
    V[5] =0.0,0.0, -l;

    // edges
    Array<Array <int> > E(12);
    for (size_t i=0; i<12; ++i) E[i].Resize(2);
    E[ 0] = 0, 2;
    E[ 1] = 0, 3;
    E[ 2] = 0, 4;
    E[ 3] = 0, 5;
    E[ 4] = 1, 2;
    E[ 5] = 1, 3;
    E[ 6] = 1, 4;
    E[ 7] = 1, 5;
    E[ 8] = 2, 4;
    E[ 9] = 2, 5;
    E[10] = 3, 4;
    E[11] = 3, 5;

    // faces
    Array<Array <int> > F(8);
    for (size_t i=0; i<8; i++) F[i].Resize(3);
    F[0] = 0, 2, 4;
    F[1] = 0, 4, 3;
    F[2] = 0, 3, 5;
    F[3] = 0, 5, 2;
    F[4] = 1, 2, 5;
    F[5] = 1, 5, 3;
    F[6] = 1, 3, 4;
    F[7] = 1, 4, 2;

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
    Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));

    // clean up
    if (!ThereisanAxis) delete Axis;
    q(0) = 1.0;
    q(1) = 0.0;
    q(2) = 0.0;
    q(3) = 0.0;
    q = q/norm(q);

    Particles[Particles.Size()-1]->Q          = q;
    Particles[Particles.Size()-1]->Props.V    = sqrt(2.0)/3.0*L*L*L;
    Particles[Particles.Size()-1]->Props.m    = rho*sqrt(2.0)/3.0*L*L*L;
    Particles[Particles.Size()-1]->I          = L*L, L*L, L*L;
    Particles[Particles.Size()-1]->I         *= Particles[Particles.Size()-1]->Props.m/10.0;
    Particles[Particles.Size()-1]->x          = X;
    Particles[Particles.Size()-1]->Ekin       = 0.0;
    Particles[Particles.Size()-1]->Erot       = 0.0;
    Particles[Particles.Size()-1]->Dmax       = L/sqrt(2.0)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index      = Particles.Size()-1;

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
    Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
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

inline void Domain::GenBox (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, double rho, bool Cohesion)
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
    AddPlane (InitialTag,   Vec3_t(Lx/2.0,0.0,0.0),  R, Cf*Lz, Cf*Ly, rho, M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-1, Vec3_t(-Lx/2.0,0.0,0.0), R, Cf*Lz, Cf*Ly, rho, 3.0*M_PI/2.0, &axis1);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-2, Vec3_t(0.0,Ly/2.0,0.0),  R, Cf*Lx, Cf*Lz, rho, 3.0*M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-3, Vec3_t(0.0,-Ly/2.0,0.0), R, Cf*Lx, Cf*Lz, rho, M_PI/2.0, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-4, Vec3_t(0.0,0.0,Lz/2.0),  R, Cf*Lx, Cf*Ly, rho);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);
    AddPlane (InitialTag-5, Vec3_t(0.0,0.0,-Lz/2.0), R, Cf*Lx, Cf*Ly, rho, M_PI, &axis0);
    Particles[Particles.Size()-1]->Initialize(Particles.Size()-1);

    // define some tolerance for comparissions
    if (Cohesion)
    {
        double tol1 = 1.0e-8;
        double tol2 = 1.0e-3;
        for (size_t i=0;i<IIndex;i++)
        {
            DEM::Particle * P1 = Particles[i];
            for (size_t j=IIndex;j<Particles.Size();j++)
            {
                DEM::Particle * P2 = Particles[j];
                for (size_t k=0;k<P1->Faces.Size();k++)
                {
                    DEM::Face * F1 = P1->Faces[k];
                    Vec3_t n1,c1;
                    F1->Normal  (n1);
                    F1->Centroid(c1);
                    DEM::Face * F2 = P2->Faces[0];
                    Vec3_t n2,c2;
                    F2->Normal  (n2);
                    F2->Centroid(c2);
                    Vec3_t n = 0.5*(n1-n2);
                    n/=norm(n);
                    if ((fabs(dot(n1,n2)+1.0)<tol1)
                       &&(fabs(dot(c2-c1,n)-2*R)<tol2))
                    {
                        BInteractons.Push(new DEM::BInteracton(P1,P2,k,1));
                        break;
                    }
                }
            }        
        }
    }
}

inline void Domain::GenBoundingBox (int InitialTag, double R, double Cf, double rho, bool Cohesion)
{
    Center();
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    GenBox(InitialTag, maxX(0)-minX(0)+2*R, maxX(1)-minX(1)+2*R, maxX(2)-minX(2)+2*R, R, Cf, rho, Cohesion);
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
        DEM::PolyhedraMP(V,F,vol,CM,It);
        DEM::Erosion(V,E,F,R);

        // add particle
        Particles.Push (new DEM::Particle(M.Cells[i]->Tag, V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
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
            double Dmax = DEM::Distance(CM,V[0])+R;
            for (size_t i=1; i<V.Size(); ++i)
            {
                if (DEM::Distance(CM,V[i])+R > Dmax) Dmax = DEM::Distance(CM,V[i])+R;
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
                if ((size_t)Neigh[i][j]>i) BInteractons.Push(new DEM::BInteracton(Particles[i+IIndex],Particles[Neigh[i][j]+IIndex],FNeigh[i][j],FNeigh[Neigh[i][j]][index]));
            }
        }
        
    }

    // info
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::AddVoroCell (int Tag, voro::voronoicell & VC, double R, double rho, bool Erode, Vec3_t nv)
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
    //VC.reset_edges();
    double vol; // volume of the polyhedron
    Vec3_t CM;  // Center of mass of the polyhedron
    Mat3_t It;  // Inertia tensor of the polyhedron
    DEM::PolyhedraMP(V,F,vol,CM,It);
    if (Erode) DEM::Erosion(V,E,F,R);
    // add particle
    Particles.Push (new DEM::Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
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
    double Dmax = DEM::Distance(CM,V[0])+R;
    for (size_t i=1; i<V.Size(); ++i)
    {
        if (DEM::Distance(CM,V[i])+R > Dmax) Dmax = DEM::Distance(CM,V[i])+R;
    }
    Particles[Particles.Size()-1]->Ekin = 0.0;
    Particles[Particles.Size()-1]->Erot = 0.0;
    Particles[Particles.Size()-1]->Dmax  = Dmax;
    Particles[Particles.Size()-1]->PropsReady = true;
    Particles[Particles.Size()-1]->Index = Particles.Size()-1;
}

inline void Domain::AddVoroPack (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho
                                 , bool Cohesion, bVec3_t Periodic,size_t Randomseed, double fraction, Vec3_t qin)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Adding Voronoi particles packing --------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    srand(Randomseed);
    const double x_min=-(nx/2.0), x_max=nx/2.0;
    const double y_min=-(ny/2.0), y_max=ny/2.0;
    const double z_min=-(nz/2.0), z_max=nz/2.0;
    voro::container con(x_min,x_max,y_min,y_max,z_min,z_max,nx,ny,nz, Periodic(0),Periodic(1),Periodic(2),8);
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

    Array<Array <size_t> > ListBpairs(n);
    size_t IIndex = Particles.Size();
    voro::c_loop_all vl(con);
    voro::voronoicell c;
    if(vl.start()) do if(con.compute_cell(c,vl)) 
    {
        {
            if (rand()<fraction*RAND_MAX)
            {
                AddVoroCell(Tag,c,R,rho,true,Vec3_t(Lx/nx,Ly/ny,Lz/nz));
                double *pp = con.p[vl.ijk]+3*vl.q;
                Vec3_t trans(Lx*pp[0]/nx,Ly*pp[1]/ny,Lz*pp[2]/nz);
                DEM::Particle * P = Particles[Particles.Size()-1];
                P->Translate(trans);
            }
        }
    } while(vl.inc());


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
            DEM::Particle * P1 = Particles[i];
            for (size_t j=i+1;j<Particles.Size();j++)
            {
                DEM::Particle * P2 = Particles[j];
                if (DEM::Distance(P1->x,P2->x)<P1->Dmax+P2->Dmax)
                {
                    for (size_t k=0;k<P1->Faces.Size();k++)
                    {
                        DEM::Face * F1 = P1->Faces[k];
                        Vec3_t n1,c1;
                        F1->Normal  (n1);
                        F1->Centroid(c1);
                        bool found = false;
                        for (size_t l=0;l<P2->Faces.Size();l++)
                        {
                            DEM::Face * F2 = P2->Faces[l];
                            Vec3_t n2,c2;
                            F2->Normal  (n2);
                            F2->Centroid(c2);
                            Vec3_t n = 0.5*(n1-n2);
                            n/=norm(n);
                            if ((fabs(dot(n1,n2)+1.0)<tol1)
                               &&(fabs(DEM::Distance(c1,*F2)-2*R)<tol2)
                               &&(fabs(DEM::Distance(c2,*F1)-2*R)<tol2))
                            {
                                BInteractons.Push(new DEM::BInteracton(P1,P2,k,l));
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

inline void Domain::AddVoroPack (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho
                                 , bool Cohesion, bool Periodic,size_t Randomseed, double fraction, Vec3_t qin)
{
    AddVoroPack(Tag,R,Lx,Ly,Lz,nx,ny,nz,rho,Cohesion,bVec3_t(Periodic,Periodic,Periodic),Randomseed,fraction,qin);
}

//Particle access methods
inline DEM::Particle * Domain::GetParticle (int Tag, bool Check)
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

inline DEM::Particle const & Domain::GetParticle (int Tag, bool Check) const
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
}

//Dynamic methods
inline void Domain::Initialize (double dt)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing particles ------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    // set flag
    if (!Initialized)
    {
        Initialized = true;
        // initialize all particles
        Voltot = 0.0;
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->Initialize(i);
            Particles[i]->InitializeVelocity(dt);
            Voltot+=Particles[i]->Props.V;
        }
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

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{

    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Finished = false;
    Nproc = TheNproc;

    // initialize particles
    Initialize (dt);
    // calc the total volume of particles (solids)
    FreePar.Resize(0);
    NoFreePar.Resize(0);
    Vs = 0.0;
    Ms = 0.0;
    MaxDmax        =  0.0;
    double MaxKn   = -1.0;
    double MaxBn   = -1.0;
    double MinDmax = -1.0;
    double MinMass = -1.0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        if (Particles[i]->IsFree())
        {
            Vs += Particles[i]->Props.V;
            Ms += Particles[i]->Props.m;
            if (Particles[i]->Dmax     > MaxDmax) MaxDmax = Particles[i]->Dmax;
            if (Particles[i]->Props.Kn > MaxKn  ||(MinDmax<0.0)) MaxKn   = Particles[i]->Props.Kn;
            if (Particles[i]->Dmax     < MinDmax||(MinDmax<0.0)) MinDmax = Particles[i]->Dmax;
            if (Particles[i]->Props.m  < MinMass||(MinMass<0.0)) MinMass = Particles[i]->Props.m;
            FreePar.Push(i);
        }
        else NoFreePar.Push(i);
    }
    for (size_t i=0; i<BInteractons.Size(); i++)
    {
        double pbn = BInteractons[i]->Bn/BInteractons[i]->L0;
        if (pbn > MaxBn||(MinDmax<0.0)) MaxBn = pbn;
    }

    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    printf("%s  Porosity                         =  %g%s\n"       ,TERM_CLR4, 1.0 - Lat[0].SolidFraction()         , TERM_RST);
    printf("%s  Total mass   of free particles   =  %g%s\n"       ,TERM_CLR4, Ms                                   , TERM_RST);
    printf("%s  Total volume of free particles   =  %g%s\n"       ,TERM_CLR4, Vs                                   , TERM_RST);
    printf("%s  Total number of particles        =  %zd%s\n"      ,TERM_CLR2, Particles.Size()                     , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                   , TERM_RST);
    printf("%s  Verlet distance                  =  %g%s\n"       ,TERM_CLR2, Alpha                                , TERM_RST);
    for (size_t i=0;i<Lat.Size();i++)
    {
    printf("%s  Tau of Lattice %zd                 =  %g%s\n"       ,TERM_CLR2, i, Lat[i].Tau                        , TERM_RST);
    }
    printf("%s  Suggested Time Step              =  %g%s\n"       ,TERM_CLR5, 0.1*sqrt(MinMass/(MaxKn+MaxBn))      , TERM_RST);
    printf("%s  Suggested Verlet distance        =  %g or %g%s\n" ,TERM_CLR5, 0.5*MinDmax, 0.25*(MinDmax + MaxDmax), TERM_RST);

    if (Alpha > MinDmax&&MinDmax>0.0)
    {
        Alpha = MinDmax;
        printf("%s  Verlet distance changed to       =  %g%s\n"   ,TERM_CLR2, Alpha                                    , TERM_RST);
    }

    if (Alpha < 0.0) throw new Fatal("Verlet distance cannot be negative");



    //std::cout << "1" << std::endl;
    // Creates pair of cells to speed up body force calculation
    for (size_t i=0;i<Lat[0].Ncells;i++)
    {
        Cell * c = Lat[0].Cells[i];
        for (size_t j=1;j<c->Nneigh;j++)
        {
            Cell * nb = Lat[0].Cells[c->Neighs[j]];
            if (nb->ID>c->ID) 
            {
                if (!c->IsSolid||!nb->IsSolid) CellPairs.Push(iVec3_t(i,nb->ID,j));
            }
        }
    }
    //std::cout << "2" << std::endl;
    
    MTD = new LBM::MtData[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = this;
        MTD[i].Dmx      = 0.0;
        MTD[i].dt       = Lat[0].dt;
    }
    //std::cout << "3" << std::endl;
    //std::cout << "4" << std::endl;
#ifdef USE_THREAD
    pthread_t thrs[Nproc];   
    
    if (Particles.Size() > 0)
    {
        for (size_t i=0; i<Particles.Size()-1; i++)
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            ListPosPairs.Push(make_pair(i,j));
        }
    }

    //GlobalResetDisplacement
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_create(&thrs[i], NULL, GlobalResetDisplacement, &MTD[i]);
    }
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_join(thrs[i], NULL);
    }
    ParCellPairs.Resize(0);
    //GlobalResetContacts1
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
                CInteractons.Push (new DEM::CInteractonSphere(Particles[n],Particles[m]));
            }
            else
            {
                CInteractons.Push (new DEM::CInteracton(Particles[n],Particles[m]));
            }
        }

        for (size_t j=0;j<MTD[i].LPC.Size();j++)
        {
            ParCellPairs.Push(MTD[i].LPC[j]);
        }
        MTD[i].LPC.Resize(0);
    }
    Interactons.Resize(0);
    //GlobalResetContacts2
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_create(&thrs[i], NULL, GlobalResetContacts2, &MTD[i]);
    }
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
    //std::cout << "1" << std::endl;
    //GlobalImprint
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_create(&thrs[i], NULL, GlobalImprint, &MTD[i]);
    }
    for (size_t i=0;i<Nproc;i++)
    {
        pthread_join(thrs[i], NULL);
    }
    //std::cout << "2" << std::endl;
#elif USE_OMP
    LinkedCell.Resize(0);
    BoundingBox(LCxmin,LCxmax);
    LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax) + iVec3_t(1,1,1);
    LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));
    //std::cout << "1" << std::endl;
    ResetDisplacements();

    //std::cout << "2" << std::endl;
    UpdateLinkedCells();

    //std::cout << "3" << std::endl;
    ResetContacts();

    //std::cout << "4" << std::endl;
    ImprintLattice(1,Nproc);    
#else

    //Connect particles and lattice
    ImprintLattice();

    // build the map of possible contacts (for the Halo)
    ResetContacts();
    
    // set the displacement of the particles to zero (for the Halo)
    ResetDisplacements();

#endif
    double tout = Time;
    while (Time < Tf)
    {
        //std::cout << Interactons.Size() << " " << CInteractons.Size() << " " << BInteractons.Size() << " " << ParCellPairs.Size() << " " << Particles.Size() << std::endl;
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time >= tout)
        {
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    WriteXDMF(fn.CStr());
                    #else
                    //WriteVTK (fn.CStr());
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            if (BInteractons.Size()>0) Clusters();
            tout += dtOut;
            idx_out++;
        }


#ifdef USE_THREAD
        //GlobalIni
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalIni, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << "1" <<std::endl;
        //GlobalImprint
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalImprint, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << "2" <<std::endl;
        //GlobalForce
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalForce, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        double maxdis = 0.0;
        //std::cout << "3" <<std::endl;
        //GlobalMove
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalMove, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
            if (maxdis<MTD[i].Dmx) maxdis = MTD[i].Dmx;
        }
        //GlobalResetDisplacement
        if (maxdis>Alpha)
        {
            //std::cout << "Remaking cells " << maxdis << std::endl;
            //GlobalResetDisplacement
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalResetDisplacement, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
            }
            ParCellPairs.Resize(0);
            //GlobalResetContacts1
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
                        CInteractons.Push (new DEM::CInteractonSphere(Particles[n],Particles[m]));
                    }
                    else
                    {
                        CInteractons.Push (new DEM::CInteracton(Particles[n],Particles[m]));
                    }
                }

                for (size_t j=0;j<MTD[i].LPC.Size();j++)
                {
                    ParCellPairs.Push(MTD[i].LPC[j]);
                }
                MTD[i].LPC.Resize(0);
            }
            //GlobalResetContacts2
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
        }
        //std::cout << "4" <<std::endl;
        if (Lat.Size()>1||fabs(Lat[0].G)>1.0e-12)
        {
            //GlobalApplyForce
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalApplyForce, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
            }
        }
        //GlobalCollide
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalCollide, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << "5" <<std::endl;
        //GlobalBounceBack
        //for (size_t i=0;i<Nproc;i++)
        //{
            //pthread_create(&thrs[i], NULL, GlobalBounceBack, &MTD[i]);
        //}
        //for (size_t i=0;i<Nproc;i++)
        //{
            //pthread_join(thrs[i], NULL);
        //}
        //GlobalStream1
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalStream1, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << "6" <<std::endl;
        //GlobalStream2
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalStream2, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << "5" <<std::endl;
#elif USE_OMP 
        //Initialize all the particles and cells
        for (size_t i=0;i<Lat.Size();i++)
        {
            Lat[i].SetZeroGamma(0,Nproc);
        }
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for(size_t i=0;i<Particles.Size();i++)
        {
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Tf;
        }
        
        //Imprint the particles into the lattice
        ImprintLattice(0,Nproc);

        //Calculate interparticle forces
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<Interactons.Size(); i++)
        {
		    if (Interactons[i]->CalcForce(dt))
            {
                WriteXDMF("error");
                std::cout << "Maximun overlap detected between particles at time " << Time << std::endl;
                sleep(1);
                throw new Fatal("Maximun overlap detected between particles");
            }
            omp_set_lock  (&Interactons[i]->P1->lck);
            Interactons[i]->P1->F += Interactons[i]->F1;
            Interactons[i]->P1->T += Interactons[i]->T1;
            omp_unset_lock(&Interactons[i]->P1->lck);
            omp_set_lock  (&Interactons[i]->P2->lck);
            Interactons[i]->P2->F += Interactons[i]->F2;
            Interactons[i]->P2->T += Interactons[i]->T2;
            omp_unset_lock(&Interactons[i]->P2->lck);
        }

        //Checking if particles have moved beyond the verlet distance
        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0;i<Nproc;i++)
        {
            MTD[i].Dmx = 0.0;
        }

        #pragma omp parallel for schedule(static) num_threads(Nproc)
        for (size_t i=0; i<Particles.Size(); i++)
        {
            //std::cout << "1" << std::endl;
		    Particles[i]->Translate(dt);
            //std::cout << "2" << std::endl;
		    Particles[i]->Rotate(dt);
            //std::cout << "3" << std::endl;
            if (Particles[i]->MaxDisplacement()>MTD[omp_get_thread_num()].Dmx) MTD[omp_get_thread_num()].Dmx = Particles[i]->MaxDisplacement();
        }

        double maxdis = 0.0;
        for (size_t i=0;i<Nproc;i++)
        {
            if (maxdis<MTD[i].Dmx) maxdis = MTD[i].Dmx;
        }
        
        if (maxdis>Alpha)
        {
            LinkedCell.Resize(0);
            BoundingBox(LCxmin,LCxmax);
            LCellDim = (LCxmax - LCxmin)/(2.0*Beta*MaxDmax) + iVec3_t(1,1,1);
            LinkedCell.Resize(LCellDim(0)*LCellDim(1)*LCellDim(2));

            ResetDisplacements();

            UpdateLinkedCells();

            ResetContacts();

        }

        //Apply molecular forces
        if (Lat.Size()>1||fabs(Lat[0].G)>1.0e-12)
        {
            bool MC = false;
            if (Lat.Size()==2)
            {
                if (fabs(Lat[0].G)<1.0e-9&&fabs(Lat[1].G)<1.0e-9) MC = true;
            }
            ApplyForce(0,Nproc,MC);
        }

        //Apply collision operator
        Collide(0,Nproc);

        //Stream the distribution functions
        for (size_t i=0;i<Lat.Size();i++)
        {
            Lat[i].Stream(0,Nproc);
        }
        
#else
        //Assigning a vlaue of zero to the particles forces and torques
        for(size_t i=0;i<Particles.Size();i++)
        {
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Tf;
        }

        //Set Gamma values of the lattice cell to zero
        for(size_t j=0;j<Lat.Size();j++)
        {
            Lat[j].SetZeroGamma();
        }
        //Connect particles and lattice
        ImprintLattice();

        //Move Particles
        for(size_t i=0;i<Interactons.Size();i++) Interactons[i]->CalcForce(dt);
        for(size_t i=0;i<Particles.Size()  ;i++) 
        {
            Particles[i]->Translate(dt);
            Particles[i]->Rotate(dt);
        }

        //Move fluid
        if (Lat.Size()>1||fabs(Lat[0].G)>1.0e-12) ApplyForce();
        Collide();
        for(size_t j=0;j<Lat.Size();j++)
        {
            Lat[j].BounceBack();
            Lat[j].Stream();
        }

        if (MaxDisplacement()>Alpha)
        {
            ResetContacts();
            ResetDisplacements();
        }
#endif

        Time += dt;
    }
    // last output
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}
}


#endif

