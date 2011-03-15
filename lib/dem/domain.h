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

#ifndef MECHSYS_DEM_DOMAIN_H
#define MECHSYS_DEM_DOMAIN_H

// Std lib
#include <cmath>
#include <stdlib.h> // for M_PI
#include <iostream>
#include <fstream>
#include <set>

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// Voro++
#include "src/voro++.cc"

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
    size_t Randomseed, double fraction, double RminFraction = 1.0);                                                              ///< General spheres
    void GenRice         (int Tag, double L, size_t N, double R, double rho, size_t Randomseed, double fraction);                ///< General rices
    void GenBox          (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf, bool Cohesion=false);            ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenOpenBox      (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf);                                 ///< Generate five walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBoundingBox  (int InitialTag, double R, double Cf,bool Cohesion=false);                                              ///< Generate o bounding box enclosing the previous included particles.
    void GenBoundingPlane(int InitialTag, double R, double Cf,bool Cohesion=false);                                              ///< Same as GenBounding but only generates one pair of planes.
    void GenFromMesh     (Mesh::Generic & M, double R, double rho, bool cohesion=false, bool MC=true, double thickness = 0.0);   ///< Generate particles from a FEM mesh generator
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz,
    double rho, bool Cohesion, bool Periodic,size_t Randomseed, double fraction, double q=0.0);                                  ///< Generate a Voronoi Packing with dimensions Li and polihedra per side ni
    // Single particle addition
    void AddSphere   (int Tag, Vec3_t const & X, double R, double rho);                                                          ///< Add sphere
    void AddCube     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddTetra    (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a tetrahedron at position X with spheroradius R, side of length L and density rho
    void AddDrill    (int Tag, Vec3_t const & X, double R, double Lt, double Ll, double rho);                                    ///< A drill made as a combination of a cube and a pyramid.
    void AddRice     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a rice at position X with spheroradius R, side of length L and density rho
    void AddPlane    (int Tag, Vec3_t const & X, double R, double Lx,double Ly, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddVoroCell (int Tag, voronoicell_neighbor & VC, double R, double rho, bool Erode);                                     ///< Add a single voronoi cell, it should be built before tough
    void AddTorus    (int Tag, Vec3_t const & X, Vec3_t & N, double Rmax, double R, double rho);                                 ///< Add a single torus at position X with a normal N, circunference Rmax and spheroradius R

    // Methods
    void SetProps          (Dict & D);                                                                          ///< Set the properties of individual grains by dictionaries
    void Initialize        (double dt=0.0);                                                                     ///< Set the particles to a initial state and asign the possible insteractions
    void Solve             (double tf, double dt, double dtOut, ptFun_t ptSetup=NULL, ptFun_t ptReport=NULL,
                            char const * FileKey=NULL, bool RenderVideo=true);                                  ///< Run simulation
    void WritePOV          (char const * FileKey);                                                              ///< Write POV file
    void WriteBPY          (char const * FileKey);                                                              ///< Write BPY (Blender) file
#ifdef USE_HDF5    
    void Save              (char const * FileKey);                                                              ///< Save the current domain
    void Load              (char const * FileKey);                                                              ///< Load the domain form a file
#endif
    void BoundingBox       (Vec3_t & minX, Vec3_t & maxX);                                                      ///< Defines the rectangular box that encloses the particles.
    void Center            (Vec3_t C = Vec3_t(0.0,0.0,0.0));                                                    ///< Centers the domain around C
    void ResetInteractons  ();                                                                                  ///< Reset the interactons
    void ResetDisplacements();                                                                                  ///< Reset the displacements
    double MaxDisplacement ();                                                                                  ///< Calculate maximun displacement
    void ResetContacts     ();                                                                                  ///< Reset the displacements
    void EnergyOutput      (size_t IdxOut, std::ostream & OutFile);                                             ///< Output of the energy variables
    void GetGSD            (Array<double> & X, Array<double> & Y, Array<double> & D, size_t NDiv=10) const;     ///< Get the Grain Size Distribution
    void Clusters          ();                                                                                  ///< Check the bounded particles in the domain and how many connected clusters are still present

    // Access methods
    Particle       * GetParticle  (int Tag, bool Check=true);       ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    Particle const & GetParticle  (int Tag, bool Check=true) const; ///< Find first particle with Tag. Check => check if there are more than one particle with tag=Tag
    void             GetParticles (int Tag, Array<Particle*> & P);  ///< Find all particles with Tag

    // Auxiliar methods
    void   LinearMomentum  (Vec3_t & L);                    ///< Return total momentum of the system
    void   AngularMomentum (Vec3_t & L);                    ///< Return total angular momentum of the system
    double CalcEnergy      (double & Ekin, double & Epot);  ///< Return total energy of the system

    // Data
    bool                                              Initialized;                 ///< System (particles and interactons) initialized ?
    bool                                              Finished;                    ///< Has the simulation finished
    Array<Particle*>                                  Particles;                   ///< All particles in domain
    Array<Interacton*>                                Interactons;                 ///< All interactons
    Array<CInteracton*>                               CInteractons;                ///< Contact interactons
    Array<BInteracton*>                               BInteractons;                ///< Cohesion interactons
    Vec3_t                                            CamPos;                      ///< Camera position for POV
    double                                            Time;                        ///< Current time
    double                                            Evis;                        ///< Energy dissipated by the viscosity of the grains
    double                                            Efric;                       ///< Energy dissipated by friction
    double                                            Wext;                        ///< Work done by external forces
    double                                            Vs;                          ///< Volume occupied by the grains
    double                                            Ms;                          ///< Total mass of the particles
    double                                            Alpha;                       ///< Verlet distance
    void *                                            UserData;                    ///< Some user data
    String                                            FileKey;                     ///< File Key for output files
    size_t                                            idx_out;                     ///< Index of output
    set<pair<Particle *, Particle *> >                Listofpairs;                 ///< List of pair of particles associated per interacton for memory optimization
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

    // MPI Data and methods
#ifdef USE_MPI
    //Data
    int My_Id;                            ///< Id of the processor
    Array<Particle*> TheirBdryParticles;  ///< Array of particles in the boundary belonging to a different domain
    Array<Particle*> MyBdryParticles;     ///< Array of particles in the boundary belonging to this domain
    map<size_t,size_t> Global2LocalID;    ///< Map that register the relation between local and global ID


    //Methods
    void UpdateBoundaries ();             ///< Method to update boundaries for parallel computing
    double GeneralMaxDisplacement ();     ///< Maxdisplacement among all domains
    void CommunicateForce ();             ///< Communicate the force variables
    void CommunicateDynamic ();           ///< Communicate dynamic variables
#endif
    
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor & Destructor

inline Domain::Domain (void * UD)
    :  Initialized(false), Time(0.0), Alpha(0.05), UserData(UD)
{
    CamPos = 1.0, 2.0, 3.0;
#ifdef USE_MPI
    My_Id  = MPI::COMM_WORLD.Get_rank(); // processor ID
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
                                 , bool Cohesion, bool Periodic,size_t Randomseed, double fraction, double qin)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Adding Voronoi particles packing --------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    srand(Randomseed);
    const double x_min=-Lx/2.0, x_max=Lx/2.0;
    const double y_min=-Ly/2.0, y_max=Ly/2.0;
    const double z_min=-Lz/2.0, z_max=Lz/2.0;
    container con(x_min,x_max,y_min,y_max,z_min,z_max,nx,ny,nz, Periodic,Periodic,Periodic,8);
    int n = 0;
    for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = x_min+(i+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*(x_max-x_min)/nx;
                double y = y_min+(j+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*(y_max-y_min)/ny;
                double z = z_min+(k+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*(z_max-z_min)/nz;
                con.put (n,x,y,z);
                n++;
            }
        }
    }

    fpoint x,y,z,px,py,pz;
    container *cp = & con;
    voropp_loop l1(cp);
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
                        AddVoroCell(Tag,c,R,rho,true);
                        Vec3_t trans(x,y,z);
                        Particle * P = Particles[Particles.Size()-1];
                        P->Translate(trans);
                        P->Props.V = c.volume();
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

// Single particle addition

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

inline void Domain::AddVoroCell (int Tag, voronoicell_neighbor & VC, double R, double rho,bool Erode)
{
    Array<Vec3_t> V(VC.p);
    Array<Array <int> > E;
    Array<int> Eaux(2);
    for(int i=0;i<VC.p;i++) 
    {
        V[i] = Vec3_t(0.5*VC.pts[3*i],0.5*VC.pts[3*i+1],0.5*VC.pts[3*i+2]);
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

inline void Domain::AddTorus (int Tag, Vec3_t const & X, Vec3_t & N, double Rmax, double R, double rho)
{
    // Normalize normal vector
    N /= norm(N);

    // Create the 2 vertices that define the torus
    Vec3_t P1 = OrthoSys::e0 - dot(OrthoSys::e0,N)*N;
    P1       /= norm(P1);
    Vec3_t P2 = cross(N,P1);
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
        ResetInteractons();
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

inline void Domain::Solve (double tf, double dt, double dtOut, ptFun_t ptSetup, ptFun_t ptReport, char const * TheFileKey, bool RenderVideo)
{
    // Assigning some domain particles especifically to the output
    FileKey.Printf("%s",TheFileKey);
    idx_out = 0;

    // initialize particles
    Initialize (dt);

    // set the displacement of the particles to zero (for the Halo)
    ResetDisplacements();

    // build the map of possible contacts (for the Halo)
    ResetContacts();

    // Define boundary and domain particles
#ifdef USE_MPI
     UpdateBoundaries();
#endif

    // calc the total volume of particles (solids)
    Vs = 0.0;
    Ms = 0.0;
    for (size_t i=0; i<Particles.Size(); i++) 
    { 
        if (Particles[i]->IsFree())
        {
            Vs += Particles[i]->Props.V;
            Ms += Particles[i]->Props.m;
        }
    }

    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    // solve
    double t0   = Time;     // initial time
    double tout = t0+dtOut; // time position for output

    // report
    Finished = false;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    // string to output energy data, if user gives the FileKey
    std::ostringstream oss_energy; 
    EnergyOutput (idx_out, oss_energy);

    // run
    while (Time<tf)
    {

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

            // initialize the coordination (number of contacts per particle) number
            Particles[i]->Cn = 0.0;

            // external work added to the system by the fixed forces Ff
            Wext += dot(Particles[i]->Ff,Particles[i]->v)*dt;
        }
#ifdef USE_MPI
        for (size_t i=0;i<TheirBdryParticles.Size();i++)
        {
            TheirBdryParticles[i]->F = TheirBdryParticles[i]->Ff;
            TheirBdryParticles[i]->T = TheirBdryParticles[i]->Tf;
        }
#endif

        // calc contact forces: collision and bonding (cohesion)
        for (size_t i=0; i<Interactons.Size(); i++)
        {
            Interactons[i]->CalcForce (dt);
        }

        // calculate the collision energy
        for (size_t i=0; i<CInteractons.Size(); i++)
        {
            Evis  += CInteractons[i]-> dEvis;
            Efric += CInteractons[i]-> dEfric;
        }

        // tell the user function to update its data
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
#ifdef USE_MPI
        //Communicate all the relevant variables between different domains
        CommunicateForce();
#endif

        // move particles
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->Rotate    (dt);
            Particles[i]->Translate (dt);
        }

        // output
        if (Time>=tout)
        {
            idx_out++;
            if (BInteractons.Size()>0) Clusters();
            if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%08d", TheFileKey, idx_out);
                if(RenderVideo) WritePOV     (fn.CStr());
                EnergyOutput (idx_out, oss_energy);
            }
            tout += dtOut;
        }

#ifdef USE_MPI
        //Communicate all the relevant variables between different domains
        CommunicateDynamic();
#endif

        // next time position
        Time += dt;


#ifdef USE_MPI
        double maxdis = GeneralMaxDisplacement();
#else 
        double maxdis = MaxDisplacement();
#endif
        // update the Halos
        if (maxdis>Alpha)
        {
#ifdef USE_MPI
            UpdateBoundaries();
#endif
            ResetDisplacements();
            ResetContacts();
        }

        //if (fabs(Time-7.15013)<1.0e-5) WriteBPY("test");
    }

    // last output
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    // save energy data
    if (TheFileKey!=NULL)
    {
        String fn;
        fn.Printf("%s_energy.res",TheFileKey);
        std::ofstream fe(fn.CStr());
        fe << oss_energy.str();
        fe.close();
    }

    // info
    double Ekin, Epot, Etot;
    Etot = CalcEnergy (Ekin, Epot);
    printf("%s  Kinematic energy   = %g%s\n",TERM_CLR4, Ekin, TERM_RST);
    printf("%s  Potential energy   = %g%s\n",TERM_CLR4, Epot, TERM_RST);
    printf("%s  Total energy       = %g%s\n",TERM_CLR2, Etot, TERM_RST);
}

inline void Domain::WritePOV (char const * FileKey)
{
#ifdef USE_MPI
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors
    String fn(FileKey);
    fn.append(".pov");
    std::ofstream of;
    if (my_id==0)
    {
        of.open(fn.CStr(), std::ios::out);
        POVHeader (of);
        POVSetCam (of, CamPos, OrthoSys::O);
    }
    else
    {
        int command;
        MPI::COMM_WORLD.Recv (&command, /*number*/1, MPI::INT, /*destination*/my_id-1, 1);
        of.open(fn.CStr(), std::ios::out | std::ios::app);
    }
    for (size_t i=0; i<Particles.Size(); i++)
    {
        if (Particles[i]->IsFree())
        {
            if(Particles[i]->IsBroken) Particles[i]->Draw(of,"Black");
            else                       Particles[i]->Draw(of,"Red");
        }
        else Particles[i]->Draw(of,"Col_Glass_Bluish");
    }
    of.close();
    int command = 0;
    if (my_id < nprocs-1) MPI::COMM_WORLD.Send (&command, /*number*/1, MPI::INT, /*destination*/my_id+1, 1);
    
#else
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


    //if (BInteractons.Size()==0)
    //{ 
        //for (size_t i=0; i<Particles.Size(); i++)
        //{
            //if (Particles[i]->IsFree())
            //{
                //Particles[i]->Draw(of,"Red");
            //}
            //else Particles[i]->Draw(of,"Col_Glass_Bluish");
        //}
        //of.close();
    //}
    //else
    //{
        //Clusters();
        //Array <String> Colors(7);
        //Colors = "Red","Blue","Yellow","Gold","Green","Blue","Orange";
        //for (size_t i=0;i<Listofclusters.Size();i++)
        //{
            //for (size_t j=0;j<Listofclusters[i].Size();j++)
            //{
                //Particles[Listofclusters[i][j]]->Draw(of,Colors[i%7].CStr());
            //}
        //}
        //for (size_t i=0;i<Particles.Size();i++)
        //{
        //}
    //}
#endif
}

inline void Domain::WriteBPY (char const * FileKey)
{
#ifdef USE_MPI
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors
    String fn(FileKey);
    fn.append(".bpy");
    std::ofstream of;
    if (my_id==0)
    {
        of.open(fn.CStr(), std::ios::out);
        BPYHeader(of);
    }
    else
    {
        int command;
        MPI::COMM_WORLD.Recv (&command, /*number*/1, MPI::INT, /*destination*/my_id-1, 1);
        of.open(fn.CStr(), std::ios::out | std::ios::app);
    }
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Draw (of,"",true);
    of.close();
    int command = 0;
    if (my_id < nprocs-1) MPI::COMM_WORLD.Send (&command, /*number*/1, MPI::INT, /*destination*/my_id+1, 1);
#else
    String fn(FileKey);
    fn.append(".bpy");
    std::ofstream of(fn.CStr(), std::ios::out);
    BPYHeader(of);
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Draw (of,"",true);
    of.close();
#endif
}

#ifdef USE_HDF5
inline void Domain::Save (char const * FileKey)
{

    // Opening the file for writing
    String fn(FileKey);
    fn.append(".hdf5");
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
   
    H5Fclose(file_id);


}

inline void Domain::Load (char const * FileKey)
{
#ifdef USE_MPI
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors
#endif

    // Opening the file for reading
    String fn(FileKey);
    fn.append(".hdf5");
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

#ifdef USE_MPI
        //Domain decomposition by the yz plane crossing the origin
        if (nprocs>1) if((my_id==0&&X[0]>0)||(my_id==1&&X[0]<0)) continue;
#endif
        //std::cout << "hi" << std::endl;




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

}
#endif

inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
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

inline void Domain::Center(Vec3_t C)
{
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    Vec3_t Transport(-0.5*(maxX+minX));
    Transport += C;
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Translate(Transport);
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
#ifdef USE_MPI
    for (size_t i=0; i<Particles.Size(); i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=0; j<TheirBdryParticles.Size(); j++)
        {
            bool pj_has_vf = !TheirBdryParticles[j]->IsFree();

            bool close = (Distance(Particles[i]->x,TheirBdryParticles[j]->x)<=Particles[i]->Dmax+TheirBdryParticles[j]->Dmax+2*Alpha);

            // if both particles have any component specified or they are far away, don't create any intereactor
            if ((pi_has_vf && pj_has_vf) || !close ) continue;

            // if both particles are spheres (just one vertex)
            if (Particles[i]->Verts.Size()==1 && TheirBdryParticles[j]->Verts.Size()==1)
            {
                CInteractons.Push (new CInteractonSphere(Particles[i],TheirBdryParticles[j]));
            }

            // normal particles
            else
            {
                CInteractons.Push (new CInteracton(Particles[i],TheirBdryParticles[j]));
            }
        }
    }
#endif
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

#ifdef USE_MPI
    for (size_t i=0; i<Particles.Size(); i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        //std::cout << Time << " " << TheirBdryParticles.Size()<<std::endl;
        for (size_t j=0; j<TheirBdryParticles.Size(); j++)
        {
            bool pj_has_vf = !TheirBdryParticles[j]->IsFree();

            bool close = (Distance(Particles[i]->x,TheirBdryParticles[j]->x)<=Particles[i]->Dmax+TheirBdryParticles[j]->Dmax+2*Alpha);
            //std::cout << Time<< " " << std::endl;
            if ((pi_has_vf && pj_has_vf) || !close) continue;
            
            // checking if the interacton exist for that pair of particles
            bool exist = false;
            for (size_t k=0; k<CInteractons.Size(); k++)
            {
                bool index_i = CInteractons[k]->P1->Index==Particles[i]->Index;
                bool index_j = CInteractons[k]->P2->Index==TheirBdryParticles[j]->Index;
                if (index_i&&index_j)
                {
                    exist = true;
                    break;
                }
            }

            // If it doesn't add it to the CInteracton array
            if (!exist&&My_Id==1)
            {
                // if both particles are spheres (just one vertex)
                if (Particles[i]->Verts.Size()==1 && TheirBdryParticles[j]->Verts.Size()==1)
                {
                    CInteractons.Push (new CInteractonSphere(Particles[i],TheirBdryParticles[j]));
                }

                // normal particles
                else
                {
                    CInteractons.Push (new CInteracton(Particles[i],TheirBdryParticles[j]));
                    std::cout << "Interacton created between "<<Particles[i]->Index<<" and: "<<TheirBdryParticles[j]->Index<< " at time:" << Time <<std::endl;
                }
            }
        }
    }
#endif

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
        Ekin += Particles[i]->Ekin + Particles[i]->Erot;
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

#ifdef USE_MPI
// Methods needed for parallelization
inline void Domain::UpdateBoundaries ()
{
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    //int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors
    size_t number_new = 0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        if ((Particles[i]->MaxX()+Alpha>0.0&&my_id==0)||(Particles[i]->MinX()-Alpha<0.0&&my_id==1))
        {
            bool exist = false;
            for (size_t j=0; j<MyBdryParticles.Size();j++)
            {
                if (MyBdryParticles[j]->Index==Particles[i]->Index)
                {
                    exist = true;
                    break;
                }
            }
            if (!exist)
            {
                MyBdryParticles.Push(Particles[i]);
                number_new++;
            }
        }
    }

    // Communicating the new boundary particles
    if (my_id==0)
    {
        MPI::COMM_WORLD.Send (&number_new, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/1, 1);
        for (size_t i=0;i<number_new;i++)
        {
            //std::cout << "im" << my_id << " " <<std::endl;
            MyBdryParticles[MyBdryParticles.Size()-number_new+i]->SendParticle(1,2);
            //std::cout << "im" << my_id << std::endl;
        }

        size_t number_new_other;
        MPI::COMM_WORLD.Recv (&number_new_other, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/1, 3);
        for (size_t i=0;i<number_new_other;i++)
        {
            // Dummy particle
            Particle  *p = new Particle();
            p->ReceiveParticle(4);
            TheirBdryParticles.Push(p);
            Global2LocalID[TheirBdryParticles[TheirBdryParticles.Size()-1]->Index] = TheirBdryParticles.Size()-1;
            std::cout << "hi im "<<my_id<<" and i receive particle"<<TheirBdryParticles[0]->Index<< "at time " << Time <<std::endl;
        }
    }

    
    if (my_id==1)
    {
        size_t number_new_other;
        MPI::COMM_WORLD.Recv (&number_new_other, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/0, 1);
        for (size_t i=0;i<number_new_other;i++)
        {
            // Dummy particle
            Particle  *p = new Particle();
            p->ReceiveParticle(2);
            TheirBdryParticles.Push(p);
            Global2LocalID[TheirBdryParticles[TheirBdryParticles.Size()-1]->Index] = TheirBdryParticles.Size()-1;
            std::cout << "hi im "<<my_id<<" and i receive particle"<<TheirBdryParticles[0]->Index<< "at time " << Time <<std::endl;
            //std::cout << "im" << my_id << std::endl;
        }

        MPI::COMM_WORLD.Send (&number_new, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/0, 3);
        for (size_t i=0;i<number_new;i++)
        {
            MyBdryParticles[MyBdryParticles.Size()-number_new+i]->SendParticle(0,4);
        }
    }

}

inline double Domain::GeneralMaxDisplacement ()
{
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    //int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors
    double localmaxdis = MaxDisplacement();
    MPI::Request req_send = MPI::COMM_WORLD.Isend(&localmaxdis,/*number*/1,MPI::DOUBLE,/*destination*/(my_id==0?1:0),my_id);
    double othermaxdis;
    MPI::Request req_recv = MPI::COMM_WORLD.Irecv(&othermaxdis,/*number*/1,MPI::DOUBLE,/*destination*/MPI::ANY_SOURCE,(my_id==0?1:0));
    req_send.Wait();
    req_recv.Wait();
    return (localmaxdis>othermaxdis?localmaxdis:othermaxdis);
}

inline void Domain::CommunicateForce()
{
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    //int nprocs = MPI::COMM_WORLD.Get_size(); // Number of processors
    if (my_id==0)
    {
        for (size_t i=0;i<MyBdryParticles.Size();i++)
        {
             MyBdryParticles[i]->ReceiveForce(2);
        }
        for (size_t i=0;i<TheirBdryParticles.Size();i++)
        {
            TheirBdryParticles[i]->SendForce(1,2);
        }
    }
    if (my_id==1)
    {
        for (size_t i=0;i<TheirBdryParticles.Size();i++)
        {
            TheirBdryParticles[i]->SendForce(0,2);
        }
        for (size_t i=0;i<MyBdryParticles.Size();i++)
        {
             MyBdryParticles[i]->ReceiveForce(2);
        }
    }
}

inline void Domain::CommunicateDynamic()
{
    int my_id  = MPI::COMM_WORLD.Get_rank(); // processor ID
    if (my_id==0)
    {
        for (size_t i=0;i<MyBdryParticles.Size();i++)
        {
             MyBdryParticles[i]->SendDynamicParticle(1,2);
        }
        for (size_t i=0;i<TheirBdryParticles.Size();i++)
        {
             TheirBdryParticles[i]->ReceiveDynamicParticle(2);
        }
    }
    if (my_id==1)
    {
        for (size_t i=0;i<TheirBdryParticles.Size();i++)
        {
            TheirBdryParticles[i]->ReceiveDynamicParticle(2);
        }
        for (size_t i=0;i<MyBdryParticles.Size();i++)
        {
            MyBdryParticles[i]->SendDynamicParticle(0,2);
        }
    }
}
#endif

}; // namespace DEM

#endif // MECHSYS_DEM_DOMAIN_H
