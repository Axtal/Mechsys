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
#include <ctime>    // for std::clock

// Voro++
#include "src/voro++.cc"

// MechSys
#include <mechsys/dem/interacton.h>
#include <mechsys/util/array.h>
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/util/maps.h>

namespace DEM
{

class Domain
{
public:
    // Constructor
    Domain();

    // Destructor
    ~Domain();

    // Particle generation
    void GenSpheres      (int Tag, double L, size_t N, double rho, char const * Type, 
                          size_t Randomseed, double fraction);                                                  ///< General spheres
    void GenBox          (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf);                ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenBoundingBox  (int InitialTag, double R, double Cf);                                                 ///< Generate o bounding box enclosing the previous included particles.
    void GenFromMesh     (int Tag, Mesh::Generic const & M, double R, double rho);                              ///< Generate particles from a FEM mesh generator
    void GenFromVoro     (int Tag, container & VC, double R, double rho, double fraction=1.0);                  ///< Generate Particles from a Voronoi container
    void AddVoroPack     (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, 
                          double rho, bool Periodic,size_t Randomseed, double fraction);                        ///< Generate a Voronoi Packing wiht dimensions Li and polihedra per side ni

    // Single particle addition
    void AddSphere   (int Tag, Vec3_t const & X, double R, double rho);                                                          ///< Add sphere
    void AddCube     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddTetra    (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a tetrahedron at position X with spheroradius R, side of length L and density rho
    void AddRice     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a rice at position X with spheroradius R, side of length L and density rho
    void AddPlane    (int Tag, Vec3_t const & X, double R, double Lx,double Ly, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddVoroCell (int Tag, voronoicell & VC, double R, double rho, bool Erode);                                              ///< Add Voronoi cell

    // Methods
    void SetBC            (Dict & D);                                                                          ///< Set the dynamic conditions of individual grains by dictionaries
    void SetProps         (Dict & D);                                                                          ///< Set the properties of individual grains by dictionaries
    void Initialize       (double dt=0.0);                                                                     ///< Set the particles to a initial state and asign the possible insteractions
    void Solve    (double tf, double dt, double dtOut, char const * FileKey, bool RenderVideo = true);         ///< Run simulation
    void WritePOV         (char const * FileKey);                                                              ///< Write POV file
    void WriteBPY         (char const * FileKey);                                                              ///< Write BPY (Blender) file
    void BoundingBox      (Vec3_t & minX, Vec3_t & maxX);                                                      ///< Defines the rectangular box that encloses the particles.
    void Center           ();                                                                                  ///< Centers the domain
    void ResetInteractons ();                                                                                  ///< Reset the interactons
    virtual void Setup    (double,double) {};                                                                  ///< Special method depends on the Setup
    void EnergyOutput     (size_t IdxOut, std::ostream & OutFile);                                             ///< Output of the energy variables
    virtual void Output   (size_t IdxOut, std::ostream & OutFile) {};                                          ///< Output current state depends on the setup

    // Auxiliar methods
    void   LinearMomentum  (Vec3_t & L);                                                                     ///< Return total momentum of the system
    void   AngularMomentum (Vec3_t & L);                                                                     ///< Return total angular momentum of the system
    double CalcEnergy      (double & Ekin, double & Epot);                                                   ///< Return total energy of the system

    // Data
    double             Time;          ///< Current time
    bool               Initialized;   ///< System (particles and interactons) initialized ?
    size_t             InitialIndex;  ///< Tag the first index of the container
    Array<Particle*>   Particles;     ///< All particles in domain
    Array<Particle*>   FreeParticles; ///< Particles without constrains
    Array<Particle*>   TParticles;    ///< Particles with translation fixed
    Array<Particle*>   RParticles;    ///< Particles with rotation fixed
    Array<Particle*>   FParticles;    ///< Particles with applied force
    Array<Interacton*> Interactons;   ///< All interactons
    Vec3_t             CamPos;        ///< Camera position for POV
    double             Evis;          ///< Energy dissipated by the viscosity of the grains
    double             Efric;         ///< Energy dissipated by friction
    double             Wext;          ///< Work done by external forces


#ifdef USE_BOOST_PYTHON
    void PyAddSphere (int Tag, BPy::tuple const & X, double R, double rho)                                                         { AddSphere (Tag,Tup2Vec3(X),R,rho); }
    void PyAddCube   (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddCube  (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddTetra  (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddTetra (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddRice   (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddRice  (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddPlane  (int Tag, BPy::tuple const & X, double R, double Lx,double Ly, double rho, double Ang, BPy::tuple const & Ax) { Vec3_t a(Tup2Vec3(Ax)); AddPlane (Tag,Tup2Vec3(X),R,Lx,Ly,rho,Ang,&a); }
    void PySetCamPos (BPy::tuple const & PyCamPos)                                                                                 { CamPos = Tup2Vec3(PyCamPos); }
#endif
};

class TriaxialDomain: public Domain
{
public:
        
    //Constructor and destructor
    TriaxialDomain ();

    // Methods for specific simulations
    void SetTxTest (Vec3_t  const & Sigf,                                                            ///< Final state of stress
                    bVec3_t const & pEps,                                                            ///< Are strain rates prescribed ?
                    Vec3_t  const & dEpsdt);                                                         ///< Prescribed values of strain rate
    void ResetEps  ();                                                                               ///< Reset strains and re-calculate initial lenght of packing
    void Output    (size_t IdxOut, std::ostream & OutFile);                                          ///< Output current state of stress and strains.
    void Setup     (double,double);                                                                  ///< For the triaxial test it will measure or set the strains and stresses

    //Data
    double  Vs;       ///< Volume occupied by the grains 
    Vec3_t  Sig;      ///< Current stress state
    Vec3_t  DSig;     ///< Total stress increment to be applied by Solve => after 
    bVec3_t pSig;     ///< Prescribed stress ?
    Vec3_t  L0;       ///< Initial length of the packing
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor & Destructor

inline Domain::Domain ()
    : Time(0.0), Initialized(false)
{
    CamPos = 1.0, 2.0, 3.0;
}

inline Domain::~Domain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    for (size_t i=0; i<Interactons.Size(); ++i) if (Interactons[i]!=NULL) delete Interactons[i];
}

// Particle generation

inline void Domain::GenSpheres (int Tag, double L, size_t N, double rho,char const * Type, size_t Randomseed, double fraction)
{
    // find radius from the edge's length
    double start = std::clock();
    std::cout << "[1;33m\n--- Generating packing of spheres -------------[0m\n";
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
            if (rand()<fraction*RAND_MAX) AddSphere (Tag,pos,R,rho);
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
                    if (rand()<fraction*RAND_MAX) AddSphere(Tag,X,R,rho);
                }
            }
        }

    }
    else throw new Fatal ("Right now there are only two possible packings available the Normal and the HCP, packing %s is not implemented yet",Type);
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    std::cout << "[1;32m    Number of particles   = " << Particles.Size() << "[0m\n";
}

inline void Domain::GenBox (int InitialTag, double Lx, double Ly, double Lz, double R, double Cf)
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

    
    InitialIndex = Particles.Size();

    // add faces of box
    Vec3_t axis0(OrthoSys::e0); // rotation of face
    Vec3_t axis1(OrthoSys::e1); // rotation of face
    AddPlane (InitialTag,   Vec3_t(Lx/2.0,0.0,0.0),  R, Cf*Lz, Cf*Ly, 0.5, M_PI/2.0, &axis1);
    AddPlane (InitialTag-1, Vec3_t(-Lx/2.0,0.0,0.0), R, Cf*Lz, Cf*Ly, 0.5, M_PI/2.0, &axis1);
    AddPlane (InitialTag-2, Vec3_t(0.0,Ly/2.0,0.0),  R, Cf*Lx, Cf*Lz, 0.5, M_PI/2.0, &axis0);
    AddPlane (InitialTag-3, Vec3_t(0.0,-Ly/2.0,0.0), R, Cf*Lx, Cf*Lz, 0.5, M_PI/2.0, &axis0);
    AddPlane (InitialTag-4, Vec3_t(0.0,0.0,Lz/2.0),  R, Cf*Lx, Cf*Ly, 0.5);
    AddPlane (InitialTag-5, Vec3_t(0.0,0.0,-Lz/2.0), R, Cf*Lx, Cf*Ly, 0.5);

    // initialize walls (required when calculating the length of packing, since will use the CG)
    for (size_t i=0; i<6; ++i) Particles[InitialIndex+i]->Initialize();

}

inline void Domain::GenBoundingBox (int InitialTag, double R, double Cf)
{
    Center();
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    GenBox(InitialTag, maxX(0)-minX(0)+2*R, maxX(1)-minX(1)+2*R, maxX(2)-minX(2)+2*R, R, Cf);
}

inline void Domain::GenFromMesh (int Tag, Mesh::Generic const & M, double R, double rho)
{
    // info
    double start = std::clock();
    std::cout << "[1;33m\n--- Generating particles from mesh -----------------------------[0m\n";

    Array <Array <int> > Empty;
    for (size_t i=0; i<M.Cells.Size(); ++i)
    {
        Array<Mesh::Vertex*> const & verts = M.Cells[i]->V;
        size_t nverts = verts.Size();

        // verts
        Array<Vec3_t> V(nverts);
        for (size_t j=0; j<nverts; ++j)
        {
            V[j] = verts[j]->C;
        }

        // edges
        size_t nedges = Mesh::NVertsToNEdges3D[nverts];
        Array<Array <int> > E(nedges);
        for (size_t j=0; j<nedges; ++j)
        {
            E[j].Push (Mesh::NVertsToEdge3D[nverts][j][0]);
            E[j].Push (Mesh::NVertsToEdge3D[nverts][j][1]);
        }

        size_t nfaces = Mesh::NVertsToNFaces3D[nverts];
        size_t nvperf = Mesh::NVertsToNVertsPerFace3D[nverts];
        Array<Array <int> > F(nfaces);
        for (size_t j=0; j<nfaces; ++j)
        {
            for (size_t k=0; k<nvperf; ++k)
            {
                // TODO: check if face is planar or not
                F[j].Push(Mesh::NVertsToFace3D[nverts][j][k]);
            }
        }
        Erosion(V,E,F,R);

        // add particle
        Particles.Push (new Particle(Tag, V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    }

    // info
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    std::cout << "[1;32m    Number of particles   = " << Particles.Size() << "[0m\n";
}

inline void Domain::GenFromVoro (int Tag, container & VC, double R, double rho, double fraction)
{
    // info
    double start = std::clock();
    std::cout << "[1;33m\n--- Generating particles from Voronoi tessellation -------------[0m\n";

    fpoint x,y,z,px,py,pz;
    container *cp = & VC;
    voropp_loop l1(cp);
    int q,s;
    voronoicell c;
    s=l1.init(VC.ax,VC.bx,VC.ay,VC.by,VC.az,VC.bz,px,py,pz);
    do 
    {
        for(q=0;q<VC.co[s];q++) 
        {
            x=VC.p[s][VC.sz*q]+px;y=VC.p[s][VC.sz*q+1]+py;z=VC.p[s][VC.sz*q+2]+pz;
            if(x>VC.ax&&x<VC.bx&&y>VC.ay&&y<VC.by&&z>VC.az&&z<VC.bz) 
            {
                if(VC.compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z)) 
                {
                    if (rand()<fraction*RAND_MAX)
                    {
                        AddVoroCell(Tag,c,R,rho,true);
                        Vec3_t trans(x,y,z);
                        Particles[Particles.Size()-1]->Translate(trans);
                    }
                }
            }
        }
    } while((s=l1.inc(px,py,pz))!=-1);

    // info
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    std::cout << "[1;32m    Number of particles   = " << Particles.Size() << "[0m\n";
}

inline void Domain::AddVoroPack (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho, bool Periodic,size_t Randomseed, double fraction)
{
    srand(Randomseed);
    const double x_min=-Lx/2.0, x_max=Lx/2.0;
    const double y_min=-Ly/2.0, y_max=Ly/2.0;
    const double z_min=-Lz/2.0, z_max=Lz/2.0;
    container con(x_min,x_max,y_min,y_max,z_min,z_max,nx,ny,nz, Periodic,Periodic,Periodic,8);
    size_t n = 0;
    double qinter = 0.0;
    for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = x_min+(i+qinter+(1-2*qinter)*double(rand())/RAND_MAX)*(x_max-x_min)/nx;
                double y = y_min+(j+qinter+(1-2*qinter)*double(rand())/RAND_MAX)*(y_max-y_min)/ny;
                double z = z_min+(k+qinter+(1-2*qinter)*double(rand())/RAND_MAX)*(z_max-z_min)/nz;
                con.put (n,x,y,z);
                n++;
            }
        }
    }
    GenFromVoro(Tag,con,R,rho,fraction);
}

// Single particle addition

inline void Domain::AddSphere (int Tag, const Vec3_t & X, double R, double rho)
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

inline void Domain::AddCube (int Tag, const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
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
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
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
    delete Axis;
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
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
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
    delete Axis;
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
    if (Axis==NULL)
    {
        Angle   = (1.0*rand())/RAND_MAX*2*M_PI;
        Axis = new Vec3_t((1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX, (1.0*rand())/RAND_MAX);
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
    delete Axis;
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
    Particles[Particles.Size()-1]->Q    = q;
    Particles[Particles.Size()-1]->I    = 1.0,1.0,1.0;
    Particles[Particles.Size()-1]->V    = Lx*Ly*2*R;
    Particles[Particles.Size()-1]->m    = rho*Lx*Ly*2*R;
    Particles[Particles.Size()-1]->x    = X;
    Particles[Particles.Size()-1]->Ekin = 0.0;
    Particles[Particles.Size()-1]->Erot = 0.0;
    Particles[Particles.Size()-1]->Dmax = sqrt(Lx*Lx+Ly*Ly)+R;
    Particles[Particles.Size()-1]->PropsReady = true;
    // clean up
    if (!ThereisanAxis) delete Axis;
}

inline void Domain::AddVoroCell (int Tag, voronoicell & VC, double R, double rho,bool Erode)
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
    if (Erode) Erosion(V,E,F,R);

    // add particle
    Particles.Push (new Particle(Tag,V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
}

// Methods

inline void Domain::SetBC (Dict & D)
{
    TParticles.Resize(0);
    RParticles.Resize(0);
    FParticles.Resize(0);
    for (size_t i =0 ; i<Particles.Size(); i++)
    {
        bool free_particle = true;
        for (size_t j=0; j<D.Keys.Size(); ++j)
        {
            int tag = D.Keys[j];
            if (tag==Particles[i]->Tag)
            {
                SDPair const & p = D(tag);
                if (p.HasKey("vx") || p.HasKey("vy") || p.HasKey("vz")) // translation specified
                {
                    TParticles.Push (Particles[i]);
                    TParticles[TParticles.Size()-1]->v = Vec3_t(p("vx"),p("vy"),p("vz"));
                    free_particle = false;
                }
                if (p.HasKey("wx") || p.HasKey("wy") || p.HasKey("wz")) // rotation specified
                {
                    RParticles.Push (Particles[i]);
                    RParticles[RParticles.Size()-1]->w = Vec3_t(p("wx"),p("wy"),p("wz"));
                    free_particle = false;
                }
                if (p.HasKey("fx") || p.HasKey("fy") || p.HasKey("fz")) // Applied force specified
                {
                    FParticles.Push (Particles[i]);
                    FParticles[FParticles.Size()-1]->Ff = Vec3_t(p("fx"),p("fy"),p("fz"));
                    free_particle = false;
                }
            }
        }
        if (free_particle) FreeParticles.Push (Particles[i]);
    }
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
                    Particles[i]->Gn = p("Gn");
                }
                if (p.HasKey("Gt"))
                {
                    Particles[i]->Gt = p("Gt");
                }
                if (p.HasKey("Kn"))
                {
                    Particles[i]->Kn = p("Kn");
                }
                if (p.HasKey("Kt"))
                {
                    Particles[i]->Kt = p("Kt");
                }
                if (p.HasKey("Mu"))
                {
                    Particles[i]->Mu = p("Mu");
                }
                if (p.HasKey("Beta"))
                {
                    Particles[i]->Beta = p("Beta");
                }
                if (p.HasKey("Eta"))
                {
                    Particles[i]->Eta = p("Eta");
                }
            }
        }
    }
}

inline void Domain::Initialize (double dt)
{
    if (!Initialized)
    {
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->Initialize();
            Particles[i]->InitializeVelocity(dt);
        }

        // info
        double start = std::clock();
        std::cout << "[1;33m\n--- Initializing particles -------------------------------------[0m\n";
        ResetInteractons();
        // set flag
        Initialized = true;

        // info
        double Ekin, Epot, Etot;
        Etot = CalcEnergy (Ekin, Epot);
        double total = std::clock() - start;
        std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
        std::cout << "[1;35m    Kinematic energy      = " << Ekin << "[0m\n";
        std::cout << "[1;35m    Potential energy      = " << Epot << "[0m\n";
        std::cout << "[1;35m    Total energy          = " << Etot << "[0m\n";
    }
    else
    {
        for (size_t i=0; i<TParticles.Size(); i++)
        {
            TParticles[i]->InitializeVelocity(dt);
        }
    }
}

inline void Domain::Solve (double tf, double dt, double dtOut, char const * FileKey, bool RenderVideo)
{
    // initialize
    if (FreeParticles.Size()==0) FreeParticles = Particles;
    Initialize (dt);

    // info
    double start = std::clock();
    std::cout << "[1;33m\n--- Solving ----------------------------------------------------[0m\n";

    // open file for walls
    String fnw;
    fnw.Printf("%s_walls.res",FileKey);
    std::ofstream fw(fnw.CStr());

    // open file for energy
    String fne;
    fne.Printf("%s_energy.res",FileKey);
    std::ofstream fe(fne.CStr());

    // solve
    double t0      = Time; // initial time
    size_t idx_out = 0;    // index of output
    double tout    = t0;  // time position for output

    //Initializing the energies
    Evis = 0.0;
    Efric = 0.0;
    Wext = 0.0;

    // run
    while (Time<tf)
    {
        // initialize forces and torques
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Tf;
            Wext += dot(Particles[i]->Ff,Particles[i]->v)*dt;
        }

        // calc force
        for (size_t i=0; i<Interactons.Size(); i++)
        {
            Interactons[i]->CalcForce (dt);
            Evis += Interactons[i]-> dEvis;
            Efric += Interactons[i]-> dEfric;
        }

        Setup(dt, tf-t0);

        // move free particles
        for (size_t i=0; i<FreeParticles.Size(); i++)
        {
            FreeParticles[i]->Rotate    (dt);
            FreeParticles[i]->Translate (dt);
        }

        // particles with translation constrained
        for (size_t i=0; i<TParticles.Size(); i++) 
        {
            TParticles[i]->F = 0.0,0.0,0.0;
            TParticles[i]->Translate(dt);
        }

        // particles with rotation constrained
        for (size_t i=0; i<RParticles.Size(); i++)
        {
            RParticles[i]->T = 0.0,0.0,0.0;
            RParticles[i]->Rotate(dt);
        }

        // particles with forces applied
        for (size_t i=0; i<FParticles.Size(); i++)
        {
            double norm_Ff = norm(FParticles[i]->Ff);
            if (norm_Ff>1.0e-7)
            {
                // set F as the projection of Ff
                Vec3_t unit_Ff = FParticles[i]->Ff/norm_Ff; // unitary vector parallel to Ff
                FParticles[i]->F = dot(FParticles[i]->F,unit_Ff)*unit_Ff;
                FParticles[i]->Translate (dt);
            }
        }


        // next time position
        Time += dt;

        // output
        if (Time>=tout)
        {
            String fn;
            fn.Printf ("%s_%08d", FileKey, idx_out);
            if (RenderVideo) WritePOV  (fn.CStr());
            EnergyOutput (idx_out, fe);
            Output (idx_out, fw);
            tout += dtOut;
            idx_out++;
        }
    }

    // close file
    fw.close();

    // info
    double Ekin, Epot, Etot;
    Etot = CalcEnergy (Ekin, Epot);
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    std::cout << "[1;35m    Kinematic energy      = " << Ekin << "[0m\n";
    std::cout << "[1;35m    Potential energy      = " << Epot << "[0m\n";
    std::cout << "[1;35m    Total energy          = " << Etot << "[0m\n";
}

inline void Domain::WritePOV (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".pov");
    std::ofstream of(fn.CStr(), std::ios::out);
    POVHeader (of);
    POVSetCam (of, CamPos, OrthoSys::O);
    for (size_t i=0; i<FreeParticles.Size(); i++) FreeParticles[i]->Draw (of,"Red");
    for (size_t i=0; i<TParticles.Size(); i++) TParticles[i]->Draw (of,"Col_Glass_Bluish");
    for (size_t i=0; i<RParticles.Size(); i++) RParticles[i]->Draw (of,"Col_Glass_Bluish");
    for (size_t i=0; i<FParticles.Size(); i++) FParticles[i]->Draw (of,"Col_Glass_Bluish");
    of.close();
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

inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    minX = Vec3_t(Particles[0]->MinX(), Particles[0]->MinY(), Particles[0]->MinZ());
    maxX = Vec3_t(Particles[0]->MaxX(), Particles[0]->MaxY(), Particles[0]->MaxZ());
    for (size_t i=1; i<Particles.Size(); i++)
    {
        if (minX(0)>Particles[i]->MinX()) minX(0) = Particles[i]->MinX();
        if (minX(1)>Particles[i]->MinY()) minX(1) = Particles[i]->MinY();
        if (minX(2)>Particles[i]->MinZ()) minX(2) = Particles[i]->MinZ();
        if (maxX(0)<Particles[i]->MaxX()) maxX(0) = Particles[i]->MaxX();
        if (maxX(1)<Particles[i]->MaxY()) maxX(1) = Particles[i]->MaxY();
        if (maxX(2)<Particles[i]->MaxZ()) maxX(2) = Particles[i]->MaxZ();
    }
}

inline void Domain::Center()
{
    Vec3_t minX,maxX;
    BoundingBox(minX,maxX);
    Vec3_t Transport(-0.5*(maxX+minX));
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Translate(Transport);
}

inline void Domain::ResetInteractons()
{
    for (size_t i=0; i<Interactons.Size(); ++i) if (Interactons[i]!=NULL) delete Interactons[i];
    Interactons.Resize(0);
    for (size_t i=0; i<FreeParticles.Size()-1; i++)
    {
        for (size_t j=i+1; j<FreeParticles.Size(); j++)
        {
            if (FreeParticles[i]->Verts.Size()==1&&FreeParticles[j]->Verts.Size()==1) Interactons.Push (new InteractonSphere(FreeParticles[i],FreeParticles[j]));
            else Interactons.Push (new Interacton(FreeParticles[i],FreeParticles[j]));
        }
    }

    for (size_t i=0; i<FreeParticles.Size(); i++)
    {
        for (size_t j=0; j<FParticles.Size(); j++)
        {
            Interactons.Push (new Interacton(FreeParticles[i],FParticles[j]));
        }

        for (size_t j=0; j<RParticles.Size(); j++)
        {
            Interactons.Push (new Interacton(FreeParticles[i],RParticles[j]));
        }

        for (size_t j=0; j<TParticles.Size(); j++)
        {
            Interactons.Push (new Interacton(FreeParticles[i],TParticles[j]));
        }
    }
}

// Auxiliar methods

inline void Domain::LinearMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        L += Particles[i]->m*Particles[i]->v;
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
        L += Particles[i]->m*cross(Particles[i]->x,Particles[i]->v)+t2;
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
    for (size_t i=0; i<Interactons.Size(); i++)
    {
        Epot += Interactons[i]->Epot;
    }

    // total energy
    return Ekin + Epot;
}

inline void Domain::EnergyOutput (size_t IdxOut, std::ostream & OF)
{

    // output triaxial test data
    // header
    if (IdxOut==0)
    {
        OF << Util::_10_6 << "Time" << Util::_8s << "Ekin" << Util::_8s << "Epot" << Util::_8s << "Evis" << Util::_8s << "Efric" << Util::_8s << "Wext" << std::endl;
    }
    double Ekin,Epot;
    CalcEnergy(Ekin,Epot);
    OF << Util::_10_6 << Time << Util::_8s << Ekin << Util::_8s << Epot << Util::_8s << Evis << Util::_8s << Efric << Util::_8s << Wext << std::endl;
    

}
// Methods for specific simulations

inline TriaxialDomain::TriaxialDomain ()
{
    Time = 0.0;
    Initialized = false;
    Sig    = 0.0,   0.0,   0.0;
    DSig   = 0.0,   0.0,   0.0;
    pSig   = false, false, false;
    CamPos = 1.0, 2.0, 3.0;
}

inline void TriaxialDomain::SetTxTest (Vec3_t const & Sigf, bVec3_t const & pEps, Vec3_t const & dEpsdt)
{
    // info
    std::cout << "[1;33m\n--- Setting up Triaxial Test -------------------------------------[0m\n";
    double start = std::clock();

    // resize arrays
    TParticles   .Resize(0);
    RParticles   .Resize(0);
    FParticles   .Resize(0);
    FreeParticles.Resize(0);

    // initialize some of the free particles
    for (size_t i=0; i<InitialIndex;     i++) FreeParticles.Push(Particles[i]);
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Initialize ();

    Vs = 0.0;
    for (size_t i=0; i<FreeParticles.Size(); i++)
    {
        Vs += FreeParticles[i]->V;
    }
    // total stress increment
    DSig = Sigf - Sig;

    // assume all strains prescribed by default
    pSig = false, false, false;

    // Eps(0) prescribed ?
    Vec3_t veloc, force;
    if (pEps(0))
    {
        double height = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0));
        veloc = 0.5*dEpsdt(0)*height, 0.0, 0.0;
        TParticles.Push(Particles[InitialIndex]);
        TParticles[TParticles.Size()-1]->v = veloc;
        TParticles[TParticles.Size()-1]->Ff = 0.0,0.0,0.0;
        TParticles.Push(Particles[InitialIndex+1]);
        TParticles[TParticles.Size()-1]->v = -veloc;
        TParticles[TParticles.Size()-1]->Ff = 0.0,0.0,0.0;
    }
    else // Sig(0) prescribed
    {
        double area = (Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1))*(Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2));
        force = Sig(0)*area, 0.0, 0.0;
        FParticles.Push(Particles[InitialIndex]);
        FParticles[FParticles.Size()-1]->Ff = force;
        FParticles.Push(Particles[InitialIndex+1]);
        FParticles[FParticles.Size()-1]->Ff = -force;
        pSig(0) = true;
    }

    // Eps(1) prescribed ?
    if (pEps(1))
    {
        double height = (Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1));
        veloc = 0.0, 0.5*dEpsdt(1)*height, 0.0;
        TParticles.Push(Particles[InitialIndex+2]);
        TParticles[TParticles.Size()-1]->v = veloc;
        TParticles[TParticles.Size()-1]->Ff = 0.0,0.0,0.0;
        TParticles.Push(Particles[InitialIndex+3]);
        TParticles[TParticles.Size()-1]->v = -veloc;
        TParticles[TParticles.Size()-1]->Ff = 0.0,0.0,0.0;
    }
    else // Sig(1) presscribed
    {
        double area = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2));
        force = 0.0, Sig(1)*area, 0.0;
        FParticles.Push(Particles[InitialIndex+2]);
        FParticles[FParticles.Size()-1]->Ff = force;
        FParticles.Push(Particles[InitialIndex+3]);
        FParticles[FParticles.Size()-1]->Ff = -force;
        pSig(1) = true;
    }

    // Eps(2) prescribed ?
    if (pEps(2))
    {
        double height = (Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2));
        veloc = 0.0, 0.0, 0.5*dEpsdt(2)*height;
        TParticles.Push(Particles[InitialIndex+4]);
        TParticles[TParticles.Size()-1]->v = veloc;
        TParticles[TParticles.Size()-1]->Ff = 0.0,0.0,0.0;
        TParticles.Push(Particles[InitialIndex+5]);
        TParticles[TParticles.Size()-1]->v = -veloc;
        TParticles[TParticles.Size()-1]->Ff = 0.0,0.0,0.0;
    }
    else // Sig(2) presscribed
    {
        double area = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1));
        force = 0.0, 0.0, Sig(2)*area;
        FParticles.Push(Particles[InitialIndex+4]);
        FParticles[FParticles.Size()-1]->Ff = force;
        FParticles.Push(Particles[InitialIndex+5]);
        FParticles[FParticles.Size()-1]->Ff = -force;
        pSig(2) = true;
    }

    // info
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
}

inline void TriaxialDomain::ResetEps ()
{
    L0(0) = Particles[InitialIndex  ]->x(0)-Particles[InitialIndex+1]->x(0);
    L0(1) = Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1);
    L0(2) = Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2);
}

inline void TriaxialDomain::Setup (double dt,double tspan)
{
        Vec3_t force;
        bool   update_sig = false;
        if (pSig(0))
        {
            double area = (Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1))*(Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2));
            force = Sig(0)*area, 0.0, 0.0;
            Particles[InitialIndex  ]->Ff =  force;
            Particles[InitialIndex+1]->Ff = -force;
            update_sig = true;
        }
        else 
        {
            double area = (Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1))*(Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2));
            Sig(0) = -0.5*(fabs(Particles[InitialIndex  ]->F(0))+fabs(Particles[InitialIndex+1]->F(0)))/area;
        }
        if (pSig(1))
        {
            double area = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2));
            force = 0.0, Sig(1)*area, 0.0;
            Particles[InitialIndex+2]->Ff =  force;
            Particles[InitialIndex+3]->Ff = -force;
            update_sig = true;
        }
        else 
        {
            double area = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2));
            Sig(1) = -0.5*(fabs(Particles[InitialIndex+2]->F(1))+fabs(Particles[InitialIndex+3]->F(1)))/area;
        }
        if (pSig(2))
        {
            double area = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1));
            force = 0.0, 0.0, Sig(2)*area;
            Particles[InitialIndex+4]->Ff =  force;
            Particles[InitialIndex+5]->Ff = -force;
            update_sig = true;
        }
        else 
        {
            double area = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1));
            Sig(2) = -0.5*(fabs(Particles[InitialIndex+4]->F(2))+fabs(Particles[InitialIndex+5]->F(2)))/area;
        }
        if (update_sig) Sig += dt*DSig/(tspan);
}

inline void TriaxialDomain::Output (size_t IdxOut, std::ostream & OF)
{

    // output triaxial test data
    // header
    if (IdxOut==0)
    {
        OF << Util::_10_6 << "Time" << Util::_8s << "sx" << Util::_8s << "sy" << Util::_8s << "sz";
        OF <<                          Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez";
        OF <<                          Util::_8s << "vr" << Util::_8s << "Nc" << Util::_8s << "Nsc" << "\n";
    }

    // stress
    OF << Util::_10_6 << Time << Util::_8s << Sig(0) << Util::_8s << Sig(1) << Util::_8s << Sig(2);

    // strain
    OF << Util::_8s << (Particles[InitialIndex  ]->x(0)-Particles[InitialIndex+1]->x(0)-L0(0))/L0(0);
    OF << Util::_8s << (Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1)-L0(1))/L0(1);
    OF << Util::_8s << (Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2)-L0(2))/L0(2);

    // void ratio
    double volumecontainer = (Particles[InitialIndex  ]->x(0)-Particles[InitialIndex+1]->x(0)-Particles[InitialIndex  ]->R+Particles[InitialIndex+1]->R)*
                             (Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1)-Particles[InitialIndex+2]->R+Particles[InitialIndex+3]->R)*
                             (Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2)-Particles[InitialIndex+4]->R+Particles[InitialIndex+5]->R);

    OF << Util::_8s << (volumecontainer-Vs)/Vs;

    // Number of contacts and number of sliding contacts
    double Nc = 0;
    double Nsc = 0;
    for (size_t i=0; i<Interactons.Size(); i++)
    {
        Nc += Interactons[i]->Nc;
        Nsc += Interactons[i]->Nsc;
    }

    OF << Util::_8s << Nc << Util::_8s << Nsc;

    OF << std::endl;

    
}
}; // namespace DEM

#endif // MECHSYS_DEM_DOMAIN_H
