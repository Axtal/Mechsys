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
#include "dem/interacton.h"
#include "util/array.h"
#include "util/util.h"
#include "util/numstreams.h"
#include "mesh/mesh.h"
#include "util/maps.h"

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
    void GenSpheres  (int Tag, double L, size_t N, double rho, char const * Type);                          ///< General spheres
    void GenBox      (int InitialTag, double Lx, double Ly, double Lz, double R, bool Triaxial, double Cf); ///< Generate six walls with successive tags. Cf is a coefficient to make walls bigger than specified in order to avoid gaps
    void GenFromMesh (int Tag, Mesh::Generic const & M, double R, double rho);                              ///< Generate particles from a FEM mesh generator
    void GenFromVoro (int Tag, container & VC, double R, double rho);                                       ///< Generate Particles from a Voronoi container
    void AddVoroPack (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, 
                      double rho, bool Periodic);                                                           ///< Generate a Voronoi Packing wiht dimensions Li and polihedra per side ni

    // Single particle addition
    void AddSphere   (int Tag, Vec3_t const & X, double R, double rho);                                                          ///< Add sphere
    void AddCube     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddTetra    (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a tetrahedron at position X with spheroradius R, side of length L and density rho
    void AddRice     (int Tag, Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL);            ///< Add a rice at position X with spheroradius R, side of length L and density rho
    void AddPlane    (int Tag, Vec3_t const & X, double R, double Lx,double Ly, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho
    void AddVoroCell (int Tag, voronoicell & VC, double R, double rho, bool Erode);                                              ///< Add Voronoi cell

    // Methods
    void SetBC      (Dict & D);                                                 ///< Set the properties of individual grains by dictionaries
    void Initialize (double dt=0.0);                                            ///< Set the particles to a initial state and asign the possible insteractions
    void Solve      (double tf, double dt, double dtOut, char const * FileKey); ///< Run simulation
    void WritePOV   (char const * FileKey);                                     ///< Write POV file
    void WriteBPY   (char const * FileKey);                                     ///< Write BPY (Blender) file

    // Methods for specific simulations
    void SetTxTest (Vec3_t  const & Sigf,    ///< Final state of stress
                    bVec3_t const & pEps,    ///< Are strain rates prescribed ?
                    Vec3_t  const & dEpsdt); ///< Prescribed values of strain rate
    void ResetEps  ();                       ///< Reset strains and re-calculate initial lenght of packing

    // Auxiliar methods
    void   LinearMomentum  (Vec3_t & L);                                                  ///< Return total momentum of the system
    void   AngularMomentum (Vec3_t & L);                                                  ///< Return total angular momentum of the system
    double CalcEnergy      (double & Ekin, double & Epot);                                ///< Return total energy of the system
    void   Output          (char const * FileKey, size_t IdxOut, std::ostream & OutFile); ///< Output current state

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

    // Data for specific simulations
    bool    IsTxTest; ///< Check if the test is a triaxial test
    Vec3_t  Sig;      ///< Current stress state
    Vec3_t  DSig;     ///< Total stress increment to be applied by Solve => after 
    bVec3_t pSig;     ///< Prescribed stress ?
    Vec3_t  L0;       ///< Initial length of the packing

#ifdef USE_BOOST_PYTHON
    void PyAddSphere (int Tag, BPy::tuple const & X, double R, double rho)                                                         { AddSphere (Tag,Tup2Vec3(X),R,rho); }
    void PyAddCube   (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddCube  (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddTetra  (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddTetra (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddRice   (int Tag, BPy::tuple const & X, double R, double L, double rho, double Ang, BPy::tuple const & Ax)            { Vec3_t a(Tup2Vec3(Ax)); AddRice  (Tag,Tup2Vec3(X),R,L,rho,Ang,&a); }
    void PyAddPlane  (int Tag, BPy::tuple const & X, double R, double Lx,double Ly, double rho, double Ang, BPy::tuple const & Ax) { Vec3_t a(Tup2Vec3(Ax)); AddPlane (Tag,Tup2Vec3(X),R,Lx,Ly,rho,Ang,&a); }
    void PySetTxTest (BPy::tuple const & Sigf, BPy::tuple const & pEps, BPy::tuple  const & dEpsdt)                                { SetTxTest (Tup2Vec3(Sigf),Tup2Vec3(pEps),Tup2Vec3(dEpsdt)); }
    void PySetCamPos (BPy::tuple const & PyCamPos)                                                                                 { CamPos = Tup2Vec3(PyCamPos); }
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor & Destructor

inline Domain::Domain ()
    : Time(0.0), Initialized(false), IsTxTest(false)
{
    Sig    = 0.0,   0.0,   0.0;
    DSig   = 0.0,   0.0,   0.0;
    pSig   = false, false, false;
    CamPos = 1.0, 2.0, 3.0;
}

inline Domain::~Domain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    for (size_t i=0; i<Interactons.Size(); ++i) if (Interactons[i]!=NULL) delete Interactons[i];
}

// Particle generation

inline void Domain::GenSpheres (int Tag, double L, size_t N, double rho,char const * Type)
{
    // find radius from the edge's length
    double start = std::clock();
    std::cout << "[1;33m\n--- Generating packing of spheres -------------[0m\n";

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
            AddSphere (Tag,pos,R,rho);
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
                    AddSphere(Tag,X,R,rho);
                }
            }
        }

    }
    else throw new Fatal ("Right now there are only two possible packings available the Normal and the HCP, packing %s is not implemented yet",Type);
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    std::cout << "[1;32m    Number of particles   = " << Particles.Size() << "[0m\n";
}

inline void Domain::GenBox (int InitialTag, double Lx, double Ly, double Lz, double R, bool Triaxial, double Cf)
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

    // set triaxial test
    if (Triaxial) InitialIndex = Particles.Size();
    IsTxTest = Triaxial;

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

    // reset strains and find initial length of packing
    ResetEps ();
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

inline void Domain::GenFromVoro (int Tag, container & VC, double R, double rho)
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
                    AddVoroCell(Tag,c,R,rho,true);
                    Vec3_t trans(x,y,z);
                    Particles[Particles.Size()-1]->Translate(trans);
                }
            }
        }
    } while((s=l1.inc(px,py,pz))!=-1);

    // info
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    std::cout << "[1;32m    Number of particles   = " << Particles.Size() << "[0m\n";
}

inline void Domain::AddVoroPack (int Tag, double R, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho, bool Periodic)
{
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
    GenFromVoro(Tag,con,R,rho);
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
    //std::cout << Particles.Size() << " " << FParticles.Size() << " " << FreeParticles.Size() << " " << RParticles.Size() << std::endl;
}

inline void Domain::Initialize (double dt)
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->Initialize (dt);
    }
    if (!Initialized)
    {
        // info
        double start = std::clock();
        std::cout << "[1;33m\n--- Initializing particles -------------------------------------[0m\n";

        for (size_t i=0; i<FreeParticles.Size()-1; i++)
        {
            // initialize interactons
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
        for (size_t i = 0;i < TParticles.Size();i++)
        {
            TParticles[i]->Initialize(dt);
        }
    }
}

inline void Domain::Solve (double tf, double dt, double dtOut, char const * FileKey)
{
    // initialize
    if (FreeParticles.Size()==0) FreeParticles = Particles;
    Initialize (dt);

    // info
    double start = std::clock();
    std::cout << "[1;33m\n--- Solving ----------------------------------------------------[0m\n";

    // open file
    String fnw;
    fnw.Printf("%s_walls.res",FileKey);
    std::ofstream fw(fnw.CStr());

    // solve
    double t0      = Time; // initial time
    size_t idx_out = 0;    // index of output
    double tout    = 0.0;  // time position for output

    // output initial state
    Output (FileKey, idx_out, fw);
    idx_out++;

    // run
    while (Time<tf)
    {
        // initialize forces and torques
        for (size_t i=0; i<Particles.Size(); i++)
        {
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Tf;
        }

        // calc force
        for (size_t i=0; i<Interactons.Size(); i++) Interactons[i]->CalcForce ();

        // move free particles
        for (size_t i=0; i<FreeParticles.Size(); i++)
        {
            FreeParticles[i]->Rotate    (dt);
            FreeParticles[i]->Translate (dt);
        }

        // particles with translation constrained
        Array<Vec3_t> FT(TParticles.Size());
        for (size_t i=0; i<TParticles.Size(); i++) 
        {
            FT[i] = TParticles[i]->F;
            TParticles[i]->F = 0.0,0.0,0.0;
            TParticles[i]->Translate(dt);
        }

        // particles with rotation constrained
        Array<Vec3_t> TT(RParticles.Size());
        for (size_t i=0; i<RParticles.Size(); i++)
        {
            TT[i] = RParticles[i]->T;
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

        // calculate force for triaxial test and update stresses
        if (IsTxTest)
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
            if (pSig(1))
            {
                double area = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2));
                force = 0.0, Sig(1)*area, 0.0;
                Particles[InitialIndex+2]->Ff =  force;
                Particles[InitialIndex+3]->Ff = -force;
                update_sig = true;
            }
            if (pSig(2))
            {
                double area = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1));
                force = 0.0, 0.0, Sig(2)*area;
                Particles[InitialIndex+4]->Ff =  force;
                Particles[InitialIndex+5]->Ff = -force;
                update_sig = true;
            }
            if (update_sig) Sig += dt*DSig/(tf-t0);
        }

        // next time position
        Time += dt;

        // output
        if (Time>=tout)
        {
            Output (FileKey, idx_out, fw);
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

// Methods for specific simulations

inline void Domain::SetTxTest (Vec3_t const & Sigf, bVec3_t const & pEps, Vec3_t const & dEpsdt)
{
    // check if it is a triaxial test or not
    if (!IsTxTest) throw new Fatal("Domain::SetTxTest: The container for the triaxial test has not be defined");

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

    // total stress increment
    DSig = Sigf - Sig;

    // assume all strains prescribed by default
    pSig = false, false, false;

    // Eps(0) prescribed ?
    Vec3_t veloc, force;
    if (pEps(0))
    {
        double height = (Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0));
        veloc = -0.5*dEpsdt(0)*height, 0.0, 0.0;
        TParticles.Push(Particles[InitialIndex]);
        TParticles[TParticles.Size()-1]->v = veloc;
        TParticles.Push(Particles[InitialIndex+1]);
        TParticles[TParticles.Size()-1]->v = -veloc;
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
        veloc = 0.0, -0.5*dEpsdt(1)*height, 0.0;
        TParticles.Push(Particles[InitialIndex+2]);
        TParticles[TParticles.Size()-1]->v = veloc;
        TParticles.Push(Particles[InitialIndex+3]);
        TParticles[TParticles.Size()-1]->v = -veloc;
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
        veloc = 0.0, 0.0, -0.5*dEpsdt(2)*height;
        TParticles.Push(Particles[InitialIndex+4]);
        TParticles[TParticles.Size()-1]->v = veloc;
        TParticles.Push(Particles[InitialIndex+5]);
        TParticles[TParticles.Size()-1]->v = -veloc;
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

inline void Domain::ResetEps ()
{
    L0(0) = Particles[InitialIndex  ]->x(0)-Particles[InitialIndex+1]->x(0);
    L0(1) = Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1);
    L0(2) = Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2);
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

inline void Domain::Output (char const * FileKey, size_t IdxOut, std::ostream & OF)
{
    String fn;
    fn.Printf ("%s_%08d", FileKey, IdxOut);
    WritePOV  (fn.CStr());

    // output triaxial test data
    if (IsTxTest)
    {
        // header
        if (IdxOut==0)
        {
            OF << Util::_10_6 << "Time" << Util::_8s << "sx" << Util::_8s << "sy" << Util::_8s << "sz";
            OF <<                          Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez" << "\n";
        }

        // stress
        double sig2 = Sig(2);
        //if (sig2==0.0) sig2 = 0.5*(fabs(FT[0](2))+fabs(FT[1](2)))/((Particles[InitialIndex]->x(0)-Particles[InitialIndex+1]->x(0))*(Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1)));
        OF << Util::_10_6 << Time << Util::_8s << Sig(0) << Util::_8s << Sig(1) << Util::_8s << sig2;

        // strain
        OF << Util::_8s << (Particles[InitialIndex  ]->x(0)-Particles[InitialIndex+1]->x(0)-L0(0))/L0(0);
        OF << Util::_8s << (Particles[InitialIndex+2]->x(1)-Particles[InitialIndex+3]->x(1)-L0(1))/L0(1);
        OF << Util::_8s << (Particles[InitialIndex+4]->x(2)-Particles[InitialIndex+5]->x(2)-L0(2))/L0(2);
        OF << std::endl;
    }
}

}; // namespace DEM

#endif // MECHSYS_DEM_DOMAIN_H
