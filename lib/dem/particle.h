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

#ifndef MECHSYS_DEM_PARTICLE_H
#define MECHSYS_DEM_PARTICLE_H

// Std lib
#include <iostream>
#ifdef USE_THREAD
    #include <thread>
#endif

// MechSys
#include <mechsys/dem/face.h>
#include <mechsys/dem/special_functions.h>
#include <mechsys/util/array.h>
#include <mechsys/numerical/montecarlo.h>
#include <mechsys/mesh/mesh.h>

// MPI
#ifdef USE_MPI
  #include <mpi.h>
#endif

struct ParticleProps
{
    double  Kn;   ///< Normal stiffness
    double  Kt;   ///< Tengential stiffness
    double  Bn;   ///< Spring constant for normal bonding
    double  Bt;   ///< Spring constant for tangential bonding
    double  Bm;   ///< Spring constant for torque bonding
    double  Gn;   ///< Normal viscous coefficient
    double  Gt;   ///< Tangential viscous coefficient
    double  Mu;   ///< Microscopic coefficient of friction
    double  eps;  ///< Maximun strain supported before breaking
    double  Beta; ///< Rolling stiffness coeffcient
    double  Eta;  ///< Plastic moment coefficient
    double  R;    ///< Spheroradius
    double  rho;  ///< Density
    double  V;    ///< Volume
    double  m;    ///< Mass
};

#ifdef USE_MPI

MPI::Datatype MPI_Part_Props_Type;

void BuildParticlePropsDataType (MPI::Datatype & MPIType)
{
    // dummy properties
    ParticleProps p;

    // blocks
    MPI::Datatype types [] = {MPI::DOUBLE};
    int           blklen[] = {15};

    // addresses
    const int nblks = 1;
    MPI::Aint origin, addr[nblks];

    origin  = MPI::Get_address (&p);
    addr[0] = MPI::Get_address (&(p.Kn));

    // displacements
    MPI::Aint disps[nblks];
    for (int i=0; i<nblks; ++i) disps[i] = addr[i] - origin;

    // create type
    MPIType = MPI::Datatype::Create_struct (nblks, blklen, disps, types); 
    MPIType.Commit();
}

#endif


class Particle
{
public:
    // Constructor
    Particle() {}
    Particle(int                         Tag,      ///< Tag of the particle
             Array<Vec3_t>       const & V,        ///< List of vertices
             Array<Array <int> > const & E,        ///< List of edges with connectivity
             Array<Array <int> > const & F,        ///< List of faces with connectivity
             Vec3_t              const & v0,       ///< Initial velocity
             Vec3_t              const & w0,       ///< Initial angular velocity
             double                      R,        ///< Spheroradius
             double                      rho=1.0); ///< Density of the material

    // Alternative constructor
    Particle (int Tag, Mesh::Generic const & M, double R, double rho=1.0);

    // Destructor
    ~Particle ();

    // Methods
    void   Initialize         (size_t i=0, size_t NCalls=5000);                                           ///< Initialize this particle
    void   InitializeVelocity (double dt = 1.0);                                              ///< Initialize this particle
    void   Rotate             (double dt);                                                    ///< Apply rotation on the particle once the total torque is found
    void   Rotate             (Quaternion_t & Q, Vec3_t & V);                                 ///< Apply rotation given by Quaternion Q at point v
    void   Translate          (double dt);                                                    ///< Apply translation once the total force is found
    void   Translate          (Vec3_t & t);                                                   ///< Apply translation by vector t
    void   ResetDisplacements ();                                                             ///< Reset the displacements for the verlet algorithm
    double MaxDisplacement    ();                                                             ///< Maximun displacement for the verlet algorithm
    void   Draw               (std::ostream & os, char const * Color="Blue", bool BPY=false); ///< Draw the particle
    void   FixVeloc           (double vx=0.0, double vy=0.0, double vz=0.0);                  ///< Fix all velocities
    bool   IsFree             () {return !vxf&&!vyf&&!vzf&&!vxf&&!vyf&&!vzf;};                ///< Ask if the particle has any constrain in its movement

// Methods to be used when MPI is used
#ifdef USE_MPI
    void SendParticle               (int n, int MsgID); ///< Send particle to process n
    void ReceiveParticle            (int MsgID);        ///< Receives a particle from any process
    void SendDynamicParticle        (int n, int MsgID); ///< Send just relevant data for force calculation
    void ReceiveDynamicParticle     (int MsgID);        ///< Receive just the dynamic relevant data
    void SendForce                  (int n, int MsgID); ///< Send force
    void ReceiveForce               (int MsgID);        ///< Receive force
#endif

#ifdef USE_THREAD
    std::mutex mtex;                 ///< to protect variables in multithreading
#endif

    // Data -- in MPI data type
    int             Tag;             ///< Tag of the particle
    size_t          Index;           ///< index of the particle in the domain
    bool            PropsReady;      ///< Are the properties calculated ready ?
    bool            IsBroken;        ///< True if the particle has at least one broken bond in cohesive simulations
    bool            vxf, vyf, vzf;   ///< Fixed components of velocity
    bool            wxf, wyf, wzf;   ///< Fixed components of angular velocity
    Vec3_t          x;               ///< Position of the center of mass
    Vec3_t          xb;              ///< Former position for the Verlet algorithm
    Vec3_t          v;               ///< Velocity
    Vec3_t          w;               ///< Angular velocity
    Vec3_t          wb;              ///< Former angular velocity for the leap frog algorithm
    Vec3_t          F;               ///< Force over the particle
    Vec3_t          Ff;              ///< Fixed Force over the particle
    Vec3_t          T;               ///< Torque over the particle
    Vec3_t          Tf;              ///< Fixed Torque over the particle
    Vec3_t          I;               ///< Vector containing the principal components of the inertia tensor
    Quaternion_t    Q;               ///< The quaternion representing the rotation
    Mat3_t          M;               ///< Momment tensor for the calculation of stress
    Mat3_t          B;               ///< Branch tensor for the study of isotropy
    double          Erot;            ///< Rotational energy of the particle
    double          Ekin;            ///< Kinetical energy of the particle
    double          Dmax;            ///< Maximal distance from the center of mass to the surface of the body
    double          Diam;            ///< Diameter of the parallelogram containing the particle
    double          Cn;              ///< Coordination number (number of contacts)

    // Data -- not in MPI data type
    Array<Vec3_t*> Verts;            ///< Vertices

    // Data -- not needed during communications
    ParticleProps       Props;       ///< Properties
    Array<Vec3_t*>      Vertso;      ///< Original postion of the Vertices
    Array<Array <int> > EdgeCon;     ///< Conectivity of Edges 
    Array<Array <int> > FaceCon;     ///< Conectivity of Faces 
    Array<Edge*>        Edges;       ///< Edges
    Array<Face*>        Faces;       ///< Faces
    Array<Torus*>       Tori;        ///< Toroidal features

    // Auxiliar methods
    void   CalcProps (size_t NCalls=5000); ///< Calculate properties: mass, center of mass, and moment of inertia
    bool   IsInside  (Vec3_t & V);         ///< Find whether the point V is inside the particle or not
    double IsInside  (double * V);         ///< Find whether the point V is inside the particle or not
    double MaxX      ();                   ///< Find Maximun X coordinate
    double MaxY      ();                   ///< Find Maximun Y coordinate
    double MaxZ      ();                   ///< Find Maximun Y coordinate
    double MinX      ();                   ///< Find Minimun X coordinate
    double MinY      ();                   ///< Find Minimun Y coordinate
    double MinZ      ();                   ///< Find Minimun Y coordinate

    // Integrants for the calc of properties
    double Vol (double * X); ///< Calculate the volume of the sample at X
    double Xc  (double * X); ///< Calculate the coordinates of the center of mass at X
    double Yc  (double * X); ///< Calculate the coordinates of the center of mass at X
    double Zc  (double * X); ///< Calculate the coordinates of the center of mass at X
    double Ixx (double * X); ///< Calculate the inertia tensor at X
    double Iyy (double * X); ///< Calculate the inertia tensor at X
    double Izz (double * X); ///< Calculate the inertia tensor at X
    double Ixy (double * X); ///< Calculate the inertia tensor at X
    double Ixz (double * X); ///< Calculate the inertia tensor at X
    double Iyz (double * X); ///< Calculate the inertia tensor at X

#ifdef USE_BOOST_PYTHON
    double PyGetFeatures (BPy::list & V, BPy::list & E, BPy::list & F) const
    {
        // vertex-ID map
        typedef std::map<Vec3_t const *,int> VertID_t;
        VertID_t vids;

        // vertices
        for (size_t i=0; i<Verts.Size(); ++i)
        {
            V.append (BPy::make_tuple((*Verts[i])(0),(*Verts[i])(1),(*Verts[i])(2)));
            vids[Verts[i]] = i;
        }

        // edges
        for (size_t i=0; i<Edges.Size(); ++i)
            E.append (BPy::make_tuple(vids[Edges[i]->X0], vids[Edges[i]->X1]));

        // faces
        for (size_t i=0; i<Faces.Size(); ++i)
        {
            /* // list of edges
            BPy::list edges;
            for (size_t j=0; j<Faces[i]->Edges.Size(); ++j)
                edges.append (BPy::make_tuple(vids[Faces[i]->Edges[j]->X0], vids[Faces[i]->Edges[j]->X1]));
            F.append (edges);
            */
            // list of vertices
            BPy::list verts;
            for (size_t j=0; j<Faces[i]->Edges.Size(); ++j)
                verts.append (vids[Faces[i]->Edges[j]->X0]);
            F.append (verts);
        }

        return Props.R;
    }
#endif
};


std::ostream & operator<< (std::ostream & os, Particle const & P)
{
    os << "Tag           = "  << P.Tag        << std::endl;
    os << "Index         = "  << P.Index      << std::endl;
    os << "PropsReady    = "  << P.PropsReady << std::endl;
    os << "IsBroken      = "  << P.IsBroken   << std::endl;
    os << "vxf, vyf, vzf = "  << P.vxf << ", " << P.vyf << ", " << P.vzf << std::endl;
    os << "wxf, wyf, wzf = "  << P.wxf << ", " << P.wyf << ", " << P.wzf << std::endl;
    os << "x             = "  << PrintVector(P.x );
    os << "xb            = "  << PrintVector(P.xb);
    os << "v             = "  << PrintVector(P.v );
    os << "w             = "  << PrintVector(P.w );
    os << "wb            = "  << PrintVector(P.wb);
    os << "F             = "  << PrintVector(P.F );
    os << "Ff            = "  << PrintVector(P.Ff);
    os << "T             = "  << PrintVector(P.T );
    os << "Tf            = "  << PrintVector(P.Tf);
    os << "I             = "  << PrintVector(P.I );
    os << "Q             = "  << P.Q << std::endl;
    os << "M             =\n" << PrintMatrix(P.M );
    os << "B             =\n" << PrintMatrix(P.B );
    os << "Erot          = "  << P.Erot << std::endl;
    os << "Ekin          = "  << P.Ekin << std::endl;
    os << "Dmax          = "  << P.Dmax << std::endl;
    os << "Diam          = "  << P.Diam << std::endl;
    os << "Cn            = "  << P.Cn   << std::endl;
    return os;
}


#ifdef USE_MPI

MPI::Datatype MPI_Particle_Type;

void BuildParticleDataType (MPI::Datatype & MPIType)
{
    // dummy particle
    Particle p;

    // blocks
    MPI::Datatype types [] = {MPI::INT, MPI::UNSIGNED_LONG, MPI::BOOL, MPI::DOUBLE};
    int           blklen[] = {1,1,8,57};

    // addresses
    const int nblks = 4;
    MPI::Aint origin, addr[nblks];

    origin  = MPI::Get_address (&p);
    addr[0] = MPI::Get_address (&(p.Tag));
    addr[1] = MPI::Get_address (&(p.Index));
    addr[2] = MPI::Get_address (&(p.PropsReady));
    addr[3] = MPI::Get_address (p.x.data());

    // displacements
    MPI::Aint disps[nblks];
    for (int i=0; i<nblks; ++i) disps[i] = addr[i] - origin;

    // create type
    MPIType = MPI::Datatype::Create_struct (nblks, blklen, disps, types); 
    MPIType.Commit();
}

#endif


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


// Constructor and destructor

inline Particle::Particle (int TheTag, Array<Vec3_t> const & V, Array<Array <int> > const & E, Array<Array <int> > const & Fa, Vec3_t const & v0, Vec3_t const & w0, double TheR, double TheRho)
    : Tag(TheTag), PropsReady(false), IsBroken(false), v(v0), w(w0)
{
    Props.Kn = 1.0e4;   
    Props.Kt = 5.0e3;   
    Props.Bn = 1.0e4;   
    Props.Bt = 5.0e3;   
    Props.Bm = 5.0e3;
    Props.Gn = 16.0;   
    Props.Gt = 8.0;   
    Props.Mu = 0.4;   
    Props.eps = 0.01;  
    Props.Beta = 0.12; 
    Props.Eta = 1.0;  
    Props.R = TheR;    
    Props.rho = TheRho;  


    vxf = false;
    vyf = false;
    vzf = false;

    wxf = false;
    wyf = false;
    wzf = false;

    F  = 0.0,0.0,0.0;
    Ff = 0.0,0.0,0.0;
    Tf = 0.0,0.0,0.0;

    EdgeCon = E;
    FaceCon = Fa;

    for (size_t i=0; i<V.Size(); i++)
    {
        Verts.Push (new Vec3_t(V[i]));
        Vertso.Push (new Vec3_t(V[i]));
    }
    for (size_t i=0; i<Fa.Size(); i++)
    {
        Array<Vec3_t*> verts(Fa[i].Size());
        for (size_t j=0; j<Fa[i].Size(); ++j) verts[j] = Verts[Fa[i][j]];
        Faces.Push (new Face(verts));
    }
    for (size_t i=0; i<E.Size(); i++) Edges.Push (new Edge((*Verts[E[i][0]]), (*Verts[E[i][1]])));
}

inline Particle::Particle (int TheTag, Mesh::Generic const & M, double TheR, double TheRho)
    : Tag(TheTag), PropsReady(false), IsBroken(false), v(Vec3_t(0.0,0.0,0.0)), w(Vec3_t(0.0,0.0,0.0))
{
    Props.Kn = 1.0e4;   
    Props.Kt = 5.0e3;   
    Props.Bn = 1.0e4;   
    Props.Bt = 5.0e3;   
    Props.Bm = 5.0e3;
    Props.Gn = 16.0;   
    Props.Gt = 8.0;   
    Props.Mu = 0.4;   
    Props.eps = 0.01;  
    Props.Beta = 0.12; 
    Props.Eta = 1.0;  
    Props.R = TheR;    
    Props.rho = TheRho;  

    Ff = 0.0,0.0,0.0;
    Tf = 0.0,0.0,0.0;

    // check if mesh is Shell
    if (!M.IsShell) throw new Fatal("Particle::Particle: Mesh must be of Shell type");

    // vertices
    size_t nv = M.Verts.Size();
    for (size_t i=0; i<nv; ++i)
    {
        Verts .Push (new Vec3_t(M.Verts[i]->C(0), M.Verts[i]->C(1), M.Verts[i]->C(2)));
        Vertso.Push (new Vec3_t(M.Verts[i]->C(0), M.Verts[i]->C(1), M.Verts[i]->C(2)));
    }

    // edges and faces
    typedef std::map<std::pair<int,int>,Edge*> Key2Edge_t;
    Key2Edge_t key2edge;        // map edge pair (v0,v1) to Edge* in Edges
    size_t nf = M.Cells.Size(); // number of faces: each cell is one face
    for (size_t i=0; i<nf; ++i)
    {
        // number of vertices per face
        size_t nvf = M.Cells[i]->V.Size();

        // edges
        size_t v0, v1;
        std::pair<int,int> keya, keyb;
        for (size_t j=0; j<nvf; ++j)
        {
            v0   = M.Cells[i]->V[j]->ID;
            v1   = M.Cells[i]->V[(j+1)%nvf]->ID;
            keya = std::make_pair(v0,v1);
            keyb = std::make_pair(v1,v0);
            Key2Edge_t::const_iterator ita = key2edge.find(keya);
            Key2Edge_t::const_iterator itb = key2edge.find(keyb);
            if (ita==key2edge.end() && itb==key2edge.end()) // new edge
            {
                Edges.Push (new Edge((*Verts[v0]), (*Verts[v1])));
                key2edge[keya] = Edges[Edges.Size()-1];
                EdgeCon.Push (Array<int>((int)v0,(int)v1)); // TODO: we may remove this
            }
        }

        // faces
        Array<Vec3_t*> verts(nvf);
        FaceCon.Push (Array<int>()); // TODO: we may remove this
        for (size_t j=0; j<nvf; ++j)
        {
            v0 = M.Cells[i]->V[j]->ID;
            verts[j] = Verts[v0];
            FaceCon[FaceCon.Size()-1].Push (v0); // TODO: we may remove this
        }
        Faces.Push (new Face(verts));
    }
}

inline Particle::~Particle()
{
    for (size_t i=0; i<Verts.Size(); ++i) delete Verts[i];
    for (size_t i=0; i<Vertso.Size(); ++i) delete Vertso[i];
    for (size_t i=0; i<Edges.Size(); ++i) delete Edges[i];
    for (size_t i=0; i<Faces.Size(); ++i) delete Faces[i];
}

// Methods

inline void Particle::Initialize (size_t i, size_t NCalls)
{
    // calc properties
    if (!PropsReady)
    {
        Index = i;
        CalcProps (NCalls);
    }
}

inline void Particle::InitializeVelocity (double dt)
{
    // initialize the particle for the Verlet algorithm
    xb = x-v*dt;
    wb = w;
}

inline void Particle::Rotate (double dt)
{
    double q0,q1,q2,q3,wx,wy,wz;
    q0 = 0.5*Q(0);
    q1 = 0.5*Q(1);
    q2 = 0.5*Q(2);
    q3 = 0.5*Q(3);

    if (wxf) T(0) = 0.0;
    if (wyf) T(1) = 0.0;
    if (wzf) T(2) = 0.0;

    Vec3_t Td;
    Td(0)=(T(0)+(I(1)-I(2))*wb(1)*wb(2))/I(0);
    Td(1)=(T(1)+(I(2)-I(0))*wb(0)*wb(2))/I(1);
    Td(2)=(T(2)+(I(0)-I(1))*wb(1)*wb(0))/I(2);
    w = wb+0.5*dt*Td;
    wx = w(0);
    wy = w(1);
    wz = w(2);
    Quaternion_t dq(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz),qm;

    wb  = wb+Td*dt;
    qm  = Q+dq*(0.5*dt);
    q0  = 0.5*qm(0);
    q1  = 0.5*qm(1);
    q2  = 0.5*qm(2);
    q3  = 0.5*qm(3);
    wx  = wb(0);
    wy  = wb(1);
    wz  = wb(2);
    dq  = Quaternion_t(-(q1*wx+q2*wy+q3*wz),q0*wx-q3*wy+q2*wz,q3*wx+q0*wy-q1*wz,-q2*wx+q1*wy+q0*wz);
    Quaternion_t Qd = (qm+dq*0.5*dt),temp;
    Conjugate (Q,temp);
    Rotate    (temp,x);
    Q  = Qd/norm(Qd);
    Rotate (Q,x);
    Erot=0.5*(I(0)*wx*wx+I(1)*wy*wy+I(2)*wz*wz);

}

inline void Particle::Rotate (Quaternion_t & Q,Vec3_t & V)
{
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size();
    for (size_t i = 0; i < nv; i++)
    {
        Vec3_t xt = *Verts[i]-V;
        Rotation(xt,Q,*Verts[i]);
        *Verts[i] += V;
    }

    for (size_t i = 0; i < ne; i++)
    {
        Edges[i]->UpdatedL();
    }

    for (size_t i = 0; i < nf; i++)
    {
        Faces[i]->UpdatedL();
    }
}

inline void Particle::Translate (double dt)
{
    if (vxf) F(0) = 0.0;
    if (vyf) F(1) = 0.0;
    if (vzf) F(2) = 0.0;
    if(isnan(norm(F))) 
    {
        throw new Fatal("Particle::Translate: The force is not a number %d(%d), try reducing the time step",Index,Tag);
    }
    Vec3_t temp,xa;
    xa    = 2*x - xb + F*(dt*dt/Props.m);
    temp  = xa - x;
    v    = 0.5*(xa - xb)/dt;
    xb   = x;
    x    = xa;
    Ekin = 0.5*Props.m*dot(v,v);

    size_t nv = Verts.Size();
    for (size_t i = 0; i < nv; i++)
    {
        *Verts[i] += temp;
    }
}

inline void Particle::Translate (Vec3_t & V)
{
    size_t nv = Verts.Size();
    for (size_t i = 0; i < nv; i++)
    {
        *Verts[i] += V;
    }
    x += V;
    xb += V;
}

inline void Particle::ResetDisplacements ()
{
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        (*Vertso[i]) = (*Verts[i]);
    }
}

inline double Particle::MaxDisplacement ()
{
    double md = 0.0;
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        double mpd = Distance((*Vertso[i]),(*Verts[i]));
        if (mpd>md) md = mpd;
    }
    return md;
}

inline void Particle::Draw (std::ostream & os, char const * Color, bool BPY)
{
    for (size_t i=0; i<Verts.Size(); ++i)
    {
        if (BPY) BPYDrawVert((*Verts[i]),os,Props.R);
        else POVDrawVert((*Verts[i]),os,Props.R,Color); 
    }
    for (size_t i=0; i<Edges.Size(); ++i)
    {
        Edges[i]->Draw(os,Props.R,Color,BPY);
    }
    for (size_t i=0; i<Faces.Size(); ++i)
    {
        Faces[i]->Draw(os,Props.R,Color,BPY);
    }
    for (size_t i=0; i<Tori.Size(); ++i)
    {
        Tori[i]->Draw(os,Props.R,Color,BPY);
    }
}

inline void Particle::FixVeloc (double vx, double vy, double vz)
{
    v   = vx, vy, vz;
    vxf = true; vyf = true; vzf = true; 
    wxf = true; wyf = true; wzf = true;
}

// MPI exclusive methods
#ifdef USE_MPI

inline void Particle::SendParticle(int n, int MsgID)
{
    // Sending properties
    MPI::COMM_WORLD.Send (&this->Props,  /*number of props.*/1, MPI_Part_Props_Type,  /*destination*/n, MsgID);
    MPI::COMM_WORLD.Send (this,  /*number of part.*/1, MPI_Particle_Type,  /*destination*/n, MsgID);

    // Sending vertices
    size_t verts_size = Verts.Size();
    MPI::COMM_WORLD.Send (&verts_size, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/n, MsgID);
    for (size_t i=0; i<verts_size; ++i)
    {
        MPI::COMM_WORLD.Send (Verts[i]->data(), /*number*/3, MPI::DOUBLE, /*destination*/n, MsgID);
    }

    // Sending edges
    size_t edge_size = Edges.Size();
    MPI::COMM_WORLD.Send (&edge_size, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/n, MsgID);
    for (size_t i=0; i<edge_size; ++i)
    {
        MPI::COMM_WORLD.Send (&EdgeCon[i][0], /*number*/2, MPI::INT, /*destination*/n, MsgID);
    }

    // Sending faces
    size_t face_size = Faces.Size();
    MPI::COMM_WORLD.Send (&face_size, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/n, MsgID);
    for (size_t i=0; i<face_size; ++i)
    {
        size_t verts_face = FaceCon[i].Size();
        MPI::COMM_WORLD.Send (&verts_face, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/n, MsgID);
        MPI::COMM_WORLD.Send (&FaceCon[i][0], /*number*/verts_face, MPI::INT, /*destination*/n, MsgID);
    }
}

inline void Particle::ReceiveParticle(int MsgID)
{
    // Receiving properties
    MPI::COMM_WORLD.Recv (&this->Props, /*number of part.*/1, MPI_Part_Props_Type, MPI::ANY_SOURCE, MsgID);
    MPI::COMM_WORLD.Recv (this, /*number of part.*/1, MPI_Particle_Type, MPI::ANY_SOURCE, MsgID);

    // Receiving vertices
    size_t verts_size;
    MPI::COMM_WORLD.Recv (&verts_size, /*number*/1, MPI::UNSIGNED_LONG, MPI::ANY_SOURCE, MsgID);
    for (size_t i=0; i<verts_size; ++i)
    {
        Verts.Push (new Vec3_t(0,0,0));
        MPI::COMM_WORLD.Recv (Verts[i]->data(), /*number*/3, MPI::DOUBLE, MPI::ANY_SOURCE, MsgID);
    }

    // Receiving edges
    size_t edge_size;
    MPI::COMM_WORLD.Recv (&edge_size, /*number*/1, MPI::UNSIGNED_LONG, MPI::ANY_SOURCE, MsgID);
    for (size_t i=0; i<edge_size; ++i)
    {
        Array<int> D(2);
        EdgeCon.Push(D);
        MPI::COMM_WORLD.Recv (&EdgeCon[i][0], /*number*/2, MPI::INT, MPI::ANY_SOURCE, MsgID);
    }
    for (size_t i=0; i<EdgeCon.Size(); i++) Edges.Push (new Edge((*Verts[EdgeCon[i][0]]), (*Verts[EdgeCon[i][1]])));

    // Receiving faces
    size_t face_size;
    MPI::COMM_WORLD.Recv (&face_size, /*number*/1, MPI::UNSIGNED_LONG, MPI::ANY_SOURCE, MsgID);
    for (size_t i=0; i<face_size; ++i)
    {
        size_t verts_face;
        MPI::COMM_WORLD.Recv (&verts_face, /*number*/1, MPI::UNSIGNED_LONG, MPI::ANY_SOURCE, MsgID);
        Array<int> D(verts_face);
        FaceCon.Push(D);
        MPI::COMM_WORLD.Recv (&FaceCon[i][0], /*number*/verts_face, MPI::INT, MPI::ANY_SOURCE, MsgID);
    }
    for (size_t i=0; i<FaceCon.Size(); i++)
    {
        Array<Vec3_t*> verts(FaceCon[i].Size());
        for (size_t j=0; j<FaceCon[i].Size(); ++j) verts[j] = Verts[FaceCon[i][j]];
        Faces.Push (new Face(verts));
    }
}

inline void Particle::SendDynamicParticle(int n, int MsgID)
{
    // Sending vertices
    size_t verts_size = Verts.Size();
    MPI::COMM_WORLD.Send (&verts_size, /*number*/1, MPI::UNSIGNED_LONG, /*destination*/n, MsgID);
    for (size_t i=0; i<verts_size; ++i)
    {
        MPI::COMM_WORLD.Send (Verts[i]->data(), /*number*/3, MPI::DOUBLE, /*destination*/n, MsgID);
    }

    // Sending vectorial properties
    MPI::COMM_WORLD.Send (x.data(), /*number*/3, MPI::DOUBLE, /*destination*/n, MsgID);
    MPI::COMM_WORLD.Send (v.data(), /*number*/3, MPI::DOUBLE, /*destination*/n, MsgID);
    MPI::COMM_WORLD.Send (w.data(), /*number*/3, MPI::DOUBLE, /*destination*/n, MsgID);
    MPI::COMM_WORLD.Send (Q.data(), /*number*/4, MPI::DOUBLE, /*destination*/n, MsgID);

    //std::cout << Index << " " << x << std::endl;
}

inline void Particle::ReceiveDynamicParticle(int MsgID)
{
    // Sending vertices
    size_t verts_size;
    MPI::COMM_WORLD.Recv (&verts_size, /*number*/1, MPI::UNSIGNED_LONG, MPI::ANY_SOURCE, MsgID);
    for (size_t i=0; i<verts_size; ++i)
    {
        MPI::COMM_WORLD.Recv (Verts[i]->data(), /*number*/3, MPI::DOUBLE, MPI::ANY_SOURCE, MsgID);
    }

    // Sending vectorial properties
    MPI::COMM_WORLD.Recv (x.data(), /*number*/3, MPI::DOUBLE, MPI::ANY_SOURCE, MsgID);
    MPI::COMM_WORLD.Recv (v.data(), /*number*/3, MPI::DOUBLE, MPI::ANY_SOURCE, MsgID);
    MPI::COMM_WORLD.Recv (w.data(), /*number*/3, MPI::DOUBLE, MPI::ANY_SOURCE, MsgID);
    MPI::COMM_WORLD.Recv (Q.data(), /*number*/4, MPI::DOUBLE, MPI::ANY_SOURCE, MsgID);


    //std::cout << Index << " " << v << std::endl;
}

inline void Particle::SendForce(int n, int MsgID)
{
    MPI::COMM_WORLD.Send (F.data(), /*number*/3, MPI::DOUBLE, /*destination*/n, MsgID);
    MPI::COMM_WORLD.Send (T.data(), /*number*/3, MPI::DOUBLE, /*destination*/n, MsgID);
}

inline void Particle::ReceiveForce(int MsgID)
{
    // Temporal force and torque
    Vec3_t Ft,Tt;
    MPI::COMM_WORLD.Recv (Ft.data(), /*number*/3, MPI::DOUBLE, MPI::ANY_SOURCE, MsgID);
    MPI::COMM_WORLD.Recv (Tt.data(), /*number*/3, MPI::DOUBLE, MPI::ANY_SOURCE, MsgID);

    //if(norm(Ft)>0.0) std::cout << Index << " " << Ft << std::endl;

    //Add this quantities to my force and torque
    F += Ft;
    T += Tt;
}
#endif

// Auxiliar methods

inline void Particle::CalcProps (size_t NCalls)
{
    if (Verts.Size()==1 && Edges.Size()==0 && Faces.Size()==0)
    {
        Props.V = (4./3.)*M_PI*Props.R*Props.R*Props.R;
        I = Vec3_t((8./15.)*M_PI*pow(Props.R,5.),(8./15.)*M_PI*pow(Props.R,5.),(8./15.)*M_PI*pow(Props.R,5.));
        x = *Verts[0];
        Q = 1,0,0,0;
        Props.m = Props.rho*Props.V;
        I*= Props.rho;
        Ekin = 0.5*Props.m*dot(v,v);
        Erot = 0.5*(I(0)*w(0)*w(0)+I(1)*w(1)*w(1)+I(2)*w(2)*w(2));
        Dmax = Props.R;
    }
    else 
    {
        Mat3_t It;
        double Xi[3] = { MinX() , MinY() , MinZ() };
        double Xs[3] = { MaxX() , MaxY() , MaxZ() };
        Numerical::MonteCarlo<Particle> MC(this, Numerical::VEGAS, NCalls);
        Props.V = MC.Integrate(&Particle::Vol, Xi,Xs);
        x(0)    = MC.Integrate(&Particle::Xc,  Xi,Xs)/Props.V;
        x(1)    = MC.Integrate(&Particle::Yc,  Xi,Xs)/Props.V;
        x(2)    = MC.Integrate(&Particle::Zc,  Xi,Xs)/Props.V;
        It(0,0) = MC.Integrate(&Particle::Ixx, Xi,Xs);
        It(1,1) = MC.Integrate(&Particle::Iyy, Xi,Xs);
        It(2,2) = MC.Integrate(&Particle::Izz, Xi,Xs);
        It(1,0) = MC.Integrate(&Particle::Ixy, Xi,Xs);
        It(2,0) = MC.Integrate(&Particle::Ixz, Xi,Xs);
        It(2,1) = MC.Integrate(&Particle::Iyz, Xi,Xs);
        It(0,1) = It(1,0);
        It(0,2) = It(2,0);
        It(1,2) = It(2,1);

        Vec3_t xp,yp,zp;
        Eig(It,I,xp,yp,zp);
        std::cout << x(0) << std::endl;
        I *= Props.rho;
        Q(0) = 0.5*sqrt(1+xp(0)+yp(1)+zp(2));
        Q(1) = (yp(2)-zp(1))/(4*Q(0));
        Q(2) = (zp(0)-xp(2))/(4*Q(0));
        Q(3) = (xp(1)-yp(0))/(4*Q(0));
        Q = Q/norm(Q);
        Rotation(w,Q,wb);
        w = wb;
        Props.m = Props.rho*Props.V;
        Ekin = 0.5*Props.m*dot(v,v);
        Erot = 0.5*(I(0)*w(0)*w(0)+I(1)*w(1)*w(1)+I(2)*w(2)*w(2));
        Dmax = Distance(x,(*Verts[0]))+Props.R;
        for (size_t i=1; i<Verts.Size(); ++i)
        {
            if (Distance(x,(*Verts[i]))+Props.R > Dmax) Dmax = Distance(x,(*Verts[i]))+Props.R;
        }
    }
    PropsReady = true;
}

inline bool Particle::IsInside (Vec3_t & V)
{
    bool inside = false;
    size_t nv = Verts.Size(),ne = Edges.Size(),nf = Faces.Size();
    if (Distance(x,V)>Dmax) return inside;
    for (size_t i = 0; i < nv; i++)
    {
        if (Distance(V,*Verts[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }

    for (size_t i = 0; i < ne; i++)
    {
        if (Distance(V,*Edges[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }
    for (size_t i = 0; i < nf; i++)
    {
        if (Distance(V,*Faces[i]) < Props.R) {
            inside = true;
            return inside;
        }
    }
    if (nf>3)
    {
        size_t k = 0;
        double Mindistance = Distance(V,*Faces[k]);
        for (size_t i = 1; i < nf; i++)
        {
            if (Distance(V,*Faces[i])<Mindistance) 
            {
                k = i;
                Mindistance = Distance(V,*Faces[k]);
            }
        }
        Vec3_t ct(0,0,0);
        Faces[k]->Centroid(ct);
        Vec3_t pro = V - ct;
        Vec3_t nor;
        Faces[k]->Normal(nor);
        if (dot(pro,nor)<0) inside =true;
    }


    return inside;
}

inline double Particle::IsInside (double * V)
{
    Vec3_t p(V);
    return static_cast<double>(IsInside(p));
}

inline double Particle::MaxX ()
{
    double result = (*Verts[0])(0)+Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(0)+Props.R > result) result = (*Verts[i])(0)+Props.R; 
    }
    return result;
}

inline double Particle::MaxY ()
{
    double result = (*Verts[0])(1)+Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(1)+Props.R > result) result = (*Verts[i])(1)+Props.R; 
    }
    return result;
}

inline double Particle::MaxZ ()
{
    double result = (*Verts[0])(2)+Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(2)+Props.R > result) result = (*Verts[i])(2)+Props.R; 
    }
    return result;
}

inline double Particle::MinX ()
{
    double result = (*Verts[0])(0)-Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(0)-Props.R < result) result = (*Verts[i])(0)-Props.R; 
    }
    return result;
}

inline double Particle::MinY ()
{
    double result = (*Verts[0])(1)-Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(1)-Props.R < result) result = (*Verts[i])(1)-Props.R; 
    }
    return result;
}

inline double Particle::MinZ ()
{
    double result = (*Verts[0])(2)-Props.R;
    for (size_t i = 1; i < Verts.Size(); i++)
    {
        if ((*Verts[i])(2)-Props.R < result) result = (*Verts[i])(2)-Props.R; 
    }
    return result;
}

// Integrants for the calc of properties

inline double Particle::Vol (double * X)
{
    return IsInside(X);
}

inline double Particle::Xc (double * X)
{
    return X[0]*IsInside(X);
}

inline double Particle::Yc (double * X)
{
    return X[1]*IsInside(X);
}

inline double Particle::Zc (double * X)
{
    return X[2]*IsInside(X);
}

inline double Particle::Ixx (double * X)
{
    return ((X[1]-x(1))*(X[1]-x(1))+(X[2]-x(2))*(X[2]-x(2)))*IsInside(X);
}

inline double Particle::Iyy (double * X)
{
    return ((X[0]-x(0))*(X[0]-x(0))+(X[2]-x(2))*(X[2]-x(2)))*IsInside(X);
}

inline double Particle::Izz (double * X)
{
    return ((X[0]-x(0))*(X[0]-x(0))+(X[1]-x(1))*(X[1]-x(1)))*IsInside(X);
}

inline double Particle::Ixy (double * X)
{
    return -(X[0]-x(0))*(X[1]-x(1))*IsInside(X);
}

inline double Particle::Ixz (double * X)
{
    return -(X[0]-x(0))*(X[2]-x(2))*IsInside(X);
}

inline double Particle::Iyz (double * X)
{
    return -(X[1]-x(1))*(X[2]-x(2))*IsInside(X);
}

#endif // MECHSYS_DEM_PARTICLE_H
