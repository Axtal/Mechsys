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
#include <ctime>    // for std::clock

// Voro++
#include "src/voro++.cc"

// MechSys
#include "dem/interacton.h"
#include "util/array.h"
#include "util/util.h"
#include "mesh/mesh.h"

class Domain
{
public:
    // Destructor
    ~Domain();

    // Particle generation
    void GenSpheres (size_t N,     ///< Number of spheres
                     double Xmin,  ///< Left boundary
                     double Xmax,  ///< Right boundary
                     double Ymin,  ///< Back boundary
                     double Ymax,  ///< Front boundary
                     double Zmin,  ///< Bottom boundary
                     double Zmax,  ///< Top boundary
                     double rho,   ///< Density of the material
                     double Rmin); ///< Minimun radius in units of the maximun radius

    void GenBox (double Xmin,       ///< Left boundary
                 double Xmax,       ///< Right boundary
                 double Ymin,       ///< Back boundary
                 double Ymax,       ///< Front boundary
                 double Zmin,       ///< Bottom boundary
                 double Zmax,       ///< Top boundary
                 double Thickness); ///< Thickness of the wall, cannot be zero

    void GenFromMesh (Mesh::Generic const & M, double R, double rho=1.0);
    void GenFromVoro (container & VC, double R, double rho=1.0);

    // Single particle addition
    void AddVoroCell (voronoicell & VC, double R, double rho=1.0);
    void AddTetra    (Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a tetrahedron at position X with spheroradius R, side of length L and density rho
    void AddRice     (Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a rice at position X with spheroradius R, side of length L and density rho
    void AddCube     (Vec3_t const & X, double R, double L, double rho, double Angle=0, Vec3_t * Axis=NULL); ///< Add a cube at position X with spheroradius R, side of length L and density rho

    // Methods
    void CopyParticle (const Particle & P);  ///< Create a new particle as a copy of particle P, it should be translated to another position.
    void Initialize   (double dt);           ///< Set the particles to a initial state and asign the possible insteractions
    void OneStep      (double dt);           ///< One simualtion step
    void Solve        (double t0, double tf, double dt, double dtOut, char const * FileKey); // Run simulation

    // Auxiliar methods
    void   LinearMomentum  (Vec3_t & L); ///< Return total momentum of the system
    void   AngularMomentum (Vec3_t & L); ///< Return total angular momentum of the system
    double TotalEnergy     ();           ///< Return total energy of the system

    // Drawing
    void WritePOV (std::ostream & os,char const *Color);
    void WriteBPY (std::ostream & os);

    // Data
    Array<Particle*>   Particles;   ///< All particles in domain
    Array<Interacton*> Interactons; ///< All interactons in domain
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Domain::~Domain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    for (size_t i=0; i<Interactons.Size(); ++i) if (Interactons[i]!=NULL) delete Interactons[i];
}

// Particle generation
    
inline void Domain::GenSpheres (size_t N, double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, double rho, double Rmin)
{
    double Lx   = Xmax-Xmin;
    double Ly   = Ymax-Ymin;
    double Lz   = Zmax-Zmin;
    double Rmax = pow(Lx*Ly*Lz/(8*N),1./3.);
    size_t nx   = Lx/(2*Rmax);
    size_t ny   = Ly/(2*Rmax);
    Array <Array <int> > Empty;
    for (size_t i=0; i<N; i++)
    {
        Array<Vec3_t> V(1);
        V[0] = Xmin + Rmax+2*Rmax*(i%nx), Ymin + Rmax+2*Rmax*((i/nx)%ny), Zmin + Rmax+2*Rmax*(i/(nx*ny));
        double R = Rmax*Rmin + (1.*rand())/RAND_MAX*Rmax*(1-Rmin);
        Particles.Push (new Particle(V,Empty,Empty, OrthoSys::O,OrthoSys::O, R,rho));
    }
}

inline void Domain::GenBox (double Xmin, double Xmax, double Ymin, double Ymax, double Zmin, double Zmax, double Thickness)
{
}

inline void Domain::GenFromMesh (Mesh::Generic const & M, double R, double rho)
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

        size_t nfaces = Mesh::NVertsToNFaces[nverts];
        size_t nvperf = Mesh::NVertsToNVertsPerFace[nverts];
        Array<Array <int> > F(nfaces);
        for (size_t j=0; j<nfaces; ++j)
        {
            for (size_t k=0; k<nvperf; ++k)
            {
                // TODO: check if face is planar or not
                F[j].Push(Mesh::NVertsToFace[nverts][j][k]);
            }
        }

        // add particle
        Particles.Push (new Particle(V,E,F,OrthoSys::O,OrthoSys::O,R,rho,true));

        // Calculate properties
        Particles[i]->CalcMassProperties();
    }

    // info
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    std::cout << "[1;32m    Number of particles   = " << Particles.Size() << "[0m\n";
}

inline void Domain::GenFromVoro (container & VC,double R,double rho)
{
	fpoint x,y,z,px,py,pz;
    container *cp = & VC;
	voropp_loop l1(cp);
	int q,s;
	voronoicell c;
	s=l1.init(VC.ax,VC.bx,VC.ay,VC.by,VC.az,VC.bz,px,py,pz);
	do {
		for(q=0;q<VC.co[s];q++) {
			x=VC.p[s][VC.sz*q]+px;y=VC.p[s][VC.sz*q+1]+py;z=VC.p[s][VC.sz*q+2]+pz;
			if(x>VC.ax&&x<VC.bx&&y>VC.ay&&y<VC.by&&z>VC.az&&z<VC.bz) {
				if(VC.compute_cell(c,l1.ip,l1.jp,l1.kp,s,q,x,y,z)) AddVoroCell(c,R,rho);
			}
		}
	} while((s=l1.inc(px,py,pz))!=-1);
}

// Single particle addition
    
inline void Domain::AddVoroCell (voronoicell & VC, double R, double rho)
{

    Array<Vec3_t> V(VC.p);
    Array<Array <int> > E;
    Array<int> Eaux(2);
    for(size_t i=0;i<VC.p;i++) 
    {
        V[i] = Vec3_t(VC.pts[3*i],VC.pts[3*i+1],VC.pts[3*i+2]);
        //ux=0.5*VC.pts[3*i]; 
        //uy=0.5*VC.pts[3*i+1]; 
        //uz=0.5*VC.pts[3*i+2];

        for(size_t j=0;j<VC.nu[i];j++) 
        {
            size_t k=VC.ed[i][j];
            if (VC.ed[i][j]<i) 
            {
                Eaux[0] = i;
                Eaux[1] = k;
                E.Push(Eaux);
            }
            //if (ed[i][j]<i) os << ux << " " << uy << " " << uz << "\n" << x+0.5*pts[3*k] << " " << y+0.5*pts[3*k+1] << " " << z+0.5*pts[3*k+2] << "\n\n\n";
        }
    }
    Array<Array <int> > F;
    Array<int> Faux;
    bool later=false;
    for(int i=0;i<VC.p;i++) 
    {
        for(int j=0;j<VC.nu[i];j++) 
        {
            int k=VC.ed[i][j];
            if (k>=0) 
            {
                //if(later) os << " ";
                //else later=true;
                //os << "(" << i;
                Faux.Push(i);
                VC.ed[i][j]=-1-k;
				int l=VC.cycle_up(VC.ed[i][VC.nu[i]+j],k);
                do 
                {
                    //os << "," << k;
                    Faux.Push(k);
                    int m=VC.ed[k][l];
                    VC.ed[k][l]=-1-m;
					l=VC.cycle_up(VC.ed[k][VC.nu[k]+l],m);
                    k=m;
                } while (k!=i);
                F.Push(Faux);
                Faux.Clear();
                //os << ")";
            }
        }
    }
    VC.reset_edges();



    // add particle
    Particles.Push (new Particle(V,E,F,OrthoSys::O,OrthoSys::O,R,rho,true));
}

inline void Domain::AddTetra (Vec3_t const & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
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
    Particles.Push (new Particle(V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    delete Axis;
}

inline void Domain::AddRice (const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
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
    Particles.Push (new Particle(V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    delete Axis;
}

inline void Domain::AddCube (const Vec3_t & X, double R, double L, double rho, double Angle, Vec3_t * Axis)
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
    Particles.Push (new Particle(V,E,F,OrthoSys::O,OrthoSys::O,R,rho));
    delete Axis;
}

// Methods

inline void Domain::CopyParticle (const Particle & P)
{

}
    
inline void Domain::Initialize (double dt)
{
    for(size_t i = 0;i < Particles.Size();i++)
    {
        //Particles[i]->CalcMassProperties(5000);
        Particles[i]->Start(dt);
    }

    for(size_t i = 0;i < Particles.Size()-1;i++)
    {
        for(size_t j = i+1;j < Particles.Size();j++)
        {
            Interactons.Push(new Interacton(Particles[i],Particles[j]));
        }
    }
}

inline void Domain::OneStep (double dt)
{
    for(size_t i = 0;i < Particles.Size();i++)
    {
        Particles[i]->StartForce();
    }
    for(size_t i = 0;i < Interactons.Size();i++)
    {
        Interactons[i]->CalcForce();
    }
    for(size_t i = 0;i < Particles.Size();i++)
    {
        Particles[i]->DynamicRotation(dt);
        Particles[i]->DynamicTranslation(dt);
    }
}

inline void Domain::Solve (double t0, double tf, double dt, double dtOut, char const * FileKey)
{
    // info
    double start = std::clock();
    std::cout << "[1;33m\n--- Solving ----------------------------------------------------[0m\n";

    size_t I=0;
    double tout = t0 + dtOut;
    for (double t=t0; t<tf; t+=dt)
    {
        OneStep(dt);
        if(t>=tout)
        {
            String fn;  
            fn.Printf("%s_%08d.pov",FileKey,I);
            std::ofstream of(fn.CStr());
            POVHeader(of);
            Vec3_t p(0,10,0);
            POVSetCam(of,p,OrthoSys::O);
            WritePOV(of,"Blue");
            of.close();
            tout+=dtOut;
            I++;
        }
    }

    // info
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
}

// Auxiliar methods
    
inline void Domain::LinearMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for(size_t i = 0;i < Particles.Size();i++)
    {
        L+=Particles[i]->m*Particles[i]->v;
    }
}

inline void Domain::AngularMomentum (Vec3_t & L)
{
    L = 0.,0.,0.;
    for(size_t i = 0;i < Particles.Size();i++)
    {
        Vec3_t t1,t2;
        t1 = Particles[i]->I(0)*Particles[i]->w(0),Particles[i]->I(1)*Particles[i]->w(1),Particles[i]->I(2)*Particles[i]->w(2);
        Rotation(t1,Particles[i]->Q,t2);
        L += Particles[i]->m*cross(Particles[i]->x,Particles[i]->v)+t2;
    }
    
}

inline double Domain::TotalEnergy ()
{
    double E = 0;
    for(size_t i = 0;i < Particles.Size();i++)
    {
        E += Particles[i]->Ekin+Particles[i]->Erot;
    }

    for(size_t i = 0;i < Interactons.Size();i++)
    {
        E += Interactons[i]->Epot;
    }
    return E;
}

// Drawing
    
inline void Domain::WritePOV (std::ostream & os, char const * Color)
{
    for(size_t i = 0;i < Particles.Size();i++) Particles[i]->Draw(os,Color);
}

inline void Domain::WriteBPY (std::ostream & os)
{   
    for(size_t i = 0;i < Particles.Size();i++) Particles[i]->Draw(os,"",true);
}

#endif // MECHSYS_DEM_DOMAIN_H
