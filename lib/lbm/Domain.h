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
#include <mechsys/lbm/Interacton.h>

using std::set;
using std::pair;
using std::make_pair;

class Domain
{
public:
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructor
    Domain (LBMethod Method, double nu, iVec3_t Ndim, double dx, double dt);

    //Methods
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
               char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                          ///< Solve the Domain dynamics
    void AddDisk(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt);        ///< Add a disk element
    void ResetContacts();
    void ResetDisplacements();
    double  MaxDisplacement();



    //Data
    Lattice                            Lat;         // Fluid Lattice
    Array <Disk *>               Particles;         // Array of Disks
    Array <Interacton *>       Interactons;         // Array of insteractons
    Array <Interacton *>      CInteractons;         // Array of valid interactons
    set<pair<Disk *, Disk *> > Listofpairs;         // List of pair of particles associated per interacton for memory optimization
    double                            Time;         // Time of the simulation
    double                              dt;         // Timestep
    double                           Alpha;         // Verlet distance
    void *                        UserData;         // User Data
    size_t                         idx_out;         // The discrete time step
};

inline Domain::Domain(LBMethod Method, double nu, iVec3_t Ndim, double dx, double Thedt)
{
    Lat  = Lattice(Method,nu,Ndim,dx,Thedt);
    Time = 0.0;
    dt   = Thedt;
    Alpha= 10.0;
}

inline void Domain::ResetDisplacements()
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->X0 = Particles[i]->X;
    }
}

inline double Domain::MaxDisplacement()
{
    double md = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        double mdp = norm(Particles[i]->X - Particles[i]->X0);
        if (mdp>md) md = mdp;
    }
    return md;
}

inline void Domain::ResetContacts()
{
	if (Particles.Size()==0) return;
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        bool pi_has_vf = !Particles[i]->IsFree();
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            bool pj_has_vf = !Particles[j]->IsFree();

            bool close = (norm(Particles[i]->X-Particles[j]->X)<=Particles[i]->R+Particles[j]->R+2*Alpha);
            if ((pi_has_vf && pj_has_vf) || !close) continue;
            
            // checking if the interacton exist for that pair of particles
            set<pair<Disk *, Disk *> >::iterator it = Listofpairs.find(make_pair(Particles[i],Particles[j]));
            if (it != Listofpairs.end())
            {
                continue;
            }
            Listofpairs.insert(make_pair(Particles[i],Particles[j]));
            CInteractons.Push (new Interacton(Particles[i],Particles[j]));
        }
    }

    Interactons.Resize(0);
    for (size_t i=0; i<CInteractons.Size(); i++)
    {
        if(CInteractons[i]->UpdateContacts(Alpha)) Interactons.Push(CInteractons[i]);
    }
}

inline void Domain::AddDisk(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt)
{
    Particles.Push(new Disk(TheTag,TheX,TheV,TheW,Therho,TheR,dt));
}

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * FileKey, bool RenderVideo, size_t Nproc)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    idx_out     = 0;
    //Connect particles and lattice
    for(size_t i=0;i<Particles.Size();i++)
    {
    	Particles[i]->ImprintDisk(Lat);
    }
    ResetContacts();
    ResetDisplacements();

    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);

        //Assigning a vlaue of zero to the particles forces and torques
        for(size_t i=0;i<Particles.Size();i++)
        {
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Ff;
        }

        //Set Gamma values of the lattice cell to zero
        Lat.SetZeroGamma();
        
        //Connect particles and lattice
        for(size_t i=0;i<Particles.Size();i++) Particles[i]->ImprintDisk(Lat);

        //Move Particles
        for(size_t i=0;i<Interactons.Size();i++) Interactons[i]->CalcForce(dt);
        for(size_t i=0;i<Particles.Size()  ;i++) Particles[i]->Translate(dt);

        //Move fluid
        if (fabs(Lat.G)>1.0e-12)
        {
            Lat.ApplyForce();
            Lat.CollideAlt();
        }
        else Lat.Collide();
        Lat.BounceBack();
        Lat.Stream();
        
        if (Time >= tout)
        {
            if (FileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%08d", FileKey, idx_out);
                Lat.WriteVTK (fn.CStr());
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }

        if (MaxDisplacement()>Alpha)
        {
            ResetContacts();
            ResetDisplacements();
        }

        Time += dt;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

#endif

