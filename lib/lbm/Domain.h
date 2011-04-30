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

// MechSys
#include <mechsys/lbm/Dem.h>

class Domain
{
public:
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructor
    Domain (LBMethod Method, double nu, iVec3_t Ndim, double dx, double dt);

    //Methods
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
               char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                     ///< Solve the Domain dynamics
    void AddDisk(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt);        ///< Add a disk element



    //Data
    Lattice              Lat;               // Fluid Lattice
    Array <Disk *> Particles;               // Array of Disks
    double              Time;               // Time of the simulation
    double                dt;               // Timestep
    void *          UserData;               // User Data
    size_t           idx_out;               // The discrete time step
};

inline Domain::Domain(LBMethod Method, double nu, iVec3_t Ndim, double dx, double Thedt)
{
    Lat  = Lattice(Method,nu,Ndim,dx,Thedt);
    Time = 0.0;
    dt   = Thedt;
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
    for(size_t i=0;i<Particles.Size();i++) Particles[i]->ImprintDisk(Lat);
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);

        //Assigning a vlaue of zero to the particles forces and torques
        for(size_t i=0;i<Particles.Size();i++)
        {
            Particles[i]->F = 0.0,0.0,0.0;
            Particles[i]->T = 0.0,0.0,0.0;
        }

        //Set Gamma values of the lattice cell to zero
        Lat.SetZeroGamma();
        
        //Connect particles and lattice
        for(size_t i=0;i<Particles.Size();i++) Particles[i]->ImprintDisk(Lat);

        //Move Particles
        for(size_t i=0;i<Particles.Size();i++) Particles[i]->Translate(dt);

        //Move fluid
        Lat.Collide();
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

        Time += dt;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

#endif

