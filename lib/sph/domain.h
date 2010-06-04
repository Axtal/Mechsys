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

#ifndef MECHSYS_SPH_DOMAIN_H
#define MECHSYS_SPH_DOMAIN_H

// Mechsys
#include <mechsys/sph/interacton.h>
#include <mechsys/dem/graph.h>

class SPHDomain
{
public:

    // Constructor
    SPHDomain(iVec3_t n, Vec3_t Xmin, Vec3_t Xmax);  ///< Constructor with a vector containing the number of divisions per length, and Xmin Xmax defining the limits of the rectangular domain to be plotted

    // Destructor
    ~SPHDomain ();

    // Methods
    void AddBox(Vec3_t const & x, size_t nx, size_t ny, size_t nz, double h, double s,  double rho0, bool Fixed);      ///< Add a box of SPHparticles
    void AddRandomBox(Vec3_t const &V, double Lx, double Ly, double Lz,
                                       size_t nx, size_t ny, size_t nz, double rho0, double R, size_t RandomSeed=100); ///< Add box of random positioned particles
    void StartAcceleration (Vec3_t const & a = Vec3_t(0.0,0.0,0.0));                                                   ///< Add a fixed acceleration
    void ComputeAcceleration (double dt);                                                                              ///< Compute the accleration due to the other particles
    void Move                (double dt);                                                                              ///< Compute the accleration due to the other particles
    void WriteBPY (char const * FileKey);                                                                              ///< Draw the entire domain in a POV file
    void WritePOV (char const * FileKey);                                                                              ///< Draw the entire domain in a blender file
    void WriteVTK (char const * FileKey);                                                                              ///< Draw the entire domain in a VTK file
    void ResetInteractons();                                                                                           ///< Reset the interacton array
    void ResetDisplacements();                                                                                         ///< Reset the particles displacement
    void ResetContacts();                                                                                              ///< Reset the possible interactons
    double MaxDisplacement();                                                                                          ///< Find max displacement of particles
    void Solve    (double tf, double dt, double dtOut, char const * TheFileKey, bool RenderVideo=true);                ///< The solving function

    // Data
    Vec3_t                  CamPos;         ///< Camera position
    Vec3_t                  Gravity;        ///< Gravity acceleration
    Array <SPHParticle*>    Particles;      ///< Array of SPH particles
    Array <SPHInteracton*>  Interactons;    ///< Array of SPH interactons
    Array <SPHInteracton*>  PInteractons;   ///< Array of SPH possible interactons
    size_t                  idx_out;        ///< Index for output pourposes
    double                  Time;           ///< The simulation Time
    double                  Alpha;          ///< Parameter for verlet lists
    iVec3_t                 N_side;         ///< Vector containing number of division per side of the domain
    Vec3_t                  DXmin;          ///< Point defining the bottom limit of the rectangle domain
    Vec3_t                  DXmax;          ///< Point defining the top limit of the rectangle domain
};

/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////

// Constructor
inline SPHDomain::SPHDomain (iVec3_t n, Vec3_t Xmin, Vec3_t Xmax)
{
    CamPos  = 1.0,2.0,3.0;
    Time    = 0.0;
    Gravity = 0.0,0.0,0.0;
    Alpha   = 0.1;
    N_side  = n;
    DXmin   = Xmin;
    DXmax   = Xmax;
}

inline SPHDomain::~SPHDomain ()
{
    for (size_t i=0; i<Particles.Size();   ++i) if (Particles  [i]!=NULL) delete Particles  [i];
    for (size_t i=0; i<Interactons.Size(); ++i) if (Interactons[i]!=NULL) delete Interactons[i];
}


// Methods
inline void SPHDomain::AddBox(Vec3_t const & V, size_t nx, size_t ny, size_t nz, double R, double s, double rho0, bool Fixed)
{
    Vec3_t C(V);
    C -= Vec3_t((nx-1)*R,(ny-1)*R,(nz-1)*R);
    for (size_t i = 0;i<nx;i++)
    for (size_t j = 0;j<ny;j++)
    for (size_t k = 0;k<nz;k++)
    {
        Vec3_t x(2*i*R,2*j*R,2*k*R);
        x +=C;
        Particles.Push(new SPHParticle(x,OrthoSys::O,rho0,s,Fixed));
    }
}

inline void SPHDomain::AddRandomBox(Vec3_t const & V, double Lx, double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho0, double R, size_t RandomSeed)
{ 
    double start = std::clock();
    std::cout << "[1;33m\n--- Generating random packing of spheres -------------[0m\n";
    const double x_min=-Lx/2.0+V(0), x_max=Lx/2.0+V(0);
    const double y_min=-Ly/2.0+V(1), y_max=Ly/2.0+V(1);
    const double z_min=-Lz/2.0+V(2), z_max=Lz/2.0+V(2);
    double qin = 0.95;
    srand(RandomSeed);
    for (size_t i=0; i<nx; i++)
    {
        for (size_t j=0; j<ny; j++)
        {
            for (size_t k=0; k<nz; k++)
            {
                double x = x_min+(i+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*(x_max-x_min)/nx;
                double y = y_min+(j+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*(y_max-y_min)/ny;
                double z = z_min+(k+0.5*qin+(1-qin)*double(rand())/RAND_MAX)*(z_max-z_min)/nz;
                Particles.Push(new SPHParticle(Vec3_t(x,y,z),OrthoSys::O,rho0,R,false));
            }
        }
    }
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
    std::cout << "[1;32m    Number of particles   = " << Particles.Size() << "[0m\n";
}

inline void SPHDomain::StartAcceleration (Vec3_t const & a)
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->a = a;
        Particles[i]->dDensity = 0.0;
    }
}

inline void SPHDomain::ComputeAcceleration (double dt)
{
    for (size_t i=0; i<PInteractons.Size(); i++) PInteractons[i]->CalcForce(dt);
}

inline void SPHDomain::Move (double dt)
{
    for (size_t i=0; i<Particles.Size(); i++) Particles[i]->Move(dt);
}

inline void SPHDomain::WritePOV (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".pov");
    std::ofstream of(fn.CStr(), std::ios::out);
    POVHeader (of);
    POVSetCam (of, CamPos, OrthoSys::O);
    for (size_t i=0; i<Particles.Size(); i++) 
    {
        if (Particles[i]->IsFree) POVDrawVert(Particles[i]->x,of,Particles[i]->h,"Blue");
        else                      POVDrawVert(Particles[i]->x,of,Particles[i]->h,"Col_Glass_Ruby");
    }
    of.close();
}

inline void SPHDomain::WriteBPY (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".bpy");
    std::ofstream of(fn.CStr(), std::ios::out);
    BPYHeader(of);
    for (size_t i=0; i<Particles.Size(); i++) BPYDrawVert(Particles[i]->x,of,Particles[i]->h);
    of.close();
}

inline void SPHDomain::WriteVTK (char const * FileKey)
{
	// Header
	std::ostringstream oss;
	oss << "# vtk DataFile Version 2.0\n";
	oss << "TimeStep = " << idx_out << "\n";
	oss << "ASCII\n";
	oss << "DATASET STRUCTURED_POINTS\n";
	oss << "DIMENSIONS "   << N_side(0)                     << " " << N_side(1)                      << " " << N_side(2)                      << "\n";
	oss << "ORIGIN "       << DXmin(0)                      << " " << DXmin(1)                       << " " << DXmin(2)                       << "\n";
	oss << "SPACING "      << (DXmax(0)-DXmin(0))/N_side(0) << " " << (DXmax(1)-DXmin(1))/N_side(1)  << " " << (DXmax(2)-DXmin(2))/N_side(2)  << "\n";
	oss << "POINT_DATA "   << N_side(0)*N_side(1)*N_side(2) << "\n";

    Array<double> rho(N_side(0)*N_side(1)*N_side(2));       //Array of density
    Array<double> Press(N_side(0)*N_side(1)*N_side(2));  //Array of pressures

    for (size_t i = 0;i < N_side(0)*N_side(1)*N_side(2); ++i)
    {
        Vec3_t x((DXmax(0)-DXmin(0))/N_side(0)*(i%N_side(0))          +DXmin(0),
                 (DXmax(1)-DXmin(1))/N_side(1)*(i/N_side(0)/N_side(2))+DXmin(1),
                 (DXmax(2)-DXmin(2))/N_side(2)*(i%N_side(2))          +DXmin(2));
        rho [i] = 0.0;
        Press [i] = 0.0;
        for (size_t j = 0;j < Particles.Size();j++)
        {
            if (Particles[j]->IsFree) 
            {
                rho[i]      += Particles[j]->Density*SPHKernel(x,Particles[j]->x,Particles[j]->h);
                Press[i]    += Pressure(Particles[j]->Density)*SPHKernel(x,Particles[j]->x,Particles[j]->h);
            }
        }
    }


	//Density field
	oss << "SCALARS Density float 1\n";
	oss << "LOOKUP_TABLE default\n";
    for (size_t i = 0;i < N_side(0)*N_side(1)*N_side(2); ++i)
    {
        oss << rho[i] << "\n";
    }

	//Density field
	oss << "SCALARS Pressure float 1\n";
	oss << "LOOKUP_TABLE default\n";
    for (size_t i = 0;i < N_side(0)*N_side(1)*N_side(2); ++i)
    {
        oss << Press[i] << "\n";
    }
    
	// Open/create file, write to file, and close file
    String fn(FileKey);
    fn.append(".vtk");
	std::ofstream  of;
	of.open (fn.CStr(), std::ios::out);
	of << oss.str();
	of.close();
}

inline void SPHDomain::ResetInteractons()
{
    // delete old interactors
    for (size_t i=0; i<Interactons.Size(); ++i)
    {
        if (Interactons[i]!=NULL) delete Interactons[i];
    }

    // new interactors
    Interactons.Resize(0);
    for (size_t i=0; i<Particles.Size()-1; i++)
    {
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            // if both particles are fixed, don't create any intereactor
            if (!Particles[i]->IsFree && !Particles[j]->IsFree) continue;
            else 
            {
                //SPHInteracton * I = new SPHInteracton(Particles[i],Particles[j]);
                //if (I->UpdateContacts(Alpha)) Interactons.Push(I);
                //else delete I;
                Interactons.Push(new SPHInteracton(Particles[i],Particles[j]));
            }
        }
    }
}

inline void SPHDomain::ResetDisplacements()
{
    for (size_t i=0; i<Particles.Size(); i++)
    {
        Particles[i]->ResetDisplacements();
    }
}

inline void SPHDomain::ResetContacts()
{
    PInteractons.Resize(0);
    for (size_t i=0; i<Interactons.Size(); i++)
    {
        if(Interactons[i]->UpdateContacts(Alpha)) PInteractons.Push(Interactons[i]);
    }
}

inline double SPHDomain::MaxDisplacement()
{
    double md = 0.0;
    for (size_t i=0; i<Particles.Size(); i++)
    {
        double mpd = Particles[i]->MaxDisplacement();
        if (mpd > md) md = mpd;
    }
    return md;
}

inline void SPHDomain::Solve (double tf, double dt, double dtOut, char const * TheFileKey, bool RenderVideo)
{
    double start = std::clock();
    std::cout << "[1;33m\n--- Solving  --------------------------------------[0m\n";
    idx_out = 0;
    double tout = Time;

    ResetInteractons();
    ResetDisplacements();
    ResetContacts();


    while (Time<tf)
    {
        // Calculate the acceleration for each particle
        StartAcceleration(Gravity);
        ComputeAcceleration(dt);

        // Move each particle
        Move(dt);

        // next time position
        Time += dt;


        // output
        if (Time>=tout)
        {
            idx_out++;
            if (TheFileKey!=NULL)
            {
                if(RenderVideo)
                {
                    String fn;
                    fn.Printf    ("%s_%08d", TheFileKey, idx_out);
                    WritePOV     (fn.CStr());
                    WriteVTK     (fn.CStr());
                }
            }
            tout += dtOut;
        }

        if (MaxDisplacement()>Alpha)
        {
            ResetDisplacements();
            ResetContacts();
        }
        
    }
    double total = std::clock() - start;
    std::cout << "[1;36m    Time elapsed          = [1;31m" <<static_cast<double>(total)/CLOCKS_PER_SEC<<" seconds[0m\n";
}

#endif // MECHSYS_SPH_DOMAIN_H
