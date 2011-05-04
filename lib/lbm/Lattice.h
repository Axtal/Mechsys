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

#ifndef MECHSYS_LBM_LATTICE_H
#define MECHSYS_LBM_LATTICE_H


// MechSys
#include <mechsys/lbm/Cell.h>
#include <mechsys/util/util.h>
#include <mechsys/util/stopwatch.h>

class Lattice
{
public:
    //typedefs
    typedef void (*ptFun_t) (Lattice & Lat, void * UserData);

    //Constructors
    Lattice () {};            //Default
    Lattice (LBMethod Method, double nu, iVec3_t Ndim, double dx, double dt);

    //Methods
    void Solve(double Tf, double dtOut, ptFun_t ptSetup=NULL, ptFun_t ptReport=NULL,
               char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);   ///< Solve the LBM equation in time
    void Stream();                                                                  ///< Stream the velocity distributions
    void Homogenize();                                                              ///< Homogenize the initial state
    void SetZeroGamma();                                                            ///< Set an initial value for the fluid/solid ratio of each cell
    void ApplyBC();                                                                 ///< apply the boundary conditions
    double Psi(double);                                                             ///< Interaction potential
    void ApplyForce();                                                              ///< Apply molecular forces
    void Collide();                                                                 ///< apply the collision operator
    void CollideAlt();                                                              ///< apply the collision operator
    void BounceBack();                                                              ///< apply interaction with solids
    void WriteVTK(char const * FileKey);                                            ///< Write the state in a VTK file
    Cell * GetCell(iVec3_t const & v);                                              ///< Get pointer to cell at v


     

    //Data
    size_t           idx_out;          // The discrete time step
    double           Time;             // The current time
    double           G;                // Interaction strength
    double           Gs;               // Interaction strength
    double           Nu;               // Real viscosity
    iVec3_t          Ndim;             // Integer dimension of the domain
    double           dx;               // grid space
    double           dt;               // time step
    double           Tau;              // Relaxation time
    double           Rhoref;           //Values for th intermolecular force
    double           Psiref;           // 
    Array<Cell *>    Cells;            // Array of pointer cells
    void *           UserData;         // User Data
};

inline Lattice::Lattice(LBMethod TheMethod, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Nu   = Thenu;
    Ndim = TheNdim;
    dx   = Thedx;
    dt   = Thedt;
    Tau  = 3.0*Nu*dt/(dx*dx) + 0.5;
    Rhoref = 200.0;
    Psiref = 4.0;
    G      = 0.0;

    Cells.Resize(Ndim[0]*Ndim[1]*Ndim[2]);
    size_t n = 0;
    for (size_t k=0;k<Ndim[2];k++)
    for (size_t j=0;j<Ndim[1];j++)
    for (size_t i=0;i<Ndim[0];i++)
    {
        Cells[n] =  new Cell(n,TheMethod,iVec3_t(i,j,k),Ndim,dx/dt,Tau);
        n++;
    }  
}

inline void Lattice::Stream()
{
    // Assign temporal distributions
    for (size_t i=0;i<Cells.Size()    ;i++)
    for (size_t j=0;j<Cells[i]->Nneigh;j++)
    {
        Cells[Cells[i]->Neighs[j]]->Ftemp[j] = Cells[i]->F[j];
    }

    // Swap the distribution values
    for (size_t i=0;i<Cells.Size()    ;i++)
    for (size_t j=0;j<Cells[i]->Nneigh;j++)
    {
        Cells[i]->F[j] = Cells[i]->Ftemp[j];
    }
}

inline void Lattice::SetZeroGamma()
{
    for (size_t i=0;i<Cells.Size();i++)
    {
        Cells[i]->Gamma  = 0.0;
        Cells[i]->BForce = Cells[i]->BForcef;
    }
}

inline double Lattice::Psi(double rho)
{
    return Psiref*exp(-Rhoref/rho);
}

inline void Lattice::ApplyForce()
{
    for (size_t i=0;i<Cells.Size();i++)
    {
        Cell * c = Cells[i];
        double psi = Psi(c->Density());
        if (fabs(c->Gamma-1.0)<1.0e-12) continue;
        for (size_t j=1;j<c->Nneigh;j++)
        {
            Cell * nb     = Cells[c->Neighs[j]];
            double nb_psi = Psi(nb->Density());
            double C      = G;
            if (nb->Gamma>0.0||nb->IsSolid)
            {
                nb_psi = 1.0;
                C      = Gs;
            }
            c->BForce    += -C*psi*c->W[j]*nb_psi*c->C[j];
        }
    }
}

inline void Lattice::Collide()
{
    double ome = 1.0/Tau;
    for (size_t i=0;i<Cells.Size()    ;i++)
    {
        Cell * c = Cells[i];
        if (c->IsSolid) continue;
        //if (fabs(c->Gamma-1.0)<1.0e-12) continue;
        Vec3_t V;
        double rho = c->VelDen(V);
        double Bn  = (c->Gamma*(Tau-0.5))/((1.0-c->Gamma)+(Tau-0.5));
        for (size_t j=0;j<c->Nneigh;j++)
        {
            double Feqn = c->Feq(j,       V,rho);
            //double Fvp  = c->Feq(j,c->VelP,rho);
            c->F[j] = c->F[j] - (1 - Bn)*ome*(c->F[j] - Feqn) + Bn*c->Omeis[j];
        }
    }
}

inline void Lattice::CollideAlt()
{
    double ome = 1.0/Tau;
    for (size_t i=0;i<Cells.Size()    ;i++)
    {
        Cell * c = Cells[i];
        if (c->IsSolid) continue;
        if (fabs(c->Gamma-1.0)<1.0e-12) continue;
        Vec3_t V;
        double rho = c->VelDen(V);
        Vec3_t DV  = V + c->BForce*Tau/rho;
        double Bn  = (c->Gamma*(Tau-0.5))/((1.0-c->Gamma)+(Tau-0.5));
        bool valid  = true;
        double omel = ome;
        double omet = ome;
        size_t num  = 0;
        while (valid)
        {
            valid = false;
            omel  = omet;
            for (size_t j=0;j<c->Nneigh;j++)
            {
                //double Feqn  = c->Feq(j,       V,rho);
                double FDeqn = c->Feq(j,      DV,rho);
                //double Fvp  = c->Feq(j, c->VelP,rho);

                //First method owen
                //c->F[j] = c->F[j] - (1 - Bn)*ome*(c->F[j] - Feqn) + Bn*c->Omeis[j] + dt*dt*(1 - Bn)*c->W[j]*dot(c->BForce,c->C[j])/(K*dx);

                //Second method LBM EOS
                //c->F[j] = c->F[j] - (1 - Bn)*(ome*(c->F[j] - Feqn) + FDeqn - Feqn) + Bn*c->Omeis[j];
                
                //Third method sukop
                c->Ftemp[j] = c->F[j] - (1 - Bn)*omel*(c->F[j] - FDeqn) + Bn*c->Omeis[j];

                //if (c->Ftemp[j]<-1.0e-6&&Bn<1.0e-6&&false)
                if (c->Ftemp[j]<-1.0e-6)
                {
                    double temp = (c->F[j] + Bn*c->Omeis[j])/((1 - Bn)*(c->F[j] - FDeqn));
                    if (temp<omet) omet = temp;
                    valid = true;
                    //std::cout << temp << " " << c->F[j] << " " << c->Ftemp[j] << " " << FDeqn << " " << Bn << " " << c->Omeis[j] <<  std::endl;
                }
            }
            num++;
            if (num>2) 
            {
                throw new Fatal("Lattice::Collide: Redefine your time step, the current value ensures unstability");
            }
        }
        for (size_t j=0;j<c->Nneigh;j++)
        {
            c->F[j] = c->Ftemp[j];
        }
    }
}

inline void Lattice::BounceBack()
{
    for (size_t i=0;i<Cells.Size()    ;i++)
    {
        if (!Cells[i]->IsSolid) continue;
        for (size_t j = 0;j<Cells[i]->Nneigh;j++) Cells[i]->Ftemp[j] = Cells[i]->F[j];
        for (size_t j = 0;j<Cells[i]->Nneigh;j++) Cells[i]->F[j]     = Cells[i]->Ftemp[Cells[i]->Op[j]];
    }
}

inline void Lattice::WriteVTK(char const * FileKey)
{
	// Header
	std::ostringstream oss;
	oss << "# vtk DataFile Version 2.0\n";
	oss << "TimeStep = " << idx_out << "\n";
	oss << "ASCII\n";
	oss << "DATASET STRUCTURED_POINTS\n";
	oss << "DIMENSIONS " << Ndim[0] << " " << Ndim[1] << " " << Ndim[2] << "\n";
	oss << "ORIGIN "     << 0       << " " << 0       << " " << 0       << "\n";
	oss << "SPACING "    << 1       << " " << 1       << " " << 1       << "\n";
	oss << "POINT_DATA " << Cells.Size()   << "\n";

	// Solid cells
	oss << "SCALARS Geom float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<Cells.Size(); ++i)
	{
		if (Cells[i]->IsSolid) oss << "1.0\n";
		else                   oss << "0.0\n";
	}

	// Density field
	oss << "SCALARS Density float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<Cells.Size(); i++)
		oss << Cells[i]->Density() << "\n";

	// Density field
	oss << "SCALARS Gamma float 1\n";
	oss << "LOOKUP_TABLE default\n";
	for (size_t i=0; i<Cells.Size(); i++)
		oss << Cells[i]->Gamma << "\n";

	oss << "VECTORS Velocity float\n";
	for (size_t i=0; i<Cells.Size(); ++i)
	{
		Vec3_t v;  Cells[i]->Velocity(v);
		oss << v(0) << " " << v(1) << " " << v(2) << "\n";
	}
	oss << "VECTORS Mass_flux float\n";
	for (size_t i=0; i<Cells.Size(); ++i)
	{
		Vec3_t v; Cells[i]->Velocity(v);
		v *= Cells[i]->Density();
		oss << v(0) << " " << v(1) << " " << v(2) << "\n";
	}

    String fn(FileKey);
    fn.append(".vtk");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline Cell * Lattice::GetCell(iVec3_t const & v)
{
    return Cells[v[0] + v[1]*Ndim[0] + v[2]*Ndim[0]*Ndim[1]];
}

inline void Lattice::Solve(double Tf, double dtOut, ptFun_t ptSetup, ptFun_t ptReport,
                           char const * FileKey, bool RenderVideo, size_t Nproc)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    idx_out     = 0;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        Collide();
        BounceBack();
        Stream();
        
        if (Time >= tout)
        {
            if (FileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%08d", FileKey, idx_out);
                WriteVTK     (fn.CStr());
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
