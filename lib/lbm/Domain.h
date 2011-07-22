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
using std::map;
using std::pair;
using std::make_pair;

class Domain
{
public:
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructor
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    Array<double>         nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Methods
#ifdef USE_HDF5
    void WriteXDMF (char const * FileKey);           ///< Write the domain data in xdmf file
#endif
    void ApplyForce ();                              ///< Apply the interaction forces and the collision operator
    void Collide ();                                 ///< Apply the interaction forces and the collision operator
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
               char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                          ///< Solve the Domain dynamics
    void AddDisk(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt);        ///< Add a disk element
    void ResetContacts();                            ///< Reset contacts for verlet method DEM
    void ResetDisplacements();                       ///< Reset the displacements for the verlet method DEM
    double  MaxDisplacement();                       ///< Give the maximun displacement of DEM particles
    void BoundingBox(Vec3_t & Xmin, Vec3_t & Xmax);  ///< Bounding box for DEM particles



    //Data
    Array<Lattice>                     Lat;         ///< Fluid Lattices
    Array <Disk *>               Particles;         ///< Array of Disks
    Array <Interacton *>       Interactons;         ///< Array of insteractons
    Array <Interacton *>      CInteractons;         ///< Array of valid interactons
    set<pair<Disk *, Disk *> > Listofpairs;         ///< List of pair of particles associated per interacton for memory optimization
    double                            Time;         ///< Time of the simulation
    double                              dt;         ///< Timestep
    double                           Alpha;         ///< Verlet distance
    double                            Gmix;         ///< Interaction constant for the mixture
    void *                        UserData;         ///< User Data
    size_t                         idx_out;         ///< The discrete time step
    
};

inline Domain::Domain(LBMethod Method, Array<double> nu, iVec3_t Ndim, double dx, double Thedt)
{
    for (size_t i=0;i<nu.Size();i++)
    {
        Lat.Push(Lattice(Method,nu[i],Ndim,dx,Thedt));
    }
    Time = 0.0;
    dt   = Thedt;
    Alpha= 10.0;
}

#ifdef USE_HDF5
inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    for (size_t j=0;j<Lat.Size();j++)
    {
        // Creating data sets
        float * Density   = new float[Lat[0].Ndim[0]*Lat[0].Ndim[1]*Lat[0].Ndim[2]];
        float * Gamma     = new float[Lat[0].Ndim[0]*Lat[0].Ndim[1]*Lat[0].Ndim[2]];
        float * Velocity  = new float[Lat[0].Ndim[0]*Lat[0].Ndim[1]*Lat[0].Ndim[2]];
        float * MassFlux  = new float[Lat[0].Ndim[0]*Lat[0].Ndim[1]*Lat[0].Ndim[2]];
        for (size_t i=0;i<Lat[j].Cells.Size();i++)
        {
            double rho;
            Vec3_t vel;
            rho = Lat[j].Cells[i]->VelDen(vel);
            Density  [i] = (float) rho;
            Gamma    [i] = (float) Lat[j].Cells[i]->IsSolid? 1.0:Lat[j].Cells[i]->Gamma;
            Velocity [i] = (float) norm(vel);
            MassFlux [i] = (float) rho*norm(vel);
        }

        //Write the data
        hsize_t dims[1];
        dims[0] = Lat[0].Ndim(0)*Lat[0].Ndim(1)*Lat[0].Ndim(2);
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density );
        dsname.Printf("Gamma_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Gamma   );
        dsname.Printf("Velocity_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velocity);
        dsname.Printf("MassFlux_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,MassFlux);

        delete [] Density ;
        delete [] Gamma   ;
        delete [] Velocity;
        delete [] MassFlux;
    }

    //Closing the file
    H5Fclose(file_id);

	// Writing xmf file
    std::ostringstream oss;
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
    oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Gamma_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Gamma_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"MassFlux_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/MassFlux_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    }
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";

    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}
#endif

inline void Domain::ApplyForce()
{
    for (size_t i=0;i<Lat.Size();i++)
    {
        for (size_t j=0;j<Lat[j].Cells.Size();j++)
        {
            Cell * c = Lat[i].Cells[j];
            if (fabs(c->Gamma-1.0)<1.0e-12) continue;
            double rho = c->Density();
            double psi = Lat[i].Psi(rho);
            for (size_t k=1;k<c->Nneigh;k++)
            {
                //if (c->Neighs[k]>j) continue;
                for (size_t l=0;l<Lat.Size();l++)
                {   
                    Cell * nb     = Lat[l].Cells[c->Neighs[k]];
                    double nb_rho = nb->Density();
                    double nb_psi = Lat[l].Psi(nb_rho);
                    if (i==l)
                    {
                        if (nb->Gamma>0.0||nb->IsSolid)
                        {
                            c ->BForce    += -Lat[i].Gs*psi*c->W[k]*nb_psi*c->C[k];
                            //nb->BForce    -= -Lat[i].Gs*psi*c->W[k]*nb_psi*c->C[k];
                        }
                        else
                        {
                            c ->BForce    += -Lat[i].G*psi*c->W[k]*nb_psi*c->C[k];
                            //nb->BForce    -= -Lat[i].G*psi*c->W[k]*nb_psi*c->C[k];
                        }
                    }
                    else
                    {
                        if (nb->Gamma>0.0||nb->IsSolid)
                        {
                            continue;
                        }
                        else 
                        {
                            c ->BForce    += -Gmix*rho*c->W[k]*nb_rho*c->C[k];
                        }
                    }
                }
            }
        }
    }
}

void Domain::Collide ()
{
    for (size_t i=0;i<Lat[0].Cells.Size()    ;i++)
    {
        Vec3_t num(0.0,0.0,0.0);
        double den = 0.0;
        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double tau = Lat[j].Tau;
            Vec3_t V;
            double rho = c->VelDen(V);
            num += V*rho/tau;
            den += rho/tau;
        }
        Vec3_t Vmix = num/den;


        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            if (c->IsSolid) continue;
            if (fabs(c->Gamma-1.0)<1.0e-12&&fabs(Lat[j].G)>1.0e-12) continue;
            double rho = c->Density();
            double Tau = Lat[j].Tau;
            Vec3_t DV  = Vmix + c->BForce*dt/rho;
            double Bn  = (c->Gamma*(Tau-0.5))/((1.0-c->Gamma)+(Tau-0.5));
            bool valid  = true;
            double alphal = 1.0;
            double alphat = 1.0;
            size_t num  = 0;
            while (valid)
            {
                valid = false;
                alphal  = alphat;
                for (size_t k=0;k<c->Nneigh;k++)
                {
                    double FDeqn = c->Feq(k,DV,rho);
                    c->Ftemp[k] = c->F[k] - alphal*((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]);
                    if (c->Ftemp[k]<0.0&&num<1)
                    {
                        double temp = fabs(c->F[k]/((1 - Bn)*(c->F[k] - FDeqn)/Tau - Bn*c->Omeis[k]));
                        if (temp<alphat) alphat = temp;
                        valid = true;
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
                c->F[j] = fabs(c->Ftemp[j]);
            }
        }
    }   
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

inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    minX = Vec3_t(Particles[0]->X(0) - Particles[0]->R, Particles[0]->X(1) - Particles[0]->R, Particles[0]->X(2) - Particles[0]->R);
    maxX = Vec3_t(Particles[0]->X(0) + Particles[0]->R, Particles[0]->X(1) + Particles[0]->R, Particles[0]->X(2) + Particles[0]->R);
    for (size_t i=1; i<Particles.Size(); i++)
    {
        if (minX(0)>(Particles[i]->X(0) - Particles[i]->R)&&Particles[i]->IsFree()) minX(0) = Particles[i]->X(0) - Particles[i]->R;
        if (minX(1)>(Particles[i]->X(1) - Particles[i]->R)&&Particles[i]->IsFree()) minX(1) = Particles[i]->X(1) - Particles[i]->R;
        if (minX(2)>(Particles[i]->X(2) - Particles[i]->R)&&Particles[i]->IsFree()) minX(2) = Particles[i]->X(2) - Particles[i]->R;
        if (maxX(0)<(Particles[i]->X(0) + Particles[i]->R)&&Particles[i]->IsFree()) maxX(0) = Particles[i]->X(0) + Particles[i]->R;
        if (maxX(1)<(Particles[i]->X(1) + Particles[i]->R)&&Particles[i]->IsFree()) maxX(1) = Particles[i]->X(1) + Particles[i]->R;
        if (maxX(2)<(Particles[i]->X(2) + Particles[i]->R)&&Particles[i]->IsFree()) maxX(2) = Particles[i]->X(2) + Particles[i]->R;
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
    for(size_t j=0;j<Lat.Size();j++)
    {
    	Particles[i]->ImprintDisk(Lat[j]);
    }
    ResetContacts();
    ResetDisplacements();

    double tout = Time;
    while (Time < Tf)
    {
        if (Time >= tout)
        {
            if (FileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", FileKey, idx_out);
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
            tout += dtOut;
            idx_out++;
        }
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);

        //Assigning a vlaue of zero to the particles forces and torques
        for(size_t i=0;i<Particles.Size();i++)
        {
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Ff;
        }

        //Set Gamma values of the lattice cell to zero
        for(size_t j=0;j<Lat.Size();j++)
        {
            Lat[j].SetZeroGamma();
        }
        
        //Connect particles and lattice
        for(size_t i=0;i<Particles.Size();i++)
        for(size_t j=0;j<Lat.Size();j++) Particles[i]->ImprintDisk(Lat[j]);

        //Move Particles
        for(size_t i=0;i<Interactons.Size();i++) Interactons[i]->CalcForce(dt);
        for(size_t i=0;i<Particles.Size()  ;i++) Particles[i]->Translate(dt);

        //Move fluid
        ApplyForce();
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

        Time += dt;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}

#endif

