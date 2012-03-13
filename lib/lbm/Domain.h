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

namespace LBM
{

class Domain
{
public:
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructors
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    Array<double>         nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Special constructor with only one component, the parameters are the same as above
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    double                nu,     ///< Viscosity for each fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step


    //Methods
#ifdef USE_HDF5
    void WriteXDMF (char const * FileKey);                                                                                            ///< Write the domain data in xdmf file
#endif
    void ApplyForce     (size_t n = 0, size_t Np = 1);                                                                                ///< Apply the interaction forces and the collision operator
    void Collide        (size_t n = 0, size_t Np = 1);                                                                                ///< Apply the interaction forces and the collision operator
    void ImprintLattice (size_t n = 0, size_t Np = 1);                                                                                ///< Imprint the DEM particles into the lattices
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                                ///< Solve the Domain dynamics
    void AddDisk  (int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt); ///< Add a disk element
    void AddSphere(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt); ///< Add a disk element
    void GenSpheresBox (int Tag, Vec3_t const & X0, Vec3_t const & X1, double R, double rho, size_t Randomseed, double fraction, double RminFraction); ///< Create an array of spheres
    void ResetContacts();                                                                                                             ///< Reset contacts for verlet method DEM
    void ResetDisplacements();                                                                                                        ///< Reset the displacements for the verlet method DEM
    double  MaxDisplacement();                                                                                                        ///< Give the maximun displacement of DEM particles
    void BoundingBox(Vec3_t & Xmin, Vec3_t & Xmax);                                                                                   ///< Bounding box for DEM particles

#ifdef USE_THREAD
    Array<pair<size_t, size_t> >      ListPosPairs;         ///< List of all possible particles pairs
#endif
    //Data
    bool                                    PrtVec;         ///< Print Vector data into the xdmf-h5 files
    Array<Lattice>                             Lat;         ///< Fluid Lattices
    Array <Particle *>                   Particles;         ///< Array of Disks
    Array <Interacton *>               Interactons;         ///< Array of insteractons
    Array <Interacton *>              CInteractons;         ///< Array of valid interactons
    Array <iVec3_t>                      CellPairs;         ///< pairs of cells
    set<pair<Particle *, Particle *> > Listofpairs;         ///< List of pair of particles associated per interacton for memory optimization
    double                                    Time;         ///< Time of the simulation
    double                                      dt;         ///< Timestep
    double                                   Alpha;         ///< Verlet distance
    double                                    Gmix;         ///< Interaction constant for the mixture
    void *                                UserData;         ///< User Data
    size_t                                 idx_out;         ///< The discrete time step
    
};

#ifdef USE_THREAD
struct MtData
{
    size_t                  ProcRank; ///< Rank of the thread
    size_t                    N_Proc; ///< Total number of threads
    LBM::Domain *                Dom; ///< Pointer to the lbm domain
    double                       Dmx; ///< Maximun displacement
    double                        dt; ///< Time step
    Array<pair<size_t,size_t> >   LC; ///< A temporal list of new contacts
    Array<size_t>                LCI; ///< A temporal array of posible Cinteractions
};

void * GlobalIni(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].SetZeroGamma(dat.ProcRank, dat.N_Proc);
    }
    Array<Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    for (size_t i=In;i<Fn;i++)
    {
        (*P)[i]->F = (*P)[i]->Ff;
        (*P)[i]->T = (*P)[i]->Tf;
    }
}

void * GlobalImprint(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    dat.Dom->ImprintLattice(dat.ProcRank, dat.N_Proc);
}

void * GlobalForce(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    Array<Interacton * > * I = &dat.Dom->Interactons;
	size_t Ni = I->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = I->Size() : Fn = (dat.ProcRank+1)*Ni;
	for (size_t i=In;i<Fn;i++)
	{
		(*I)[i]->CalcForce(dat.dt);
	}
}

void * GlobalMove(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.Dmx = 0.0;
	for (size_t i=In;i<Fn;i++)
	{
		(*P)[i]->Translate(dat.dt);
        if (norm((*P)[i]->X-(*P)[i]->X0)>dat.Dmx) dat.Dmx = norm((*P)[i]->X-(*P)[i]->X0);
	}
}

void * GlobalApplyForce (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    dat.Dom->ApplyForce(dat.ProcRank, dat.N_Proc);
}

void * GlobalCollide (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    dat.Dom->Collide(dat.ProcRank, dat.N_Proc);
}

void * GlobalBounceBack (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].BounceBack(dat.ProcRank, dat.N_Proc);
    }
}

void * GlobalStream (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].Stream(dat.ProcRank, dat.N_Proc);
    }
}

void * GlobalStream1 (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].Stream1(dat.ProcRank, dat.N_Proc);
    }
}

void * GlobalStream2 (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    for (size_t i=0;i<dat.Dom->Lat.Size();i++)
    {
        dat.Dom->Lat[i].Stream2(dat.ProcRank, dat.N_Proc);
    }
}

void * GlobalResetDisplacement(void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
    Array<Particle * > * P = &dat.Dom->Particles;
	size_t Ni = P->Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = P->Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.Dmx = 0.0;
	for (size_t i=In;i<Fn;i++)
    {
        (*P)[i]->X0 = (*P)[i]->X;
    }
}

void * GlobalResetContacts1 (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
	size_t Ni = dat.Dom->ListPosPairs.Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = dat.Dom->ListPosPairs.Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LC.Resize(0);

    for (size_t n=In;n<Fn;n++)
    {
        size_t i = dat.Dom->ListPosPairs[n].first;
        size_t j = dat.Dom->ListPosPairs[n].second;
        bool pi_has_vf = !dat.Dom->Particles[i]->IsFree();
        bool pj_has_vf = !dat.Dom->Particles[j]->IsFree();

        bool close = (norm(dat.Dom->Particles[i]->X-dat.Dom->Particles[j]->X)<=dat.Dom->Particles[i]->R+dat.Dom->Particles[j]->R+2*dat.Dom->Alpha);
        if ((pi_has_vf && pj_has_vf) || !close) continue;
        
        // checking if the interacton exist for that pair of particles
        set<pair<Particle *, Particle *> >::iterator it = dat.Dom->Listofpairs.find(make_pair(dat.Dom->Particles[i],dat.Dom->Particles[j]));
        if (it != dat.Dom->Listofpairs.end())
        {
            continue;
        }
        dat.LC.Push(make_pair(i,j));
    }
}

void * GlobalResetContacts2 (void * Data)
{
    LBM::MtData & dat = (*static_cast<LBM::MtData *>(Data));
	size_t Ni = dat.Dom->CInteractons.Size()/dat.N_Proc;
    size_t In = dat.ProcRank*Ni;
    size_t Fn;
    dat.ProcRank == dat.N_Proc-1 ? Fn = dat.Dom->CInteractons.Size() : Fn = (dat.ProcRank+1)*Ni;
    dat.LCI.Resize(0);
    for (size_t n=In;n<Fn;n++)
    {
        if(dat.Dom->CInteractons[n]->UpdateContacts(dat.Dom->Alpha)) dat.LCI.Push(n);
    }
}

#endif

inline Domain::Domain(LBMethod Method, Array<double> nu, iVec3_t Ndim, double dx, double Thedt)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    if (nu.Size()==0) throw new Fatal("LBM::Domain: Declare at leat one fluid please");
    if (Ndim(2) >1&&Method==D2Q9)  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (Ndim(2)==1&&Method==D3Q15) throw new Fatal("LBM::Domain: Ndim(2) is greater than 1. Either change the method to D2Q9 or increse the z-dimension");
    for (size_t i=0;i<nu.Size();i++)
    {
        Lat.Push(Lattice(Method,nu[i],Ndim,dx,Thedt));
    }
    Time   = 0.0;
    dt     = Thedt;
    Alpha  = 10.0;
    PrtVec = false;
    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Size()*Lat[0].Cells.Size(),TERM_RST);
}

inline Domain::Domain(LBMethod Method, double nu, iVec3_t Ndim, double dx, double Thedt)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    Lat.Push(Lattice(Method,nu,Ndim,dx,Thedt));
    if (Ndim(2) >1&&Method==D2Q9)  throw new Fatal("LBM::Domain: D2Q9 scheme does not allow for a third dimension, please set Ndim(2)=1 or change to D3Q15");
    if (Ndim(2)==1&&Method==D3Q15) throw new Fatal("LBM::Domain: Ndim(2) is greater than 1. Either change the method to D2Q9 or increse the z-dimension");
    Time   = 0.0;
    dt     = Thedt;
    Alpha  = 10.0;
    PrtVec = false;
    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Size()*Lat[0].Cells.Size(),TERM_RST);
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
        //float * Velocity  = new float[Lat[0].Ndim[0]*Lat[0].Ndim[1]*Lat[0].Ndim[2]];
        //float * MassFlux  = new float[Lat[0].Ndim[0]*Lat[0].Ndim[1]*Lat[0].Ndim[2]];
        float * Vvec      = new float[3*Lat[0].Ndim[0]*Lat[0].Ndim[1]*Lat[0].Ndim[2]];
        for (size_t i=0;i<Lat[j].Cells.Size();i++)
        {
            double rho = Lat[j].Cells[i]->Rho;
            Vec3_t vel = Lat[j].Cells[i]->Vel;
            Density  [i] = (float) rho;
            Gamma    [i] = (float) Lat[j].Cells[i]->IsSolid? 1.0:Lat[j].Cells[i]->Gamma;
            //Velocity [i] = (float) norm(vel);
            //MassFlux [i] = (float) rho*norm(vel);
            Vvec[3*i  ]  = (float) vel(0);
            Vvec[3*i+1]  = (float) vel(1);
            Vvec[3*i+2]  = (float) vel(2);
        }

        //Write the data
        hsize_t dims[1];
        dims[0] = Lat[0].Ndim(0)*Lat[0].Ndim(1)*Lat[0].Ndim(2);
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density );
        dsname.Printf("Gamma_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Gamma   );
        //dsname.Printf("Velocity_%d",j);
        //H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velocity);
        //dsname.Printf("MassFlux_%d",j);
        //H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,MassFlux);
        dims[0] = 3*Lat[0].Ndim(0)*Lat[0].Ndim(1)*Lat[0].Ndim(2);
        dsname.Printf("Velocity_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vvec    );


        delete [] Density ;
        delete [] Gamma   ;
        //delete [] Velocity;
        //delete [] MassFlux;
        delete [] Vvec    ;
    }

    //Creating data sets
    float * Radius = new float[  Particles.Size()];
    float * Posvec = new float[3*Particles.Size()];
    float * Velvec = new float[3*Particles.Size()];
    int   * Tags   = new int  [  Particles.Size()];

    //Writing particle data
    for (size_t i=0;i<Particles.Size();i++)
    {
        Radius[i]     = (float) Particles[i]->R;
        Posvec[3*i  ] = (float) Particles[i]->X(0);
        Posvec[3*i+1] = (float) Particles[i]->X(1);
        Posvec[3*i+2] = (float) Particles[i]->X(2);
        Velvec[3*i  ] = (float) Particles[i]->V(0);
        Velvec[3*i+1] = (float) Particles[i]->V(1);
        Velvec[3*i+2] = (float) Particles[i]->V(2);
        Tags  [i]     = (int)   Particles[i]->Tag;
    }

    hsize_t dims[1];
    dims[0] = 3*Particles.Size();
    String dsname;
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("Velocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dims[0] = Particles.Size();
    dsname.Printf("Radius");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
    dsname.Printf("Tag");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tags  );


    delete [] Radius;
    delete [] Posvec;
    delete [] Velvec;
    delete [] Tags  ;
    
    //Closing the file
    H5Fclose(file_id);

	// Writing xmf file
    std::ostringstream oss;

    if (Lat[0].Ndim(2)==1)
    {
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
        //oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        //oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        //oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        //oss << "       </DataItem>\n";
        //oss << "     </Attribute>\n";
        //oss << "     <Attribute Name=\"MassFlux_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        //oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        //oss << "        " << fn.CStr() <<":/MassFlux_" << j << "\n";
        //oss << "       </DataItem>\n";
        //oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "   </Grid>\n";
        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    else
    {
        oss << "<?xml version=\"1.0\" ?>\n";
        oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
        oss << "<Xdmf Version=\"2.0\">\n";
        oss << " <Domain>\n";
        oss << "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Lat[0].Ndim(2) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(0) << "\"/>\n";
        oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
        oss << "       </DataItem>\n";
        oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 1.0 1.0 1.0\n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        for (size_t j=0;j<Lat.Size();j++)
        {
        oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Gamma_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Gamma_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        //oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        //oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        //oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        //oss << "       </DataItem>\n";
        //oss << "     </Attribute>\n";
        //oss << "     <Attribute Name=\"MassFlux_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        //oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        //oss << "        " << fn.CStr() <<":/MassFlux_" << j << "\n";
        //oss << "       </DataItem>\n";
        //oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Lat[0].Ndim(0) << " " << Lat[0].Ndim(1) << " " << Lat[0].Ndim(2) << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        }
        oss << "   </Grid>\n";

        oss << "   <Grid Name=\"mesh2\" GridType=\"Uniform\">\n";
        oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
        oss << "     <Geometry GeometryType=\"XYZ\">\n";
        oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
        oss << "        " << fn.CStr() <<":/Position \n";
        oss << "       </DataItem>\n";
        oss << "     </Geometry>\n";
        oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Radius \n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Tag \n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
        oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
        oss << "        " << fn.CStr() <<":/Velocity\n";
        oss << "       </DataItem>\n";
        oss << "     </Attribute>\n";
        oss << "   </Grid>\n";

        oss << " </Domain>\n";
        oss << "</Xdmf>\n";
    }
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}
#endif

inline void Domain::ApplyForce(size_t n, size_t Np)
{
    size_t Ni = CellPairs.Size()/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = CellPairs.Size() : Fn = (n+1)*Ni;

    for (size_t i=In;i<Fn;i++)
    {
        size_t ind1 = CellPairs[i](0);
        size_t ind2 = CellPairs[i](1);
        size_t vec  = CellPairs[i](2);
        for (size_t j=0;j<Lat.Size();j++)
        {
            bool solid = false;
            Cell * c = Lat[j].Cells[ind1];
            double psi;
            if (c->IsSolid||c->Gamma>0.0)
            {
                psi   = 1.0;
                solid = true;
            }
            else fabs(Lat[j].G)>1.0e-12 ? psi = Lat[j].Psi(c->Rho) : psi = c->Rho;
            for (size_t k=0;k<Lat.Size();k++)
            {
                Cell * nb = Lat[k].Cells[ind2];
                double nb_psi;
                if (nb->IsSolid||nb->Gamma>0.0)
                {
                    nb_psi = 1.0;
                    solid  = true;
                }
                else fabs(Lat[j].G)>1.0e-12 ? nb_psi = Lat[k].Psi(nb->Rho) : nb_psi = nb->Rho;
                double G;
                solid ? G = Lat[j].Gs*2.0*ReducedValue(nb->Gs,c->Gs) : G = Lat[j].G; 
                Vec3_t BF(OrthoSys::O);
                if (j==k)
                {
                    BF += -G*psi*nb_psi*c->W[vec]*c->C[vec];
                }
                else if(!solid)
                {
                    BF += -Gmix*c->Rho*nb->Rho*c->W[vec]*c->C[vec];
                }
#ifdef USE_THREAD
                pthread_mutex_lock(&c ->lck);
                pthread_mutex_lock(&nb->lck);
#endif
                c ->BForce += BF;
                nb->BForce -= BF;
#ifdef USE_THREAD
                pthread_mutex_unlock(&c ->lck);
                pthread_mutex_unlock(&nb->lck);
#endif
            }
        }
    }
}

void Domain::Collide (size_t n, size_t Np)
{
	size_t Ni = Lat[0].Cells.Size()/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Lat[0].Cells.Size() : Fn = (n+1)*Ni;
    for (size_t i=In;i<Fn;i++)
    {
        Vec3_t num(0.0,0.0,0.0);
        double den = 0.0;
        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double tau = Lat[j].Tau;
            num += c->Vel*c->Rho/tau;
            den += c->Rho/tau;
        }
        Vec3_t Vmix = num/den;


        for (size_t j=0;j<Lat.Size();j++)
        {
            Cell * c = Lat[j].Cells[i];
            double rho = c->Rho;
            if (c->IsSolid||rho<1.0e-12) continue;
            if (fabs(c->Gamma-1.0)<1.0e-12&&fabs(Lat[j].G)>1.0e-12) continue;
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
                    throw new Fatal("Domain::Collide: Redefine your time step, the current value ensures unstability");
                }
            }
            for (size_t k=0;k<c->Nneigh;k++)
            {
                if (isnan(c->Ftemp[k]))
                {
                    c->Gamma = 2.0;
                    #ifdef USE_HDF5
                    WriteXDMF("error");
                    #endif
                    std::cout << c->Density() << " " << c->BForce << " " << num << " " << alphat << " " << c->Index << " " << c->IsSolid << " " << j << " " << k << std::endl;
                    throw new Fatal("Domain::Collide: Body force gives nan value, check parameters");
                }
                c->F[k] = fabs(c->Ftemp[k]);
            }
        }
    }   
}

void Domain::ImprintLattice (size_t n,size_t Np)
{
    
	size_t Ni = Particles.Size()/Np;
    size_t In = n*Ni;
    size_t Fn;
    n == Np-1 ? Fn = Particles.Size() : Fn = (n+1)*Ni;
    //std::cout << "Im proccess = " << n << std::endl;
    // 2D imprint
    if (Lat[0].Ndim(2)==1)
    {
        for (size_t i = In;i<Fn;i++)
        {
            LBM::Particle * Pa = Particles[i];
            for (size_t n=std::max(0.0,double(Pa->X(0)-Pa->R-Lat[0].dx)/Lat[0].dx);n<=std::min(double(Lat[0].Ndim(0)-1),double(Pa->X(0)+Pa->R+Lat[0].dx)/Lat[0].dx);n++)
            for (size_t m=std::max(0.0,double(Pa->X(1)-Pa->R-Lat[0].dx)/Lat[0].dx);m<=std::min(double(Lat[0].Ndim(1)-1),double(Pa->X(1)+Pa->R+Lat[0].dx)/Lat[0].dx);m++)
            {
                Cell  * cell = Lat[0].GetCell(iVec3_t(n,m,0));
                double x     = Lat[0].dx*cell->Index(0);
                double y     = Lat[0].dx*cell->Index(1);
                double z     = Lat[0].dx*cell->Index(2);
                Vec3_t  C(x,y,z);
                Array<Vec3_t> P(4);

                P[0] = C - 0.5*Lat[0].dx*OrthoSys::e0 - 0.5*Lat[0].dx*OrthoSys::e1;
                P[1] = C + 0.5*Lat[0].dx*OrthoSys::e0 - 0.5*Lat[0].dx*OrthoSys::e1;
                P[2] = C + 0.5*Lat[0].dx*OrthoSys::e0 + 0.5*Lat[0].dx*OrthoSys::e1;
                P[3] = C - 0.5*Lat[0].dx*OrthoSys::e0 + 0.5*Lat[0].dx*OrthoSys::e1;
                double dmin = 2*Pa->R;
                double dmax = 0.0;
                for (size_t j=0;j<P.Size();j++)
                {
                    double dist = norm(P[j] - Pa->X);
                    if (dmin>dist) dmin = dist;
                    if (dmax<dist) dmax = dist;
                }
                if (dmin > Pa->R + Lat[0].dx) continue;

                double len = 0.0;

                if (dmax < Pa->R)
                {
                    len = 4.0;
                }
                else
                {
                    for (size_t j=0;j<4;j++)
                    {
                        Vec3_t D = P[(j+1)%4] - P[j];
                        double a = dot(D,D);
                        double b = 2*dot(P[j]-Pa->X,D);
                        double c = dot(P[j]-Pa->X,P[j]-Pa->X) - Pa->R*Pa->R;
                        if (b*b-4*a*c>0.0)
                        {
                            double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
                            double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
                            if (ta>1.0&&tb>1.0) continue;
                            if (ta<0.0&&tb<0.0) continue;
                            if (ta<0.0) ta = 0.0;
                            if (tb>1.0) tb = 1.0;
                            len += norm((tb-ta)*D);
                        }
                    }
                }

                for (size_t j=0;j<Lat.Size();j++)
                {
                    cell = Lat[j].GetCell(iVec3_t(n,m,0));
                    cell->Gamma   = std::max(len/(4.0*Lat[0].dx),cell->Gamma);
                    if (fabs(cell->Gamma-1.0)<1.0e-12&&fabs(Lat[0].G)>1.0e-12) continue;
                    Vec3_t B      = C - Pa->X;
                    Vec3_t VelP   = Pa->V + cross(Pa->W,B);
                    double rho = cell->Rho;
                    double Bn  = (cell->Gamma*(cell->Tau-0.5))/((1.0-cell->Gamma)+(cell->Tau-0.5));
                    for (size_t k=0;k<cell->Nneigh;k++)
                    {
                        double Fvpp    = cell->Feq(cell->Op[k],VelP,rho);
                        double Fvp     = cell->Feq(k          ,VelP,rho);
                        cell->Omeis[k] = cell->F[cell->Op[k]] - Fvpp - (cell->F[k] - Fvp);
                        Vec3_t Flbm    = -Bn*cell->Omeis[k]*cell->C[k];
                        Pa->F          += Flbm;
                        Pa->T          += cross(B,Flbm);
                    }
                }
            }
        }
    }

    //3D imprint
    else
    {
        for (size_t i = In;i<Fn;i++)
        {
            LBM::Particle * Pa = Particles[i];
            for (size_t n=std::max(0.0,double(Pa->X(0)-Pa->R-Lat[0].dx)/Lat[0].dx);n<=std::min(double(Lat[0].Ndim(0)-1),double(Pa->X(0)+Pa->R+Lat[0].dx)/Lat[0].dx);n++)
            for (size_t m=std::max(0.0,double(Pa->X(1)-Pa->R-Lat[0].dx)/Lat[0].dx);m<=std::min(double(Lat[0].Ndim(1)-1),double(Pa->X(1)+Pa->R+Lat[0].dx)/Lat[0].dx);m++)
            for (size_t l=std::max(0.0,double(Pa->X(2)-Pa->R-Lat[0].dx)/Lat[0].dx);l<=std::min(double(Lat[0].Ndim(2)-1),double(Pa->X(2)+Pa->R+Lat[0].dx)/Lat[0].dx);l++)
            {
                Cell  * cell = Lat[0].GetCell(iVec3_t(n,m,l));
                double x     = Lat[0].dx*cell->Index(0);
                double y     = Lat[0].dx*cell->Index(1);
                double z     = Lat[0].dx*cell->Index(2);
                Vec3_t  C(x,y,z);
                Array<Vec3_t> P(8);

                P[0] = C - 0.5*Lat[0].dx*OrthoSys::e0 - 0.5*Lat[0].dx*OrthoSys::e1 + 0.5*Lat[0].dx*OrthoSys::e2; 
                P[1] = C + 0.5*Lat[0].dx*OrthoSys::e0 - 0.5*Lat[0].dx*OrthoSys::e1 + 0.5*Lat[0].dx*OrthoSys::e2;
                P[2] = C + 0.5*Lat[0].dx*OrthoSys::e0 + 0.5*Lat[0].dx*OrthoSys::e1 + 0.5*Lat[0].dx*OrthoSys::e2;
                P[3] = C - 0.5*Lat[0].dx*OrthoSys::e0 + 0.5*Lat[0].dx*OrthoSys::e1 + 0.5*Lat[0].dx*OrthoSys::e2;
                P[4] = C - 0.5*Lat[0].dx*OrthoSys::e0 - 0.5*Lat[0].dx*OrthoSys::e1 - 0.5*Lat[0].dx*OrthoSys::e2; 
                P[5] = C + 0.5*Lat[0].dx*OrthoSys::e0 - 0.5*Lat[0].dx*OrthoSys::e1 - 0.5*Lat[0].dx*OrthoSys::e2;
                P[6] = C + 0.5*Lat[0].dx*OrthoSys::e0 + 0.5*Lat[0].dx*OrthoSys::e1 - 0.5*Lat[0].dx*OrthoSys::e2;
                P[7] = C - 0.5*Lat[0].dx*OrthoSys::e0 + 0.5*Lat[0].dx*OrthoSys::e1 - 0.5*Lat[0].dx*OrthoSys::e2;

                double dmin = 2*Pa->R;
                double dmax = 0.0;
                for (size_t j=0;j<P.Size();j++)
                {
                    double dist = norm(P[j] - Pa->X);
                    if (dmin>dist) dmin = dist;
                    if (dmax<dist) dmax = dist;
                }
                if (dmin > Pa->R + Lat[0].dx) continue;
                
                double len = 0.0;
                
                if (dmax < Pa->R)
                {
                    len = 12.0*Lat[0].dx;
                }
                else
                {
                    for (size_t j=0;j<4;j++)
                    {
                        Vec3_t D;
                        double a; 
                        double b; 
                        double c; 
                        D = P[(j+1)%4] - P[j];
                        a = dot(D,D);
                        b = 2*dot(P[j]-Pa->X,D);
                        c = dot(P[j]-Pa->X,P[j]-Pa->X) - Pa->R*Pa->R;
                        if (b*b-4*a*c>0.0)
                        {
                            double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
                            double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
                            if (ta>1.0&&tb>1.0) continue;
                            if (ta<0.0&&tb<0.0) continue;
                            if (ta<0.0) ta = 0.0;
                            if (tb>1.0) tb = 1.0;
                            len += norm((tb-ta)*D);
                        }
                        D = P[(j+1)%4 + 4] - P[j + 4];
                        a = dot(D,D);
                        b = 2*dot(P[j + 4]-Pa->X,D);
                        c = dot(P[j + 4]-Pa->X,P[j + 4]-Pa->X) - Pa->R*Pa->R;
                        if (b*b-4*a*c>0.0)
                        {
                            double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
                            double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
                            if (ta>1.0&&tb>1.0) continue;
                            if (ta<0.0&&tb<0.0) continue;
                            if (ta<0.0) ta = 0.0;
                            if (tb>1.0) tb = 1.0;
                            len += norm((tb-ta)*D);
                        }
                        D = P[j+4] - P[j];
                        a = dot(D,D);
                        b = 2*dot(P[j]-Pa->X,D);
                        c = dot(P[j]-Pa->X,P[j]-Pa->X) - Pa->R*Pa->R;
                        if (b*b-4*a*c>0.0)
                        {
                            double ta = (-b - sqrt(b*b-4*a*c))/(2*a);
                            double tb = (-b + sqrt(b*b-4*a*c))/(2*a);
                            if (ta>1.0&&tb>1.0) continue;
                            if (ta<0.0&&tb<0.0) continue;
                            if (ta<0.0) ta = 0.0;
                            if (tb>1.0) tb = 1.0;
                            len += norm((tb-ta)*D);
                        }
                    }
                }

                for (size_t j=0;j<Lat.Size();j++)
                {
                    cell = Lat[j].GetCell(iVec3_t(n,m,l));
                    cell->Gamma   = std::max(len/(12.0*Lat[0].dx),cell->Gamma);
                    if (fabs(cell->Gamma-1.0)<1.0e-12&&fabs(Lat[0].G)>1.0e-12) 
                    {
                        continue;
                    }
                    Vec3_t B      = C - Pa->X;
                    Vec3_t VelP   = Pa->V + cross(Pa->W,B);
                    double rho = cell->Rho;
                    double Bn  = (cell->Gamma*(cell->Tau-0.5))/((1.0-cell->Gamma)+(cell->Tau-0.5));
                    for (size_t k=0;k<cell->Nneigh;k++)
                    {
                        double Fvpp    = cell->Feq(cell->Op[k],VelP,rho);
                        double Fvp     = cell->Feq(k          ,VelP,rho);
                        cell->Omeis[k] = cell->F[cell->Op[k]] - Fvpp - (cell->F[k] - Fvp);
                        Vec3_t Flbm    = -Bn*cell->Omeis[k]*cell->C[k];
                        Pa->F          += Flbm;
                        Pa->T          += cross(B,Flbm);
                    }
                }
            }
        }
    }
    //std::cout << "Im finished = " << n << std::endl;
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
            set<pair<Particle *, Particle *> >::iterator it = Listofpairs.find(make_pair(Particles[i],Particles[j]));
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

inline void Domain::AddSphere(int TheTag, Vec3_t const & TheX, Vec3_t const & TheV, Vec3_t const & TheW, double Therho, double TheR, double dt)
{
    Particles.Push(new Sphere(TheTag,TheX,TheV,TheW,Therho,TheR,dt));
}

inline void Domain::GenSpheresBox (int Tag, Vec3_t const & X0, Vec3_t const & X1, double R, double rho, size_t Randomseed, double fraction, double RminFraction)
{
    // find radius from the edge's length
    Util::Stopwatch stopwatch;
    printf("\n%s--- Generating packing of spheres -----------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    srand(Randomseed);
    //size_t nx = 0.5*(X1(0)-X0(0))/R;
    //size_t ny = 0.5*(X1(1)-X0(1))/R;
    //size_t nz = 0.5*(X1(2)-X0(2))/R;
    //std::cout << nx << " " << ny << " " << nz << std::endl;

    //for (size_t i=0; i<nx; ++i)
    //{
        //for (size_t j=0; j<ny; ++j)
        //{
            //for (size_t k=0; k<nz; ++k)
            //{
                //Vec3_t pos((i+0.5)*2.0*R,(j+0.5)*2.0*R,(k+0.5)*2.0*R);
                //pos += X0;
                //if (rand()<fraction*RAND_MAX) AddSphere (Tag,pos,OrthoSys::O,OrthoSys::O,rho,R*RminFraction+(1.0*rand())/RAND_MAX*(R-R*RminFraction),Lat[0].dt);
            //}
        //}
    //}
   



    size_t nx = 0.5*(X1(0)-X0(0))/R-1;
    size_t ny = int((X1(1)-X0(1))/(sqrt(3.0)*R));
    size_t nz = int((X1(2)-X0(2))/(sqrt(8.0/3.0)*R));
    for (size_t k = 0; k < nz; k++)
    {
        for (size_t j = 0; j < ny; j++)
        {
            Vec3_t X;
            if (k%2==0) X = Vec3_t(-R,R,2*R+k*sqrt(8.0/3.0)*R) + X0;
            else X = Vec3_t(0.0,R+sqrt(1.0/3.0)*R,2*R+k*sqrt(8.0/3.0)*R) + X0;
            if (j%2==0) X += Vec3_t(R,j*sqrt(3.0)*R,0.0);
            else X += Vec3_t(0.0,j*sqrt(3.0)*R,0.0);
            for (size_t i = 0; i < nx; i++)
            {
                X += Vec3_t(2*R,0.0,0.0);
                if (rand()<fraction*RAND_MAX) AddSphere (Tag,X,OrthoSys::O,OrthoSys::O,rho,R*RminFraction+(1.0*rand())/RAND_MAX*(R-R*RminFraction),Lat[0].dt);
            }
        }
    }
    
    printf("%s  Num of particles   = %zd%s\n",TERM_CLR2,Particles.Size(),TERM_RST);
}

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * FileKey, bool RenderVideo, size_t Nproc)
{
    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    idx_out     = 0;
    //Connect particles and lattice
    //
    
    ImprintLattice();
    ResetContacts();
    ResetDisplacements();
    
    // Creates pair of cells to speed up body force calculation
    for (size_t i=0;i<Lat[0].Cells.Size();i++)
    {
        Cell * c = Lat[0].Cells[i];
        for (size_t j=1;j<c->Nneigh;j++)
        {
            Cell * nb = Lat[0].Cells[c->Neighs[j]];
            if (nb->ID>c->ID) 
            {
                if (!c->IsSolid||!nb->IsSolid) CellPairs.Push(iVec3_t(i,nb->ID,j));
            }
        }
    }
#ifdef USE_THREAD
    LBM::MtData MTD[Nproc];
    for (size_t i=0;i<Nproc;i++)
    {
        MTD[i].N_Proc   = Nproc;
        MTD[i].ProcRank = i;
        MTD[i].Dom      = this;
        MTD[i].Dmx      = 0.0;
        MTD[i].dt       = Lat[0].dt;
    }
    pthread_t thrs[Nproc];   
    
    if (Particles.Size() > 0)
    {
        for (size_t i=0; i<Particles.Size()-1; i++)
        for (size_t j=i+1; j<Particles.Size(); j++)
        {
            ListPosPairs.Push(make_pair(i,j));
        }
    }
#endif
    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
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


#ifdef USE_THREAD
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalIni, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << "1" <<std::endl;
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalImprint, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << "2" <<std::endl;
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalForce, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        double maxdis = 0.0;
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalMove, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
            if (maxdis<MTD[i].Dmx) maxdis = MTD[i].Dmx;
        }
        //std::cout << "3" <<std::endl;
        if (maxdis>Alpha)
        {
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalResetDisplacement, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalResetContacts1, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
                for (size_t j=0;j<MTD[i].LC.Size();j++)
                {
                    size_t n = MTD[i].LC[j].first;
                    size_t m = MTD[i].LC[j].second;
                    Listofpairs.insert(make_pair(Particles[n],Particles[m]));
                    CInteractons.Push (new Interacton(Particles[n],Particles[m]));
                }
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalResetContacts2, &MTD[i]);
            }
            Interactons.Resize(0);
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
                for (size_t j=0;j<MTD[i].LCI.Size();j++)
                {
                    Interactons.Push(CInteractons[MTD[i].LCI[j]]);
                }
            }
        }
        //std::cout << "4" <<std::endl;
        if (Lat.Size()>1||fabs(Lat[0].G)>1.0e-12)
        {
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_create(&thrs[i], NULL, GlobalApplyForce, &MTD[i]);
            }
            for (size_t i=0;i<Nproc;i++)
            {
                pthread_join(thrs[i], NULL);
            }
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalCollide, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalBounceBack, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalStream1, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_create(&thrs[i], NULL, GlobalStream2, &MTD[i]);
        }
        for (size_t i=0;i<Nproc;i++)
        {
            pthread_join(thrs[i], NULL);
        }
        //std::cout << "5" <<std::endl;
#else 
        //Assigning a vlaue of zero to the particles forces and torques
        for(size_t i=0;i<Particles.Size();i++)
        {
            Particles[i]->F = Particles[i]->Ff;
            Particles[i]->T = Particles[i]->Tf;
        }

        //Set Gamma values of the lattice cell to zero
        for(size_t j=0;j<Lat.Size();j++)
        {
            Lat[j].SetZeroGamma();
        }
        //Connect particles and lattice
        ImprintLattice();

        //Move Particles
        for(size_t i=0;i<Interactons.Size();i++) Interactons[i]->CalcForce(dt);
        for(size_t i=0;i<Particles.Size()  ;i++) Particles[i]->Translate(dt);

        //Move fluid
        if (Lat.Size()>1||fabs(Lat[0].G)>1.0e-12) ApplyForce();
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
#endif

        Time += dt;
    }
    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}
}


#endif

