/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2014 Sergio Galindo                                    *
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


#ifndef MECHSYS_ADLBM_DOMAIN_H
#define MECHSYS_ADLBM_DOMAIN_H

// STD
#include <map>
#include <vector>
#include <utility>
#include <set>

// MechSys
#include <mechsys/adlbm/Lattice.h>

using std::set;
using std::map;
using std::pair;
using std::make_pair;

namespace ADLBM
{

class Domain
{
public:
    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Constructors
    Domain (
    double                Thenu,
    double                Thedif,
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Methods
#ifdef USE_HDF5
    void WriteXDMF         (char const * FileKey);  ///< Write the domain data in xdmf file
#endif

    void Initialize     (double dt=0.0);                                                                                              ///< Set the particles to a initial state and asign the possible insteractions
    void Collide        (size_t Np = 1);                                                                                ///< Apply the interaction forces and the collision operator
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);                                                                ///< Solve the Domain dynamics

    //Data
    bool                                         Initialized;         ///< System (particles and interactons) initialized ?
    bool                                              PrtVec;         ///< Print Vector data into the xdmf-h5 files
    bool                                            Finished;         ///< Has the simulation finished
    String                                           FileKey;         ///< File Key for output files
    Lattice                                              Lat;         ///< Fluid Lattices
    double                                              Time;         ///< Time of the simulation
    double                                                dt;         ///< Timestep
    void *                                          UserData;         ///< User Data
    size_t                                           idx_out;         ///< The discrete time step
    size_t                                              Step;         ///< The space step to reduce the size of the h5 file for visualization
    size_t                                             Nproc;         ///< Number of cores used for the simulation
};

inline Domain::Domain(double Thenu, double Thedif, iVec3_t Ndim, double Thedx, double Thedt)
{
    Initialized = false;
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing LBM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    Lat = Lattice(Thenu, Thedif, Ndim,Thedx,Thedt);
    Time   = 0.0;
    dt     = Thedt;
    Step   = 1;
    PrtVec = true;
    printf("%s  Num of cells   = %zd%s\n",TERM_CLR2,Lat.Ncells,TERM_RST);
}


#ifdef USE_HDF5

inline void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    size_t  Nx = Lat.Ndim[0]/Step;
    size_t  Ny = Lat.Ndim[1]/Step;
    size_t  Nz = Lat.Ndim[2]/Step;
    // Creating data sets
    float * Rho       = new float[  Nx*Ny*Nz];
    float * Con       = new float[  Nx*Ny*Nz];
    float * Vel       = new float[3*Nx*Ny*Nz];

    size_t i=0;
    for (size_t m=0;m<Lat.Ndim(2);m+=Step)
    for (size_t l=0;l<Lat.Ndim(1);l+=Step)
    for (size_t n=0;n<Lat.Ndim(0);n+=Step)
    {
        double rho    = 0.0;
        double con    = 0.0;
        Vec3_t vel    = OrthoSys::O;

        for (size_t ni=0;ni<Step;ni++)
        for (size_t li=0;li<Step;li++)
        for (size_t mi=0;mi<Step;mi++)
        {
            rho      += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Rho;
            con      += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Con;
            vel (0)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel[0];
            vel (1)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel[1];
            vel (2)  += Lat.GetCell(iVec3_t(n+ni,l+li,m+mi))->Vel[2];
        }
        rho  /= Step*Step*Step;
        con  /= Step*Step*Step;
        vel  /= Step*Step*Step;
        Rho [i]      = (float) rho;
        Con [i]      = (float) con;
        Vel [3*i  ]  = (float) vel (0);
        Vel [3*i+1]  = (float) vel (1);
        Vel [3*i+2]  = (float) vel (2);
        i++;
    } 
    //Write the data
    hsize_t dims[1];
    dims[0] = Nx*Ny*Nz;
    String dsname;
    dsname.Printf("Density");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Rho );
    dsname.Printf("Concentration");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Con );
    if (PrtVec)
    {
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vel );
    }
    dims[0] = 1;
    int N[1];
    N[0] = Nx;
    dsname.Printf("Nx");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
    dims[0] = 1;
    N[0] = Ny;
    dsname.Printf("Ny");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
    dims[0] = 1;
    N[0] = Nz;
    dsname.Printf("Nz");
    H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);

    delete [] Rho     ;
    delete [] Con     ;
    delete [] Vel     ;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"EMLBM_Mesh\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << Step*Lat.dx << " " << Step*Lat.dx  << " " << Step*Lat.dx  << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Density" << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Density" << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Concentration" << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Concentration" << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    if (PrtVec)
    {
    oss << "     <Attribute Name=\"Velocity" << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity" << "\n";
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

void Domain::Collide (size_t Np)
{
#ifdef USE_OMP
    #pragma omp parallel for schedule (static) num_threads(Np)
#endif
    for (size_t i=0;i<Lat.Ncells;i++)
    {
        Cell * c = Lat.Cells[i];
        for (size_t k=0;k<c->Nneigh;k++)
        {
            if (!c->IsSolid)
            {
                c->Ftemp[k] = c->F[k] - (c->F[k] - c->Feq(k))/Lat.Tau;
                c->Gtemp[k] = c->G[k] - (c->G[k] - c->Geq(k))/Lat.Tauc;
            }
            else
            {
                c->Ftemp[k] = c->F[Cell::Op[k]];
                c->Gtemp[k] = c->G[Cell::Op[k]];
            }
        }
        for (size_t k=0;k<c->Nneigh;k++)
        {
            c->F[k] = c->Ftemp[k];
            c->G[k] = c->Gtemp[k];
        }
    }   
}

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{

    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Finished = false;

    // info
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1   , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                   , TERM_RST);

    Nproc = TheNproc;

    //for (size_t j=0;j<Lat.Size();j++)
    //{
        //for (size_t i=0;i<Lat[j].Ncells;i++)
        //{
            //Lat[j].Cells[i]->Initialize();
            //Lat[j].Cells[i]->CalcProp();
        //}
    //}



    double tout = Time;
    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time >= tout)
        {
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
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


#ifdef USE_OMP 
        Collide(Nproc);
        Lat.Stream1  (Nproc);
        Lat.Stream2  (Nproc);
        Lat.CalcProps(Nproc);
#endif

        Time += dt;
    }
    // last output
    Finished = true;
    if (ptReport!=NULL) (*ptReport) ((*this), UserData);

    printf("%s  Final CPU time       = %s\n",TERM_CLR2, TERM_RST);
}
}


#endif

