/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raúl D. D. Farfan             *
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

// Std Lib
#include <iostream>

// Google
#include <google/dense_hash_map>
#include <google/dense_hash_set>

// MechSys
#include <mechsys/linalg/matvec.h>
#include <mechsys/util/maps.h>
#include <mechsys/mesh/structured.h>
#include "paragrid3d.h"

#define MACH_EPS 1.0e-16

using std::cout;
using std::endl;

class Particle
{
public:
    void Start () { F = 0.0,0.0,0.0; }
    void Move  (double dt, bool Verlet=true)
    {
        if (Verlet)
        {
            Vec3_t tmp(X);
            X  = 2.0*X - Xp + (F/m)*dt*dt;
            V  = (X - Xp)/(2.0*dt);
            Xp = tmp;
        }
        else
        {
            V += F*dt/m;
            X += V*dt;
        }
    }

    int        Id;       ///< Global Id of particle
    Vec3_t     X,V,F,Xp; ///< Position, velocity, force, previous position (Verlet)
    double     m,R;      ///< Mass and radius
    bool       Active;   ///< Is particle Active, or deleted from simulation ?
    int        CellTag;  ///< Tag of the cell with this particle
    Array<int> Cells;    ///< Cells touched by this particle
};

inline void CalcForce (Particle & a, Particle & b)
{
    double Kn    = 1000.0;
    Vec3_t nor   = b.X - a.X;
    double dist  = norm(nor);
    double delta = a.R + b.R - dist;
    if (delta > 0.0)
    {
        Vec3_t F = Kn*delta*nor/dist;
        a.F -= F;
        b.F += F;
    }
}

inline void Output (String const & FKey, Array<Particle*> const & Parts, int StpOut)
{
    Table tab;
    tab.SetZero ("id xc yc zc ra vx vy vz",Parts.Size());
    for (size_t i=0; i<Parts.Size(); ++i)
    {
        if (Parts[i]->Active)
        {
            tab("id",i) = Parts[i]->Id;
            tab("xc",i) = Parts[i]->X(0);
            tab("yc",i) = Parts[i]->X(1);
            tab("zc",i) = Parts[i]->X(2);
            tab("ra",i) = Parts[i]->R;
            tab("vx",i) = Parts[i]->V(0);
            tab("vy",i) = Parts[i]->V(1);
            tab("vz",i) = Parts[i]->V(2);
        }
    }
    String buf;
    buf.Printf("%s_%08d.res",FKey.CStr(),StpOut);
    tab.Write(buf.CStr());
}

typedef std::map<int,Array<Particle*> > Cell2Part_t;
typedef std::map<Particle*,Array<int> > Part2Cells_t;
typedef std::map<Particle*,Particle*>   Neighbours_t;

inline void Allocate (Array<double> const & data, Array<Particle*> & parts)
{
    /*
    // push particle into parts
    for (size_t i=0; i<data.Size(); i+=8)
    {
        // check if particle belongs to any box in this processor
        Vec3_t x(data[i+1],data[i+2],data[i+3]);
        double r(data[i+4]);
        //int n = FindBoxNumber (N, L, x, r); // need to allocate particle
        //if (mesh.Cells[n]->PartID!=my_id) continue;

        // allocate particle
        int id = static_cast<int>(data[i]);
        Vec3_t v(data[i+5],data[i+6],data[i+7]);
        Vec3_t xp = x - dt*v; // need to receive xp as well ?
        parts.Push (new Particle());
        parts.Last()->Active = true;
        parts.Last()->Dummy  = (mesh.Cells[n]->PartID!=my_id);
        parts.Last()->Id     = id;
        parts.Last()->X      = x;
        parts.Last()->Xp     = xp;
        parts.Last()->V      = v;
        parts.Last()->m      = mass;
        parts.Last()->R      = r;

        printf("Proc # %d: particle Id = %d arrived and is Dummy = %d\n",my_id,id,parts.Last()->Dummy);

        // update maps
        Array<int> nums(8); // this array needs to be resized externally
        FindBoxNumbers (N, L, parts.Last()->X, parts.Last()->R, nums);
        for (size_t j=0; j<8; ++j)
        {
            Box2Part_t::iterator it = box2part.find(nums[j]);
            if (it==box2part.end()) box2part[nums[j]].Push (parts.Last());
            parts.Last()->Boxes.Push (nums[j]);
        }
    }
    */
}

inline void Transmission ()
{
    /*
    // pack data to send
    Array<double> data; // id,xc,yc,zc,ra,vx,vy,vz
    for (size_t i=0; i<bry_cells.Size(); ++i)
    {
        int n = bry_cells[i];
        Box2Part_t::const_iterator it = box2part.find(n);
        if (it!=box2part.end()) // box has particles
        {
            Array<Particle*> const & parts = it->second; // particles in the cell
            for (size_t j=0; j<parts.Size(); ++j)
            {
                data.Push (parts[j]->Id);
                data.Push (parts[j]->X(0));
                data.Push (parts[j]->X(1));
                data.Push (parts[j]->X(2));
                data.Push (parts[j]->R);
                data.Push (parts[j]->V(0));
                data.Push (parts[j]->V(1));
                data.Push (parts[j]->V(2));
                //printf("Proc # %d, part Id = %d\n",my_id,parts[j]->Id);
            }
        }
    }

    // broadcast
    for (int i=0; i<nprocs; ++i)
    {
        if (i!=my_id)
        {
            MPI::Request req_send = MPI::COMM_WORLD.Isend (data.GetPtr(), data.Size(), MPI::DOUBLE, i, 1000);
            req_send.Wait ();
        }
    }

    // receive messages from everyone
    MPI::Status status;
    for (int i=0; i<nprocs; ++i)
    {
        if (i!=my_id)
        {
            // get message
            MPI::COMM_WORLD.Probe (MPI::ANY_SOURCE, 1000, status);
            int source = status.Get_source();
            int count  = status.Get_count(MPI::DOUBLE);
            data.Resize (count);
            MPI::COMM_WORLD.Recv (data.GetPtr(), count, MPI::DOUBLE, source, 1000);

        }
    }
    */
}

int main(int argc, char **argv) try
{
    // initialize MPI
    MECHSYS_CATCH_PARALLEL = true;
    MECHSYS_MPI_INIT
    int my_id  = MPI::COMM_WORLD.Get_rank();
    int nprocs = MPI::COMM_WORLD.Get_size();

    // filekey
    String fkey;
    fkey.Printf("dem_test2_proc_%d",my_id);

    // input
    Array<int> N(10, 10, 1); // number of cells/boxes along each side of grid
    double dt = 0.001;       // timestep
    if (argc>1) N[0]  = atoi(argv[1]);
    if (argc>2) N[1]  = atoi(argv[2]);
    if (argc>3) N[2]  = atoi(argv[3]);
    if (argc>4) dt    = atof(argv[4]);
    if (N[0]<1) throw new Fatal("nx must be greater than 0");
    if (N[1]<1) throw new Fatal("ny must be greater than 0");
    if (N[2]<1) throw new Fatal("nz must be greater than 0");

    // limits of grid
    Array<double> L(6);
    //     0    1      2    3      4    5
    //   xmi  xma    ymi  yma    zmi  zma
    L =  -2.,  2.,   -2.,  2.,    0., 0.1;

    // grid
    ParaGrid3D grid(N, L, fkey.CStr());
    
    // read data
    Table tab;
    tab.Read ("parts1.dat");
    Array<double> const & xc = tab("xc");
    Array<double> const & yc = tab("yc");
    Array<double> const & zc = tab("zc");
    Array<double> const & ra = tab("ra");
    Array<double> const & vx = tab("vx");
    Array<double> const & vy = tab("vy");
    Array<double> const & vz = tab("vz");

    //std::map<int,Particle*> Id2Part; // global Id to particle

    // allocate particles
    Array<Particle*> parts;
    double Ekin0 = 0.0;
    double mass  = 1.0;
    for (size_t id=0; id<xc.Size(); ++id)
    {
        Vec3_t x(xc[id],yc[id],zc[id]);
        int key = grid.FindCellKey (x, ra[id]);
        int tag = grid.Key2Tag     (key);
        if (tag!=ParaGrid3D::Outer_t) // if it's not outside
        {
            Vec3_t v(vx[id],vy[id],vz[id]);
            Vec3_t xp = x - dt*v;
            parts.Push (new Particle());
            Particle & p = (*parts.Last());
            p.Active  = true;
            p.CellTag = tag;
            p.Id      = id;
            p.X       = x;
            p.Xp      = xp;
            p.V       = v;
            p.m       = mass;
            p.R       = ra[id];
            Ekin0    += 0.5*mass*dot(p.V,p.V);
            //Id2Part[id] = &p;
        }
        //else Id2Part[id] = NULL;
    }

    // first output
    int stp_out = 0;
    Output (fkey, parts, stp_out);
    stp_out++;

    
    Array<double> new_parts; // id,xc,yc,zc,ra,vx,vy,vz


    // solve
    double tf    = 1.0;
    double dtout = 0.1;
    double tout  = 0.1;
    for (double t=0.0; t<tf; t+=dt)
    {
        // initialize particles and find map: cell => particles in/crossed cell
        Cell2Part_t cell2part;
        for (size_t i=0; i<parts.Size(); ++i)
        {
            parts[i]->Start       ();
            parts[i]->Cells.Clear ();
            int keys[8];
            grid.FindCellsKeys (parts[i]->X, parts[i]->R, keys);
            for (size_t j=0; j<8; ++j)
            {
                cell2part[keys[j]].XPush (parts[i]);
                parts[i]->Cells.XPush    (keys[j]);
            }
        }

        /*
        for (Box2Part_t::const_iterator it=box2part.begin(); it!=box2part.end(); ++it)
        {
            cout << it->first << ": ";
            for (size_t k=0; k<it->second.Size(); ++k) cout << it->second[k]->Id << " ";
            cout << endl;
        }
        cout << endl;
        */

        // find possible contacts
        Neighbours_t neighs;
        for (size_t i=0; i<parts.Size(); ++i)
        {
            // particle and cells touched by particle
            Particle         & pa    = (*parts[i]);
            Array<int> const & cells = parts[i]->Cells;

            // loop over the cells touched by this particle
            for (size_t k=0; k<cells.Size(); ++k)
            {
                // cell key => particles
                Cell2Part_t::const_iterator it = cell2part.find(cells[k]);

                // if there are particles touching cell touched by particle
                if (it!=cell2part.end())
                {
                    // all particles touching particle's cells
                    Array<Particle*> const & cell_parts = it->second;

                    // if there are more than one particle (itself) touching this cell
                    if (cell_parts.Size()<1) throw new Fatal("__internal_error__: All cells should have at least one particle");
                    if (cell_parts.Size()>1)
                    {
                        // loop over particles touching cell touched by particle
                        for (size_t j=0; j<cell_parts.Size(); ++j)
                        {
                            // neighbour particle
                            Particle & pb = (*cell_parts[j]);

                            // if pb is not pa
                            if (pb.Id != pa.Id)
                            {
                                // halo overlapping
                                double del = pa.R + pb.R - norm(pb.X-pa.X);

                                // there is overlapping
                                if (del>MACH_EPS)
                                {
                                    if (pa.Id < pb.Id) neighs[&pa] = &pb;
                                    else               neighs[&pb] = &pa;
                                }
                            }
                        }
                    }
                }
                else throw new Fatal("__internal_error__: There should be at least one particle in cell touched by particle; the particle itself");
            }
        }

        // update forces
        for (Neighbours_t::const_iterator r=neighs.begin(); r!=neighs.end(); ++r)
        {
            Particle & pa = (*r->first);
            Particle & pb = (*r->second);
            CalcForce (pa,pb);
        }

        // move particles
        for (size_t i=0; i<parts.Size(); ++i)
        {
            if (parts[i]->CellTag==ParaGrid3D::Inner_t || 
                parts[i]->CellTag==ParaGrid3D::BryIn_t) parts[i]->Move(dt);
        }

        // output
        if (t>=tout)
        {
            Output (fkey, parts, stp_out);
            tout += dtout;
            stp_out++;
        }
    }

    // last output
    Output (fkey, parts, stp_out);

    // energy
    double Ekin1 = 0.0;
    for (size_t i=0; i<parts.Size(); ++i)
    {
        Particle & p = (*parts[i]);
        Ekin1 += 0.5*mass*dot(p.V,p.V);
    }
    printf("\nProc # %d, Ekin (before) = %16.8e\n", my_id, Ekin0);
    printf("Proc # %d, Ekin (after ) = %16.8e\n\n", my_id, Ekin1);

    // write control file
    String buf;
    buf.Printf("%s_control.res",fkey.CStr());
    std::ofstream of(buf.CStr(), std::ios::out);
    of << "fkey  " << fkey      << "\n";
    of << "nout  " << stp_out-1 << "\n";
    of << "nx    " << N[0]      << "\n";
    of << "ny    " << N[1]      << "\n";
    of << "nz    " << N[2]      << "\n";
    of << "lxmi  " << L[0]      << "\n";
    of << "lxma  " << L[1]      << "\n";
    of << "lymi  " << L[2]      << "\n";
    of << "lyma  " << L[3]      << "\n";
    of << "lzmi  " << L[4]      << "\n";
    of << "lzma  " << L[5]      << "\n";
    of.close();

    // end
    MPI::Finalize();
    return 0;
}
MECHSYS_CATCH
