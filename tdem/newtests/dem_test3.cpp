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
#include <mechsys/mesh/paragrid3d.h>

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
    CellType   CType;    ///< The type of the cell in which this particle is located
    Array<int> Cells;    ///< Cells touched by this particle
};

typedef std::map<int,Particle*>         Id2Part_t;
typedef std::map<int,Array<Particle*> > Cell2Part_t;
typedef std::map<Particle*,Particle*>   Neighbours_t;

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
        //if (Parts[i]->CType!=Outer_t)
        if (Parts[i]->CType==Inner_t || Parts[i]->CType==BryIn_t)
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

inline void PackToData (Particle const * P, Array<double> & Data)
{
    Data.Push (P->Id+0.1);
    Data.Push (P->X(0));
    Data.Push (P->X(1));
    Data.Push (P->X(2));
    Data.Push (P->V(0));
    Data.Push (P->V(1));
    Data.Push (P->V(2));
    Data.Push (P->Xp(0));
    Data.Push (P->Xp(1));
    Data.Push (P->Xp(2));
    Data.Push (P->m);
    Data.Push (P->R);
}

int main(int argc, char **argv) try
{
    // initialize MPI
    MECHSYS_CATCH_PARALLEL = true;
    MECHSYS_MPI_INIT
    int my_id  = MPI::COMM_WORLD.Get_rank();

    // filekey
    String fkey;
    fkey.Printf("dem_test3_proc_%d",my_id);

    // limits of grid
    Array<double> L(6);
    //     0    1      2    3      4    5
    //   xmi  xma    ymi  yma    zmi  zma
    L =  -2.,  2.,   -2.,  2.,    0., 0.1;

    // input
    Array<int> N(10, 10, 1); // number of cells/boxes along each side of grid
    double dt = 0.001;       // timestep
    if (argc> 1) N[0]  = atoi(argv[ 1]);
    if (argc> 2) N[1]  = atoi(argv[ 2]);
    if (argc> 3) N[2]  = atoi(argv[ 3]);
    if (argc> 4) L[0]  = atoi(argv[ 4]);
    if (argc> 5) L[1]  = atoi(argv[ 5]);
    if (argc> 6) L[2]  = atoi(argv[ 6]);
    if (argc> 7) L[3]  = atoi(argv[ 7]);
    if (argc> 8) L[4]  = atoi(argv[ 8]);
    if (argc> 9) L[5]  = atoi(argv[ 9]);
    if (argc>10) dt    = atof(argv[10]);
    if (N[0]<1) throw new Fatal("nx must be greater than 0");
    if (N[1]<1) throw new Fatal("ny must be greater than 0");
    if (N[2]<1) throw new Fatal("nz must be greater than 0");

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

    // allocate particles
    Array<Particle*> parts;   // particles in this processor
    Id2Part_t        id2part; // global Id to particle
    double Ekin0 = 0.0;
    double mass  = 1.0;
    for (size_t id=0; id<xc.Size(); ++id)
    {
        Vec3_t x(xc[id],yc[id],zc[id]);
        int      cell = grid.FindCell  (x);
        CellType type = grid.Cell2Type (cell);
        if (type!=Outer_t) // if it's not outside
        {
            Vec3_t v(vx[id],vy[id],vz[id]);
            Vec3_t xp = x - dt*v;
            parts.Push (new Particle());
            Particle & p = (*parts.Last());
            p.CType  = type;
            p.Id     = id;
            p.X      = x;
            p.Xp     = xp;
            p.V      = v;
            p.m      = mass;
            p.R      = ra[id];
            Ekin0   += 0.5*mass*dot(p.V,p.V);
            id2part[id] = &p;
        }
    }

    // first output
    int stp_out = 0;
    Output (fkey, parts, stp_out);
    stp_out++;

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
            if (parts[i]->CType==Outer_t) continue;
            parts[i]->Start       ();
            parts[i]->Cells.Clear ();
            int cells[8];
            grid.FindCells (parts[i]->X, parts[i]->R, cells);
            for (size_t j=0; j<8; ++j)
            {
                cell2part[cells[j]].XPush (parts[i]);
                parts[i]->Cells.XPush     (cells[j]);
            }
        }

        // find possible contacts
        Neighbours_t neighs;
        for (size_t i=0; i<parts.Size(); ++i)
        {
            // skip outer particles
            if (parts[i]->CType==Outer_t) continue;

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
        Array<double> data;
        for (size_t i=0; i<parts.Size(); ++i)
        {
            CellType type = parts[i]->CType;
            if (type==Outer_t) continue;
            if (type==Inner_t || type==BryIn_t) parts[i]->Move(dt);
            if (type==BryIn_t) PackToData (parts[i], data);
            parts[i]->CType = grid.Cell2Type (grid.FindCell(parts[i]->X));
        }

        // broadcast
        int my_id  = MPI::COMM_WORLD.Get_rank();
        int nprocs = MPI::COMM_WORLD.Get_size();
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

                // unpack data
                for (size_t j=0; j<data.Size(); j+=12)
                {
                    // id and position
                    Vec3_t x(data[j+1],data[j+2],data[j+3]);
                    int      id   = static_cast<int>(data[j]);
                    int      cell = grid.FindCell  (x);
                    CellType type = grid.Cell2Type (cell);
                    if (type!=Outer_t)
                    {
                        bool has_part = id2part.count(id);
                        Particle * p;
                        if (has_part) p = id2part[id];
                        else
                        {
                            p = new Particle();
                            parts.Push (p);
                            id2part[id] = p;
                        }
                        p->Id = id;
                        p->X  = x;
                        p->V  = data[j+4], data[j+5], data[j+6];
                        p->Xp = data[j+7], data[j+8], data[j+9];
                        p->m  = data[j+10];
                        p->R  = data[j+11];
                    }
                }
            }
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
