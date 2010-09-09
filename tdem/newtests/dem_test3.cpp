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
typedef std::map<std::pair<Particle*,Particle*>,double> NeighDist_t;

inline void CalcForce (Particle & a, Particle & b)
{
    double Kn    = 1000.0;
    Vec3_t nor   = b.X - a.X;
    double dist  = norm(nor);
    double delta = a.R + b.R - dist;
    if (delta > MACH_EPS)
    {
        Vec3_t F = Kn*delta*nor/dist;
        a.F -= F;
        b.F += F;
    }
}

inline void Output (String const & FKey, Array<Particle*> const & Parts, int StpOut)
{
    // header
    Array<String> keys("id", "xc", "yc", "zc", "ra", "vx", "vy", "vz", "ct");
    std::ostringstream oss;
    oss << Util::_6 << keys[0];
    for (size_t i=1; i<keys.Size()-1; ++i) { oss << Util::_8s << keys[i]; }
    oss << Util::_6 << keys.Last() << "\n";

    // values
    for (size_t i=0; i<Parts.Size(); ++i)
    {
        oss << Util::_6  << Parts[i]->Id;
        oss << Util::_8s << Parts[i]->X(0);
        oss << Util::_8s << Parts[i]->X(1);
        oss << Util::_8s << Parts[i]->X(2);
        oss << Util::_8s << Parts[i]->R;
        oss << Util::_8s << Parts[i]->V(0);
        oss << Util::_8s << Parts[i]->V(1);
        oss << Util::_8s << Parts[i]->V(2);
        oss << Util::_6  << Parts[i]->CType;
        oss << "\n";
    }

    // open file and save data
    String buf;   buf.Printf("%s_%08d.res",FKey.CStr(),StpOut);
    std::ofstream of(buf.CStr(), std::ios::out);
    of << oss.str();
    of.close();
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
    // mpirun -np 3 ./dem_test3 10 10 1 -1 1 -1 2
    // mpirun -np 5 ./dem_test3 10 10 1 -0.2 1.2 -0.2 1.2 
    
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
    if (argc> 4) L[0]  = atof(argv[ 4]);
    if (argc> 5) L[1]  = atof(argv[ 5]);
    if (argc> 6) L[2]  = atof(argv[ 6]);
    if (argc> 7) L[3]  = atof(argv[ 7]);
    if (argc> 8) L[4]  = atof(argv[ 8]);
    if (argc> 9) L[5]  = atof(argv[ 9]);
    if (argc>10) dt    = atof(argv[10]);
    bool calc_N = false;
    if (N[0]<1) calc_N = true;
    if (N[1]<1) calc_N = true;
    if (N[2]<1) calc_N = true;
    
    // read data
    Table tab;
    tab.Read ("parts4.dat");
    Array<double> const & xc = tab("xc");
    Array<double> const & yc = tab("yc");
    Array<double> const & zc = tab("zc");
    Array<double> const & ra = tab("ra");
    Array<double> const & vx = tab("vx");
    Array<double> const & vy = tab("vy");
    Array<double> const & vz = tab("vz");

    /*
    Array<double> xc;
    Array<double> yc;
    Array<double> zc;
    Array<double> ra;
    Array<double> vx;
    Array<double> vy;
    Array<double> vz;
    // Generating random packing of large number of particles
    size_t Np  = 400; // Number of particles
    size_t i   = 0;
    double rad = 0.0494;
    double v0  = 0.05;
    while (i<Np)
    {
        double xmi = -1.5;
        double xma =  1.5;
        double ymi = -1.5;
        double yma =  1.5;
        Vec3_t x(xmi+(xma-xmi)/10 + (1.0*rand())/RAND_MAX*(xma-(xma-xmi)/10-xmi-(xma-xmi)/10),
                 ymi+(yma-ymi)/10 + (1.0*rand())/RAND_MAX*(yma-(yma-ymi)/10-ymi-(yma-ymi)/10),
                 0.05);
        bool valid = true;
        for (size_t j=0;j<xc.Size();j++)
        {
            Vec3_t xp(xc[j],yc[j],zc[j]);
            if (norm(x-xp)<2*rad)
            {
                valid = false;
                break;
            }
        }
        if (valid)
        {
            xc.Push(x(0));
            yc.Push(x(1));
            zc.Push(x(2));
            ra.Push(rad);
            vx.Push(v0*((1.0*rand())/RAND_MAX-0.5));
            vy.Push(v0*((1.0*rand())/RAND_MAX-0.5));
            vz.Push(0.0);
            i++;
        }
    }
    */

    // find N
    if (calc_N)
    {
        double maxR = 0.0;
        for (size_t i=0; i<ra.Size(); ++i) if (ra[i]>maxR) maxR = ra[i];
        N[0] = static_cast<int>((L[1]-L[0])/(1.01*2.0*maxR));
        N[1] = static_cast<int>((L[3]-L[2])/(1.01*2.0*maxR));
        N[2] = static_cast<int>((L[5]-L[4])/(1.01*2.0*maxR));
        if (N[0]==0) N[0] = 1;
        if (N[1]==0) N[1] = 1;
        if (N[2]==0) N[2] = 1;
        if (my_id==0) printf("maxR = %g, nx,ny,nz = %d,%d,%d\n",maxR,N[0],N[1],N[2]);
    }

    // grid
    ParaGrid3D grid(N, L, fkey.CStr());

    /*
    // neighbour partitions
    cout << "Proc # " << my_id << ": neighbours = ";
    for (size_t i=0; i<grid.NNeighParts(); ++i) cout << grid.NeighPart(i) << " ";
    cout << endl;
    */

    // allocate particles
    Array<Particle*> parts;   // particles in this processor
    Id2Part_t        id2part; // global Id to particle
    double Ekin0 = 0.0;
    double mass  = 1.0;
    for (int id=0; id<static_cast<int>(xc.Size()); ++id)
    {
        double d = 2.0*ra[id];
        if (d>grid.D[0] || d>grid.D[1] || d>grid.D[2]) throw new Fatal("Particles must have diameter smaller than the smallest cell in grid. diam=%g is greater than %g, %g, or %g",d,grid.D[0],grid.D[1],grid.D[2]);
        Vec3_t x(xc[id],yc[id],zc[id]);
        int      cell = grid.FindCell (x);
        CellType type = grid.Id2Type[cell];
        if (type!=Outer_t) // if it's not outside
        {
            //if (my_id==0) printf("Proc # 0 is adding particle # %d\n",id);
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
    //printf("Proc # %d has %zd particles\n",my_id,parts.Size());

    // first output
    int stp_out = 0;
    Output (fkey, parts, stp_out);
    stp_out++;

    // solve
    double tf    = 5.0;
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
        NeighDist_t neighs;
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

                            // if pb is not pa and these two particles are not in my outer boundary (the other processor will calculate this force for me)
                            if ((pb.Id!=pa.Id) && (!(pa.CType==BryOut_t && pb.CType==BryOut_t)))
                            {
                                // halo overlapping
                                double del = pa.R + pb.R - norm(pb.X-pa.X);

                                // there is overlapping
                                if (del>MACH_EPS)
                                {
                                    if (pa.Id < pb.Id) neighs[std::make_pair(&pa,&pb)] = del;
                                    else               neighs[std::make_pair(&pb,&pa)] = del;
                                }
                            }
                        }
                    }
                }
                else throw new Fatal("__internal_error__: There should be at least one particle in cell touched by particle; the particle itself");
            }
        }

        // update forces
        for (NeighDist_t::const_iterator r=neighs.begin(); r!=neighs.end(); ++r)
        //for (size_t i=0;   i<parts.Size()-1; ++i)
        //for (size_t j=i+1; j<parts.Size();   ++j)
        {
            //Particle & pa = *parts[i];
            //Particle & pb = *parts[j];
            Particle & pa = (*r->first.first);
            Particle & pb = (*r->first.second);
            CalcForce (pa,pb);
        }

        // move particles
        std::map<int,Array<double> > proc2data; // map: proc # => packed data to be sent
        for (size_t i=0; i<parts.Size(); ++i)
        {
            int      cell = grid.FindCell (parts[i]->X);
            CellType type = parts[i]->CType;
            if (type==Outer_t) continue;
            if (type==Inner_t || type==BryIn_t) parts[i]->Move(dt);
            if (type==BryIn_t)
            {
                for (size_t j=0; j<grid.Id2NeighProcs[cell].Size(); ++j)
                {
                    Array<double> & data_snd = proc2data[grid.Id2NeighProcs[cell][j]];
                    //printf("Proc # %d will send particle # %d to processor to # %d\n",my_id,parts[i]->Id,grid.Id2NeighProcs[cell][j]);
                    PackToData (parts[i], data_snd);
                }
            }
            parts[i]->CType = grid.Id2Type[grid.FindCell(parts[i]->X)];
        }

        // synchronize
        MPI::COMM_WORLD.Barrier(); // TODO: check if this is really needed

        // post messages with data_sent
        Array<MPI::Request> req(grid.NeighProcs.Size());
        for (size_t i=0; i<grid.NeighProcs.Size(); ++i)
        {
            //printf("Proc # %d is sending to # %d\n",my_id,grid.NeighProcs[i]);
            Array<double> const & data_snd = proc2data[grid.NeighProcs[i]];
            req[i] = MPI::COMM_WORLD.Isend (data_snd.GetPtr(), data_snd.Size(), MPI::DOUBLE, grid.NeighProcs[i], 1000);
        }

        // receive messages
        Array<double> data;
        for (size_t i=0; i<grid.NeighProcs.Size(); ++i)
        {
            // get size of message
            MPI::Status status;
            MPI::COMM_WORLD.Probe (MPI::ANY_SOURCE, 1000, status);
            int source = status.Get_source();
            int count  = status.Get_count(MPI::DOUBLE);
            //printf("Proc # %d is receiving from # %d\n",my_id,source);
            data.Resize (count);
            MPI::COMM_WORLD.Recv (data.GetPtr(), count, MPI::DOUBLE, source, 1000);

            // unpack data
            for (size_t j=0; j<data.Size(); j+=12)
            {
                // id and position
                Vec3_t x(data[j+1],data[j+2],data[j+3]);
                int      id   = static_cast<int>(data[j]);
                int      cell = grid.FindCell (x);
                CellType type = grid.Id2Type[cell];
                bool has_part = id2part.count(id);
                if (type==Outer_t)
                { 
                    // particle moved outside, flag this
                    if (has_part) id2part[id]->CType = type;
                    continue;
                }
                else
                {
                    Particle * p;
                    if (has_part) p = id2part[id];
                    else
                    {
                        p = new Particle();
                        parts.Push (p);
                        id2part[id] = p;
                        //printf("Proc # %d got a new particle Id=%d from %d\n",my_id,id,source);
                    }
                    p->CType = type; // the CType of existent particle may be changed after movement in another processor => needs to be updated as well
                    p->Id = id;
                    p->X  = x;
                    p->V  = data[j+4], data[j+5], data[j+6];
                    p->Xp = data[j+7], data[j+8], data[j+9];
                    p->m  = data[j+10];
                    p->R  = data[j+11];
                    //printf("Proc # %d, stp=%d, t=%g: just set particle Id=%d sent by Proc # %d\n",my_id,stp_out,t,id,source);
                }
            }
        }

        // wait for all sent messages to arrive at their destinations. Needs to be after Recv
        for (size_t i=0; i<req.Size(); ++i) req[i].Wait();

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
    //printf("\nProc # %d, Ekin (before) = %16.8e\n", my_id, Ekin0);
    //printf("Proc # %d, Ekin (after ) = %16.8e\n\n", my_id, Ekin1);

    // write control file
    if (my_id==0)
    {
        String buf;
        buf.Printf("dem_test3_control.res");
        std::ofstream of(buf.CStr(), std::ios::out);
        of << "fkey  " << "dem_test3" << "\n";
        of << "nout  " << stp_out+1   << "\n";
        of << "nx    " << N[0]        << "\n";
        of << "ny    " << N[1]        << "\n";
        of << "nz    " << N[2]        << "\n";
        of << "lxmi  " << L[0]        << "\n";
        of << "lxma  " << L[1]        << "\n";
        of << "lymi  " << L[2]        << "\n";
        of << "lyma  " << L[3]        << "\n";
        of << "lzmi  " << L[4]        << "\n";
        of << "lzma  " << L[5]        << "\n";
        of << "dt    " << dt          << "\n";
        of << "proc  " << N[0]*N[1]*N[2] << "   ";
        for (int i=0; i<N[0]*N[1]*N[2]; ++i) { of<<grid.Id2Proc[i]<<" "; } of<<"\n";
        of.close();
        cout << "File <" << buf << "> written" << endl;
    }

    // end
    MPI::Finalize();
    return 0;
}
MECHSYS_CATCH
