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

    bool       Active;   ///< Is particle Active, or deleted from simulation ?
    int        Id;       ///< Global Id of particle
    Vec3_t     X,V,F,Xp; ///< Position, velocity, force, previous position (Verlet)
    double     m,R;      ///< Mass and radius
    Array<int> Boxes;    ///< Boxes that this particle touches
};

inline void FindBoxNumbers (Array<int> const & N, Array<double> const & L, Vec3_t const & X, double R, Array<int> & Nums)
{
    // coordinates of corners
    double xa=X(0)-R;   double xb=X(0)+R;
    double ya=X(1)-R;   double yb=X(1)+R;
    double za=X(2)-R;   double zb=X(2)+R;

    // check
    if (xa<L[0] || xb>L[1] ||
        ya<L[2] || yb>L[3] ||
        za<L[4] || zb>L[5]) throw new Fatal("FindBoxNumbers: Corner of cube centered at (%g,%g,%g) with inner radius=%g is outside limits of domain (N=(%d,%d,%d), L=[x:(%g,%g) y:(%g,%g) z:(%g,%g)])", X(0),X(1),X(2),R, N[0],N[1],N[2], L[0],L[1],L[2],L[3],L[4],L[5]);

    // corners
    double C[8][3] = {{xa,  ya,  za},
                      {xb,  ya,  za},
                      {xa,  yb,  za},
                      {xb,  yb,  za},
                      {xa,  ya,  zb},
                      {xb,  ya,  zb},
                      {xa,  yb,  zb},
                      {xb,  yb,  zb}};

    // box sizes
    double D[3] = {(L[1]-L[0])/static_cast<double>(N[0]),
                   (L[3]-L[2])/static_cast<double>(N[1]),
                   (L[5]-L[4])/static_cast<double>(N[2])};

    // find box numbers
    for (size_t i=0; i<8; ++i)
    {
        // box number touched by corner
        int I = static_cast<int>((C[i][0]-L[0])/D[0]);
        int J = static_cast<int>((C[i][1]-L[2])/D[1]);
        int K = static_cast<int>((C[i][2]-L[4])/D[2]);
        Nums[i] = I + J*(N[0]+1) + K*(N[0]+1)*(N[1]+1);

        // check
        if (Nums[i]<0) throw new Fatal("Particle::Boxes:: __internal_error: box id (%d) cannot be negative",Nums[i]);
    }
}

inline int FindBoxNumber (Array<int> const & N, Array<double> const & L, Vec3_t const & X, double R)
{
    // box sizes
    double D[3] = {(L[1]-L[0])/static_cast<double>(N[0]),
                   (L[3]-L[2])/static_cast<double>(N[1]),
                   (L[5]-L[4])/static_cast<double>(N[2])};

    // box touched by (X,R)
    int I = static_cast<int>((X[0]-L[0])/D[0]);
    int J = static_cast<int>((X[1]-L[2])/D[1]);
    int K = static_cast<int>((X[2]-L[4])/D[2]);
    int n = I + J*(N[0]+1) + K*(N[0]+1)*(N[1]+1);

    // check
    if (n<0) throw new Fatal("FindBoxNumber:: __internal_error: box id (%d) cannot be negative",n);
    return n;
}

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

typedef google::dense_hash_map<int,Array<Particle*> > Box2Part_t;
typedef google::dense_hash_map<Particle*,Particle*>   Neighbours_t;
//typedef std::map<int,Array<Particle*> > Box2Part_t;
//typedef std::map<Particle*,Particle*>   Neighbours_t;

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

    // generate grid
    Array<Mesh::Block> blks(1);
    blks[0].Set (3, -1, 8,
            0., L[0],L[2],L[4],
            0., L[1],L[2],L[4],
            0., L[1],L[3],L[4],
            0., L[0],L[3],L[4],
            0., L[0],L[2],L[5],
            0., L[1],L[2],L[5],
            0., L[1],L[3],L[5],
            0., L[0],L[3],L[5], 0.,0.,0.,0.,0.,0.);
    blks[0].SetNx (N[0]);
    blks[0].SetNy (N[1]);
    blks[0].SetNz (N[2]);
    Mesh::Structured mesh(3);
    mesh.Generate   (blks);
    mesh.PartDomain (nprocs);
    mesh.WriteVTU   (fkey.CStr());

    // find what cells are at the boundaries of this domain
    Array<int> bry_cells;
    //std::map<int, Array<int> > cell2bryprocs;
    for (size_t i=0; i<mesh.Cells.Size(); ++i)
    {
        if (mesh.Cells[i]->PartID==my_id)
        {
            bool is_bry = false;
            for (size_t j=0; j<mesh.Cells[i]->V.Size(); ++j)
            {
                if (mesh.Cells[i]->V[j]->PartIDs.Size()>1)
                {
                    is_bry = true;
                    break;
                }
            }
            if (is_bry) bry_cells.Push(i);
        }
    }
    
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
        int n = FindBoxNumber (N, L, x, ra[id]);
        if (mesh.Cells[n]->PartID==my_id)
        {
            Vec3_t v(vx[id],vy[id],vz[id]);
            Vec3_t xp = x - dt*v;
            parts.Push (new Particle());
            Particle & p = (*parts.Last());
            p.Active = true;
            p.Id     = id;
            p.X      = x;
            p.Xp     = xp;
            p.V      = v;
            p.m      = mass;
            p.R      = ra[id];
            Ekin0   += 0.5*mass*dot(p.V,p.V);
            //Id2Part[id] = &p;
        }
        //else Id2Part[id] = NULL;
    }

    // solve
    double tf      = 1.0;
    double dtout   = 0.1;
    double tout    = 0.1;
    int    stp_out = 0;
    for (double t=0.0; t<tf; t+=dt)
    {
        // maps: box2part, part2box
        Box2Part_t box2part;          // map: box => particles in/crossed box
        box2part.set_empty_key (-1);
        //std::map<Particle*,Array<int> > part2box;
        for (size_t i=0; i<parts.Size(); ++i)
        {
            if (!parts[i]->Active) continue; // skip non-active particles

            int n = FindBoxNumber (N, L, parts[i]->X, parts[i]->R);
            if (mesh.Cells[n]->PartID==my_id)
            {
                parts[i]->Start       ();
                parts[i]->Boxes.Clear ();

                Array<int> nums(8); // this array needs to be resized externally
                FindBoxNumbers (N, L, parts[i]->X, parts[i]->R, nums);
                //bool found = false; // TODO: use nums instead of n to find whether this particle belongs to this processor or no
                for (size_t j=0; j<8; ++j)
                {
                    Box2Part_t::iterator it = box2part.find(nums[j]);
                    if (it==box2part.end()) box2part[nums[j]].Push (parts[i]);
                    parts[i]->Boxes.Push (nums[j]);
                }
            }
            else parts[i]->Active = false; // remove this particle from this processor
        }

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

                // push particle into parts
                for (size_t i=0; i<data.Size(); i+=8)
                {
                    // check if particle belongs to any box in this processor
                    Vec3_t x(data[i+1],data[i+2],data[i+3]);
                    double r(data[i+4]);
                    int n = FindBoxNumber (N, L, x, r);
                    if (mesh.Cells[n]->PartID==my_id) // need to allocate particle
                    {
                        // allocate particle
                        int id = static_cast<int>(data[i]);
                        Vec3_t v(data[i+5],data[i+6],data[i+7]);
                        Vec3_t xp = x - dt*v; // need to receive xp as well ?
                        parts.Push (new Particle());
                        parts.Last()->Id = id;
                        parts.Last()->X  = x;
                        parts.Last()->Xp = xp;
                        parts.Last()->V  = v;
                        parts.Last()->m  = mass;
                        parts.Last()->R  = r;

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
                }
            }
        }

        // find possible contacts
        Neighbours_t neighs;
        neighs.set_empty_key(NULL);
        for (size_t i=0; i<parts.Size(); ++i)
        {
            if (!parts[i]->Active) continue; // skip the non-active particles

            // particle
            Particle & pa = (*parts[i]);

            // loop particle's boxes
            for (size_t k=0; k<pa.Boxes.Size(); ++k)
            {
                // box => particles map
                Box2Part_t::const_iterator it = box2part.find(pa.Boxes[k]);

                // if there are particles in particle's box
                if (it!=box2part.end())
                {
                    // particles in particle's box
                    Array<Particle*> const & parts_in_box = it->second;

                    // if there are more than one particle (itself) in this box
                    if (parts_in_box.Size()<1) throw new Fatal("__internal_error__: All boxes should have at least one particle");
                    if (parts_in_box.Size()>1)
                    {
                        // loop particles in particle's box
                        for (size_t j=0; j<parts_in_box.Size(); ++j)
                        {
                            // neighbour particle
                            Particle & pb = (*parts_in_box[j]);

                            //cout << "Comparing: " << pa.Id << " with " << pb.Id << "  boxid =" << (*boxid) << endl;

                            // if pb is not pa
                            if (pb.Id != pa.Id)
                            {
                                // halo overlapping
                                double del = pa.R + pb.R - norm(pb.X-pa.X);

                                // there is overlapping
                                if (del>MACH_EPS)
                                {
                                    //cout << "Neighbours: " << pa.Id << " => " << pb.Id << ", del = " << del << endl;
                                    if (pa.Id < pb.Id)
                                    {
                                        neighs[&pa] = &pb;
                                    }
                                    else
                                    {
                                        neighs[&pb] = &pa;
                                    }
                                }
                            }
                        }
                    }
                }
                else throw new Fatal("__internal_error__: There should be at least one particle in particle's box: the particle itself");
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
            if (parts[i]->Active) parts[i]->Move(dt);
        }

        // output
        if (t>=tout)
        {
            Table tab;
            tab.SetZero ("id xc yc zc ra vx vy vz",parts.Size());
            for (size_t i=0; i<parts.Size(); ++i)
            {
                if (parts[i]->Active)
                {
                    tab("id",i) = parts[i]->Id;
                    tab("xc",i) = parts[i]->X(0);
                    tab("yc",i) = parts[i]->X(1);
                    tab("zc",i) = parts[i]->X(2);
                    tab("ra",i) = parts[i]->R;
                    tab("vx",i) = parts[i]->V(0);
                    tab("vy",i) = parts[i]->V(1);
                    tab("vz",i) = parts[i]->V(2);
                }
            }
            String buf;
            buf.Printf("%s_%08d.res",fkey.CStr(),stp_out);
            tab.Write(buf.CStr());
            tout += dtout;
            stp_out++;
        }
    }

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
