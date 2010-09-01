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
    // static
    static int    N[3];   ///< Domain: N:{nx,ny,nz}
    static double L[6];   ///< Domain: L:{Xmin,Xmax, Ymin,Ymax, Zmin,Zmax}
    static double D[3];   ///< Domain: D:{dx,dy,dz}
    static bool   LimSet; ///< Limits (N,L,D) set ?
    static void   SetLimits (Array<int> const & ArrN, Array<double> const & ArrL);

    Particle ()
    {
        Boxes.set_empty_key (-1);
    }

    void Move (double dt)
    {
        bool verlet=false;
        if (verlet)
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

    void Start()
    {
        F = 0.0,0.0,0.0;
    }

    int                          Id;
    Vec3_t                       X,V,F,Xp;
    google::dense_hash_set<int>  Boxes; ///< Boxes that this particle touches
    double                       m,R;

    double Active;
};

int    Particle::N[3]   = {0,0,0};
double Particle::L[6]   = {0,0, 0,0, 0,0};
double Particle::D[3]   = {0,0,0};
bool   Particle::LimSet = false;

inline void Particle::SetLimits (Array<int> const & ArrN, Array<double> const & ArrL)
{
    N[0]=ArrN[0];  N[1]=ArrN[1];  N[2]=ArrN[2];
    L[0]=ArrL[0];  L[1]=ArrL[1];  L[2]=ArrL[2];  L[3]=ArrL[3];  L[4]=ArrL[4];  L[5]=ArrL[5];

    D[0] = (L[1]-L[0])/(N[0]-1.0);
    D[1] = (L[3]-L[2])/(N[1]-1.0);
    D[2] = (L[5]-L[4])/(N[2]-1.0);

    LimSet = true;
}


inline void FindCellNumbers (Array<int> const & N, Array<double> const & L, Vec3_t const & X, double R, Array<int> & CellNums)
{
    // coordinates of corners
    double xa=X(0)-R;   double xb=X(0)+R;
    double ya=X(1)-R;   double yb=X(1)+R;
    double za=X(2)-R;   double zb=X(2)+R;

    // check
    double C[8][3];
    // corners
    C[0][0]=xa;  C[0][1]=ya;  C[0][2]=za;
    C[1][0]=xb;  C[1][1]=ya;  C[1][2]=za;
    C[2][0]=xa;  C[2][1]=yb;  C[2][2]=za;
    C[3][0]=xb;  C[3][1]=yb;  C[3][2]=za;
    C[4][0]=xa;  C[4][1]=ya;  C[4][2]=zb;
    C[5][0]=xb;  C[5][1]=ya;  C[5][2]=zb;
    C[6][0]=xa;  C[6][1]=yb;  C[6][2]=zb;
    C[7][0]=xb;  C[7][1]=yb;  C[7][2]=zb;

    double D[3];
    D[0] = (L[1]-L[0])/(N[0]-1.0);
    D[1] = (L[3]-L[2])/(N[1]-1.0);
    D[2] = (L[5]-L[4])/(N[2]-1.0);

    // set
    //cout << "Particle: " << Geo.Id << endl << " boxes: ";
    for (size_t i=0; i<8; ++i)
    {
        //cout << (C[i][0]-L[0])/D[0] << "   ";// << " " << (C[i][1]-L[2]) << " " << (C[i][2]-L[4]) << "   ";

        // box number (corner)
        int I = static_cast<int>((C[i][0]-L[0])/D[0]);
        int J = static_cast<int>((C[i][1]-L[2])/D[1]);
        int K = static_cast<int>((C[i][2]-L[4])/D[2]);
        int n = I + J*N[0] + K*N[0]*N[1];

        CellNums[i] = n;

        // set maps
        //cout << n << " ";
        if (n<0) throw new Fatal("Particle::Boxes:: __internal_error: box id (%d) cannot be negative",n);
    }
    //cout << endl << endl;
}

inline int FindCellNumber (Array<int> const & N, Array<double> const & L, Vec3_t const & X, double R)
{
    double D[3];
    D[0] = (L[1]-L[0])/(N[0]-1.0);
    D[1] = (L[3]-L[2])/(N[1]-1.0);
    D[2] = (L[5]-L[4])/(N[2]-1.0);

    int I = static_cast<int>((X[0]-L[0])/D[0]);
    int J = static_cast<int>((X[1]-L[2])/D[1]);
    int K = static_cast<int>((X[2]-L[4])/D[2]);
    int n = I + J*N[0] + K*N[0]*N[1];

    if (n<0) throw new Fatal("Particle::Boxes:: __internal_error: box id (%d) cannot be negative",n);

    return n;
}

inline void CalcForce (Particle & a, Particle & b)
{
    double Kn = 1000.0;
    Vec3_t nor = b.X - a.X;
    double dist = norm(nor);
    double delta = a.R + b.R - dist;
    if (delta > 0.0)
    {
        Vec3_t F = Kn*delta*nor/dist;
        a.F -= F;
        b.F += F;
    }
}

typedef google::dense_hash_set<int>                   BoxSet_t;
typedef google::dense_hash_map<int,Array<Particle*> > Box2Part_t;
typedef google::dense_hash_map<Particle*,Particle*>   Neighbours_t;

int main(int argc, char **argv) try
{
    MECHSYS_CATCH_PARALLEL = true;
    MECHSYS_MPI_INIT
    int my_id  = MPI::COMM_WORLD.Get_rank();
    int nprocs = MPI::COMM_WORLD.Get_size();

    // number:  nx ny nz
    Array<int> N(11, 11, 2);
    //Array<int> N(2, 2, 2);
    double dt = 0.001;
    if (argc>1) N[0]  = atoi(argv[1]);
    if (argc>2) N[1]  = atoi(argv[2]);
    if (argc>3) N[2]  = atoi(argv[3]);
    if (argc>4) dt    = atof(argv[4]);
    if (N[0]<2) throw new Fatal("nx must be greater than 1");
    if (N[1]<2) throw new Fatal("ny must be greater than 1");
    if (N[2]<2) throw new Fatal("nz must be greater than 1");

    // limits
    Array<double> L(6);
    //     0   1     2   3     4   5
    //   xmi xma   ymi yma   zmi zma
    //L =    0,  1.01,    0,  1.01,    0,0.101;
    L =  -2,  2,    -2,  2,    0,0.1;

    // set particle with limits
    Particle::SetLimits (N,L);

    String fkey;
    fkey.Printf("dem_test2_proc_%d",my_id);
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
    blks[0].SetNx (N[0]-1);
    blks[0].SetNy (N[1]-1);
    blks[0].SetNz (N[2]-1);
    Mesh::Structured mesh(3);
    mesh.Generate   (blks);
    mesh.PartDomain (nprocs);
    mesh.WriteVTU   (fkey.CStr());

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
    Array<double> const & X = tab("Xc");
    Array<double> const & Y = tab("Yc");
    Array<double> const & Z = tab("Zc");
    Array<double> const & R = tab("R");
    Array<Particle*> Part;

    //std::map<int,Particle*> Id2Part; // global Id to particle

    double Ekin = 0.0;
    for (size_t id=0; id<X.Size(); ++id)
    {
        double vx = static_cast<double>(rand())/static_cast<double>(RAND_MAX)-0.5;
        //double vy = static_cast<double>(rand())/static_cast<double>(RAND_MAX)-0.5;
        //double vx = 0.0;
        double vy = 0.0;
        double vz = 0.0;
        Vec3_t v(vx,vy,vz);
        Vec3_t x(X[id],Y[id],Z[id]);
        Vec3_t xp = x-dt*v;

        int n = FindCellNumber (N, L, x, R[id]);
        if (mesh.Cells[n]->PartID==my_id)
        {
            Part.Push (new Particle());
            Part.Last()->Id = id;
            Part.Last()->X  = x;
            Part.Last()->Xp = xp;
            Part.Last()->V  = v;
            Part.Last()->m  = 1.0;
            Part.Last()->R  = R[id];
            Ekin += 0.5*Part.Last()->m*dot(Part.Last()->V,Part.Last()->V);
            //Id2Part[id] = Part.Last();
        }
        //else Id2Part[id] = NULL;
    }
    printf("\nProc # %d, Ekin (before) = %16.8e\n", my_id,Ekin);

    double Tf      = 1.0;
    double dtout   = 0.1;
    double tout    = 0.1;
    int    stp_out = 0;
    for (double t = 0.0 ; t<Tf;t+=dt)
    {
        Box2Part_t box2part; // map: box => particles in/crossed box
        box2part.set_empty_key (-1);
        //std::map<Particle*,Array<int> > part2box;
        for (size_t i=0; i<Part.Size(); ++i)
        {
            if (!Part[i]->Active) continue; // skip non-active particles

            int n = FindCellNumber (N, L, Part[i]->X, Part[i]->R);
            if (mesh.Cells[n]->PartID==my_id)
            {
                Part[i]->Start    ();
                Part[i]->Boxes.clear ();

                Array<int> nums(8); // this array needs to be resize externally
                FindCellNumbers (N,L, Part[i]->X, Part[i]->R, nums);
                //bool found = false; // TODO: use nums instead of n to find whether this particle belongs to this processor or no
                for (size_t j=0; j<8; ++j)
                {
                    Box2Part_t::iterator it = box2part.find(nums[j]);
                    if (it==box2part.end()) box2part[nums[j]].Push (Part[i]);
                    Part[i]->Boxes.insert (nums[j]);
                }
            }
            else Part[i]->Active = false; // remove this particle from this processor
        }

        //Array<double> x,y,z,r;
        Array<double> xyzrid;
        for (size_t i=0; i<bry_cells.Size(); ++i)
        {
            int n = bry_cells[i];
            Box2Part_t::const_iterator it = box2part.find(n);
            if (it!=box2part.end()) // box has particles
            {
                Array<Particle*> const & parts = it->second; // particles in the cell
                for (size_t j=0; j<parts.Size(); ++j)
                {
                    xyzrid.Push(parts[j]->X(0));
                    xyzrid.Push(parts[j]->X(1));
                    xyzrid.Push(parts[j]->X(2));
                    xyzrid.Push(parts[j]->R);
                    xyzrid.Push(parts[j]->Id);
                    //printf("Proc # %d, part Id = %d\n",my_id,parts[j]->Id);
                }
            }
        }

        // broadcast
        for (int i=0; i<nprocs; ++i)
        {
            if (i!=my_id)
            {
                MPI::Request req_send = MPI::COMM_WORLD.Isend (xyzrid.GetPtr(), xyzrid.Size(), MPI::DOUBLE, i, 1000);
                req_send.Wait ();
            }
        }

        // receive messages from everyone
        MPI::Status status;
        for (int i=0; i<nprocs; ++i)
        {
            if (i!=my_id)
            {
                MPI::COMM_WORLD.Probe (MPI::ANY_SOURCE, 1000, status);
                int source = status.Get_source();
                int count  = status.Get_count(MPI::DOUBLE);
                xyzrid.Resize (count);
                MPI::COMM_WORLD.Recv (xyzrid.GetPtr(), count, MPI::DOUBLE, source, 1000);
            }
        }

        // push particle into Part
        for (size_t i=0; i<xyzrid.Size(); i+=5)
        {
            Vec3_t x(xyzrid[i],xyzrid[i+1],xyzrid[i+2]);
            double r(xyzrid[i+3]);
            int id = static_cast<int>(xyzrid[i+4]);
            int n  = FindCellNumber (N, L, x, r);
            if (mesh.Cells[n]->PartID==my_id) // need to allocate particle
            {
                Part.Push (new Particle());
                Part.Last()->Id = id;
                Part.Last()->X  = x;
                Part.Last()->Xp = x; // need to receive xp as well
                Part.Last()->V  = 0.0; // need to receive v as well
                Part.Last()->m  = 1.0;
                Part.Last()->R  = r;

                // update maps
                Array<int> nums(8); // this array needs to be resize externally
                FindCellNumbers (N,L, Part.Last()->X, Part.Last()->R, nums);
                for (size_t j=0; j<8; ++j)
                {
                    Box2Part_t::iterator it = box2part.find(nums[j]);
                    if (it==box2part.end()) box2part[nums[j]].Push (Part.Last());
                    Part.Last()->Boxes.insert (nums[j]);
                }
            }
        }

        Neighbours_t neighs;
        neighs.set_empty_key(NULL);
        for (size_t i=0; i<Part.Size(); ++i)
        {
            if (!Part[i]->Active) continue; // skip the non-active particles

            // particle
            Particle & pa = (*Part[i]);

            // loop particle's boxes
            for (BoxSet_t::const_iterator boxid=pa.Boxes.begin(); boxid!=pa.Boxes.end(); ++boxid)
            {
                // box => particles map
                Box2Part_t::const_iterator it = box2part.find((*boxid));

                // if there are particles in particle's box
                if (it!=box2part.end())
                {
                    // particles in particle's box
                    Array<Particle*> const & Part_in_box = it->second;

                    // if there are more than one particle (itself) in this box
                    if (Part_in_box.Size()<1) throw new Fatal("__internal_error__: All boxes should have at least one particle");
                    if (Part_in_box.Size()>1)
                    {
                        // loop particles in particle's box
                        //for (PartSet_t::const_iterator neigh=Part_in_box.begin(); neigh!=Part_in_box.end(); ++neigh)
                        for (size_t j=0; j<Part_in_box.Size(); ++j)
                        {
                            // neighbour particle
                            Particle & pb = (*Part_in_box[j]);//(*(*neigh));

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

        for (Neighbours_t::const_iterator r=neighs.begin(); r!=neighs.end(); ++r)
        {
            Particle & pa = (*r->first);
            Particle & pb = (*r->second);
            CalcForce (pa,pb);
        }

        for (size_t i=0; i<Part.Size(); ++i)
        {
            if (Part[i]->Active) Part[i]->Move(dt);
        }

        if (t>=tout)
        {
            Table tab;
            tab.SetZero ("Id Xc Yc Zc R",Part.Size());
            for (size_t i=0; i<Part.Size(); ++i)
            {
                if (Part[i]->Active)
                {
                    tab("Id",i) = Part[i]->Id;
                    tab("Xc",i) = Part[i]->X(0);
                    tab("Yc",i) = Part[i]->X(1);
                    tab("Zc",i) = Part[i]->X(2);
                    tab("R", i) = Part[i]->R;
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
    Ekin = 0.0;
    for (size_t i=0; i<Part.Size(); ++i)
    {
        Ekin += 0.5*Part[i]->m*dot(Part[i]->V,Part[i]->V);
    }
    printf("Proc # %d, Ekin (after ) = %16.8e\n\n", my_id,Ekin);

    // end
    MPI::Finalize();
    return 0;
}
MECHSYS_CATCH
