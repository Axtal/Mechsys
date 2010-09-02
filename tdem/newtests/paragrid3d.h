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

#ifndef MECHSYS_PARAGRID3D_H
#define MECHSYS_PARAGRID3D_H

// MechSys
#include <mechsys/util/maps.h>
#include <mechsys/linalg/matvec.h>
#include <mechsys/mesh/structured.h>

class ParaGrid3D
{
public:
    // enum
    enum CellType { Inner_t=-1, BryIn_t=-2, BryOut_t=-3, Outer_t=-4 };

    // Constructor & Destructor
     ParaGrid3D (Array<int> const & N, Array<double> & L, char const * FileKey=NULL);
    ~ParaGrid3D () { if (_mesh!=NULL) delete _mesh; }

    // Methods
    int  FindCellKey   (Vec3_t const & X, double R) const;
    void FindCellsKeys (Vec3_t const & X, double R, int Keys[8]) const;
    int  Key2Tag       (int CellKey) const { return _key2tag[CellKey]; }

    // Data
    Array<int>    N;
    Array<double> L;
    Vec3_t        D;

private:
    Mesh::Structured * _mesh;
    Array<int>         _key2tag;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline ParaGrid3D::ParaGrid3D (Array<int> const & TheN, Array<double> & TheL, char const * FileKey)
    : N(TheN), L(TheL)
{
    // cell size
    D = (L[1]-L[0])/static_cast<double>(N[0]),
        (L[3]-L[2])/static_cast<double>(N[1]),
        (L[5]-L[4])/static_cast<double>(N[2]);

    // proc id
    int my_id  = MPI::COMM_WORLD.Get_rank();
    int nprocs = MPI::COMM_WORLD.Get_size();

    // generate grid
    Array<Mesh::Block> blks(1);
    blks[0].Set (3, Outer_t, 8,
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
    _mesh = new Mesh::Structured(3);
    _mesh->Generate   (blks);
    _mesh->PartDomain (nprocs);
    if (FileKey!=NULL) _mesh->WriteVTU (FileKey);

    // cell key to tag
    _key2tag.Resize (_mesh->Verts.Size());

    // find what cells are at the boundaries of this domain
    Array<Mesh::Cell*> cells_bry_out;
    //std::map<int, Array<int> > cell2bryprocs;
    for (size_t i=0; i<_mesh->Cells.Size(); ++i)
    {
        _mesh->Cells[i]->Tag = Outer_t; // outer cell
        if (_mesh->Cells[i]->PartID==my_id)
        {
            _mesh->Cells[i]->Tag = Inner_t; // inner cell
            for (size_t j=0; j<_mesh->Cells[i]->V.Size(); ++j)
            {
                if (_mesh->Cells[i]->V[j]->PartIDs.Size()>1) // is a boundary cell
                {
                    _mesh->Cells[i]->Tag = BryIn_t; // inner boundary cell
                    Array<Mesh::Share> const & sha = _mesh->Cells[i]->V[j]->Shares;
                    for (size_t k=0; k<sha.Size(); ++k)
                    {
                        if (sha[k].C!=_mesh->Cells[i] && sha[k].C->PartID!=my_id) cells_bry_out.XPush (sha[k].C);
                    }
                }
            }
        }
    }
    for (size_t i=0; i<cells_bry_out.Size(); ++i) cells_bry_out[i]->Tag = BryOut_t; // outer boundary cell
}

inline int ParaGrid3D::FindCellKey (Vec3_t const & X, double R) const
{
    // cell touched by (X,R)
    int I = static_cast<int>((X(0)-L[0])/D(0));
    int J = static_cast<int>((X(1)-L[2])/D(1));
    int K = static_cast<int>((X(2)-L[4])/D(2));
    int n = I + J*(N[0]+1) + K*(N[0]+1)*(N[1]+1);

    // check
    if (n<0) throw new Fatal("ParaGrid3D::FindCellId:: __internal_error: cell id (%d) cannot be negative",n);
    return n;
}

inline void ParaGrid3D::FindCellsKeys (Vec3_t const & X, double R, int Keys[8]) const
{
    // coordinates of corners
    double xa=X(0)-R;   double xb=X(0)+R;
    double ya=X(1)-R;   double yb=X(1)+R;
    double za=X(2)-R;   double zb=X(2)+R;

    // check
    if (xa<L[0] || xb>L[1] ||
        ya<L[2] || yb>L[3] ||
        za<L[4] || zb>L[5]) throw new Fatal("ParaGrid3D::FindCellsIds: Corner of cube centered at (%g,%g,%g) with inner radius=%g is outside limits of domain (N=(%d,%d,%d), L=[x:(%g,%g) y:(%g,%g) z:(%g,%g)])", X(0),X(1),X(2),R, N[0],N[1],N[2], L[0],L[1],L[2],L[3],L[4],L[5]);

    // corners
    double C[8][3] = {{xa,  ya,  za},
                      {xb,  ya,  za},
                      {xa,  yb,  za},
                      {xb,  yb,  za},
                      {xa,  ya,  zb},
                      {xb,  ya,  zb},
                      {xa,  yb,  zb},
                      {xb,  yb,  zb}};

    // find cells ids
    for (size_t i=0; i<8; ++i)
    {
        // cell touched by corner
        int I = static_cast<int>((C[i][0]-L[0])/D(0));
        int J = static_cast<int>((C[i][1]-L[2])/D(1));
        int K = static_cast<int>((C[i][2]-L[4])/D(2));
        Keys[i] = I + J*(N[0]+1) + K*(N[0]+1)*(N[1]+1);

        // check
        if (Keys[i]<0) throw new Fatal("ParaGrid3D::FindCellsIds:: __internal_error: cell id (%d) cannot be negative",Keys[i]);
    }
}



/*
struct Box
{
    int ID;
    int Tag;
    int PartID;
};

class Grid3D
{
public:
    Grid3D (Array<int> const & TheN, Array<double> & TheL)
        : N(TheN), L(TheL)
    {
        Array<int> Xadj(0, true);
        Boxes.Resize (N[0]*N[1]*N[2]);
        size_t m = 0;
        for (size_t k=0; k<N[2]; ++k)
        for (size_t j=0; j<N[1]; ++j)
        for (size_t i=0; i<N[0]; ++i)
        {
            int n = i + j*(N[0]+1) + k*(N[0]+1)*(N[1]+1);
            Boxes[m].ID = n;
            int back   = (i-1) +  j   *(N[0]+1) +  k   *(N[0]+1)*(N[1]+1);
            int front  = (i+1) +  j   *(N[0]+1) +  k   *(N[0]+1)*(N[1]+1);
            int left   =  i    + (j-1)*(N[0]+1) +  k   *(N[0]+1)*(N[1]+1);
            int right  =  i    + (j+1)*(N[0]+1) +  k   *(N[0]+1)*(N[1]+1);
            int bottom =  i    +  j   *(N[0]+1) + (k-1)*(N[0]+1)*(N[1]+1);
            int top    =  i    +  j   *(N[0]+1) + (k+1)*(N[0]+1)*(N[1]+1);
            m++;
        }
    }
    Array<int>    N;
    Array<double> L;
    Array<Box>    Boxes;
};
*/

#endif // MECHSYS_PARAGRID3D_H
