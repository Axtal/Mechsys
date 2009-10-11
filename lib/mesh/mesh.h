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

#ifndef MECHSYS_MESH_H
#define MECHSYS_MESH_H

// STL
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>   // for DBL_EPSILON
#include <cstdarg>  // for va_list, va_start, va_end
#include <cstdlib>  // for atoi
#include <map>

// Boost
#include <boost/tuple/tuple_io.hpp>
#include <boost/tuple/tuple_comparison.hpp>

// MechSys
#include "util/array.h"
#include "util/fatal.h"
#include "util/util.h"
#include "util/numstreams.h"
#include "linalg/matvec.h"
#include "vtkcelltype.h"

namespace Mesh
{

#ifndef VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE_DEFINED
  #define VTU_NEWLINE(I,K,N,KMAX,OF) if (K>KMAX) { OF<<(I<N-1?"\n        ":"\n"); K=0; } else if (I==N-1) { OF<<"\n"; }
#endif

#define NOEDGE    {-1,-1,-1}
#define NOFACE    {-1,-1,-1,-1}
#define NOEDGES   {NOEDGE,NOEDGE,NOEDGE,NOEDGE}
#define NOEDGES3D {NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE}
#define NOFACES   {NOFACE,NOFACE,NOFACE,NOFACE,NOFACE,NOFACE}

size_t MaxNVerts2D               = 8;
size_t NVertsToNEdges2D[]        = {0,0,0,3,4,0,3,0,4};
size_t NVertsToNVertsPerEdge2D[] = {0,0,0,2,2,0,3,0,3};
int    NVertsToEdge2D[][4/*edges at most*/][3/*verts per edge at most*/]=
{
    NOEDGES,NOEDGES,NOEDGES,               // 0,1,2 verts
    {{0,1,-1},{1,2,-1},{2,0,-1},NOEDGE},   // 3 verts => TRIANGLE
    {{0,1,-1},{1,2,-1},{2,3,-1},{3,0,-1}}, // 4 verts => QUAD
    NOEDGES,                               // 5 verts
    {{0,1,3},{1,2,4},{2,0,5},NOEDGE},      // 6 verts => O2 TRIANGLE
    NOEDGES,                               // 7 verts
    {{0,1,4},{1,2,5},{2,3,6},{3,0,7}}      // 8 verts => O2 QUAD
};

size_t MaxNVerts3D               = 20;
size_t NVertsToNFaces3D[]        = {0,0,0,0, 4, 0,0,0, 6, 0, 4, 0,0,0,0,0,0,0,0,0, 6};
size_t NVertsToNVertsPerFace3D[] = {0,0,0,0, 3, 0,0,0, 4, 0, 6, 0,0,0,0,0,0,0,0,0, 8};
int    NVertsToFace3D[][6/*faces at most*/][4/*verts per face at most*/]=
{
    NOFACES,NOFACES,NOFACES,NOFACES,                                         //  0,1,2,3 verts
    {{0,3,2,-1},{0,1,3,-1},{0,2,1,-1},{1,2,3,-1},NOFACE,NOFACE},             //  4 verts => TETRA
    NOFACES,NOFACES,NOFACES,                                                 //  5,6,7 verts
    {{0,4,7,3},{1,2,6,5},{0,1,5,4},{2,3,7,6},{0,3,2,1},{4,5,6,7}},           //  8 verts => HEX
    NOFACES,                                                                 //  9 verts
    {{0,3,2,-1},{0,1,3,-1},{0,2,1,-1},{1,2,3,-1},NOFACE,NOFACE},             // 10 verts => O2 TETRA
    NOFACES,NOFACES,NOFACES,NOFACES,NOFACES,NOFACES,NOFACES,NOFACES,NOFACES, // 11,12,13,14,15,16,17,18,19 verts
    {{0,4,7,3},{1,2,6,5},{0,1,5,4},{2,3,7,6},{0,3,2,1},{4,5,6,7}}            // 20 verts => O2 HEX
};

size_t NVertsToNEdges3D[] = {0,0,0,0, 6, 0,0,0, 12, 0, 6, 0,0,0,0,0,0,0,0,0, 12};
int    NVertsToEdge3D[][12/*edges at most*/][3/*verts per edge at most*/]=
{
    NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,                                                                       //  0,1,2,3 verts
    {{0,1,-1},{1,2,-1},{2,0,-1},{0,3,-1},{1,3,-1},{2,3,-1},NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE},             //  4 verts => TETRA
    NOEDGES3D,NOEDGES3D,NOEDGES3D,                                                                                 //  5,6,7 verts
    {{0,1,-1},{1,2,-1},{2,3,-1},{3,0,-1},{4,5,-1},{5,6,-1},{6,7,-1},{7,4,-1},{0,4,-1},{1,5,-1},{2,6,-1},{3,7,-1}}, //  8 verts => HEX
    NOEDGES3D,                                                                                                     //  9 verts
    {{0,1,4},{1,2,5},{2,0,6},{0,3,7},{1,3,8},{2,3,9},NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE,NOEDGE},                   // 10 verts => O2 TETRA
    NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,NOEDGES3D,                     // 11,12,13,14,15,16,17,18,19 verts
    {{0,1,8},{1,2,9},{2,3,10},{3,0,11},{4,5,12},{5,6,13},{6,7,14},{7,4,15},{0,4,16},{1,5,17},{2,6,18},{3,7,19}}    // 20 verts => O2 HEX
};

#define BRYKEY(num_verts,idx_cell,idx_bry)                                      \
    int vert_a, vert_b, vert_c=-1;                                              \
    if (NDim==2)                                                                \
    {                                                                           \
        vert_a = Cells[idx_cell]->V[NVertsToEdge2D[num_verts][idx_bry][0]]->ID; \
        vert_b = Cells[idx_cell]->V[NVertsToEdge2D[num_verts][idx_bry][1]]->ID; \
        Util::Sort (vert_a,vert_b);                                             \
    }                                                                           \
    else                                                                        \
    {                                                                           \
        vert_a = Cells[idx_cell]->V[NVertsToFace3D[num_verts][idx_bry][0]]->ID; \
        vert_b = Cells[idx_cell]->V[NVertsToFace3D[num_verts][idx_bry][1]]->ID; \
        vert_c = Cells[idx_cell]->V[NVertsToFace3D[num_verts][idx_bry][2]]->ID; \
        Util::Sort (vert_a,vert_b,vert_c);                                      \
    }                                                                           \
    BryKey_t brykey(vert_a,vert_b,vert_c);

#undef NOEDGE
#undef NOFACE
#undef NOEDGES
#undef NOFACES

struct Cell; ///< Forward declaration due to the following definitions

struct Share
{
    Cell * C; ///< The cell
    int    N; ///< Local node index. Example: 2D=>0,1,2,3, 3D=>0,1,2,3,4,5,6,7
};

struct Vertex
{
    size_t       ID;     ///< ID
    int          Tag;    ///< Tag
    Vec3_t       C;      ///< X, Y, and Z coordinates
    Array<Share> Shares; ///< IDs of cells sharing this vertex
};

typedef std::map<int,int>                      BryTag_t;   ///< Map: edge/face ID => edge/face tag
typedef boost::tuple<int,int,int>              BryKey_t;   ///< Edge/face key = (left_node,right_node) for edges or (node0,node1,node2) for faces
typedef std::pair<int,Cell*>                   NeighDat_t; ///< Pair: local edge/face ID, neighbour cell
typedef std::map<BryKey_t, NeighDat_t>         Neighs_t;   ///< Map: edge/face key => neighbour data
typedef std::map<BryKey_t, Array<NeighDat_t> > BryCell_t;  ///< Map: edge/face key => cells sharing this boundary (edge/face)

struct Cell
{
    size_t         ID;      ///< ID
    int            Tag;     ///< Cell tag. Required for setting up of attributes, for example.
    Array<Vertex*> V;       ///< Connectivity
    BryTag_t       BryTags; ///< Boundary (edge/face) tags: map iEdgeFace => Tag
    Neighs_t       Neighs;  ///< Neighbours information: map edge/face key => neighbour data (local edge/face ID, pointer to neighbour Cell)
};

class Generic
{
public:
    // Constructor
    Generic (int TheNDim) : NDim(TheNDim) {}

    // Destructor
    virtual ~Generic () { Erase(); }

    // Set methods
    void ReadMesh   (char const * FileKey);                               ///< (.mesh) Erase old mesh and read mesh from python file
    void SetSize    (size_t NumVerts, size_t NumCells);                   ///< Erase old mesh and set number of vertices
    void SetVert    (int iVert, int Tag, double X, double Y, double Z=0); ///< Set vertex
    void SetCell    (int iCell, int Tag, size_t NVerts, ...);             ///< Set element ... => connectivity
    void SetBryTag  (int iCell, int iEdgeFace, int Tag);                  ///< Set element's edge or face tag
    void FindNeigh  ();                                                   ///< Find neighbours of each cell
    void GenO2Verts ();                                                   ///< Generate O2 (mid) vertices
    void Erase      ();                                                   ///< Erase current mesh (deallocate memory)

    // Method to extend mesh
    void AddLinCells (size_t NItems, ...); ///< Set linear cells given edge tags (edg_tag,edg_tag,...) or pair of vertices (v0,v1,new_elem_tag, v0,v1,new_elem_tag, ...) or both
    void BreakLin    (int VertexIdOrTag, Vec3_t const * dX=NULL, bool DoConnect=false, int ConnectionTag=0);

    // Methods
    void WriteVTU (char const * FileKey, int VolSurfOrBoth=0) const; ///< (.vtu) Write output file for ParaView. Vol=0, Surf=1, Both=2
    void WriteMPY (char const * FileKey, bool OnlyMesh=true,
                   char const * Extra=NULL) const;                   ///< (.mpy) Write Python script that calls mesh_drawing.py

    // Auxiliar methods
    void ThrowError (std::istringstream & iss, char const * Message) const; ///< Used in ReadMesh

    // Data
    int            NDim;      ///< Space dimension
    Array<Vertex*> Verts;     ///< Vertices
    Array<Cell*>   Cells;     ///< Cells
    Array<Vertex*> TgdVerts;  ///< Tagged Vertices
    Array<Cell*>   TgdCells;  ///< Tagged Cells
    BryCell_t      Bry2Cells; ///< map: bry (edge/face ids) => neighbours cells
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


std::ostream & operator<< (std::ostream & os, BryKey_t const & BK)
{
    os << "(" << boost::get<0>(BK) << "," << boost::get<1>(BK) << "," << boost::get<2>(BK) << ")";
    return os;
}

std::ostream & operator<< (std::ostream & os, Generic const & M)
{
    // verts
    os << "V=[";
    for (size_t i=0; i<M.Verts.Size(); ++i)
    {
        os << "[" << Util::_4  << M.Verts[i]->ID;
        os << "," << Util::_4  << M.Verts[i]->Tag;
        os << "," << Util::_8s << M.Verts[i]->C(0);
        os << "," << Util::_8s << M.Verts[i]->C(1);  if (M.NDim==3)
        os << "," << Util::_8s << M.Verts[i]->C(2);
        if (i==M.Verts.Size()-1) os << "]]\n";
        else                     os << "],\n   ";
    }

    // cells
    os << "\nC=[";
    for (size_t i=0; i<M.Cells.Size(); ++i)
    {
        os << "[" << Util::_4 << M.Cells[i]->ID;
        os << "," << Util::_4 << M.Cells[i]->Tag;
        os << ", [";
        for (size_t j=0; j<M.Cells[i]->V.Size(); ++j)
        {
            os << M.Cells[i]->V[j]->ID;
            if (j!=M.Cells[i]->V.Size()-1) os << ",";
        }
        os << "], {";
        size_t k = 0;
        for (BryTag_t::const_iterator p=M.Cells[i]->BryTags.begin(); p!=M.Cells[i]->BryTags.end(); ++p)
        {
            os << p->first << ":" << p->second;
            if (k!=M.Cells[i]->BryTags.size()-1) os << ",";
            k++;
        }
        os << "}";
        os << ", {";
        for (Neighs_t::const_iterator p=M.Cells[i]->Neighs.begin(); p!=M.Cells[i]->Neighs.end(); ++p)
        {
            if (p!=M.Cells[i]->Neighs.begin()) os << ",";
            os << "(" << boost::get<0>(p->first);
            os << "," << boost::get<1>(p->first);
            os << "," << boost::get<2>(p->first);
            os << "):";
            os << "[" << p->second.first << "," << p->second.second->ID << "]";
        }
        os << "}";
        if (i==M.Cells.Size()-1) os << "]]\n";
        else                     os << "],\n   ";
    }
    return os;
}

inline void Generic::ThrowError (std::istringstream & iss, char const * Message) const
{
    size_t num = 1+iss.tellg();
    String str(iss.str());
    str.resize (num);
    size_t pos = str.find(" ");
    while (pos!=String::npos)
    {
        str.replace (pos,1,"");
        pos = str.find(" ",pos+1);
    }
    throw new Fatal("Generic::ReadMesh: Mesh file format invalid\n    %s\n    %s",Message,str.CStr());
}

inline void Generic::ReadMesh (char const * FileKey)
{
    String fn(FileKey); fn.append(".mesh");
    std::fstream fil(fn.CStr(), std::ios::in);
    if (!fil.is_open()) throw new Fatal("Generic::ReadMesh: Could not open file < %s >",fn.CStr());

    String line, str, comma, colon, buf;
    size_t vid,  cid,  bid;  // vertex id,  cell id,  boundary id
    int    vtag, ctag, btag; // vertex tag, cell tag, boundary tag
    double x,y,z=0.0;
    bool   reading_verts = false;
    bool   reading_cells = false;
    bool   verts_read    = false;
    while (!fil.eof())
    {
        // read line
        std::getline (fil,line);

        // add spaces around: = [ ] ,
        size_t pos;
        pos=line.find("="); while (pos!=String::npos) { line.replace(pos,1," = "); pos=line.find("=",pos+3); }
        pos=line.find("["); while (pos!=String::npos) { line.replace(pos,1," [ "); pos=line.find("[",pos+3); }
        pos=line.find("]"); while (pos!=String::npos) { line.replace(pos,1," ] "); pos=line.find("]",pos+3); }
        pos=line.find(","); while (pos!=String::npos) { line.replace(pos,1," , "); pos=line.find(",",pos+3); }
        pos=line.find("{"); while (pos!=String::npos) { line.replace(pos,1," { "); pos=line.find("{",pos+3); }
        pos=line.find("}"); while (pos!=String::npos) { line.replace(pos,1," } "); pos=line.find("}",pos+3); }
        pos=line.find(":"); while (pos!=String::npos) { line.replace(pos,1," : "); pos=line.find(":",pos+3); }

        // parse
        std::istringstream iss(line);
        while (iss>>str)
        {
            if (str=="V")
            {
                iss>>str>>str>>str; //  = [ [
                reading_cells = false;
                reading_verts = true;
            }
            else if (str=="C")
            {
                iss>>str>>str>>str; //  = [ [
                reading_verts = false;
                reading_cells = true;
            }
            if (reading_verts)
            {
                if (str!="[") ThrowError (iss, "Verts: opening bracket '[' lacking");
                while (iss>>vid)
                {
                    // read vertex data
                    iss >> comma >> vtag >> comma >> x >> comma >> y;
                    if (NDim==3) iss >> comma >> z;
                    iss >> str;
                    if (str!="]") ThrowError (iss, "Verts: closing brackets ']' lacking");
                    iss >> str;
                    if (str=="]") { reading_verts = false; verts_read = true; }
                    else if (str!=",") ThrowError (iss, "Verts: last vertex: ',' or ']' lacking");

                    // allocate vertex
                    Verts.Push (new Vertex);
                    if (vid!=Verts.Size()-1) throw new Fatal("Generic::ReadMesh: Verts IDs must be sequential. Vertex ID (%d) must be equal to (%d)",vid,Verts.Size()-1);
                    Verts[vid]->ID  = vid;
                    Verts[vid]->Tag = vtag;
                    Verts[vid]->C   = x, y, z;

                    break;
                }
            }
            if (reading_cells)
            {
                if (!verts_read) throw new Fatal("Generic::ReadMesh: List with vertices data must come before list with cells data");
                if (str!="[") ThrowError (iss, "Cells: opening brackets '[' lacking");
                while (iss>>cid)
                {
                    // read cell data
                    iss >> comma >> ctag >> comma >> str;
                    if (str!="[") ThrowError (iss, "Cells: opening brackets of connectivity '[' lacking");

                    // allocate cell
                    Cells.Push (new Cell);
                    if (cid!=Cells.Size()-1) throw new Fatal("Generic::ReadMesh: Cells IDs must be sequential. Vertex ID (%d) must be equal to (%d)",cid,Cells.Size()-1);
                    Cells[cid]->ID  = cid;
                    Cells[cid]->Tag = ctag;

                    // read connectivity
                    while (iss>>vid)
                    {
                        Cells[cid]->V.Push (Verts[vid]);
                        iss >> str;
                        if (str=="]") break;
                        else if (str!=",") ThrowError (iss, "Cells: closing brackets of connectivity: ',' or ']' lacking");
                    }
                    iss >> comma;
                    iss >> str;
                    if (str!="{") ThrowError (iss, "Cells: opening braces of bry tags '{' lacking");

                    // check connectivity
                    if (NDim==2 && Cells[cid]->V.Size()>=3)
                    {
                        Vec3_t p0;  p0 = Cells[cid]->V[0]->C - Cells[cid]->V[1]->C;
                        Vec3_t p1;  p1 = Cells[cid]->V[0]->C - Cells[cid]->V[2]->C;
                        Vec3_t p2 = blitz::cross(p0,p1);
                        if (fabs(p2(0))>1.0e-10 || fabs(p2(1)>1.0e-10)) throw new Fatal("Generic::ReadMesh: All vertices of cells must be on the x-y plane");
                        if (p2(2)<0.0) throw new Fatal("Generic::ReadMesh: Numbering of vertices is incorrect (must be counter-clockwise)");
                    }

                    // read boundary tags
                    while (iss>>str)
                    {
                        if (str=="}") break;
                        else bid = atoi(str.CStr());
                        iss >> colon >> btag >> str;
                        Cells[cid]->BryTags[bid] = btag;
                        if (TgdCells.Find(Cells[cid])<0) TgdCells.Push (Cells[cid]);
                        if (str=="}") break;
                        else if (str!=",") ThrowError (iss, "Cells: closing braces of bry tags: ',' or '}' lacking");
                    }
                    iss >> str;
                    if (str!="]") ThrowError (iss, "Cells: closing brackets ']' lacking");
                    iss >> str;
                    if (str=="]") reading_cells = false;
                    else if (str!=",") ThrowError (iss, "Cells: last cell: ',' or ']' lacking");

                    break;
                }
            }
        }
    }
}

inline void Generic::SetSize (size_t NumVerts, size_t NumCells)
{
    // erase previous mesh
    Erase ();

    // nodes
    Verts.Resize    (NumVerts);
    Verts.SetValues (NULL);
    TgdVerts.Resize (0);

    // elements
    Cells.Resize    (NumCells);
    Cells.SetValues (NULL);
    TgdCells.Resize (0);
}

inline void Generic::SetVert (int i, int Tag, double X, double Y, double Z)
{
    // set Verts
    if (Verts[i]==NULL) Verts[i] = new Vertex;
    Verts[i]->ID  = i;
    Verts[i]->Tag = Tag;
    Verts[i]->C   = X, Y, Z;

    // set TgdVerts
    if (Tag<0) TgdVerts.Push (Verts[i]);
}

inline void Generic::SetCell (int i, int Tag, size_t NVerts, ...)
{
    // check
    if (NDim==2) { if (NVertsToVTKCell2D[NVerts]<0) throw new Fatal("Generic::SetCell: In two-dimensions (NDim=2), Cell=%d with Tag=%d has invalid NVerts=%d",i,Tag,NVerts); }
    else         { if (NVertsToVTKCell3D[NVerts]<0) throw new Fatal("Generic::SetCell: In three-dimensions (NDim=3), Cell=%d with Tag=%d has invalid NVerts=%d",i,Tag,NVerts); }

    // set elems
    if (Cells[i]==NULL) Cells[i] = new Cell;
    Cells[i]->ID  = i;
    Cells[i]->Tag = Tag;
    Cells[i]->V.Resize (NVerts);

    // set connectivity
    double    sum_z = 0.0;
    va_list   arg_list;
    va_start (arg_list, NVerts);
    for (size_t j=0; j<NVerts; ++j)
    {
        int ivert = va_arg(arg_list,int);
        if (Verts[ivert]==NULL) throw new Fatal("Generic::SetCell: Vert=%d of Cell %d, %d could not be found. Vertices must be set (with SetVert) before calling this method",ivert,i,Tag);
        Cells[i]->V[j] = Verts[ivert];
        Share sha = {Cells[i],j};
        Verts[ivert]->Shares.Push (sha);
        if (NDim==3) sum_z += Verts[ivert]->C(2);
    }
    va_end (arg_list);

    // check connectivity
    if (NDim==2 && Cells[i]->V.Size()>=3)
    {
        Vec3_t p0;  p0 = Cells[i]->V[0]->C - Cells[i]->V[1]->C;
        Vec3_t p1;  p1 = Cells[i]->V[0]->C - Cells[i]->V[2]->C;
        Vec3_t p2 = blitz::cross(p0,p1);
        if (fabs(p2(0))>1.0e-10 || fabs(p2(1)>1.0e-10)) throw new Fatal("Generic::SetCell: In 2D, all vertices of cells must lie on the x-y plane");
        if (p2(2)<0.0) throw new Fatal("Generic::SetCell: Numbering of vertices is incorrect (must be counter-clockwise)");
    }

    // check
    if (NDim==3 and fabs(sum_z)<1.0e-10) throw new Fatal("Generic::SetCell: In 3D, the sum of all z coordinates of a Cell (%d, %d) must not be zero",i,Tag);
}

inline void Generic::AddLinCells (size_t NItems, ...)
{
    if (NDim==3) throw new Fatal("Generic::AddLinCells: Method not available for 3D yet");

    va_list   arg_list;
    va_start (arg_list, NItems);
    for (size_t i=0; i<NItems; ++i)
    {
        int etag_or_v0 = va_arg(arg_list,int);
        if (etag_or_v0<0) // tag given
        {
            bool found = false;
            for (size_t j=0; j<TgdCells.Size(); ++j)
            {
                BryTag_t const & eftags = TgdCells[j]->BryTags;
                for (BryTag_t::const_iterator p=eftags.begin(); p!=eftags.end(); ++p)
                {
                    if (etag_or_v0==p->second)
                    {
                        int    eid    = p->first;
                        size_t nv     = TgdCells[j]->V.Size();
                        int    ivert0 = TgdCells[j]->V[NVertsToEdge2D[nv][eid][0]]->ID;
                        int    ivert1 = TgdCells[j]->V[NVertsToEdge2D[nv][eid][1]]->ID;
                        int    ivert2 = -1;
                        if (NVertsToNVertsPerEdge2D[nv]>2) ivert2 = TgdCells[j]->V[NVertsToEdge2D[nv][eid][2]]->ID;
                        if (ivert2<0)
                        {
                            Cells.Push (NULL);
                            SetCell (Cells.Size()-1, etag_or_v0, 2, ivert0, ivert1);
                        }
                        else
                        {
                            Cells.Push (NULL);  SetCell (Cells.Size()-1, etag_or_v0, 2, ivert0, ivert2);
                            Cells.Push (NULL);  SetCell (Cells.Size()-1, etag_or_v0, 2, ivert2, ivert1);
                        }
                        found  = true;
                        break;
                    }
                }
            }
            if (!found) throw new Fatal("Generic::SetLinCells: Could not find cell with edge with tag = %d",etag_or_v0);
        }
        else // vertices IDs given
        {
            int ivert1 = va_arg(arg_list,int);
            int tag    = va_arg(arg_list,int);
            Cells.Push (NULL);
            SetCell (Cells.Size()-1, tag, 2, etag_or_v0, ivert1);
        }
    }
    va_end (arg_list);
}

inline void Generic::BreakLin (int VertIdOrTag, Vec3_t const * dX, bool DoConnect, int ConnectionTag)
{
    // find vertex
    Vertex * vert = NULL;
    if (VertIdOrTag<0) // vertex tag given
    {
        bool found = false;
        for (size_t i=0; i<TgdVerts.Size(); ++i)
        {
            if (VertIdOrTag==TgdVerts[i]->Tag)
            {
                vert  = Verts[TgdVerts[i]->ID];
                found = true;
                break;
            }
        }
        if (!found) throw new Fatal("Mesh::AddPin: Could not find vertex with tag = %d", VertIdOrTag);
    }
    else vert = Verts[VertIdOrTag];

    // break lins
    Array<size_t> lins_to_disconnect;
    size_t idx_lin = 0;
    size_t nshares = vert->Shares.Size();
    for (size_t i=0; i<nshares; ++i)
    {
        Cell * cell = vert->Shares[i].C;
        if (cell->V.Size()==2) // cell sharing this vertex is a lin2 (linear cell with 2 vertices)
        {
            if (idx_lin>0) // do this for the second lin2 sharing this vertex on
            {
                // add new vertex
                Verts.Push (NULL);
                if (dX==NULL) SetVert (Verts.Size()-1, vert->Tag, vert->C(0),          vert->C(1),          vert->C(2));
                else          SetVert (Verts.Size()-1, vert->Tag, vert->C(0)+(*dX)(0), vert->C(1)+(*dX)(1), vert->C(2)+(*dX)(2));

                // set connectivity of lin2
                int idx = vert->Shares[i].N; // local vertex ID
                cell->V[idx] = Verts[Verts.Size()-1];
                Share sha = {cell,idx};
                Verts[Verts.Size()-1]->Shares.Push (sha);

                // set this lin2 to be disconnected from vert
                lins_to_disconnect.Push (i);

                // add new lin2
                if (DoConnect)
                {
                    AddLinCells (/*NItems*/1, /*v0*/vert->ID, /*v1*/Verts.Size()-1, ConnectionTag);
                    //Cells.Push (NULL);
                }
            }
            idx_lin++;
        }
    }
    if (idx_lin<2) throw new Fatal("Mesh::AddPin: Vertex %d (%d) must be connected to at least two linear elements",vert->ID,vert->Tag);

    // disconnect lins from vert
    Array<Share> old_shares(vert->Shares);
    vert->Shares.Resize (vert->Shares.Size()-lins_to_disconnect.Size());
    size_t k = 0;
    for (size_t i=0; i<old_shares.Size(); ++i)
    {
        if (lins_to_disconnect.Find(i)<0)
        {
            vert->Shares[k].C = old_shares[i].C;
            vert->Shares[k].N = old_shares[i].N;
            k++;
        }
    }
}

inline void Generic::SetBryTag (int i, int iEdgeFace, int Tag)
{
    if (Tag<0)
    {
        if (Cells[i]==NULL) throw new Fatal("Generic::SetBryTag: This method must be called after SetCell which allocates Cells");
        Cells[i]->BryTags[iEdgeFace] = Tag;
        if (TgdCells.Find(Cells[i])<0) TgdCells.Push (Cells[i]);
    }
}

inline void Generic::FindNeigh ()
{
    // build Bry2Cells map
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        size_t nverts = Cells[i]->V.Size();
        size_t nbrys  = (NDim==2 ? NVertsToNEdges2D[nverts] : NVertsToNFaces3D[nverts]);
        for (size_t j=0; j<nbrys; ++j)
        {
            BRYKEY(nverts,i,j)
            NeighDat_t neigh_dat(j,Cells[i]);
            Bry2Cells[brykey].Push (neigh_dat);
        }
    }

    // set neighbours information in cells
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        size_t nverts = Cells[i]->V.Size();
        size_t nbrys  = (NDim==2 ? NVertsToNEdges2D[nverts] : NVertsToNFaces3D[nverts]);
        for (size_t j=0; j<nbrys; ++j)
        {
            BRYKEY(nverts,i,j)
            Array<NeighDat_t> const & neigh_dat = Bry2Cells[brykey];
            for (size_t k=0; k<neigh_dat.Size(); ++k)
            {
                if (neigh_dat[k].second!=Cells[i])
                {
                    NeighDat_t dat(j, neigh_dat[k].second);
                    Cells[i]->Neighs[brykey] = dat;
                    //std::cout << Cells[i]->ID << "  " << neigh_dat[k].second->ID << "  " << brykey << "\n";
                }
            }
        }
        //std::cout << "\n";
    }
}

inline void Generic::GenO2Verts ()
{
    // generate neighbours map
    if (Bry2Cells.empty()) FindNeigh ();

    //for (BryCell_t::const_iterator p=Bry2Cells.begin(); p!=Bry2Cells.end(); ++p) std::cout << p->first << std::endl;

    // create vertices
    if (NDim==2)
    {
        for (BryCell_t::const_iterator p=Bry2Cells.begin(); p!=Bry2Cells.end(); ++p) // loop over edges
        {
            // add vertex
            Vec3_t mid;
            mid = 0.5*(Verts[boost::get<0>(p->first)]->C + Verts[boost::get<1>(p->first)]->C);
            //std::cout << boost::get<0>(p->first) << ", " << boost::get<1>(p->first) << "   :   " << mid << "   ";
            Verts.Push (new Vertex);
            Verts[Verts.Size()-1]->ID  = Verts.Size()-1;
            Verts[Verts.Size()-1]->Tag = 0;
            Verts[Verts.Size()-1]->C   = mid;

            // set connectivity in cells
            for (size_t i=0; i<p->second.Size(); ++i) // loop over cells sharing this edge
            {
                int idx_edge = p->second[i].first;
                Cell *  cell = p->second[i].second;
                size_t    nv = cell->V.Size();
                int idx_vert = -1;
                if (NVertsToVTKCell2D[nv]==VTK_TRIANGLE || NVertsToVTKCell2D[nv]==VTK_QUAD)
                {
                    // allocate space for nv vertices
                    cell->V.PushN (NULL,nv);
                    idx_vert = NVertsToEdge2D[2*nv][idx_edge][2];
                }
                else if (NVertsToVTKCell2D[nv]==VTK_QUADRATIC_TRIANGLE || NVertsToVTKCell2D[nv]==VTK_QUADRATIC_QUAD)
                {
                    // OK. Space already allocated
                    idx_vert = NVertsToEdge2D[nv][idx_edge][2];
                }
                else throw new Fatal("GenO2Verts::GenO2Verts: NDim==2D. Cell must be of type VTK_TRIANGLE (tri3) or VTK_QUAD (quad4) in order to generate O2 vertices. Number of vertices = %d is invalid",nv);

                // set pointer to vertex
                cell->V[idx_vert] = Verts[Verts.Size()-1];
            }
        }
    }
    else throw new Fatal("Generic::GenO2Verts: Method not available for 3D yet");
}

inline void Generic::Erase ()
{
    for (size_t i=0; i<Verts.Size(); ++i) if (Verts[i]!=NULL) delete Verts[i]; // it is only necessary to delete nodes in Verts array
    for (size_t i=0; i<Cells.Size(); ++i) if (Cells[i]!=NULL) delete Cells[i]; // it is only necessary to delete elems in Cells array
    if (Verts   .Size()>0) Verts   .Resize(0);
    if (Cells   .Size()>0) Cells   .Resize(0);
    if (TgdVerts.Size()>0) TgdVerts.Resize(0);
    if (TgdCells.Size()>0) TgdCells.Resize(0);
}

inline void Generic::WriteVTU (char const * FileKey, int VolSurfOrBoth) const
{
    // Vol=0, Surf=1, Both=2

    // boundary cells (for plotting bry tags)
    Array<Cell*> bcells;
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<Cells.Size(); ++i)
    {
        int nverts = Cells[i]->V.Size();
        for (BryTag_t::const_iterator p=Cells[i]->BryTags.begin(); p!=Cells[i]->BryTags.end(); ++p)
        {
            int ibry = p->first;
            int btag = p->second;
            if (btag<0)
            {
                int ibcell  = bcells.Size();
                int nbverts = (NDim==3 ? NVertsToNVertsPerFace3D[nverts] : 2);
                bcells.Push (new Cell);
                bcells[ibcell]->ID  = Cells.Size()+ibcell;
                bcells[ibcell]->Tag = btag;
                for (int j=0; j<nbverts; ++j)
                {
                    bcells[ibcell]->V.Push (new Vertex);
                    bcells[ibcell]->V[j]->ID = (NDim==3 ? Cells[i]->V[NVertsToFace3D[nverts][ibry][j]]->ID : Cells[i]->V[NVertsToEdge2D[nverts][ibry][j]]->ID);
                }
            }
        }
    }

    // shares data
    size_t max_nshares = 0;
    for (size_t i=0; i<Verts.Size(); ++i) if (Verts[i]->Shares.Size()>max_nshares) max_nshares = Verts[i]->Shares.Size();
    if (max_nshares<1) throw new Fatal("Mesh::WriteVTU: Max number of shares is wrong == 0");

    // data
    String fn(FileKey); fn.append(".vtu");
    std::ostringstream oss;
    size_t nn = Verts.Size();                           // number of Nodes
    size_t nc = (VolSurfOrBoth!=1 ? Cells.Size() : 0);  // number of Cells
    size_t nb = bcells.Size();                          // number of boundaries (faces or edges)

    // constants
    size_t          nimax = 40;        // number of integers in a line
    size_t          nfmax =  6;        // number of floats in a line
    Util::NumStream nsflo = Util::_8s; // number format for floats

    // header
    oss << "<?xml version=\"1.0\"?>\n";
    oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    oss << "  <UnstructuredGrid>\n";
    oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << nc+nb << "\">\n";

    // nodes: coordinates
    oss << "      <Points>\n";
    oss << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    size_t k = 0; oss << "        ";
    for (size_t i=0; i<nn; ++i)
    {
        oss << "  " << nsflo <<          Verts[i]->C(0) << " ";
        oss <<         nsflo <<          Verts[i]->C(1) << " ";
        oss <<         nsflo << (NDim==3?Verts[i]->C(2):0.0);
        k++;
        VTU_NEWLINE (i,k,nn,nfmax/3-1,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Points>\n";

    // elements: connectivity, offsets, types
    oss << "      <Cells>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    if (VolSurfOrBoth!=1)
    for (size_t i=0; i<nc; ++i)
    {
        oss << "  ";
        for (size_t j=0; j<Cells[i]->V.Size(); ++j) oss << Cells[i]->V[j]->ID << " ";
        k++;
        VTU_NEWLINE (i,k,nc,nimax/Cells[i]->V.Size(),oss);
    }
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<nb; ++i)
    {
        oss << "  ";
        for (size_t j=0; j<bcells[i]->V.Size(); ++j) oss << bcells[i]->V[j]->ID << " ";
        k++;
        VTU_NEWLINE (i,k,nb,nimax/bcells[i]->V.Size(),oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    size_t offset = 0;
    if (VolSurfOrBoth!=1)
    for (size_t i=0; i<nc; ++i)
    {
        offset += Cells[i]->V.Size();
        oss << (k==0?"  ":" ") << offset;
        k++;
        VTU_NEWLINE (i,k,nc,nimax,oss);
    }
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<nb; ++i)
    {
        offset += bcells[i]->V.Size();
        oss << (k==0?"  ":" ") << offset;
        k++;
        VTU_NEWLINE (i,k,nb,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    if (VolSurfOrBoth!=1)
    for (size_t i=0; i<nc; ++i)
    {
        if (NDim==2) oss << (k==0?"  ":" ") << NVertsToVTKCell2D[Cells[i]->V.Size()];
        else         oss << (k==0?"  ":" ") << NVertsToVTKCell3D[Cells[i]->V.Size()];
        k++;
        VTU_NEWLINE (i,k,nc,nimax,oss);
    }
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<nb; ++i)
    {
        oss << (k==0?"  ":" ") << NVertsToVTKCell2D[bcells[i]->V.Size()];
        k++;
        VTU_NEWLINE (i,k,nb,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Cells>\n";

    // data -- nodes
    oss << "      <PointData Scalars=\"TheScalars\">\n";
    oss << "        <DataArray type=\"Int32\" Name=\"" << "tag" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<nn; ++i)
    {
        oss << (k==0?"  ":" ") << Verts[i]->Tag;
        k++;
        VTU_NEWLINE (i,k,nn,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"" << "shares" << "\" NumberOfComponents=\""<< max_nshares <<"\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<nn; ++i)
    {
        oss << "  ";
        for (size_t j=0; j<max_nshares; ++j)
        {
            if (j<Verts[i]->Shares.Size()) oss << Verts[i]->Shares[j].C->ID << " ";
            else                           oss << -1 << " ";
        }
        k++;
        VTU_NEWLINE (i,k,nn,nimax/max_nshares-1,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </PointData>\n";

    // data -- elements
    oss << "      <CellData Scalars=\"TheScalars\">\n";
    oss << "        <DataArray type=\"Int32\" Name=\"" << "tag" << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    if (VolSurfOrBoth!=1)
    for (size_t i=0; i<nc; ++i)
    {
        oss << (k==0?"  ":" ") << Cells[i]->Tag;
        k++;
        VTU_NEWLINE (i,k,nc,nimax,oss);
    }
    if (VolSurfOrBoth>0)
    for (size_t i=0; i<nb; ++i)
    {
        oss << (k==0?"  ":" ") << bcells[i]->Tag;
        k++;
        VTU_NEWLINE (i,k,nb,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </CellData>\n";

    // Bottom
    oss << "    </Piece>\n";
    oss << "  </UnstructuredGrid>\n";
    oss << "</VTKFile>" << std::endl;

    // Write to file
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

inline void Generic::WriteMPY (char const * FileKey, bool OnlyMesh, char const * Extra) const
{
    String fn(FileKey); fn.append(".mpy");
    std::ostringstream oss;
    oss << "from mesh_drawing import *\n\n";
    oss << (*this) << "\n";
    oss << "d = Drawing(V,C)\n";
    if (OnlyMesh) oss << "d.draw_mesh(with_tags=False,with_ids=False)\n";
    else          oss << "d.draw_mesh()\n";
    if (Extra!=NULL) oss << Extra;
    oss << "d.show()\n";
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

}; // namespace Mesh

#endif // MECHSYS_MESH_H
