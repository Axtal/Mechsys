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

#ifndef MECHSYS_DOMAIN_H
#define MECHSYS_DOMAIN_H

// Std Lib
#include <iostream>
#include <sstream>  // for istringstream, ostringstream
#include <cstdarg>  // for va_list, va_start, va_end
#include <map>

// MechSys
#include <mechsys/geomtype.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/models/model.h>
#include <mechsys/fem/element.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/draw.h>

namespace FEM
{

inline double Multiplier (double t) { return 1.0; }

class Domain
{
public:
    // typedefs
    typedef std::map<int,pCalcM> MDatabase_t; ///< Map tag to M function pointer
    typedef std::map<int,Model*> Models_t;    ///< Map tag to model pointer

    // Constructor
    Domain (Mesh::Generic const & Mesh,  ///< The mesh
            Dict const          & Prps,  ///< Element properties
            Dict const          & Mdls,  ///< Model names and parameters
            Dict const          & Inis); ///< Initial values

    // Destructor
    ~Domain ();

    // Methods
    void SetBCs       (Dict const & BCs);
    void ClrBCs       ();
    void Gravity      ();                                                   ///< Apply gravity
    void SetUVals     (SDPair const & UVals);                               ///< Set U values
    void SetOutNods   (char const * FileKey, Array<int> const & IDsOrTags); ///< Set nodes for output
    void SetOutEles   (char const * FileKey, Array<int> const & IDsOrTags); ///< Set elements for output
    void OutResults   (double Time, Vec_t const & F_int) const;             ///< Do output results
    void PrintResults (char const * NF="%15.6g", int IdxIP=-1) const;       ///< IdxIP < 0 => Centroid
    bool CheckError   (Table const & NodSol, Table const & EleSol, 
                       SDPair const & NodTol, SDPair const & EleTol) const; ///< At nodes and centroid
    bool CheckErrorIP (Table const & EleSol, SDPair const & EleTol) const;  ///< At integration points
    void WriteMPY     (char const * FileKey, double SFCoef=1.0) const;      ///< SFCoef: Scale-factor coefficient
    void WriteVTU     (char const * FileKey) const;                         ///< Write file for ParaView

    // Data
    Mesh::Generic const & Msh;     ///< The mesh
    Dict          const & Prps;    ///< Element properties
    Dict          const & Inis;    ///< Initial values
    int                   NDim;    ///< Space dimension
    double                gAccel;  ///< Gravity acceleration
    Models_t              Mdls;    ///< Models
    Array<Node*>          Nods;    ///< Nodes
    Array<Element*>       Eles;    ///< Elements
    Array<size_t>         OutNods; ///< ID of nodes for which output (results) is generated
    Array<size_t>         OutEles; ///< ID of elements for which output (results) is generated
    Array<std::ofstream*> FilNods; ///< Files with results at selected nodes (OutNods)
    Array<std::ofstream*> FilEles; ///< Files with results at selected elements (OutEles)
    Array<size_t>         Beams;   ///< Subset of elements of type Beam
    NodBCs_t              pU;      ///< Nodes with prescribed U
    NodBCs_t              pF;      ///< Nodes with prescribed F
    MDatabase_t           MFuncs;  ///< Database of pointers to M functions

#ifdef USE_BOOST_PYTHON
    void PySetOutNods (BPy::str const & FileKey, BPy::list const & IDsOrTags) { SetOutNods (BPy::extract<char const *>(FileKey)(), Array<int>(IDsOrTags)); }
    void PySetOutEles (BPy::str const & FileKey, BPy::list const & IDsOrTags) { SetOutEles (BPy::extract<char const *>(FileKey)(), Array<int>(IDsOrTags)); }
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Domain::Domain (Mesh::Generic const & TheMesh, Dict const & ThePrps, Dict const & TheMdls, Dict const & TheInis)
    : Msh(TheMesh), Prps(ThePrps), Inis(TheInis), NDim(TheMesh.NDim), gAccel(9.81)
{
    // allocate models
    for (size_t i=0; i<TheMdls.Keys.Size(); ++i)
    {
        int tag = TheMdls.Keys[i];
        if (!Prps.HasKey(tag)) throw new Fatal("Domain::Domain: Prps and Mdls dictionaries must have the same tags");
        if (TheMdls(tag).HasKey("name"))
        {
            String model_name;
            MODEL.Val2Key (TheMdls(tag)("name"), model_name);
            Mdls[tag] = AllocModel (model_name, NDim, TheMdls(tag));
        }
        else throw new Fatal("Domain::Domain: Dictionary of models must have keyword 'name' defining the name of the model");
    }

    // set nodes from mesh
    for (size_t i=0; i<Msh.Verts.Size(); ++i) Nods.Push (new Node((*Msh.Verts[i])));

    // set elements from mesh
    for (size_t i=0; i<Msh.Cells.Size(); ++i)
    {
        int tag = Msh.Cells[i]->Tag;
        if (!Prps.HasKey(tag)) throw new Fatal("Domain::SetMesh: Dictionary of element properties does not have tag=%d",tag);
        if (Prps(tag).HasKey("prob"))
        {
            // problem name
            String prob_name;
            PROB.Val2Key (Prps(tag)("prob"), prob_name);

            // model
            Model const * mdl = NULL;
            Models_t::const_iterator m = Mdls.find(tag);
            if (m!=Mdls.end()) mdl = m->second;

            // inis
            if (Inis.HasKey(tag)) Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, Prps(tag), Inis(tag), Nods));
            else                  Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, Prps(tag), SDPair() , Nods));

            // set array of Beams
            if (prob_name=="Beam") Beams.Push (Eles.Size()-1);
        }
        else throw new Fatal("Domain::SetMesh: Dictionary of properties must have keyword 'prob' defining the type of element corresponding to a specific problem");
    }
}

inline Domain::~Domain()
{
    for (Domain::Models_t::iterator p=Mdls.begin(); p!=Mdls.end(); ++p) delete p->second;
    for (size_t i=0; i<Nods   .Size(); ++i) if (Nods   [i]!=NULL) delete Nods   [i];
    for (size_t i=0; i<Eles   .Size(); ++i) if (Eles   [i]!=NULL) delete Eles   [i];
    for (size_t i=0; i<FilNods.Size(); ++i) if (FilNods[i]!=NULL) { FilNods[i]->close(); delete FilNods[i]; }
    for (size_t i=0; i<FilEles.Size(); ++i) if (FilEles[i]!=NULL) { FilEles[i]->close(); delete FilEles[i]; }
}

inline void Domain::SetBCs (Dict const & BCs)
{
    // clear previous BCs
    ClrBCs ();

    // map tag/key => keys set
    std::map<int,bool> keys_set;
    for (size_t i=0; i<BCs.Keys.Size(); ++i) keys_set[BCs.Keys[i]] = false;

    // elements
    for (size_t i=0; i<Msh.TgdCells.Size(); ++i)
    {
        Mesh::BryTag_t const & eftags = Msh.TgdCells[i]->BryTags;
        for (Mesh::BryTag_t::const_iterator p=eftags.begin(); p!=eftags.end(); ++p)
        {
            int idx_edge_or_face = p->first;
            int tag              = p->second;
            if (tag<0)
            {
                size_t eid = Msh.TgdCells[i]->ID;
                if (BCs.HasKey(tag))
                {
                    pCalcM calcm = &Multiplier;
                    if (BCs(tag).HasKey("mfunc")) // callback specified
                    {
                        MDatabase_t::const_iterator p = MFuncs.find(tag);
                        if (p!=MFuncs.end()) calcm = p->second;
                        else throw new Fatal("Domain::SetBCs: Multiplier function with tag=%d was not found in MFuncs database",tag);
                    }
                    Eles[eid]->SetBCs (idx_edge_or_face, BCs(tag), pF, pU, calcm);
                    keys_set[tag] = true;
                }
            }
        }
    }

    // nodes
    for (size_t i=0; i<Msh.TgdVerts.Size(); ++i)
    {
        int tag = Msh.TgdVerts[i]->Tag;
        if (tag<0)
        {
            size_t nid = Msh.TgdVerts[i]->ID;
            if (BCs.HasKey(tag))
            {
                pCalcM calcm = &Multiplier;
                if (BCs(tag).HasKey("mfunc")) // callback specified
                {
                    MDatabase_t::const_iterator p = MFuncs.find(tag);
                    if (p!=MFuncs.end()) calcm = p->second;
                    else throw new Fatal("Domain::SetBCs: Multiplier function with tag=%d was not found in MFuncs database",tag);
                }
                SDPair const & bcs = BCs(tag);
                for (StrDbl_t::const_iterator p=bcs.begin(); p!=bcs.end(); ++p)
                {
                    if      (Nods[nid]->UMap.HasKey(p->first)) pU[Nods[nid]].first[Nods[nid]->UMap(p->first)]  = p->second;
                    else if (Nods[nid]->FMap.HasKey(p->first))
                    {
                        pF[Nods[nid]].first[Nods[nid]->FMap(p->first)] += p->second;
                        pF[Nods[nid]].second = calcm;
                    }
                }
                keys_set[tag] = true;
            }
        }
    }

    // check if all keys were set
    for (std::map<int,bool>::const_iterator p=keys_set.begin(); p!=keys_set.end(); ++p)
    {
        if (p->second==false)
        {
            // try elements (property) tags => Beams for example
            bool found = false;
            for (size_t i=0; i<Msh.Cells.Size(); ++i)
            {
                size_t eid = Msh.Cells[i]->ID;
                int    tag = p->first;
                if (Msh.Cells[i]->Tag==tag)
                {
                    pCalcM calcm = &Multiplier;
                    if (BCs(tag).HasKey("mfunc")) // callback specified
                    {
                        MDatabase_t::const_iterator p = MFuncs.find(tag);
                        if (p!=MFuncs.end()) calcm = p->second;
                        else throw new Fatal("Domain::SetBCs: Multiplier function with tag=%d was not found in MFuncs database",tag);
                    }
                    Eles[eid]->SetBCs (/*ignored*/0, BCs(tag), pF, pU, calcm);
                    found = true;
                }
            }
            if (!found) 
            {
                std::ostringstream oss;
                oss << BCs(p->first);
                throw new Fatal("Domain::SetBCs: Keys=%s of BCs dictionary were not set. Probably neither Verts or Cells have tag=%d",oss.str().c_str(),p->first);
            }
        }
    }
}

inline void Domain::ClrBCs ()
{
    for (size_t i=0; i<Eles.Size(); ++i) Eles[i]->ClrBCs ();
    pF.clear();
    pU.clear();
}

inline void Domain::Gravity ()
{
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        int tag = Eles[i]->Cell.Tag;
        pCalcM calcm = &Multiplier;
        if (Prps(tag).HasKey("mfunc"))
        {
            MDatabase_t::const_iterator p = MFuncs.find(tag);
            if (p!=MFuncs.end()) calcm = p->second;
            else throw new Fatal("Domain::Gravity: Multiplier function with tag=%d was not found in MFuncs database",tag);
        }
        Eles[i]->Gravity (pF, calcm, gAccel);
    }
}

inline void Domain::SetUVals (SDPair const & UVals)
{
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        for (StrDbl_t::const_iterator p=UVals.begin(); p!=UVals.end(); ++p)
        {
            Nods[i]->U[Nods[i]->UMap(p->first)] = p->second;
        }
    }
}

inline void Domain::SetOutNods (char const * FNKey, Array<int> const & IDsOrTags)
{
    // find ids
    Array<int> ids;
    for (size_t i=0; i<IDsOrTags.Size(); ++i)
    {
        int nod_or_tag = IDsOrTags[i];
        int nod = 0;
        if (nod_or_tag<0) // tag
        {
            bool found = false;
            for (size_t j=0; j<Msh.TgdVerts.Size(); ++j)
            {
                if (nod_or_tag==Msh.TgdVerts[j]->Tag)
                {
                    nod   = Msh.TgdVerts[j]->ID;
                    found = true;
                    ids.Push (nod);
                }
            }
            if (!found) throw new Fatal("Domain::SetOutNods: Could not find any node with tag = %d", nod_or_tag);
        }
        else ids.Push (nod_or_tag);
    }

    // set files
    for (size_t i=0; i<ids.Size(); ++i)
    {
        int nod = ids[i];
        if (OutNods.Find(nod)<0)
        {
            String buf; buf.Printf("%s_nod_%d_%d.res",FNKey,nod,Nods[nod]->Vert.Tag);
            std::ofstream * of = new std::ofstream (buf.CStr(),std::ios::out);
            OutNods.Push (nod);
            FilNods.Push (of);
            (*of) << Util::_8s << "Time";
            (*of) << Util::_8s << "x";
            (*of) << Util::_8s << "y";  if (NDim==3)
            (*of) << Util::_8s << "z";
            for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*of) << Util::_8s << Nods[nod]->UMap.Keys[j];
            for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*of) << Util::_8s << Nods[nod]->FMap.Keys[j];
            for (size_t j=0; j<Nods[nod]->nDOF(); ++j)
            {
                buf.Printf ("%s_int",Nods[nod]->FMap.Keys[j].CStr());
                (*of) << Util::_8s << buf;
            }
            (*of) << "\n";
        }
    }
}

inline void Domain::SetOutEles (char const * FNKey, Array<int> const & IDsOrTags)
{
    // find ids
    Array<int> ids;
    for (size_t i=0; i<IDsOrTags.Size(); ++i)
    {
        int ele_or_tag = IDsOrTags[i];
        int ele = 0;
        if (ele_or_tag<0) // tag
        {
            bool found = false;
            for (size_t j=0; j<Msh.Cells.Size(); ++j)
            {
                if (ele_or_tag==Msh.Cells[j]->Tag)
                {
                    ele   = Msh.Cells[j]->ID;
                    found = true;
                    ids.Push (ele);
                }
            }
            if (!found) throw new Fatal("Domain::SetOutEles: Could not find any element with tag = %d", ele_or_tag);
        }
        else ids.Push (ele_or_tag);
    }

    // set files
    for (size_t i=0; i<ids.Size(); ++i)
    {
        int ele = ids[i];
        if (OutEles.Find(ele)<0)
        {
            String buf; buf.Printf("%s_ele_%d_%d.res",FNKey,ele,Eles[ele]->Cell.Tag);
            std::ofstream * of = new std::ofstream (buf.CStr(),std::ios::out);
            OutEles.Push (ele);
            FilEles.Push (of);
            SDPair dat;
            Eles[ele]->GetState (dat);
            (*of) << Util::_8s << "Time";
            (*of) << Util::_8s << "x";
            (*of) << Util::_8s << "y";  if (NDim==3)
            (*of) << Util::_8s << "z";
            for (size_t j=0; j<dat.Keys.Size(); ++j) (*of) << Util::_8s << dat.Keys[j];
            (*of) << "\n";
        }
    }
}

inline void Domain::OutResults (double Time, Vec_t const & F_int) const
{
    // nodes
    for (size_t i=0; i<OutNods.Size(); ++i)
    {
        size_t nod = OutNods[i];
        (*FilNods[i]) << Util::_8s << Time;
        (*FilNods[i]) << Util::_8s << Nods[nod]->Vert.C(0);
        (*FilNods[i]) << Util::_8s << Nods[nod]->Vert.C(1);  if (NDim==3)
        (*FilNods[i]) << Util::_8s << Nods[nod]->Vert.C(2);
        for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << Nods[nod]->U[j];
        for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << Nods[nod]->F[j];
        for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << F_int(Nods[nod]->EQ[j]);
        (*FilNods[i]) << "\n";
    }

    // elements
    for (size_t i=0; i<OutEles.Size(); ++i)
    {
        size_t ele = OutEles[i];
        (*FilEles[i]) << Util::_8s << Time;
        SDPair dat;
        Vec_t  Xct;
        Eles[ele]->GetState (dat);
        Eles[ele]->Centroid (Xct);
        (*FilEles[i]) << Util::_8s << Xct(0);
        (*FilEles[i]) << Util::_8s << Xct(1);  if (NDim==3)
        (*FilEles[i]) << Util::_8s << Xct(2);
        for (size_t j=0; j<dat.Keys.Size(); ++j) (*FilEles[i]) << Util::_8s << dat(dat.Keys[j]);
        (*FilEles[i]) << "\n";
    }
}

inline void Domain::PrintResults (char const * NF, int IdxIP) const
{
    std::cout << "\n[1;37m--- Results ------------------------------------------------------------------\n";

    // nodes: ukeys, fkeys
    size_t neq = 0;
    Array<String> ukeys, fkeys;
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        for (size_t j=0; j<Nods[i]->nDOF(); ++j)
        {
            String ukey = Nods[i]->UMap.Keys[j];
            String fkey = Nods[i]->FMap.Keys[j];
            if (ukeys.Find(ukey)<0) ukeys.Push(ukey);
            if (fkeys.Find(fkey)<0) fkeys.Push(fkey);
            neq++;
        }
    }

    // number format for text
    String nf(NF);
    size_t pos;
    pos=nf.find("g"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("g",pos+1); }
    pos=nf.find("f"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("f",pos+1); }
    pos=nf.find("e"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("e",pos+1); }

    // nodes: header
    String buf;
    std::cout << Util::_6 << "Node";
    for (size_t i=0; i<ukeys.Size(); ++i) { buf.Printf(nf, ukeys[i].CStr());  std::cout<<buf; }
    for (size_t i=0; i<fkeys.Size(); ++i) { buf.Printf(nf, fkeys[i].CStr());  std::cout<<buf; }
    for (size_t i=0; i<ukeys.Size(); ++i)
    {
        buf.Printf ("R%s",ukeys[i].CStr());
        String tmp;
        tmp.Printf (nf,buf.CStr());
        std::cout << tmp;
    }
    std::cout << "[0m\n";

    // nodes: data
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        std::cout << Util::_6 << Nods[i]->Vert.ID;
        for (size_t j=0; j<ukeys.Size(); ++j)
        {
            if (Nods[i]->UMap.HasKey(ukeys[j]))
            {
                size_t idx = Nods[i]->UMap(ukeys[j]); // idx of DOF
                buf.Printf(NF, Nods[i]->U[idx]);
            }
            else buf.Printf(nf, "---");
            std::cout << buf;
        }
        for (size_t j=0; j<fkeys.Size(); ++j)
        {
            if (Nods[i]->FMap.HasKey(fkeys[j]))
            {
                size_t idx = Nods[i]->FMap(fkeys[j]); // idx of DOF
                buf.Printf(NF, Nods[i]->F[idx]);
            }
            else buf.Printf(nf, "---");
            std::cout << buf;
        }
        for (size_t j=0; j<fkeys.Size(); ++j)
        {
            //if (Nods[i]->FMap.HasKey(fkeys[j]))
            //{
                //size_t idx = Nods[i]->FMap(fkeys[j]); // idx of DOF
                //buf.Printf(NF, Nods[i]->F[idx] - Nods[i]->Fa(idx,1.0));
            //}
            //else
                buf.Printf(nf, "---");
            std::cout << buf;
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    // elems: keys
    Array<String> keys;
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        SDPair dat;
        Eles[i]->GetState (dat);
        for (size_t j=0; j<dat.Keys.Size(); ++j)
        {
            if (keys.Find(dat.Keys[j])<0) keys.Push (dat.Keys[j]);
        }
    }

    // elems: header
    std::cout << "[1;37m" << Util::_6 << "Elem";
    buf.Printf(nf, "x");  std::cout << buf;
    buf.Printf(nf, "y");  std::cout << buf;  if (NDim==3) {
    buf.Printf(nf, "z");  std::cout << buf; }
    for (size_t i=0; i<keys.Size(); ++i) { buf.Printf(nf, keys[i].CStr());  std::cout<<buf; }
    std::cout << "[0m\n";

    // elems: data
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        std::cout << Util::_6 << Eles[i]->Cell.ID;
        Vec_t  X;
        SDPair dat;
        Eles[i]->GetState (dat, IdxIP);
        if (IdxIP<0) Eles[i]->Centroid   (X);
        else         Eles[i]->CoordsOfIP (IdxIP, X);
        buf.Printf(NF, X(0)); std::cout << buf;
        buf.Printf(NF, X(1)); std::cout << buf;  if (NDim==3) {
        buf.Printf(NF, X(2)); std::cout << buf; }
        for (size_t j=0; j<keys.Size(); ++j)
        {
            if (dat.HasKey(keys[j])) buf.Printf(NF, dat(keys[j]));
            else                     buf.Printf(nf, "---");
            std::cout << buf;
        }
        std::cout << "\n";
    }
}

inline bool Domain::CheckError (Table const & NodSol, Table const & EleSol, SDPair const & NodTol, SDPair const & EleTol) const
{
    /*
        NodBCs_t::const_iterator p = pF.find(Nods[nod]);
        if (p==pF.end()) { for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << 0.0; }
        else
        {
            for (size_t j=0; j<Nods[nod]->nDOF(); ++j)
            {
                IntDbl_t::const_iterator q = p->second.find(j);
                if (q==p->second.end()) (*FilNods[i]) << Util::_8s << 0.0;
                else                    (*FilNods[i]) << Util::_8s << Nods[nod]->F[j] - Nods[nod]->Fa(j,Time);
            }
        }
    */

    // header
    std::cout << "\n[1;37m--- Error Summary --- nodes and centroid -------------------------------------\n";
    std::cout << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm" << "[0m\n";

    // results
    bool error = false;

    // nodes
    for (size_t i=0; i<NodSol.Keys.Size(); ++i)
    {
        // calc error
        String key = NodSol.Keys[i];
        Array<double> err(NodSol.NRows);
        if (key[0]=='R') // reaction
        {
            String ukey;
            for (size_t j=1; j<key.size(); ++j) ukey.Printf("%s%c",ukey.CStr(),key[j]);
            //for (size_t j=0; j<Nods.Size(); ++j)
            //{
                //size_t idx = Nods[j]->UMap(ukey); // idx of DOF
                //err[j] = fabs(Nods[j]->F[idx] - Nods[j]->Fa(idx,1.0) - NodSol(key,j));
            //}
        }
        else
        {
            for (size_t j=0; j<Nods.Size(); ++j)
            {
                if (Nods[j]->UMap.HasKey(key)) err[j] = fabs(Nods[j]->U[Nods[j]->UMap(key)] - NodSol(key,j));
                else                           err[j] = fabs(Nods[j]->F[Nods[j]->FMap(key)] - NodSol(key,j));
            }
        }

        // summary
        double max_err = err[err.Max()];
        double tol     = NodTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err[err.Min()] << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? "[1;31m" : "[1;32m") << Util::_8s<<max_err << "[0m" << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    // elements
    for (size_t i=0; i<EleSol.Keys.Size(); ++i)
    {
        // calc error
        String key = EleSol.Keys[i];
        Array<double> err(EleSol.NRows);
        for (size_t j=0; j<Eles.Size(); ++j)
        {
            SDPair dat;
            Eles[j]->GetState (dat);
            err [j] = fabs(dat(key) - EleSol(key,j));
        }

        // summary
        double max_err = err[err.Max()];
        double tol     = EleTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err[err.Min()] << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? "[1;31m" : "[1;32m") << Util::_8s<<max_err << "[0m" << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    return error;
}

inline bool Domain::CheckErrorIP (Table const & EleSol, SDPair const & EleTol) const
{
    // header
    std::cout << "\n[1;37m--- Error Summary --- integration points -------------------------------------\n";
    std::cout << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm" << "[0m\n";

    // results
    bool error = false;

    // elements
    for (size_t i=0; i<EleSol.Keys.Size(); ++i)
    {
        // calc error
        String key = EleSol.Keys[i];
        Array<double> err(EleSol.NRows);
        for (size_t j=0; j<Eles.Size(); ++j)
        {
            if (Eles[j]->GE==NULL) throw new Fatal("Domain::CheckError: This method works only when GE (geometry element) is not NULL");
            Array<SDPair> res;
            Eles[j]->GetState (res);
            for (size_t k=0; k<Eles[j]->GE->NIP; ++k)
            {
                size_t row = k + j*Eles[j]->GE->NIP;
                err [row] = fabs(res[k](key) - EleSol(key,row));
            }
        }

        // summary
        double max_err = err[err.Max()];
        double tol     = EleTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err[err.Min()] << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? "[1;31m" : "[1;32m") << Util::_8s<<max_err << "[0m" << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    std::cout << "\n";

    return error;
}

inline void Domain::WriteMPY (char const * FNKey, double SFCoef) const
{
    // bounding box
    Vec3_t min, max, del;
    min = Nods[0]->Vert.C(0), Nods[0]->Vert.C(1), Nods[0]->Vert.C(2);
    max = min;
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        if (Nods[i]->Vert.C(0)<min(0)) min(0) = Nods[i]->Vert.C(0);
        if (Nods[i]->Vert.C(1)<min(1)) min(1) = Nods[i]->Vert.C(1);
        if (Nods[i]->Vert.C(2)<min(2)) min(2) = Nods[i]->Vert.C(2);
        if (Nods[i]->Vert.C(0)>max(0)) max(0) = Nods[i]->Vert.C(0);
        if (Nods[i]->Vert.C(1)>max(1)) max(1) = Nods[i]->Vert.C(1);
        if (Nods[i]->Vert.C(2)>max(2)) max(2) = Nods[i]->Vert.C(2);
    }
    del  = max - min;
    double max_dist = Norm(del);

    /*
    // max element diagonal
    Vec3_t min, max, del;
    double diag, max_diag = 0.0;
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        min = Eles[i]->Con[0]->Vert.C(0), Eles[i]->Con[0]->Vert.C(1), Eles[i]->Con[0]->Vert.C(2);
        max = min;
        for (size_t j=0; j<Eles[i]->Con.Size(); ++j)
        {
            if (Eles[i]->Con[j]->Vert.C(0)<min(0)) min(0) = Eles[i]->Con[j]->Vert.C(0);
            if (Eles[i]->Con[j]->Vert.C(1)<min(1)) min(1) = Eles[i]->Con[j]->Vert.C(1);
            if (Eles[i]->Con[j]->Vert.C(2)<min(2)) min(2) = Eles[i]->Con[j]->Vert.C(2);
            if (Eles[i]->Con[j]->Vert.C(0)>max(0)) max(0) = Eles[i]->Con[j]->Vert.C(0);
            if (Eles[i]->Con[j]->Vert.C(1)>max(1)) max(1) = Eles[i]->Con[j]->Vert.C(1);
            if (Eles[i]->Con[j]->Vert.C(2)>max(2)) max(2) = Eles[i]->Con[j]->Vert.C(2);
        }
        del  = max - min;
        diag = Norm(del);
        if (diag>max_diag) max_diag = diag;
    }
    */

    String fn(FNKey);  fn.append(".mpy");
    std::ofstream of(fn.CStr(), std::ios::out);
    MPL::Header   (of);
    for (size_t i=0; i<Eles.Size(); ++i) Eles[i]->Draw (of, SFCoef*max_dist);
    MPL::AddPatch (of);
    MPL::Show     (of);
    of.close      ();
}

inline void Domain::WriteVTU (char const * FNKey) const
{
    // data
    String fn(FNKey); fn.append(".vtu");
    std::ostringstream oss;
    size_t nn = Nods.Size(); // number of nodes
    size_t ne = Eles.Size(); // number of elements

    // constants
    size_t          nimax = 40;        // number of integers in a line
    size_t          nfmax =  6;        // number of floats in a line
    Util::NumStream nsflo = Util::_8s; // number format for floats

    // check available data
    bool nod_displacements = false;
    bool nod_velocities    = false;
    bool ele_velocities    = false;
    bool ele_gradients     = false;
    Array<String> nod_keys;
    Array<String> ele_keys;
    for (size_t i=0; i<nn; ++i)
    {
        for (size_t j=0; j<Nods[i]->UMap.Keys.Size(); ++j)
        {
            if (nod_keys.Find(Nods[i]->UMap.Keys[j])<0) nod_keys.Push (Nods[i]->UMap.Keys[j]);
        }
        if (Nods[i]->UMap.HasKey("ux") || Nods[i]->UMap.HasKey("uy") || Nods[i]->UMap.HasKey("uz")) nod_displacements = true;
        if (Nods[i]->UMap.HasKey("vx") || Nods[i]->UMap.HasKey("vy") || Nods[i]->UMap.HasKey("vz")) nod_velocities    = true;
    }
    for (size_t i=0; i<ne; ++i)
    {
        SDPair dat;
        Eles[i]->GetState (dat);
        for (size_t j=0; j<dat.Keys.Size(); ++j)
        {
            if (ele_keys.Find(dat.Keys[j])<0) ele_keys.Push (dat.Keys[j]);
        }
        if (dat.HasKey("vx") || dat.HasKey("vy") || dat.HasKey("vz")) ele_velocities = true;
        if (dat.HasKey("gx") || dat.HasKey("gy") || dat.HasKey("gz")) ele_gradients  = true;
        if (dat.HasKey("sx") || dat.HasKey("sy") || dat.HasKey("sz"))
        {
            if (ele_keys.Find("pcam")<0) ele_keys.Push ("pcam");
            if (ele_keys.Find("qcam")<0) ele_keys.Push ("qcam");
        }
    }
    //std::cout << "nod_keys          = " << nod_keys          << std::endl;
    //std::cout << "ele_keys          = " << ele_keys          << std::endl;
    //std::cout << "nod_displacements = " << nod_displacements << std::endl;
    //std::cout << "nod_velocities    = " << nod_velocities    << std::endl;
    //std::cout << "ele_velocities    = " << ele_velocities    << std::endl;

    // extrapolate data
    Array<String> nod_keys_extrap; // keys of extrapolated results
    Array<SDPair> nod_results(nn); // results at each node
    for (size_t i=0; i<ne; ++i)
    {
        Array<SDPair> loc_res; // local results: size==number of nodes in element
        Eles[i]->StateAtNodes (loc_res);
        for (size_t j=0; j<Eles[i]->Con.Size(); ++j) 
        {
            size_t vid = Eles[i]->Con[j]->Vert.ID;
            for (StrDbl_t::const_iterator it=loc_res[j].begin(); it!=loc_res[j].end(); ++it)
            {
                String key("count_"+it->first);
                if (nod_results[vid].HasKey(it->first))
                {
                    nod_results[vid][it->first] += it->second;
                    nod_results[vid][key]       += 1.0;
                }
                else
                {
                    nod_results[vid].Set(it->first.CStr(), it->second);
                    nod_results[vid].Set(key.CStr(),       1.0);
                }
                if (nod_keys_extrap.Find(it->first)<0) nod_keys_extrap.Push (it->first);
            }
        }
    }

    // header
    oss << "<?xml version=\"1.0\"?>\n";
    oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    oss << "  <UnstructuredGrid>\n";
    oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << ne << "\">\n";

    // nodes: coordinates
    oss << "      <Points>\n";
    oss << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    size_t k = 0; oss << "        ";
    for (size_t i=0; i<nn; ++i)
    {
        oss << "  " << nsflo <<          Nods[i]->Vert.C(0) << " ";
        oss <<         nsflo <<          Nods[i]->Vert.C(1) << " ";
        oss <<         nsflo << (NDim==3?Nods[i]->Vert.C(2):0.0);
        k++;
        VTU_NEWLINE (i,k,nn,nfmax/3-1,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Points>\n";

    // elements: connectivity, offsets, types
    oss << "      <Cells>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<ne; ++i)
    {
        oss << "  ";
        for (size_t j=0; j<Eles[i]->Con.Size(); ++j) oss << Eles[i]->Con[j]->Vert.ID << " ";
        k++;
        VTU_NEWLINE (i,k,ne,nimax/Eles[i]->Con.Size(),oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    size_t offset = 0;
    for (size_t i=0; i<ne; ++i)
    {
        offset += Eles[i]->Con.Size();
        oss << (k==0?"  ":" ") << offset;
        k++;
        VTU_NEWLINE (i,k,ne,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<ne; ++i)
    {
        if (NDim==2) oss << (k==0?"  ":" ") << NVertsToVTKCell2D[Eles[i]->Con.Size()];
        else         oss << (k==0?"  ":" ") << NVertsToVTKCell3D[Eles[i]->Con.Size()];
        k++;
        VTU_NEWLINE (i,k,ne,nimax,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Cells>\n";

    // data -- nodes
    oss << "      <PointData Scalars=\"point_scalars\">\n";
    for (size_t i=0; i<nod_keys.Size(); ++i)
    {
        if (!(nod_keys[i]=="ux" || nod_keys[i]=="uy" || nod_keys[i]=="uz"))
        {
            oss << "        <DataArray type=\"Float32\" Name=\"" << nod_keys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                oss << (k==0?"  ":" ") << Nods[j]->U[Nods[j]->UMap(nod_keys[i])];
                k++;
                VTU_NEWLINE (j,k,nn,nfmax-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }
    if (nod_displacements)
    {
        oss << "        <DataArray type=\"Float32\" Name=\"" << "u" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<nn; ++i)
        {
            oss << "  " << nsflo <<          Nods[i]->U[Nods[i]->UMap("ux")] << " ";
            oss <<         nsflo <<          Nods[i]->U[Nods[i]->UMap("uy")] << " ";
            oss <<         nsflo << (NDim==3?Nods[i]->U[Nods[i]->UMap("uz")]:0.0);
            k++;
            VTU_NEWLINE (i,k,nn,nfmax/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }
    for (size_t i=0; i<nod_keys_extrap.Size(); ++i)
    {
        oss << "        <DataArray type=\"Float32\" Name=\"" << nod_keys_extrap[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t j=0; j<nn; ++j)
        {
            String ckey("count_"+nod_keys_extrap[i]);
            oss << (k==0?"  ":" ") << nod_results[j](nod_keys_extrap[i])/nod_results[j](ckey);
            k++;
            VTU_NEWLINE (j,k,nn,nfmax-1,oss);
        }
        oss << "        </DataArray>\n";
    }
    oss << "      </PointData>\n";

    // data -- elements
    oss << "      <CellData Scalars=\"cell_scalars\">\n";
    for (size_t i=0; i<ele_keys.Size(); ++i)
    {
        if (!(ele_keys[i]=="vx" || ele_keys[i]=="vy" || ele_keys[i]=="vz") &&
            !(ele_keys[i]=="gx" || ele_keys[i]=="gy" || ele_keys[i]=="gz"))
        {
            oss << "        <DataArray type=\"Float32\" Name=\"" << ele_keys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<ne; ++j)
            {
                SDPair dat;
                Eles[j]->GetState (dat);
                if (ele_keys[i]=="pcam" && dat.HasKey("sx"))
                {
                    double pcam = -(dat("sx")+dat("sy")+dat("sz"))/3.0;
                    oss << (k==0?"  ":" ") << pcam;
                }
                else if (ele_keys[i]=="qcam" && dat.HasKey("sx"))
                {
                    double m    = (NDim==3 ? pow(dat("syz"),2.0)+pow(dat("szx"),2.0) : 0.0);
                    double qcam = sqrt(pow(dat("sx")-dat("sy"),2.0) + pow(dat("sy")-dat("sz"),2.0) + pow(dat("sz")-dat("sx"),2.0) + 3.0*(pow(dat("sxy"),2.0)+m))/sqrt(2.0);
                    oss << (k==0?"  ":" ") << qcam;
                }
                else
                {
                    oss << (k==0?"  ":" ") << (dat.HasKey(ele_keys[i]) ? dat(ele_keys[i]) : 0.0);
                }
                k++;
                VTU_NEWLINE (j,k,ne,nfmax-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }
    if (ele_velocities)
    {
        oss << "        <DataArray type=\"Float32\" Name=\"" << "v" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<ne; ++i)
        {
            SDPair dat;
            Eles[i]->GetState (dat);
            oss << "  " << nsflo <<          dat("vx") << " ";
            oss <<         nsflo <<          dat("vy") << " ";
            oss <<         nsflo << (NDim==3?dat("vz"):0.0);
            k++;
            VTU_NEWLINE (i,k,ne,nfmax/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }
    if (ele_gradients)
    {
        oss << "        <DataArray type=\"Float32\" Name=\"" << "g" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<ne; ++i)
        {
            SDPair dat;
            Eles[i]->GetState (dat);
            oss << "  " << nsflo <<          dat("gx") << " ";
            oss <<         nsflo <<          dat("gy") << " ";
            oss <<         nsflo << (NDim==3?dat("gz"):0.0);
            k++;
            VTU_NEWLINE (i,k,ne,nfmax/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }
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

std::ostream & operator<< (std::ostream & os, Domain const & D)
{
    os << "\n[1;37m--- Models -------------------------------------------------------------------[0m\n";
    for (Domain::Models_t::const_iterator p=D.Mdls.begin(); p!=D.Mdls.end(); ++p)
        os << p->first << " " << (*p->second) << std::endl;

    os << "\n[1;37m--- Elements properties ------------------------------------------------------[0m\n";
    os << D.Prps << std::endl;

    os << "\n[1;37m--- Initial values -----------------------------------------------------------[0m\n";
    os << D.Inis << std::endl;

    os << "\n[1;37m--- Nodes --------------------------------------------------------------------[0m\n";
    for (size_t i=0; i<D.Nods.Size(); ++i) os << (*D.Nods[i]) << "\n";

    os << "\n[1;37m--- Elements -----------------------------------------------------------------[0m\n";
    for (size_t i=0; i<D.Eles.Size(); ++i) os << (*D.Eles[i]) << "\n";

    os << "\n[1;37m--- Boundary Conditions ------------------------------------------------------[0m\n";
    os << "  Nodes with prescribed F:\n";
    for (NodBCs_t::const_iterator p=D.pF.begin(); p!=D.pF.end(); ++p)
    {
        os << Util::_6 << p->first->Vert.ID << "  ";
        for (IntDbl_t::const_iterator q=p->second.first.begin(); q!=p->second.first.end(); ++q)
        {
            size_t idof = q->first;
            if (q!=p->second.first.begin()) os << "   ";
            os << "D" << p->first->FMap.Keys[idof] << "=" << Util::_10_4 << q->second;
        }
        os << "\n";
    }
    os << "\n  Nodes with prescribed U:\n";
    for (NodBCs_t::const_iterator p=D.pU.begin(); p!=D.pU.end(); ++p)
    {
        os << Util::_6 << p->first->Vert.ID << "  ";
        for (IntDbl_t::const_iterator q=p->second.first.begin(); q!=p->second.first.end(); ++q)
        {
            size_t idof = q->first;
            if (q!=p->second.first.begin()) os << ", ";
            os << "D" << p->first->UMap.Keys[idof] << "=" << Util::_10_6 << q->second;
        }
        os << "\n";
    }

    return os;
}

}; // namespace FEM

#endif // MECHSYS_DOMAIN_H
