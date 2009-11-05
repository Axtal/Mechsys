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
#include "geomtype.h"
#include "util/maps.h"
#include "util/fatal.h"
#include "util/numstreams.h"
#include "models/model.h"
#include "fem/element.h"
#include "mesh/mesh.h"
#include "draw.h"

namespace FEM
{

class Domain
{
public:
    typedef void (*pCalcF) (double t, double & fx, double & fy, double & fz); ///< Callback function
    typedef std::map<int,pCalcF>   FDatabase_t; ///< Map tag to F function pointer
    typedef std::map<int,Model*>   Models_t;    ///< Map tag to model pointer

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
    void SetUVals     (SDPair const & UVals);
    void SetOutNods   (char const * FileKey, Array<int> const & IDsOrTags);
    void SetOutEles   (char const * FileKey, Array<int> const & IDsOrTags);
    void OutResults   (double Time, Vec_t const & F_int) const;
    void PrintResults (char const * NF="%15.6g", int IdxIP=-1) const;       ///< IdxIP < 0 => Centroid
    bool CheckError   (Table const & NodSol, Table const & EleSol, 
                       SDPair const & NodTol, SDPair const & EleTol) const; ///< At nodes and centroid
    bool CheckErrorIP (Table const & EleSol, SDPair const & EleTol) const;  ///< At integration points
    void WriteMPY     (char const * FileKey, double SFCoef=1.0) const;      ///< SFCoef: Scale-factor coefficient
    void WriteVTU     (char const * FileKey) const;

    // Data
    Mesh::Generic const & Msh;     ///< The mesh
    Dict          const & Prps;    ///< Element properties
    Dict          const & Inis;    ///< Initial values
    int                   NDim;    ///< Space dimension
    Models_t              Mdls;    ///< Models
    Array<Node*>          Nods;    ///< Nodes
    Array<Element*>       Eles;    ///< Elements
    Array<size_t>         OutNods; ///< ID of nodes for which output (results) is generated
    Array<size_t>         OutEles; ///< ID of elements for which output (results) is generated
    Array<std::ofstream*> FilNods; ///< Files with results at selected nodes (OutNods)
    Array<std::ofstream*> FilEles; ///< Files with results at selected elements (OutEles)
    FDatabase_t           FFuncs;  ///< Database of pointers to F functions
    Array<Node*>          NodsF;   ///< Nodes with F specified through callback function
    Array<pCalcF>         CalcF;   ///< Array with the F callbacks corresponding to FNodes
    Array<size_t>         Beams;   ///< Subset of elements of type Beam

#ifdef USE_BOOST_PYTHON
    void PySetOutNods (BPy::str const & FileKey, BPy::list const & IDsOrTags) { SetOutNods (BPy::extract<char const *>(FileKey)(), Array<int>(IDsOrTags)); }
    void PySetOutEles (BPy::str const & FileKey, BPy::list const & IDsOrTags) { SetOutEles (BPy::extract<char const *>(FileKey)(), Array<int>(IDsOrTags)); }
#endif
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Domain::Domain (Mesh::Generic const & TheMesh, Dict const & ThePrps, Dict const & TheMdls, Dict const & TheInis)
    : Msh(TheMesh), Prps(ThePrps), Inis(TheInis), NDim(TheMesh.NDim)
{
    // check
    //if (Prps.Keys.Size()!=TheMdls.Keys.Size()) throw new Fatal("Domain::Domain: Prps and Mdls dictionaries must have the same number of tags");
    //if (Prps.Keys.Size()!=TheInis.Keys.Size()) throw new Fatal("Domain::Domain: Prps and Inis dictionaries must have the same number of tags");

    // allocate models
    for (size_t i=0; i<TheMdls.Keys.Size(); ++i)
    {
        int tag = TheMdls.Keys[i];
        if (!Prps.HasKey(tag)) throw new Fatal("Domain::Domain: Prps and Mdls dictionaries must have the same tags");
        if (!Inis.HasKey(tag)) throw new Fatal("Domain::Domain: Inis and Mdls dictionaries must have the same tags");
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

    // map key => keys set
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
                if (BCs.HasKey(tag)) Eles[eid]->SetBCs (idx_edge_or_face, BCs(tag));
                else
                {
                    //std::cout << "Domain::SetBCs: BCs dictionary does not have tag=" << tag << " for edge/face\n";
                    /*
                    std::ostringstream oss;
                    oss << (*Eles[i]);
                    throw new Fatal("Domain::SetBCs: BCs dictionary does not have tag=%d for edge/face. Elem: %s", tag, oss.str().c_str());
                    */
                }
                keys_set[tag] = true;
            }
        }
    }

    // nodes (must be after elements)
    for (size_t i=0; i<Msh.TgdVerts.Size(); ++i)
    {
        int tag = Msh.TgdVerts[i]->Tag;
        if (tag<0)
        {
            size_t nid = Msh.TgdVerts[i]->ID;
            if (BCs.HasKey(tag))
            {
                if (BCs(tag).HasKey("ffunc")) // callback specified
                {
                    NodsF.Push (Nods[nid]);
                    FDatabase_t::const_iterator p = FFuncs.find(tag);
                    if (p!=FFuncs.end()) CalcF.Push (p->second);
                    else throw new Fatal("Domain::SetBCs: Callback function with tag=%d was not found in FFuncs database",tag);
                }
                else Nods[nid]->SetBCs (BCs(tag));
            }
            else
            {
                //std::cout << "Domain::SetBCs: BCs dictionary does not have tag=" << tag << " for node\n";
                /*
                std::ostringstream oss;
                oss << (*Nods[i]);
                throw new Fatal("Domain::SetBCs: BCs dictionary does not have tag=%d for node. Node: %s", tag, oss.str().c_str());
                */
            }
            keys_set[tag] = true;
        }
    }

    // check if all keys were set
    for (std::map<int,bool>::const_iterator p=keys_set.begin(); p!=keys_set.end(); ++p)
    {
        if (p->second==false)
        {
            // try elements (property) tags
            bool found = false;
            for (size_t i=0; i<Msh.Cells.Size(); ++i)
            {
                size_t eid = Msh.Cells[i]->ID;
                if (Msh.Cells[i]->Tag==p->first)
                {
                    Eles[eid]->SetBCs (/*ignored*/0, BCs(p->first));
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
    for (size_t i=0; i<Nods.Size(); ++i) Nods[i]->ClrBCs ();
    for (size_t i=0; i<Eles.Size(); ++i) Eles[i]->ClrBCs ();
    NodsF.Resize (0);
    CalcF.Resize (0);
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
            for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*of) << Util::_8s << Nods[nod]->UMap.Keys[j];
            for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*of) << Util::_8s << Nods[nod]->FMap.Keys[j];
            for (size_t j=0; j<Nods[nod]->nDOF(); ++j)
            {
                buf.Printf ("R%s",Nods[nod]->UMap.Keys[j].CStr());
                (*of) << Util::_8s << buf;
            }
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
        for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << Nods[nod]->U[j];
        for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << Nods[nod]->F[j];
        for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << Nods[nod]->F[j] - Nods[nod]->Fa(j,Time);
        for (size_t j=0; j<Nods[nod]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << F_int(Nods[nod]->EQ[j]);
        (*FilNods[i]) << "\n";
    }

    // elements
    for (size_t i=0; i<OutEles.Size(); ++i)
    {
        size_t ele = OutEles[i];
        (*FilEles[i]) << Util::_8s << Time;
        SDPair dat;
        Eles[ele]->GetState (dat);
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
            if (Nods[i]->FMap.HasKey(fkeys[j]))
            {
                size_t idx = Nods[i]->FMap(fkeys[j]); // idx of DOF
                buf.Printf(NF, Nods[i]->F[idx] - Nods[i]->Fa(idx,1.0));
            }
            else buf.Printf(nf, "---");
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
            for (size_t j=0; j<Nods.Size(); ++j)
            {
                size_t idx = Nods[j]->UMap(ukey); // idx of DOF
                err[j] = fabs(Nods[j]->F[idx] - Nods[j]->Fa(idx,1.0) - NodSol(key,j));
            }
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
    if (Nods[0]->U.Size()==1)
    {
        String key = Nods[0]->UMap.Keys[0];
        oss << "      <PointData Scalars=\"TheScalars\">\n";
        oss << "        <DataArray type=\"Float32\" Name=\"" << key << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<nn; ++i)
        {
            oss << (k==0?"  ":" ") << Nods[i]->Vert.Tag;
            k++;
            VTU_NEWLINE (i,k,nn,nimax,oss);
        }
        oss << "        </DataArray>\n";
        oss << "      </PointData>\n";
    }
    else if (Nods[0]->U.Size()==(size_t)NDim)
    {
        oss << "      <PointData Vectors=\"TheVectors\">\n";
        oss << "        <DataArray type=\"Float32\" Name=\"" << "U" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
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
        oss << "      </PointData>\n";
    }

    // data -- elements
    oss << "      <CellData Scalars=\"TheScalars\">\n";
    SDPair dat;
    Eles[0]->GetState (dat);
    for (size_t i=0; i<dat.Keys.Size(); ++i)
    {
        oss << "        <DataArray type=\"Float32\" Name=\"" << dat.Keys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t j=0; j<ne; ++j)
        {
            Eles[j]->GetState (dat);
            oss << (k==0?"  ":" ") << dat(dat.Keys[i]);
            k++;
            VTU_NEWLINE (j,k,nn,nfmax/3-1,oss);
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

    return os;
}

}; // namespace FEM

#endif // MECHSYS_DOMAIN_H
