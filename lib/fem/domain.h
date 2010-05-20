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

typedef std::map<Node*,double> InclSupport_t; ///< Inclined support type. Maps Node ==> alpha

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
    void Deactivate   (int EleTag);                                         ///< Deactivate all elements with EleTag
    void SetUVals     (SDPair const & UVals);                               ///< Set U values
    void SetOutNods   (char const * FileKey, Array<int> const & IDsOrTags); ///< Set nodes for output
    void SetOutEles   (char const * FileKey, Array<int> const & IDsOrTags); ///< Set elements for output
    void OutResults   (double Time, Vec_t const & F_int) const;             ///< Do output results
    void PrintResults (char const * NF="%15.6g") const;                     ///< Print results
    bool CheckError   (Table const & NodSol, SDPair const & NodTol) const;  ///< At nodes
    bool CheckError   (Table const & NodSol, Table const & EleSol, 
                       SDPair const & NodTol, SDPair const & EleTol) const; ///< At nodes and centroid
    bool CheckErrorIP (Table const & EleSol, SDPair const & EleTol) const;  ///< At integration points
    void WriteMPY     (char const * FileKey, double SFCoef=1.0,
                       char const * Extra=NULL) const;                      ///< SFCoef: Scale-factor coefficient
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
    Array<String>         DisplKeys;   ///< Displacement keys
    InclSupport_t         InclSupport; ///< Inclined support

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

    // displacement keys
    DisplKeys.Resize (NDim);
    DisplKeys[0] = "ux";
    DisplKeys[1] = "uy";  if (NDim==3)
    DisplKeys[2] = "uz";
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

    // set maps with bry info
    typedef std::pair<Element*,int> eleside_t;       // (element,side) pair
    typedef std::pair<int,SDPair>   tagbcs_t;        // (tag,bcs data) pair
    std::map<eleside_t,tagbcs_t> eleside2tag_to_ubc; // map: (ele,side) ==> U bry cond
    std::map<eleside_t,tagbcs_t> eleside2tag_to_fbc; // map: (ele,side) ==> F bry cond
    std::map<Node*,tagbcs_t>     nod2tag_to_ubc;     // map: node ==> U bry cond
    std::map<Node*,tagbcs_t>     nod2tag_to_fbc;     // map: node ==> F bry cond
    for (size_t i=0; i<BCs.Keys.Size(); ++i)
    {
        int bc_tag = BCs.Keys[i];
        SDPair const & bcs = BCs(bc_tag);
        bool found = false;

        // find elements with edges set with this bc_tag
        for (size_t j=0; j<Msh.TgdCells.Size(); ++j)
        {
            size_t eid = Msh.TgdCells[j]->ID;
            if (Eles[eid]->Active)
            {
                Mesh::BryTag_t const & eftags = Msh.TgdCells[j]->BryTags;
                for (Mesh::BryTag_t::const_iterator p=eftags.begin(); p!=eftags.end(); ++p)
                {
                    int idx_side = p->first;
                    int side_tag = p->second;
                    if (side_tag==bc_tag) // found
                    {
                        found = true;
                        eleside_t es(Eles[eid],idx_side);
                        for (size_t k=0; k<bcs.Keys.Size(); ++k) // we have to split Ubcs from Fbcs
                        {
                            if (Eles[eid]->UKeys.Find(bcs.Keys[k])>=0)
                            {
                                eleside2tag_to_ubc[es].first = bc_tag;
                                eleside2tag_to_ubc[es].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else
                            {
                                eleside2tag_to_fbc[es].first = bc_tag;
                                eleside2tag_to_fbc[es].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                        }
                    }
                }
            }
        }

        // find nodes with tags equal to this bc_tag
        for (size_t j=0; j<Msh.TgdVerts.Size(); ++j)
        {
            size_t nid = Msh.TgdVerts[j]->ID;
            if (Nods[nid]->NShares>0) // active
            {
                if (Msh.TgdVerts[j]->Tag==bc_tag) // found
                {
                    found = true;
                    for (size_t k=0; k<bcs.Keys.Size(); ++k) // we have to split Ubcs from Fbcs
                    {
                        if (bcs.Keys[k]=="inclsupport")
                        {
                            nod2tag_to_ubc[Nods[nid]].first  = bc_tag;
                            nod2tag_to_ubc[Nods[nid]].second = bcs;
                            break;
                        }
                        else
                        {
                            if (Nods[nid]->UMap.HasKey(bcs.Keys[k]))
                            {
                                nod2tag_to_ubc[Nods[nid]].first = bc_tag;
                                nod2tag_to_ubc[Nods[nid]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else
                            {
                                nod2tag_to_fbc[Nods[nid]].first = bc_tag;
                                nod2tag_to_fbc[Nods[nid]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                        }
                    }
                }
            }
        }

        if (!found) 
        {
            // find elements with tags equal to this bc_tag. Special cases: qn of beams, s (source term), cbx (cetrifugal body force), ...
            for (size_t j=0; j<Msh.Cells.Size(); ++j)
            {
                size_t eid = Msh.Cells[j]->ID;
                if (Eles[eid]->Active)
                {
                    if (Msh.Cells[j]->Tag==bc_tag) // found
                    {
                        found = true;
                        int idx_side = 0; // irrelevant
                        eleside_t es(Eles[eid],idx_side);
                        eleside2tag_to_fbc[es].first  = bc_tag;
                        eleside2tag_to_fbc[es].second = bcs;
                    }
                }
            }
            if (!found) throw new Fatal("FEM::Domain::SetBCs: Could not find any edge/face of element or node (or line element) with boundary Tag=%d",bc_tag);
        }
    }

    // set F bcs at sides (edges/faces) of elements (or at the element itself => Line elements/beams)
    for (std::map<eleside_t,tagbcs_t>::iterator p=eleside2tag_to_fbc.begin(); p!=eleside2tag_to_fbc.end(); ++p)
    {
        Element      * ele      = p->first.first;
        int            idx_side = p->first.second;
        int            bc_tag   = p->second.first;
        SDPair const & bcs      = p->second.second;
        pCalcM         calcm    = &Multiplier;

        if (bcs.HasKey("mfunc")) // callback specified
        {
            MDatabase_t::const_iterator q = MFuncs.find(bc_tag);
            if (q!=MFuncs.end()) calcm = q->second;
            else throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary Tag=%d was not found in MFuncs database",bc_tag);
        }

        ele->SetBCs (idx_side, bcs, pF, pU, calcm);
    }

    // set F bcs at nodes
    for (std::map<Node*,tagbcs_t>::iterator p=nod2tag_to_fbc.begin(); p!=nod2tag_to_fbc.end(); ++p)
    {
        Node         * nod    = p->first;
        int            bc_tag = p->second.first;
        SDPair const & bcs    = p->second.second;
        pCalcM         calcm  = &Multiplier;

        if (bcs.HasKey("mfunc")) // callback specified
        {
            MDatabase_t::const_iterator q = MFuncs.find(bc_tag);
            if (q!=MFuncs.end()) calcm = q->second;
            else throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary Tag=%d was not found in MFuncs database",bc_tag);
        }

        for (StrDbl_t::const_iterator q=bcs.begin(); q!=bcs.end(); ++q)
        {
            pF[nod].first[nod->FMap(q->first)] += q->second;
            pF[nod].second = calcm;
        }
    }

    // set U bcs at sides (edges/faces) of elements
    for (std::map<eleside_t,tagbcs_t>::iterator p=eleside2tag_to_ubc.begin(); p!=eleside2tag_to_ubc.end(); ++p)
    {
        Element      * ele      = p->first.first;
        int            idx_side = p->first.second;
        SDPair const & bcs      = p->second.second;

        ele->SetBCs (idx_side, bcs, pF, pU, NULL);
    }

    // set U bcs at nodes
    for (std::map<Node*,tagbcs_t>::iterator p=nod2tag_to_ubc.begin(); p!=nod2tag_to_ubc.end(); ++p)
    {
        Node         * nod = p->first;
        SDPair const & bcs = p->second.second;

        if (bcs.HasKey("inclsupport"))
        {
            if (NDim!=2) throw new Fatal("FEM::Domain::SetBCs: Inclined support is not implemented for 3D problems yet");
            InclSupport[nod] = bcs("alpha");
        }
        else
        {
            for (StrDbl_t::const_iterator q=bcs.begin(); q!=bcs.end(); ++q)
                pU[nod].first[nod->UMap(q->first)] = q->second;
        }
    }
}

inline void Domain::ClrBCs ()
{
    for (size_t i=0; i<Eles.Size(); ++i) Eles[i]->ClrBCs ();
    pF.clear();
    pU.clear();
    InclSupport.clear();
}

inline void Domain::Gravity ()
{
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        if (Eles[i]->Active)
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
}

inline void Domain::Deactivate (int EleTag)
{
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        int tag = Eles[i]->Cell.Tag;
        if (tag==EleTag)
        {
            pCalcM calcm = &Multiplier;
            if (Prps(tag).HasKey("mfunc"))
            {
                MDatabase_t::const_iterator p = MFuncs.find(tag);
                if (p!=MFuncs.end()) calcm = p->second;
                else throw new Fatal("Domain::Deactivate: Multiplier function with tag=%d was not found in MFuncs database",tag);
            }
            Eles[i]->Deactivate (pF, calcm, gAccel, pU);
        }
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
            Array<String> keys;
            Eles[ele]->StateKeys (keys);
            (*of) << Util::_8s << "Time";
            (*of) << Util::_8s << "x";
            (*of) << Util::_8s << "y";  if (NDim==3)
            (*of) << Util::_8s << "z";
            for (size_t j=0; j<keys.Size(); ++j) (*of) << Util::_8s << keys[j];
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
        Eles[ele]->StateAtCt (dat);
        Eles[ele]->Centroid  (Xct);
        (*FilEles[i]) << Util::_8s << Xct(0);
        (*FilEles[i]) << Util::_8s << Xct(1);  if (NDim==3)
        (*FilEles[i]) << Util::_8s << Xct(2);
        for (size_t j=0; j<dat.Keys.Size(); ++j) (*FilEles[i]) << Util::_8s << dat(dat.Keys[j]);
        (*FilEles[i]) << "\n";
    }
}

inline void Domain::PrintResults (char const * NF) const
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
        Eles[i]->StateAtCt (dat);
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
        Eles[i]->StateAtCt (dat);
        Eles[i]->Centroid  (X);
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

inline bool Domain::CheckError (Table const & NodSol, SDPair const & NodTol) const
{
    // header
    std::cout << "\n[1;37m--- Error Summary --- nodes --------------------------------------------------\n";
    std::cout << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm" << "[0m\n";

    // results
    bool error = false;

    // nodes
    for (size_t i=0; i<NodSol.Keys.Size(); ++i)
    {
        // calc error
        String key = NodSol.Keys[i];
        Array<double> err(NodSol.NRows);
        for (size_t j=0; j<Nods.Size(); ++j)
        {
            if (Nods[j]->UMap.HasKey(key)) err[j] = fabs(Nods[j]->U[Nods[j]->UMap(key)] - NodSol(key,j));
            else                           err[j] = fabs(Nods[j]->F[Nods[j]->FMap(key)] - NodSol(key,j));
        }

        // summary
        double max_err = err[err.Max()];
        double tol     = NodTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err[err.Min()] << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? "[1;31m" : "[1;32m") << Util::_8s<<max_err << "[0m" << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    return error;
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
            Eles[j]->StateAtCt (dat);
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
            Eles[j]->StateAtIPs (res);
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

inline void Domain::WriteMPY (char const * FNKey, double SFCoef, char const * Extra) const
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
    if (Extra!=NULL) of << Extra;
    MPL::SaveFig  (FNKey, of);
    of.close      ();
}

inline void Domain::WriteVTU (char const * FNKey) const
{
    // data
    size_t nn     = Nods.Size(); // number of nodes
    size_t ne     = Eles.Size(); // number of elements
    size_t ne_act = 0;           // number of active elements

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
        if (Nods[i]->NShares>0)
        {
            for (size_t j=0; j<Nods[i]->UMap.Keys.Size(); ++j)
            {
                if (nod_keys.Find(Nods[i]->UMap.Keys[j])<0) nod_keys.Push (Nods[i]->UMap.Keys[j]);
            }
            if (Nods[i]->UMap.HasKey("ux") || Nods[i]->UMap.HasKey("uy") || Nods[i]->UMap.HasKey("uz")) nod_displacements = true;
            if (Nods[i]->UMap.HasKey("vx") || Nods[i]->UMap.HasKey("vy") || Nods[i]->UMap.HasKey("vz")) nod_velocities    = true;
        }
    }
    for (size_t i=0; i<ne; ++i)
    {
        if (Eles[i]->Active)
        {
            Array<String> keys;
            Eles[i]->StateKeys (keys);
            for (size_t j=0; j<keys.Size(); ++j)
            {
                if (ele_keys.Find(keys[j])<0) ele_keys.Push (keys[j]);
            }
            if (keys.Find("vx")>=0 || keys.Find("vy")>=0 || keys.Find("vz")>=0) ele_velocities = true;
            if (keys.Find("gx")>=0 || keys.Find("gy")>=0 || keys.Find("gz")>=0) ele_gradients  = true;
            ne_act++;
        }
    }

    // extrapolate data
    size_t nreskeys = ele_keys.Size(); // number of results keys
    Mat_t nod_results (nn,nreskeys);   // results at nodes
    Mat_t nod_count   (nn,nreskeys);   // count how many times a variable was added to this node
    set_to_zero (nod_results);
    set_to_zero (nod_count);
    for (size_t i=0; i<ne; ++i)
    {
        if (Eles[i]->Active)
        {
            Array<SDPair> loc_res; // local results: size==number of nodes in element
            Eles[i]->StateAtNodes (loc_res);
            for (size_t j=0; j<nreskeys; ++j)
            {
                if (loc_res[0].HasKey(ele_keys[j]))
                {
                    for (size_t k=0; k<Eles[i]->Con.Size(); ++k)
                    {
                        size_t vid = Eles[i]->Con[k]->Vert.ID;
                        nod_results(vid, j) += loc_res[k](ele_keys[j]);
                        nod_count  (vid, j) += 1.0;
                    }
                }
            }
        }
    }

    // header
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\"?>\n";
    oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    oss << "  <UnstructuredGrid>\n";
    oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << ne_act << "\">\n";

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
        if (Eles[i]->Active)
        {
            oss << "  ";
            for (size_t j=0; j<Eles[i]->Con.Size(); ++j) oss << Eles[i]->Con[j]->Vert.ID << " ";
            k++;
            VTU_NEWLINE (i,k,ne,nimax/Eles[i]->Con.Size(),oss);
        }
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    size_t offset = 0;
    for (size_t i=0; i<ne; ++i)
    {
        if (Eles[i]->Active)
        {
            offset += Eles[i]->Con.Size();
            oss << (k==0?"  ":" ") << offset;
            k++;
            VTU_NEWLINE (i,k,ne,nimax,oss);
        }
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<ne; ++i)
    {
        if (Eles[i]->Active)
        {
            if (NDim==2) oss << (k==0?"  ":" ") << NVertsToVTKCell2D[Eles[i]->Con.Size()];
            else         oss << (k==0?"  ":" ") << NVertsToVTKCell3D[Eles[i]->Con.Size()];
            k++;
            VTU_NEWLINE (i,k,ne,nimax,oss);
        }
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
                double val = (Nods[j]->NShares>0 ? Nods[j]->U[Nods[j]->UMap(nod_keys[i])] : 0.0);
                oss << (k==0?"  ":" ") << val;
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
            if (Nods[i]->NShares>0)
            {
                oss << "  " << nsflo <<          Nods[i]->U[Nods[i]->UMap("ux")] << " ";
                oss <<         nsflo <<          Nods[i]->U[Nods[i]->UMap("uy")] << " ";
                oss <<         nsflo << (NDim==3?Nods[i]->U[Nods[i]->UMap("uz")]:0.0);
            }
            else
            {
                oss << "  " << nsflo << 0.0 << " ";
                oss <<         nsflo << 0.0 << " ";
                oss <<         nsflo << 0.0;
            }
            k++;
            VTU_NEWLINE (i,k,nn,nfmax/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }
    for (size_t i=0; i<ele_keys.Size(); ++i)
    {
        oss << "        <DataArray type=\"Float32\" Name=\"" << ele_keys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t j=0; j<nn; ++j)
        {
            double den = (nod_count(j,i)>0.0 ? nod_count(j,i) : 1.0);
            oss << (k==0?"  ":" ") << nod_results(j,i)/den;
            k++;
            VTU_NEWLINE (j,k,nn,nfmax-1,oss);
        }
        oss << "        </DataArray>\n";
    }
    oss << "      </PointData>\n";

    // Bottom
    oss << "    </Piece>\n";
    oss << "  </UnstructuredGrid>\n";
    oss << "</VTKFile>" << std::endl;

    // Write to file
    String fn(FNKey); fn.append(".vtu");
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
