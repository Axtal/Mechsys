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
#include <mechsys/util/stopwatch.h>
#include <mechsys/models/model.h>
#include <mechsys/fem/element.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/draw.h>

namespace FEM
{

typedef std::map<Node*,double> InclSupport_t; ///< Inclined support type. Maps Node ==> alpha
typedef std::map<Node*,Array<double> > Res_t;

inline double Multiplier (double t) { return 1.0; }

class Domain
{
public:
    // static
    static bool PARA; ///< Parallel code ?

    // typedefs
    typedef std::map<int,pCalcM> MDatabase_t; ///< Map tag to M function pointer
    typedef std::map<int,Model*> Models_t;    ///< Map tag to model pointer

    // Constructor
    Domain (Mesh::Generic const & Msh,        ///< The mesh
            Dict const          & Prps,       ///< Element properties
            Dict const          & Mdls,       ///< Model names and parameters
            Dict const          & Inis,       ///< Initial values
            char const          * FKey=NULL,  ///< File key to be used during output of results
            Array<int>    const * OutV=NULL,  ///< IDs or Tags of vertices to generate output
            Array<int>    const * OutC=NULL); ///< IDs or Tags of cells to generate output

    // Destructor
    ~Domain ();

    // Methods
    void SetBCs       (Dict const & BCs);
    void ClrBCs       ();
    void Gravity      ();                                                        ///< Apply gravity
    void Deactivate   (int EleTag);                                              ///< Deactivate all elements with EleTag
    void SetUVals     (SDPair const & UVals);                                    ///< Set U values
    void OutResults   (double Time, Vec_t const & F_int) const;                  ///< Do output results
    void PrintResults (char const * NF="%15.6g", bool WithElems=true) const;     ///< Print results
    bool CheckError   (Table const & NodSol, SDPair const & NodTol) const;       ///< At nodes
    bool CheckError   (Table const & NodSol, Table const & EleSol, 
                       SDPair const & NodTol, SDPair const & EleTol) const;      ///< At nodes and centroid
    bool CheckErrorIP (Table const & EleSol, SDPair const & EleTol) const;       ///< At integration points
    void WriteMPY     (char const * FileKey, double SFCoef=1.0, bool PNG=false,
                       char const * Extra=NULL) const;                           ///< SFCoef: Scale-factor coefficient
    void AvailableData();                                                        ///< Check available data and resize results matrices
    void NodalResults () const;                                                  ///< Extrapolate results from element to nodes
    void WriteVTU     (char const * FileKey) const;                              ///< Write file for ParaView

    // Data
    Dict          const & Prps;        ///< Element properties
    Dict          const & Inis;        ///< Initial values
    int                   NDim;        ///< Space dimension
    double                gAccel;      ///< Gravity acceleration
    Models_t              Mdls;        ///< Models
    Array<Node*>          Nods;        ///< (Allocated memory) Nodes
    Array<Element*>       Eles;        ///< (Allocated memory) Elements
    Array<Node*>          TgdNods;     ///< Tagged Nodes (at boundaries)
    Array<Element*>       TgdEles;     ///< Tagged Elements (at boundaries)
    Array<Node*>          OutNods;     ///< Nodes for which output (results) is generated
    Array<Element*>       OutEles;     ///< Elements for which output (results) is generated
    Array<std::ofstream*> FilNods;     ///< Files with results at selected nodes (OutNods)
    Array<std::ofstream*> FilEles;     ///< Files with results at selected elements (OutEles)
    Array<size_t>         Beams;       ///< Subset of elements of type Beam
    NodBCs_t              pU;          ///< Nodes with prescribed U. The values are the Delta over the previous stage
    NodBCs_t              pF;          ///< Nodes with prescribed F. The values are the Delta over the previous stage
    MDatabase_t           MFuncs;      ///< Database of pointers to M functions
    Array<String>         DisplKeys;   ///< Displacement keys
    InclSupport_t         InclSupport; ///< Inclined support
#ifdef HAS_MPI
    Array<Node*>          InterNodes;  ///< Nodes on the inferface between partitions
#endif

    // Nodal results
    size_t        NActNods;    ///< Number of active nodes
    size_t        NActEles;    ///< Number of active elements
    bool          HasDisps;    ///< Has displacements (ux, ...) ?
    bool          HasVeloc;    ///< Has velocities (vx, ...) ?
    Array<String> NodKeys;     ///< All node keys (ux, uy, T, ...)
    Array<String> EleKeys;     ///< All element keys (sx, sy, vx, vy, ...)
    mutable Res_t NodResults;  ///< Extrapolated nodal results
    mutable Res_t NodResCount; ///< Count how many times a variable (key) was added to a node

#ifdef USE_BOOST_PYTHON
    void PySetOutNods (BPy::str const & FileKey, BPy::list const & IDsOrTags) { SetOutNods (BPy::extract<char const *>(FileKey)(), Array<int>(IDsOrTags)); }
    void PySetOutEles (BPy::str const & FileKey, BPy::list const & IDsOrTags) { SetOutEles (BPy::extract<char const *>(FileKey)(), Array<int>(IDsOrTags)); }
#endif
};

bool Domain::PARA = false;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


std::ostream & operator<< (std::ostream & os, Domain const & D)
{
    String buf;
    buf.Printf("\n%s--- Models -------------------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (Domain::Models_t::const_iterator p=D.Mdls.begin(); p!=D.Mdls.end(); ++p)
        os << p->first << " " << (*p->second) << std::endl;

    buf.Printf("\n%s--- Elements properties ------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    os << D.Prps << std::endl;

    buf.Printf("\n%s--- Initial values -----------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    os << D.Inis << std::endl;

    buf.Printf("\n%s--- Nodes --------------------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (size_t i=0; i<D.Nods.Size(); ++i) os << (*D.Nods[i]) << "\n";

    buf.Printf("\n%s--- Elements -----------------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (size_t i=0; i<D.Eles.Size(); ++i) os << (*D.Eles[i]) << "\n";

    buf.Printf("\n%s--- Boundary Conditions ------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
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

inline Domain::Domain (Mesh::Generic const & Msh, Dict const & ThePrps, Dict const & TheMdls, Dict const & TheInis, char const * FKey, Array<int> const * OutV, Array<int> const * OutC)
    : Prps(ThePrps), Inis(TheInis), NDim(Msh.NDim), gAccel(9.81)
{
    // info
    Util::Stopwatch stopwatch(/*only_root*/PARA);
    bool root = true;
#ifdef HAS_MPI
    if (PARA && MPI::COMM_WORLD.Get_rank()!=0) root = false;
#endif
    if (root) printf("\n%s--- Domain --- allocating nodes and elements ---------------------------------------%s\n",TERM_CLR1,TERM_RST);

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
    std::map<int,Node*> VertID2Node;
    for (size_t i=0; i<Msh.Verts.Size(); ++i)
    {
#ifdef HAS_MPI
        if (PARA)
        {
            // skip vertices that aren't in this partition
            if (Msh.Verts[i]->PartIDs.Find(MPI::COMM_WORLD.Get_rank())<0) continue;
        }
#endif
        // new node
        Nods.Push (new Node((*Msh.Verts[i])));
        if (Msh.Verts[i]->Tag<0) TgdNods.Push (Nods.Last());
        VertID2Node[Msh.Verts[i]->ID] = Nods.Last();

#ifdef HAS_MPI
        if (PARA)
        {
            // add DOFs to interface nodes (between partitions)
            if (Msh.Verts[i]->PartIDs.Size()>1) // this is a shared node (between partitions)
            {
                InterNodes.Push (Nods.Last());
                for (size_t k=0; k<Msh.Verts[i]->Shares.Size(); ++k)
                {
                    int elem_tag = Msh.Verts[i]->Shares[k].C->Tag;
                    String prob_name_ND;
                    PROB.Val2Key (Prps(elem_tag)("prob"), prob_name_ND);
                    prob_name_ND.Printf("%s%dD", prob_name_ND.CStr(), NDim);
                    Nods.Last()->AddDOF (ElementVarKeys[prob_name_ND].first.CStr(), ElementVarKeys[prob_name_ND].second.CStr());
                }
            }
        }
#endif
        // set list of nodes for output
        if (FKey!=NULL && OutV!=NULL)
        {
            if (OutV->Find(Msh.Verts[i]->ID)>=0 || OutV->Find(Msh.Verts[i]->Tag)>=0)
            {
                String buf; buf.Printf("%s_nod_%d_%d.res", FKey, Msh.Verts[i]->ID, Msh.Verts[i]->Tag);
                std::ofstream * of = new std::ofstream (buf.CStr(),std::ios::out);
                OutNods.Push (Nods.Last());
                FilNods.Push (of);
            }
        }
    }

    // set elements from mesh
    for (size_t i=0; i<Msh.Cells.Size(); ++i)
    {
#ifdef HAS_MPI
        if (PARA)
        {
            // skip elements that aren't in this partition
            if (Msh.Cells[i]->PartID!=MPI::COMM_WORLD.Get_rank()) continue;
        }
#endif
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

            // connectivity
            Array<Node*> nodes(Msh.Cells[i]->V.Size());
            for (size_t k=0; k<nodes.Size(); ++k) nodes[k] = VertID2Node[Msh.Cells[i]->V[k]->ID];

            // allocate element
            if (Inis.HasKey(tag)) Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, Prps(tag), Inis(tag), nodes));
            else                  Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, Prps(tag), SDPair() , nodes));

            // set array of Beams
            if (prob_name=="Beam") Beams.Push (Eles.Size()-1);

            // tagged elements
            if (Msh.Cells[i]->BryTags.size()>0 || prob_name=="Beam") TgdEles.Push (Eles.Last());

            // set list of elements for output
            if (FKey!=NULL && OutC!=NULL)
            {
                if (OutC->Find(Msh.Cells[i]->ID)>=0 || OutC->Find(Msh.Cells[i]->Tag)>=0)
                {
                    String buf; buf.Printf("%s_ele_%d_%d.res", FKey, Msh.Cells[i]->ID, Msh.Cells[i]->Tag);
                    std::ofstream * of = new std::ofstream (buf.CStr(),std::ios::out);
                    OutEles.Push (Eles.Last());
                    FilEles.Push (of);
                }
            }
        }
        else throw new Fatal("Domain::SetMesh: Dictionary of properties must have keyword 'prob' defining the type of element corresponding to a specific problem");
    }

    // divide U values in nodes by the number of times elements added to a U var
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        if (Nods[i]->NShares>0) // active
        {
            for (size_t j=0; j<Nods[i]->U.Size(); ++j)
            {
                if (Nods[i]->NaddU[j]>0.0) Nods[i]->U[j] /= Nods[i]->NaddU[j];
            }
        }
    }

    // displacement keys
    DisplKeys.Resize (NDim);
    DisplKeys[0] = "ux";
    DisplKeys[1] = "uy";  if (NDim==3)
    DisplKeys[2] = "uz";

    // check available data
    AvailableData ();
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
        for (size_t j=0; j<TgdEles.Size(); ++j)
        {
            if (TgdEles[j]->Active)
            {
                Mesh::BryTag_t const & eftags = TgdEles[j]->Cell.BryTags;
                for (Mesh::BryTag_t::const_iterator p=eftags.begin(); p!=eftags.end(); ++p)
                {
                    int idx_side = p->first;
                    int side_tag = p->second;
                    if (side_tag==bc_tag) // found
                    {
                        found = true;
                        eleside_t es(TgdEles[j],idx_side);
                        for (size_t k=0; k<bcs.Keys.Size(); ++k) // we have to split Ubcs from Fbcs
                        {
                            if (TgdEles[j]->UKeys.Find(bcs.Keys[k])>=0)
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
        for (size_t j=0; j<TgdNods.Size(); ++j)
        {
            if (TgdNods[j]->NShares>0) // active
            {
                if (TgdNods[j]->Vert.Tag==bc_tag) // found
                {
                    found = true;
                    for (size_t k=0; k<bcs.Keys.Size(); ++k) // we have to split Ubcs from Fbcs
                    {
                        if (bcs.Keys[k]=="inclsupport")
                        {
                            nod2tag_to_ubc[TgdNods[j]].first  = bc_tag;
                            nod2tag_to_ubc[TgdNods[j]].second = bcs;
                            break;
                        }
                        else
                        {
                            if (TgdNods[j]->UMap.HasKey(bcs.Keys[k]))
                            {
                                nod2tag_to_ubc[TgdNods[j]].first = bc_tag;
                                nod2tag_to_ubc[TgdNods[j]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else
                            {
                                nod2tag_to_fbc[TgdNods[j]].first = bc_tag;
                                nod2tag_to_fbc[TgdNods[j]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                        }
                    }
                }
            }
        }

        if (!found) 
        {
            // find elements with tags equal to this bc_tag. Special cases: qn of beams, s (source term), cbx (cetrifugal body force), ...
            for (size_t j=0; j<Eles.Size(); ++j)
            {
                if (Eles[j]->Active)
                {
                    if (Eles[j]->Cell.Tag==bc_tag) // found
                    {
                        // problem name
                        String prob_name_ND;
                        PROB.Val2Key (Prps(Eles[j]->Cell.Tag)("prob"), prob_name_ND);
                        prob_name_ND.Printf("%s%dD", prob_name_ND.CStr(), NDim);

                        // check if bc key is not U or F
                        for (size_t k=0; k<bcs.Keys.Size(); ++k)
                        {
                            if (Util::HasKey(ElementVarKeys[prob_name_ND].first, bcs.Keys[k])) throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' with tag==%d cannot be applied to the element (%d,%d) itself", bcs.Keys[k].CStr(), bc_tag, Eles[j]->Cell.ID, Eles[j]->Cell.Tag);
                            if (Util::HasKey(ElementVarKeys[prob_name_ND].second,bcs.Keys[k])) throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' with tag==%d cannot be applied to the element (%d,%d) itself", bcs.Keys[k].CStr(), bc_tag, Eles[j]->Cell.ID, Eles[j]->Cell.Tag);
                        }

                        // map bcs
                        found = true;
                        int idx_side = 0; // irrelevant
                        eleside_t es(Eles[j],idx_side);
                        eleside2tag_to_fbc[es].first  = bc_tag;
                        eleside2tag_to_fbc[es].second = bcs;
                    }
                }
            }
            if (!PARA && !found) throw new Fatal("FEM::Domain::SetBCs: Could not find any edge/face of element or node (or line element) with boundary Tag=%d",bc_tag);
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
            if (q->first!="mfunc")
            {
                pF[nod].first[nod->FMap(q->first)] += q->second;
                pF[nod].second = calcm;
            }
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

    // check available data
    AvailableData ();

    // write header to output files
    // TODO: this should be done at each stage, since new data may become available
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

inline void Domain::OutResults (double Time, Vec_t const & F_int) const
{
    // nodes
    for (size_t i=0; i<OutNods.Size(); ++i)
    {
        if (Time<1.0e-14)
        {
            String buf;
            (*FilNods[i]) << Util::_8s << "Time";
            (*FilNods[i]) << Util::_8s << "x";
            (*FilNods[i]) << Util::_8s << "y";  if (NDim==3)
            (*FilNods[i]) << Util::_8s << "z";
            for (size_t j=0; j<OutNods[i]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->UMap.Keys[j];
            for (size_t j=0; j<OutNods[i]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->FMap.Keys[j];
            for (size_t j=0; j<OutNods[i]->nDOF(); ++j)
            {
                buf.Printf ("%s_int", OutNods[i]->FMap.Keys[j].CStr());
                (*FilNods[i]) << Util::_8s << buf;
            }
            (*FilNods[i]) << "\n";
        }
        (*FilNods[i]) << Util::_8s << Time;
        (*FilNods[i]) << Util::_8s << OutNods[i]->Vert.C(0);
        (*FilNods[i]) << Util::_8s << OutNods[i]->Vert.C(1);  if (NDim==3)
        (*FilNods[i]) << Util::_8s << OutNods[i]->Vert.C(2);
        for (size_t j=0; j<OutNods[i]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->U[j];
        for (size_t j=0; j<OutNods[i]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->F[j];
        for (size_t j=0; j<OutNods[i]->nDOF(); ++j) (*FilNods[i]) << Util::_8s << F_int(OutNods[i]->EQ[j]);
        (*FilNods[i]) << "\n";
    }

    // elements
    for (size_t i=0; i<OutEles.Size(); ++i)
    {
        if (Time<1.0e-14)
        {
            Array<String> keys;
            OutEles[i]->StateKeys (keys);
            (*FilEles[i]) << Util::_8s << "Time";
            (*FilEles[i]) << Util::_8s << "x";
            (*FilEles[i]) << Util::_8s << "y";  if (NDim==3)
            (*FilEles[i]) << Util::_8s << "z";
            for (size_t j=0; j<keys.Size(); ++j) (*FilEles[i]) << Util::_8s << keys[j];
            (*FilEles[i]) << "\n";
        }
        (*FilEles[i]) << Util::_8s << Time;
        SDPair dat;
        Vec_t  Xct;
        OutEles[i]->StateAtCt (dat);
        OutEles[i]->Centroid  (Xct);
        (*FilEles[i]) << Util::_8s << Xct(0);
        (*FilEles[i]) << Util::_8s << Xct(1);  if (NDim==3)
        (*FilEles[i]) << Util::_8s << Xct(2);
        for (size_t j=0; j<dat.Keys.Size(); ++j) (*FilEles[i]) << Util::_8s << dat(dat.Keys[j]);
        (*FilEles[i]) << "\n";
    }
}

inline void Domain::PrintResults (char const * NF, bool WithElems) const
{
    printf("\n%s--- Results ------------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

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
    std::cout << TERM_CLR2 << Util::_6 << "Node";
    for (size_t i=0; i<ukeys.Size(); ++i) { buf.Printf(nf, ukeys[i].CStr());  std::cout<<buf; }
    for (size_t i=0; i<fkeys.Size(); ++i) { buf.Printf(nf, fkeys[i].CStr());  std::cout<<buf; }
    for (size_t i=0; i<ukeys.Size(); ++i)
    {
        buf.Printf ("R%s",ukeys[i].CStr());
        String tmp;
        tmp.Printf (nf,buf.CStr());
        std::cout << tmp;
    }
    printf("%s\n",TERM_RST);

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
    if (WithElems)
    {
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
        std::cout << TERM_CLR2 << Util::_6 << "Elem";
        buf.Printf(nf, "x");  std::cout << buf;
        buf.Printf(nf, "y");  std::cout << buf;  if (NDim==3) {
        buf.Printf(nf, "z");  std::cout << buf; }
        for (size_t i=0; i<keys.Size(); ++i) { buf.Printf(nf, keys[i].CStr());  std::cout<<buf; }
        std::cout << TERM_RST << std::endl;

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
}

inline bool Domain::CheckError (Table const & NodSol, SDPair const & NodTol) const
{
    // header
    printf("\n%s--- Error Summary --- nodes --------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    std::cout << TERM_CLR2 << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm";
    printf("%s\n",TERM_RST);

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
        double max_err = err.TheMax();
        double tol     = NodTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
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
    printf("\n%s--- Error Summary --- nodes and centroid -------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    std::cout << TERM_CLR2 << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm";
    printf("%s\n",TERM_RST);

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
        double max_err = err.TheMax();
        double tol     = NodTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
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
        double max_err = err.TheMax();
        double tol     = EleTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    return error;
}

inline bool Domain::CheckErrorIP (Table const & EleSol, SDPair const & EleTol) const
{
    // header
    printf("\n%s--- Error Summary --- integration points -------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    std::cout << TERM_CLR2 << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm";
    printf("%s\n",TERM_RST);

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
        double max_err = err.TheMax();
        double tol     = EleTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    std::cout << "\n";

    return error;
}

inline void Domain::WriteMPY (char const * FNKey, double SFCoef, bool PNG, char const * Extra) const
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
    if (PNG) MPL::SaveFig (FNKey, of);
    else     MPL::Show    (of);
    of.close ();
}

inline void Domain::AvailableData ()
{
    // check available nodal data
    NActNods = 0;
    HasDisps = false;
    String key;
    NodKeys.Resize (0);
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        if (Nods[i]->NShares>0) // active node
        {
            for (size_t j=0; j<Nods[i]->UMap.Keys.Size(); ++j)
            {
                key = Nods[i]->UMap.Keys[j];
                if (NodKeys.Find(key)<0) NodKeys.Push (key);
                if (key=="ux") HasDisps = true;
            }
            NActNods++;
        }
    }

    // check available element data
    NActEles = 0;
    HasVeloc = false;
    Array<String> keys;
    EleKeys.Resize (0);
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        if (Eles[i]->Active) // active element
        {
            Eles[i]->StateKeys (keys);
            for (size_t j=0; j<keys.Size(); ++j)
            {
                if (EleKeys.Find(keys[j])<0) EleKeys.Push (keys[j]);
                if (keys[j]=="vx") HasVeloc = true;
            }
            NActEles++;
        }
    }

    // resize matrices
    //NodResults .change_dim (Nods.Size(), EleKeys.Size()); // results at nodes
    //NodResCount.change_dim (Nods.Size(), EleKeys.Size()); // times a var was added to a node
}

inline void Domain::NodalResults () const
{
    // extrapolate data
    for (Res_t::iterator p=NodResults .begin(); p!=NodResults .end(); ++p) p->second.SetValues(0.0);
    for (Res_t::iterator p=NodResCount.begin(); p!=NodResCount.end(); ++p) p->second.SetValues(0.0);
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        if (Eles[i]->Active) // active element
        {
            Array<SDPair> loc_res; // local results: size==number of nodes in element
            Eles[i]->StateAtNodes (loc_res);
            for (size_t j=0; j<EleKeys.Size(); ++j)
            {
                if (loc_res[0].HasKey(EleKeys[j]))
                {
                    for (size_t k=0; k<Eles[i]->Con.Size(); ++k)
                    {
                        NodResults [Eles[i]->Con[k]][j] += loc_res[k](EleKeys[j]);
                        NodResCount[Eles[i]->Con[k]][j] += 1.0;
                    }
                }
            }
        }
    }
}

inline void Domain::WriteVTU (char const * FNKey) const
{
    // extrapolate results
    //NodalResults ();

    // data
    size_t nn = Nods.Size(); // number of nodes
    size_t ne = Eles.Size(); // number of elements

    // header
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\"?>\n";
    oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    oss << "  <UnstructuredGrid>\n";
    oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << NActEles << "\">\n";

    // nodes: coordinates
    oss << "      <Points>\n";
    oss << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    size_t k = 0; oss << "        ";
    std::map<int,int> VertID2LocID;
    for (size_t i=0; i<nn; ++i)
    {
        VertID2LocID[Nods[i]->Vert.ID] = i;
        oss << "  " << Util::_8s <<          Nods[i]->Vert.C(0) << " ";
        oss <<         Util::_8s <<          Nods[i]->Vert.C(1) << " ";
        oss <<         Util::_8s << (NDim==3?Nods[i]->Vert.C(2):0.0);
        k++;
        VTU_NEWLINE (i,k,nn,6/3-1,oss);
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
            size_t nne = (NDim==2 ? NVertsToVTKNVerts2D[Eles[i]->Con.Size()] : NVertsToVTKNVerts3D[Eles[i]->Con.Size()]);
            oss << "  ";
            for (size_t j=0; j<nne; ++j) oss << VertID2LocID[Eles[i]->Con[j]->Vert.ID] << " ";
            k++;
            VTU_NEWLINE (i,k,ne,40/Eles[i]->Con.Size(),oss);
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
            size_t nne = (NDim==2 ? NVertsToVTKNVerts2D[Eles[i]->Con.Size()] : NVertsToVTKNVerts3D[Eles[i]->Con.Size()]);
            offset += nne;
            oss << (k==0?"  ":" ") << offset;
            k++;
            VTU_NEWLINE (i,k,ne,40,oss);
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
            VTU_NEWLINE (i,k,ne,40,oss);
        }
    }
    oss << "        </DataArray>\n";
    oss << "      </Cells>\n";

    // data -- nodes
    oss << "      <PointData Scalars=\"point_scalars\">\n";
    for (size_t i=0; i<NodKeys.Size(); ++i)
    {
        if (!(NodKeys[i]=="ux" || NodKeys[i]=="uy" || NodKeys[i]=="uz"))
        {
            oss << "        <DataArray type=\"Float32\" Name=\"" << NodKeys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                double val = (Nods[j]->NShares>0 ? Nods[j]->U[Nods[j]->UMap(NodKeys[i])] : 0.0);
                oss << (k==0?"  ":" ") << val;
                k++;
                VTU_NEWLINE (j,k,nn,6-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }

    // data -- nodes -- displacements
    if (HasDisps)
    {
        oss << "        <DataArray type=\"Float32\" Name=\"" << "u" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<nn; ++i)
        {
            if (Nods[i]->NShares>0)
            {
                oss << "  " << Util::_8s <<          Nods[i]->U[Nods[i]->UMap("ux")] << " ";
                oss <<         Util::_8s <<          Nods[i]->U[Nods[i]->UMap("uy")] << " ";
                oss <<         Util::_8s << (NDim==3?Nods[i]->U[Nods[i]->UMap("uz")]:0.0);
            }
            else
            {
                oss << "  " << Util::_8s << 0.0 << " ";
                oss <<         Util::_8s << 0.0 << " ";
                oss <<         Util::_8s << 0.0;
            }
            k++;
            VTU_NEWLINE (i,k,nn,6/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }

    // data -- nodes -- extrapolated values
    /*
    for (size_t i=0; i<EleKeys.Size(); ++i)
    {
        if (!(EleKeys[i]=="vx" || EleKeys[i]=="vy" || EleKeys[i]=="vz"))
        {
            oss << "        <DataArray type=\"Float32\" Name=\"" << EleKeys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                //Nods[j]
                //double den = (NodResCount(j,i)>0.0 ? NodResCount(j,i) : 1.0);
                //oss << (k==0?"  ":" ") << NodResults(j,i)/den;
                k++;
                VTU_NEWLINE (j,k,nn,6-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }
    */

    // data -- nodes -- velocities
    if (HasVeloc)
    {
        size_t idx_vx=-1, idx_vy=-1, idx_vz=-1;
        for (size_t i=0; i<EleKeys.Size(); ++i)
        {
            if (EleKeys[i]=="vx") idx_vx = i;
            if (EleKeys[i]=="vy") idx_vy = i;
            if (EleKeys[i]=="vz") idx_vz = i;
        }
        oss << "        <DataArray type=\"Float32\" Name=\"" << "v" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t j=0; j<nn; ++j)
        {
            if (Nods[j]->NShares>0)
            {
                //oss << "  " << Util::_8s <<            NodResults(j,idx_vx) / (NodResCount(j,idx_vx)>0.0 ? NodResCount(j,idx_vx) : 1.0) << " ";
                //oss <<         Util::_8s <<            NodResults(j,idx_vy) / (NodResCount(j,idx_vy)>0.0 ? NodResCount(j,idx_vy) : 1.0) << " ";
                //oss <<         Util::_8s << (NDim==3 ? NodResults(j,idx_vz) / (NodResCount(j,idx_vz)>0.0 ? NodResCount(j,idx_vz) : 1.0) : 0.0);
            }
            else
            {
                oss << "  " << Util::_8s << 0.0 << " ";
                oss <<         Util::_8s << 0.0 << " ";
                oss <<         Util::_8s << 0.0;
            }
            k++;
            VTU_NEWLINE (j,k,nn,6/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }

    // data -- nodes -- end
    oss << "      </PointData>\n";

    // bottom
    oss << "    </Piece>\n";
    oss << "  </UnstructuredGrid>\n";
    oss << "</VTKFile>" << std::endl;

    // write to file
    String fn(FNKey); fn.append(".vtu");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();
}

}; // namespace FEM

#endif // MECHSYS_DOMAIN_H
