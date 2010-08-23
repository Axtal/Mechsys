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

typedef std::map<Node*,double>          InclSupport_t; ///< Inclined support type. Maps Node ==> alpha
typedef std::map<Node*,Array<double> >  Res_t;         ///< Maps node to results at nodes
typedef std::map<int, Array<Element*> > Tag2Eles_t;    ///< Maps tag to elements

inline double Multiplier (double t) { return 1.0; }

class Domain
{
public:
    // static
    static bool PARA;     ///< Parallel code ?
    static bool WithInfo; ///< Print information ?

    // enum
    enum BryTagType_t { None_t, Element_t, Border_t, Node_t }; ///< type of boundary condition tag

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
    void SetBCs         (Dict const & BCs);                                                                                ///< Set boundary conditions
    void PrintResults   (char const * NF="%15.6e", bool OnlySummary=false, bool WithElems=true, double Tol=1.0e-10) const; ///< Print results (Tol:tolerance to ignore zeros)
    void WriteMPY       (char const * FileKey, double SFCoef=1.0, bool PNG=false, char const * Extra=NULL) const;          ///< SFCoef: Scale-factor coefficient
    void WriteVTU       (char const * FileKey) const;                                                                      ///< Write file for ParaView
    bool CheckErrorNods (Table const & NodSol, SDPair const & NodTol) const;                                               ///< Check error at nodes
    bool CheckErrorEles (Table const & EleSol, SDPair const & EleTol) const;                                               ///< Check error at the centre of elements
    bool CheckErrorIP   (Table const & EleSol, SDPair const & EleTol) const;                                               ///< Check error at integration points of elements

    // Internal methods
    void OutResults    (double Time, Vec_t const & F_int) const;              ///< Do output results
    void CalcReactions (Table & NodesReactions, SDPair & SumReactions) const; ///< Calculate reactions
    void AvailableData ();                                                    ///< Check available data and resize results matrices
    void NodalResults  () const;                                              ///< Extrapolate results from element to nodes

    // Data
    Dict          const & Prps;        ///< Element properties
    Dict          const & Inis;        ///< Initial values
    int                   NDim;        ///< Space dimension
    Models_t              Mdls;        ///< Models
    Array<Node*>          Nods;        ///< (Allocated memory) Nodes
    Array<Element*>       Eles;        ///< (Allocated memory) Elements
    Array<Node*>          TgdNods;     ///< Tagged Nodes (at boundaries)
    Array<Element*>       TgdEles;     ///< Tagged edges or faces of Elements (at boundaries)
    Array<Node*>          OutNods;     ///< Nodes for which output (results) is generated
    Array<Element*>       OutEles;     ///< Elements for which output (results) is generated
    Array<std::ofstream*> FilNods;     ///< Files with results at selected nodes (OutNods)
    Array<std::ofstream*> FilEles;     ///< Files with results at selected elements (OutEles)
    Array<Element*>       Beams;       ///< Subset of elements of type Beam
    NodBCs_t              pU;          ///< Nodes with prescribed U. The values are the Delta over the previous stage
    NodBCs_t              pF;          ///< Nodes with prescribed F. The values are the Delta over the previous stage
    NodBCs_t              pFaccum;     ///< Accumulated pF (used to calculated reactions)
    MDatabase_t           MFuncs;      ///< Database of pointers to M functions
    Array<String>         DisplKeys;   ///< Displacement keys
    InclSupport_t         InclSupport; ///< Inclined support
    Array<Node*>          InterNodes;  ///< Nodes on the inferface between partitions (if PARA==true)
    Tag2Eles_t            Tag2Eles;    ///< Map tag to elements. Useful to activate/deactivate layers

    // Nodal results
    size_t        NActNods;    ///< Number of active nodes
    size_t        NActEles;    ///< Number of active elements
    bool          HasDisps;    ///< Has displacements (ux, ...) ?
    bool          HasVeloc;    ///< Has velocities (vx, ...) ?
    Array<String> AllUKeys;    ///< All U node keys (ux, uy, H, ...)
    Array<String> AllFKeys;    ///< All F node keys (fx, fy, Q, ...)
    Array<String> EleKeys;     ///< All element keys (sx, sy, vx, vy, ...)
    mutable Res_t NodResults;  ///< Extrapolated nodal results
    mutable Res_t NodResCount; ///< Count how many times a variable (key) was added to a node
};

bool Domain::PARA     = false;
bool Domain::WithInfo = true;


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

    buf.Printf("\n%s--- Elements by tags ---------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (Tag2Eles_t::const_iterator it=D.Tag2Eles.begin(); it!=D.Tag2Eles.end(); ++it)
    {
        os << "  Tag = " << it->first << " : Elements = ";
        for (size_t i=0; i<it->second.Size(); ++i) os << it->second[i]->Cell.ID << " ";
        os << "\n";
    }

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
    os << "\n  Accumulated prescribed F at nodes:\n";
    for (NodBCs_t::const_iterator p=D.pFaccum.begin(); p!=D.pFaccum.end(); ++p)
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

    return os;
}

// Constructor and destructor

inline Domain::Domain (Mesh::Generic const & Msh, Dict const & ThePrps, Dict const & TheMdls, Dict const & TheInis, char const * FKey, Array<int> const * OutV, Array<int> const * OutC)
    : Prps(ThePrps), Inis(TheInis), NDim(Msh.NDim)
{
    // info
#ifdef HAS_MPI
    if (PARA && MPI::COMM_WORLD.Get_rank()!=0) WithInfo = false;
#endif
    Util::Stopwatch stopwatch(/*activated*/WithInfo);
    if (WithInfo) printf("\n%s--- Domain --- allocating nodes and elements ---------------------------------------%s\n",TERM_CLR1,TERM_RST);

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

        // add DOFs to node
        for (size_t k=0; k<Msh.Verts[i]->Shares.Size(); ++k)
        {
            std::pair<String,String> const & varkeys = CellTag2VarKeys (Prps, NDim, Msh.Verts[i]->Shares[k].C->Tag);
            Nods.Last()->AddDOF (varkeys.first.CStr(), varkeys.second.CStr());
        }

#ifdef HAS_MPI
        if (PARA)
        {
            // set subset of shared nodes (between partitions)
            if (Msh.Verts[i]->PartIDs.Size()>1) InterNodes.Push (Nods.Last());
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

    // check size of Nods
    if (Nods.Size()==0) throw new Fatal("FEM::Domain::Domain: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   The array of Nodes is empty.");

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
            for (size_t k=0; k<nodes.Size(); ++k)
            {
                std::map<int,Node*>::const_iterator it = VertID2Node.find(Msh.Cells[i]->V[k]->ID);
                if (it==VertID2Node.end()) throw new Fatal("FEM::Domain::Domain:: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   Could not find Node # %d necessary for the definition of Element # %d.",Msh.Cells[i]->V[k]->ID,Msh.Cells[i]->ID);
                nodes[k] = it->second;
            }

            // allocate element
            if (Inis.HasKey(tag)) Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, Prps(tag), Inis(tag), nodes));
            else                  Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, Prps(tag), SDPair() , nodes));

            // set array of Beams
            if (prob_name=="Beam") Beams.Push (Eles.Last());

            // elements with tagged borders (edges/faces)
            if (Msh.Cells[i]->BryTags.size()>0 || prob_name=="Beam") TgdEles.Push (Eles.Last());

            // map tag to element
            Tag2Eles[tag].Push (Eles.Last());

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

    // check size of Eles
    if (Eles.Size()==0) throw new Fatal("FEM::Domain::Domain: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   The array of Elements is empty.");

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

// Methods

inline void Domain::SetBCs (Dict const & BCs)
{
    // info
#ifdef HAS_MPI
    if (PARA && MPI::COMM_WORLD.Get_rank()!=0) WithInfo = false;
#endif
    Util::Stopwatch stopwatch(/*activated*/WithInfo);
    if (WithInfo) printf("\n%s--- Domain --- setting BCs and activating/deactivating elements --------------------%s\n",TERM_CLR1,TERM_RST);

    // clear previous BCs
    for (size_t i=0; i<Eles.Size(); ++i) Eles[i]->ClrBCs ();
    pF.clear();
    pU.clear();
    InclSupport.clear();

    // map used to track whether BC was already set or not
    std::map<int,BryTagType_t> bc_set; // maps: btag => type of BC: none, element, border, node
    for (Dict_t::const_iterator it=BCs.begin(); it!=BCs.end(); ++it) bc_set[it->first] = None_t;

    // 'BCs' at the element itself (ex: 'activate', 's':source, 'qn' at beams, etc.).
    // This must come before the setting up of 'real' boundary conditions, since, for example,
    // some nodes may become activated or deactivated
    for (Dict_t::const_iterator dict_it=BCs.begin(); dict_it!=BCs.end(); ++dict_it)
    {
        // auxiliar variables
        int            btag = dict_it->first;  // the element/boundary tag
        SDPair const & bcs  = dict_it->second; // the boundary conditions

        // check if tag corresponds to element tag (and not boundary tag)
        Tag2Eles_t::const_iterator it = Tag2Eles.find(btag);
        if (it!=Tag2Eles.end()) // found elements with this btag
        {
            // flag BC set
            bc_set[btag] = Element_t;

            // multiplier
            pCalcM calcm = &Multiplier;
            if (bcs.HasKey("mfunc")) // callback specified
            {
                MDatabase_t::const_iterator q = MFuncs.find(btag);
                if (q!=MFuncs.end()) calcm = q->second;
                else throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary Tag=%d was not found in MFuncs database",btag);
            }

            // set BCs: ex: 'activate', 's', 'qn', etc.
            for (size_t i=0; i<it->second.Size(); ++i) // for each element
            {
                // check if bc key is not U and is not F (it cannot be, for example: 'ux', 'fx', etc., since these are applied to the edges/faces of the elements only)
                Element * ele = it->second[i];
                std::pair<String,String> const & varkeys = CellTag2VarKeys (Prps, NDim, ele->Cell.Tag);
                for (size_t k=0; k<bcs.Keys.Size(); ++k)
                {
                    if (Util::HasKey(varkeys.first, bcs.Keys[k])) throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' with tag==%d cannot be applied to the element (%d,%d) itself", bcs.Keys[k].CStr(), btag, ele->Cell.ID, ele->Cell.Tag);
                    if (Util::HasKey(varkeys.second,bcs.Keys[k])) throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' with tag==%d cannot be applied to the element (%d,%d) itself", bcs.Keys[k].CStr(), btag, ele->Cell.ID, ele->Cell.Tag);
                }

                // set BCs
                ele->SetBCs (/*idx_side(ignored)*/0, bcs, pF, pU, calcm);
            }
        }
    }

    // BCs at the edges/faces of elements or at nodes
    typedef std::pair<Element*,int> eleside_t;          // (element,side) pair
    typedef std::pair<int,SDPair>   tagbcs_t;           // (tag,bcs data) pair
    std::map<eleside_t,tagbcs_t>    eleside2tag_to_ubc; // maps: (ele,side) ==> U bry cond: (tag,Ubcs)
    std::map<eleside_t,tagbcs_t>    eleside2tag_to_fbc; // maps: (ele,side) ==> F bry cond: (tag,Fbcs)
    std::map<Node*,tagbcs_t>        nod2tag_to_ubc;     // map: node ==> U bry cond: (tag,Ubcs)
    std::map<Node*,tagbcs_t>        nod2tag_to_fbc;     // map: node ==> F bry cond: (tag,Fbcs)
    for (Dict_t::const_iterator dict_it=BCs.begin(); dict_it!=BCs.end(); ++dict_it)
    {
        // auxiliar variables
        int            btag  = dict_it->first;  // the element/boundary tag
        SDPair const & bcs   = dict_it->second; // the boundary conditions

        // find elements with edges equal to this btag
        for (size_t j=0; j<TgdEles.Size(); ++j) // for each element with tagged border
        {
            if (TgdEles[j]->Active) // active element
            {
                Mesh::BryTag_t           const & eftags  = TgdEles[j]->Cell.BryTags;
                std::pair<String,String> const & varkeys = CellTag2VarKeys (Prps, NDim, TgdEles[j]->Cell.Tag);
                for (Mesh::BryTag_t::const_iterator p=eftags.begin(); p!=eftags.end(); ++p)
                {
                    int idx_side = p->first;
                    int side_tag = p->second;
                    if (side_tag==btag) // found
                    {
                        // flag BC set
                        std::map<int,BryTagType_t>::const_iterator it = bc_set.find(btag);
                        if (it->second==None_t) bc_set[btag] = Border_t;
                        else if (it->second!=Border_t) throw new Fatal("FEM::Domain::SetBCs: Boundary tag==%d cannot be applied to %s and borders (edges/faces) at the same time",btag,(it->second==Element_t?"elements":"nodes"));

                        // split Ubcs from Fbcs
                        eleside_t es(TgdEles[j],idx_side); // (edge/face,side) pair
                        for (size_t k=0; k<bcs.Keys.Size(); ++k)
                        {
                            if (Util::HasKey(varkeys.first, bcs.Keys[k])) // is U key
                            {
                                eleside2tag_to_ubc[es].first = btag;
                                eleside2tag_to_ubc[es].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else // is F key or 'qn', 'qt', etc.
                            {
                                eleside2tag_to_fbc[es].first = btag;
                                eleside2tag_to_fbc[es].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                        }
                    }
                }
            }
        }

        // find nodes with tags equal to this btag
        for (size_t j=0; j<TgdNods.Size(); ++j)
        {
            if (TgdNods[j]->NShares>0) // active node
            {
                if (TgdNods[j]->Vert.Tag==btag) // found
                {
                    // flag BC set
                    std::map<int,BryTagType_t>::const_iterator it = bc_set.find(btag);
                    if (it->second==None_t) bc_set[btag] = Node_t;
                    else if (it->second!=Node_t) throw new Fatal("FEM::Domain::SetBCs: Boundary tag==%d cannot be applied to %s and nodes at the same time",btag,(it->second==Element_t?"elements":"borders (edges/faces)"));

                    // split U bcs from F bcs
                    for (size_t k=0; k<bcs.Keys.Size(); ++k)
                    {
                        if (bcs.Keys[k]=="inclsupport") // inclined support
                        {
                            nod2tag_to_ubc[TgdNods[j]].first  = btag;
                            nod2tag_to_ubc[TgdNods[j]].second = bcs;
                            break;
                        }
                        else
                        {
                            if (TgdNods[j]->UMap.HasKey(bcs.Keys[k])) // is U key
                            {
                                nod2tag_to_ubc[TgdNods[j]].first = btag;
                                nod2tag_to_ubc[TgdNods[j]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else if (TgdNods[j]->FMap.HasKey(bcs.Keys[k])) // is F key
                            {
                                nod2tag_to_fbc[TgdNods[j]].first = btag;
                                nod2tag_to_fbc[TgdNods[j]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k]));
                            }
                            else throw new Fatal("FEM::Domain::SetBCs: BC==%s with tag==%d cannot be specified to Node # %d", bcs.Keys[k].CStr(), btag, TgdNods[j]->Vert.ID);
                        }
                    }
                }
            }
        }
    }

    // check if all tags were set
    if (!PARA) // only in serial problems, since this processor may not know the tags of other elements/nodes in other processors
    {
        for (std::map<int,BryTagType_t>::const_iterator it=bc_set.begin(); it!=bc_set.end(); ++it)
        {
            if (it->second==None_t) throw new Fatal("FEM::Domain::SetBCs: Could not find any element, edge/face, or node with boundary tag == %d",it->first);
        }
    }

    // set F bcs at sides (edges/faces) of elements
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

    // add pF to pFaccum
    for (NodBCs_t::const_iterator p=pF.begin(); p!=pF.end(); ++p) // for each node with prescribed F
    {
        Node const * nod = p->first;
        for (IntDbl_t::const_iterator q=p->second.first.begin(); q!=p->second.first.end(); ++q)
        {
            size_t idof = q->first;
            pFaccum[nod].first[idof] += q->second;
        }
    }

    // check available data
    AvailableData ();

    // write header to output files
    // TODO: this should be done at each stage, since new data may become available
}

inline void Domain::PrintResults (char const * NF, bool OnlySummary, bool WithElems, double Tol) const
{
    printf("\n%s--- Results ------------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    // number format for text
    String nf(NF);
    size_t pos;
    pos=nf.find("g"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("g",pos+1); }
    pos=nf.find("f"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("f",pos+1); }
    pos=nf.find("e"); while (pos!=String::npos) { nf.replace(pos,1,"s"); pos=nf.find("e",pos+1); }

    if (!OnlySummary)
    {
        // nodes: header
        String buf;
        std::cout << TERM_CLR2 << Util::_6 << "Node";
        for (size_t i=0; i<AllUKeys.Size(); ++i) { buf.Printf(nf, AllUKeys[i].CStr());  std::cout<<buf; }
        for (size_t i=0; i<AllFKeys.Size(); ++i) { buf.Printf(nf, AllFKeys[i].CStr());  std::cout<<buf; }
        printf("%s\n",TERM_RST);

        // nodes: data
        for (size_t i=0; i<Nods.Size(); ++i)
        {
            if (Nods[i]->NShares>0)
            {
                std::cout << Util::_6 << Nods[i]->Vert.ID;
                for (size_t j=0; j<AllUKeys.Size(); ++j)
                {
                    if (Nods[i]->UMap.HasKey(AllUKeys[j]))
                    {
                        size_t idx = Nods[i]->UMap(AllUKeys[j]); // idx of DOF
                        buf.Printf(NF, Nods[i]->U[idx]);
                    }
                    else buf.Printf(nf, "   ");
                    std::cout << buf;
                }
                for (size_t j=0; j<AllFKeys.Size(); ++j)
                {
                    if (Nods[i]->FMap.HasKey(AllFKeys[j]))
                    {
                        size_t idx = Nods[i]->FMap(AllFKeys[j]); // idx of DOF
                        buf.Printf(NF, Nods[i]->F[idx]);
                    }
                    else buf.Printf(nf, "   ");
                    std::cout << buf;
                }
                printf("\n");
            }
        }
        printf("\n");
    }

    // elems: keys
    if (WithElems && !OnlySummary)
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
        String buf;
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

    // reactions
    String buf;
    Table  reac;
    SDPair sum;
    CalcReactions (reac, sum);
    printf("\n   Reactions\n%s %6s",TERM_CLR2,"Node");
    for (size_t i=1; i<reac.Keys.Size(); ++i)
    {
        String tmp("R");
        tmp.append (reac.Keys[i]);
        buf.Printf (nf, tmp.CStr());
        std::cout<<buf; 
    }
    printf("%s\n",TERM_RST);
    for (size_t i=0; i<reac.NRows; ++i)
    {
        printf(" %6d",static_cast<int>(reac("Node",i)));
        for (size_t j=1; j<reac.Keys.Size(); ++j) printf(NF,reac(reac.Keys[j],i));
        printf("\n");
    }
    printf("   -------------------------------------------------------------------\n");
    printf(" %6s","Sum");
    for (size_t j=0; j<sum.Keys.Size(); ++j)
    {
        double val = sum(sum.Keys[j]);
        printf(NF,(fabs(val)>Tol?val:0.0));
    }
    printf("\n");
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
    for (size_t i=0; i<AllUKeys.Size(); ++i)
    {
        if (!(AllUKeys[i]=="ux" || AllUKeys[i]=="uy" || AllUKeys[i]=="uz"))
        {
            oss << "        <DataArray type=\"Float32\" Name=\"" << AllUKeys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                double val = 0.0;
                if (Nods[j]->NShares>0 && Nods[j]->UMap.HasKey(AllUKeys[i]))
                {
                    val = Nods[j]->U[Nods[j]->UMap(AllUKeys[i])];
                }
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

inline bool Domain::CheckErrorNods (Table const & NodSol, SDPair const & NodTol) const
{
    // check
    if (NodSol.NRows!=Nods.Size()) throw new Fatal("Domain::CheckErrorNods: Number of rows (%d) in table NodSol with the nodes solution must be equal to the number of nodes (%d)",NodSol.NRows,Nods.Size());

    // header
    printf("\n%s--- Error Summary --- nodes --------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    std::cout << TERM_CLR2 << Util::_4<< "Key" << Util::_8s<<"Min" << Util::_8s<<"Mean" << Util::_8s<<"Max" << Util::_8s<<"Norm";
    printf("%s\n",TERM_RST);

    // results
    bool error = false;

    // nodes
    Table  reac; // reactions
    SDPair sumR; // sum of reactions
    bool reactions_calculated = false;
    for (size_t i=0; i<NodSol.Keys.Size(); ++i)
    {
        // calc error
        String key = NodSol.Keys[i];
        Array<double> err(NodSol.NRows);
        if (key[0]=='R') // reaction
        {
            if (!reactions_calculated)
            {
                CalcReactions (reac, sumR);
                reactions_calculated = true;
            }
            String ukey;
            for (size_t j=1; j<key.size(); ++j) ukey.Printf("%s%c",ukey.CStr(),key[j]);
            for (size_t j=0; j<Nods.Size(); ++j)
            {
                long pos = reac("Node").Find (Nods[j]->Vert.ID);
                if (pos<0) err[j] = 0.0;
                else       err[j] = fabs(reac(ukey,pos) - NodSol(key,j));
            }
        }
        else
        {
            for (size_t j=0; j<Nods.Size(); ++j)
            {
                if (Nods[j]->NShares>0)
                {
                    if (Nods[j]->UMap.HasKey(key)) err[j] = fabs(Nods[j]->U[Nods[j]->UMap(key)] - NodSol(key,j));
                    else                           err[j] = fabs(Nods[j]->F[Nods[j]->FMap(key)] - NodSol(key,j));
                }
                else err[j] = 0.0;
            }
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

inline bool Domain::CheckErrorEles (Table const & EleSol, SDPair const & EleTol) const
{
    // check
    if (EleSol.NRows!=Eles.Size()) throw new Fatal("Domain::CheckErrorEles: Number of rows (%d) in table EleSol with the elements solution must be equal to the number of elements (%d)",EleSol.NRows,Eles.Size());

    // header
    printf("\n%s--- Error Summary --- centre of elements -------------------------------------------%s\n",TERM_CLR1,TERM_RST);
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

// Internal methods

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

inline void Domain::CalcReactions (Table & NodesReactions, SDPair & SumReactions) const
{
    // resize tables
    String buf("Node");
    for (size_t i=0; i<AllUKeys.Size(); ++i)
    {
        buf.Printf ("%s %s",buf.CStr(),AllUKeys[i].CStr());
        SumReactions.Set (AllUKeys[i].CStr(), 0.0);
    }
    NodesReactions.SetZero (buf.CStr(), pU.size());

    // find reactions
    size_t k = 0;
    for (NodBCs_t::const_iterator p=pU.begin(); p!=pU.end(); ++p) // for each node with prescribed U
    {
        Node const *             nod  = p->first;
        NodBCs_t::const_iterator itpf = pFaccum.find(nod);
        NodesReactions.SetVal ("Node", k, nod->Vert.ID);
        for (IntDbl_t::const_iterator q=p->second.first.begin(); q!=p->second.first.end(); ++q)
        {
            size_t idof = q->first;
            double reac = nod->F[idof];
            if (itpf!=pFaccum.end()) // has prescribed F
            {
                IntDbl_t::const_iterator it = itpf->second.first.find(idof);
                if (it!=itpf->second.first.end()) reac -= it->second; // subtract prescribed F
            }
            NodesReactions.SetVal (nod->UMap.Keys[idof], k, reac);
            SumReactions[nod->UMap.Keys[idof]] += reac;
        }
        k++;
    }
}

inline void Domain::AvailableData ()
{
    // check available nodal data
    NActNods = 0;
    HasDisps = false;
    AllUKeys.Resize (0);
    AllFKeys.Resize (0);
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        if (Nods[i]->NShares>0) // active node
        {
            for (size_t j=0; j<Nods[i]->UMap.Keys.Size(); ++j)
            {
                if (Nods[i]->UMap.Keys[j]=="ux") HasDisps = true;
                AllUKeys.XPush (Nods[i]->UMap.Keys[j]);
            }
            for (size_t j=0; j<Nods[i]->FMap.Keys.Size(); ++j)
            {
                AllFKeys.XPush (Nods[i]->FMap.Keys[j]);
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
                if (keys[j]=="vx") HasVeloc = true;
                EleKeys.XPush (keys[j]);
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

}; // namespace FEM

#endif // MECHSYS_DOMAIN_H
