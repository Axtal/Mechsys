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

// Hdf5
#ifdef USE_HDF5
  #include <hdf5.h>
  #include <hdf5_hl.h>
#endif

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

#include <mechsys/fem/beam.h> ///< Needed by WriteMPY: TODO: move WriteMPY to a Python script, then remove this.

namespace FEM
{

// boundary condition function
struct BCfunc {
    String        Name;  // ex: zero, load, myfunction1, etc.
    String        Type;  // ex: cte, rmp
    Array<String> Prms;  // ex: c, m, ta, tb, tc, A, om
    Array<double> Vals;  // ex: 1, 2, 3,  4,  5,  6, 7
    Array<String> Units; // ex: m kPa -   s   s   -  -
};

// face boundary condition
struct FaceBc {
    int           Tag;   // ex: -11
    Array<String> Keys;  // ex: qn, pw, ux, uy, uz, wwx, wwy, wwz
    Array<String> Funcs; // ex: zero, load, myfunction1, etc.
};

// node boundary condition
struct NodeBc {
    int           Tag;   // ex: -100
    Array<String> Keys;  // ex: pw, ux, uy, uz, wwx, wwy, wwz
    Array<String> Funcs; // ex: zero, load, myfunction1, etc.
};

//void AddFaceBc(Array<FaceBc> & Fbcs, int tag, 

typedef std::map<Node*,Array<double> >  Res_t;      ///< Maps node to results at nodes
typedef std::map<int, Array<Element*> > Tag2Eles_t; ///< Maps tag to elements

typedef void (*pVTUfunc) (Vec3_t const & X, double Time, SDPair & KeysVals); ///< Callback to be called in WriteVTU

class Domain
{
public:
    // static
    static bool PARA;     ///< Parallel code ?
    static bool WithInfo; ///< Print information ?

    // enum
    enum BryTagType_t { None_t, Element_t, Border_t, Node_t }; ///< type of boundary condition tag

    // typedefs
    typedef std::map<String,BCFuncs*> MDatabase_t; ///< Map 'tag_varkey' to M function pointer (ex.: '-10_ux')
    typedef std::map<int,Model*>      Models_t;    ///< Map tag to model pointer

    // Constructor
    Domain (Mesh::Generic const & Msh,        ///< The mesh
            Dict const          & Prps,       ///< Element properties
            Dict const          & Mdls,       ///< Model names and parameters
            Dict const          & Inis,       ///< Initial values
            char const          * FNKey=NULL, ///< Filename key to be used during the output of results
            Array<int>    const * OutV=NULL,  ///< IDs or Tags of vertices to generate output
            pVTUfunc              pVTU=NULL); ///< Pointer to function to be called in WriteVTU (ex.: solution values)

    // Destructor
    ~Domain ();

    // Methods
    void SetBCsNew      (Array<BCfunc> const & Bfuns, Array<FaceBc> const & Fbcs, Array<NodeBc> const & Nbcs);
    void SetBCs         (Dict const & BCs);                                                                                ///< Set boundary conditions
    void NewNodsClearU  ();                                                                                                ///< Clear U values of StgNewNods (new nodes just activated by SetBCs)
    void PrintResults   (char const * NF="%15.6e", bool OnlySummary=false, bool WithElems=true, double Tol=1.0e-10) const; ///< Print results (Tol:tolerance to ignore zeros)
    void WriteMPY       (char const * FileKey, MPyPrms const & Prms)      const; // TODO: WriteHDF5 file with results and then create a post-processing script (in Python?) to generate .mpy files ///< Write MeshPython file .mpy (run with: python FileKey.mpy, needs PyLab)
    void WriteVTU       (char const * FileKey, bool DoExtrapolation=true) const; // TODO: WriteHDF5 file with results and then create a post-processing script (in Python?) to generate .vtu files ///< Write file for ParaView
    bool CheckErrorNods (Table const & NodSol, SDPair const & NodTol)     const;                                               ///< Check error at nodes
    bool CheckErrorEles (Table const & EleSol, SDPair const & EleTol)     const;                                               ///< Check error at the centre of elements
    bool CheckErrorIPs  (Table const & EleSol, SDPair const & EleTol)     const;                                               ///< Check error at integration points of elements
    void PrintBCs       (std::ostream & os, double tf=1.0)                const; ///< Print boundary conditions
    void SaveState      (char const * FileKey)                            const; ///< Save Nodes and Elements' state to an HDF5 file
    void LoadState      (char const * FileKey);                                  ///< Load Nodes and Elements' state from an HDF5 file
    void SetGeostatic   (SDPair const & Data);                                   ///< Set geostatic state

    // Internal methods
    void NodalResults  (bool OnlyOutNods=false) const;                        ///< Calculate extrapolated values at elments IPs' to nodes
    void OutResults    (char const * VTUFName=NULL, bool HeaderOnly=false);   ///< Do output results
    void CalcReactions (Table & NodesReactions, SDPair & SumReactions) const; ///< Calculate reactions

    // Data
    double                Time;        ///< Current time (t)
    size_t                IdxOut;      ///< Index for output of Time files (VTU)
    Dict          const & Prps;        ///< Element properties
    Dict          const & Inis;        ///< Initial values
    int                   NDim;        ///< Space dimension
    Models_t              Mdls;        ///< Models
    Models_t              XMdls;       ///< Extra Models
    pVTUfunc              VTUfunc;     ///< Pointer to function to be called in WriteVTU
    std::map<int,Node*>   VertID2Node; ///< Map vertex ID to Node
    Array<Node*>          Nods;        ///< (Allocated memory) Nodes
    Array<Element*>       Eles;        ///< (Allocated memory) Elements
    Array<Node*>          TgdNods;     ///< Tagged Nodes (at boundaries)
    Array<Element*>       TgdEles;     ///< Tagged edges or faces of Elements (at boundaries)
    Array<Node*>          OutNods;     ///< Nodes for which output (results) is generated
    Array<Element*>       OutEles;     ///< Elements corresponding to Nodes for which output is to be generated (needed during extrapolation)
    Array<std::ofstream*> FilNods;     ///< Files with results at selected nodes (OutNods)
    Array<Element*>       Beams;       ///< Subset of elements of type Beam
    MDatabase_t           MFuncs;      ///< Database of pointers to M functions
    Array<Node*>          InterNodes;  ///< Nodes on the inferface between partitions (if PARA==true)
    Tag2Eles_t            Tag2Eles;    ///< Map tag to elements. Useful to activate/deactivate layers
    Array<Node*>          ActNods;     ///< Subset of active nodes (set once per call to SetBCs)
    Array<Element*>       ActEles;     ///< Subset of active elements (set once per call to SetBCs)
    Array<Node*>          NodsWithPF;  ///< Subset of nodes with prescribed F
    Array<Node*>          NodsWithPU;  ///< Subset of nodes with prescribed U
    Array<Element*>       ElesWithBC;  ///< Subset of elements with prescribed boundary conditions (ex.: beams, flowelem with convection, etc.)
    Array<Node*>          NodsPins;    ///< Subset of nodes that are pins
    Array<Node*>          NodsIncSup;  ///< Subset of nodes with inclined supports
    Array<Node*>          StgNewNods;  ///< Nodes just activated in the new stage set by SetBCs

    // Nodal results
    bool          HasDisps;        ///< Has displacements (ux, ...) ?
    bool          HasVeloc;        ///< Has velocities (vx, ...) ?
    bool          HasWDisch;       ///< Has specific water discharge (qwx, ...) ?
    Array<String> AllUKeys;        ///< All U node keys (ux, uy, hh, ...)
    Array<String> AllFKeys;        ///< All F node keys (fx, fy, qq, ...)
    Array<String> AllUKeysBCF;     ///< All U BC callback functions (UX, UY, HH, ...) : upper case
    Array<String> AllFKeysBCF;     ///< All F BC callback functions (FX, FY, QQ, ...) : upper case
    Array<String> AllExtraKeys;    ///< All extra BC keys such as 'qn', 'qt', 's', 'activate', 'grav'
    Array<String> AllExtraKeysBCF; ///< All extra BC keys such as 'qn', 'qt', 's', 'activate', 'grav'
    Array<String> AllEKeys;        ///< All element keys (sx, sy, vx, vy, ...)
    mutable Res_t NodResults;      ///< Extrapolated nodal results
    mutable Res_t NodResCount;     ///< Count how many times a variable (key) was added to a node
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

    buf.Printf("\n%s--- E(x)tra Models -----------------------------------------------------------%s\n",TERM_CLR3,TERM_RST);
    os << buf;
    for (Domain::Models_t::const_iterator p=D.XMdls.begin(); p!=D.XMdls.end(); ++p)
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
    return os;
}

inline void Domain::PrintBCs (std::ostream & os, double tf) const
{
    Array<double> T(2);
    T = 0.0, tf;
    for (size_t k=0; k<T.Size();    ++k)
    {
        String buf;
        os    << "\n  -----------------------------------------------------------------------------\n";
        buf.Printf("  Boundary conditions at t = %g\n",T[k]);
        os << buf;
        for (size_t i=0; i<Nods.Size(); ++i)
        {
            Node const & nod = (*Nods[i]);
            if (nod.NPU()>0 || nod.NPF()>0 || nod.HasIncSup())
            {
                os << Util::_4 << nod.Vert.ID << " ";
                if (nod.NPU()>0)
                {
                    os << TERM_CLR4 << "PU:{" << TERM_RST;
                    for (size_t i=0; i<nod.NPU(); ++i)
                    {
                        os << nod.PUKey(i) << "=" << nod.PU(i,T[k]);
                        if (i!=nod.NPU()-1) os << " ";
                    }
                    os << TERM_CLR4 << "} " << TERM_RST;
                }
                if (nod.NPF()>0)
                {
                    os << TERM_CLR5 << "PF:{" << TERM_RST;
                    for (size_t i=0; i<nod.NPF(); ++i)
                    {
                        os << nod.PFKey(i) << "=" << nod.PF(i,T[k]);
                        if (i!=nod.NPF()-1) os << " ";
                    }
                    os << TERM_CLR5 << "} " << TERM_RST;
                }
                if (nod.HasIncSup()) os << TERM_CLR1 << "INCSUP" << TERM_RST;
                os << std::endl;
            }
        }
    }
}

// Constructor and destructor

inline Domain::Domain (Mesh::Generic const & Msh, Dict const & ThePrps, Dict const & TheMdls, Dict const & TheInis, char const * FNKey, Array<int> const * OutV, pVTUfunc pVTU)
    : Time(0.0), IdxOut(0), Prps(ThePrps), Inis(TheInis), NDim(Msh.NDim), VTUfunc(pVTU), HasDisps(false), HasVeloc(false), HasWDisch(false)
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
            Mdls[tag] = AllocModel (model_name, NDim, TheMdls(tag), NULL);
        }
        else throw new Fatal("Domain::Domain: Dictionary of models must have keyword 'name' defining the name of the model");
        if (TheMdls(tag).HasKey("xname"))
        {
            String model_name;
            MODEL.Val2Key (TheMdls(tag)("xname"), model_name);
            XMdls[tag] = AllocModel (model_name, NDim, TheMdls(tag), Mdls[tag]);
        }
    }

    // set nodes from mesh
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

        // add DOFs to node and set BCs keys, BCs callbacks keys, extra keys, and extra keys callbacks' keys
        for (size_t k=0; k<Msh.Verts[i]->Shares.Size(); ++k)
        {
            // find variable keys
            int cell_tag = Msh.Verts[i]->Shares[k].C->Tag;
            String prob, prob_nd;
            PROB.Val2Key   (Prps(cell_tag)("prob"), prob);
            prob_nd.Printf ("%s%dD", prob.CStr(), NDim);
            ElementVarKeys_t  ::const_iterator it = ElementVarKeys  .find (prob_nd);
            ElementExtraKeys_t::const_iterator jt = ElementExtraKeys.find (prob_nd);
            if (it==ElementVarKeys  .end()) throw new Fatal("Domain::Domain: Could not find %s in ElementVarKeys map",prob_nd.CStr());
            if (jt==ElementExtraKeys.end()) throw new Fatal("Domain::Domain: Could not find %s in ElementExtraKeys map",prob_nd.CStr());

            // add DOF to node
            Nods.Last()->AddDOF (it->second.first.CStr(), it->second.second.CStr());

            // collect U keys
            Array<String> ukeys;
            Util::Keys2Array (it->second.first,  ukeys);
            for (size_t m=0; m<ukeys.Size(); ++m)
            {
                String upukey(ukeys[m]);
                upukey     .ToUpper ();
                AllUKeys   .XPush   (ukeys[m]);
                AllUKeysBCF.XPush   (upukey);
                if (ukeys[m]=="ux") HasDisps = true;
            }

            // collect F keys
            Array<String> fkeys;
            Util::Keys2Array (it->second.second, fkeys);
            for (size_t m=0; m<fkeys.Size(); ++m)
            {
                String upfkey(fkeys[m]);
                upfkey     .ToUpper ();
                AllFKeys   .XPush   (fkeys[m]);
                AllFKeysBCF.XPush   (upfkey);
            }

            // collect extra keys
            for (size_t m=0; m<jt->second.Size(); ++m)
            {
                String upextkey(jt->second[m]);
                upextkey       .ToUpper ();
                AllExtraKeys   .XPush   (jt->second[m]);
                AllExtraKeysBCF.XPush   (upextkey);
            }
        }

        /*
        std::cout << "AllUKeys        = " << AllUKeys        << std::endl;
        std::cout << "AllFKeys        = " << AllFKeys        << std::endl;
        std::cout << "AllUKeysBCF     = " << AllUKeysBCF     << std::endl;
        std::cout << "AllFKeysBCF     = " << AllFKeysBCF     << std::endl;
        std::cout << "AllExtraKeys    = " << AllExtraKeys    << std::endl;
        std::cout << "AllExtraKeysBCF = " << AllExtraKeysBCF << std::endl;
        */

#ifdef HAS_MPI
        if (PARA)
        {
            // set subset of shared nodes (between partitions)
            if (Msh.Verts[i]->PartIDs.Size()>1) InterNodes.Push (Nods.Last());
        }
#endif
        // set list of nodes for output
        if (FNKey!=NULL && OutV!=NULL)
        {
            bool has_tag = (Msh.Verts[i]->Tag<0 ? OutV->Has(Msh.Verts[i]->Tag) : false);
            if (OutV->Has(Msh.Verts[i]->ID) || has_tag)
            {
                String buf; buf.Printf("%s_nod_%d.res", FNKey, Msh.Verts[i]->ID);
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

            // extra model
            Model const * xmdl = NULL;
            Models_t::const_iterator xm = XMdls.find(tag);
            if (xm!=XMdls.end()) xmdl = xm->second;

            // connectivity
            bool has_outnode = false; // does this new element have any node belonging to the output array ?
            Array<Node*> nodes(Msh.Cells[i]->V.Size());
            for (size_t k=0; k<nodes.Size(); ++k)
            {
                std::map<int,Node*>::const_iterator it = VertID2Node.find(Msh.Cells[i]->V[k]->ID);
                if (it==VertID2Node.end()) throw new Fatal("FEM::Domain::Domain:: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   Could not find Node # %d necessary for the definition of Element # %d.",Msh.Cells[i]->V[k]->ID,Msh.Cells[i]->ID);
                nodes[k] = it->second;

                // check if node of element belongs to output array
                if (OutNods.Size()>0) { if (OutNods.Has(nodes[k])) has_outnode = true; }
            }

            // allocate element
            if (Inis.HasKey(tag)) Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, xmdl, Prps(tag), Inis(tag), nodes));
            else                  Eles.Push (AllocElement(prob_name, NDim, (*Msh.Cells[i]), mdl, xmdl, Prps(tag), SDPair() , nodes));

            // set subset of active elements
            if (Eles.Last()->Active) ActEles.Push (Eles.Last());

            // set array of Beams
            if (prob_name=="Beam") Beams.Push (Eles.Last());

            // elements with tagged borders (edges/faces)
            if (Msh.Cells[i]->BryTags.size()>0 || prob_name=="Beam") TgdEles.Push (Eles.Last());

            // map tag to element
            Tag2Eles[tag].Push (Eles.Last());

            // element belongs to output array
            if (has_outnode) OutEles.Push (Eles.Last());

            // collect all element keys
            Array<String> keys;
            Eles.Last()->StateKeys (keys);
            for (size_t j=0; j<keys.Size(); ++j)
            {
                if (keys[j]=="vx")  HasVeloc  = true;
                if (keys[j]=="qwx") HasWDisch = true;
                AllEKeys.XPush (keys[j]); // push only if item is not in array yet
            }
        }
        else throw new Fatal("Domain::SetMesh: Dictionary of properties must have keyword 'prob' defining the type of element corresponding to a specific problem");
    }

    // check size of Eles
    if (Eles.Size()==0) throw new Fatal("FEM::Domain::Domain: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   The array of Elements is empty.");

    // set pin nodes
    if (Msh.Pins.size()>0)
    {
        for (Mesh::Pin_t::const_iterator p=Msh.Pins.begin(); p!=Msh.Pins.end(); ++p)
        {
            std::map<int,Node*>::const_iterator it0 = VertID2Node.find (p->first->ID);
            if (it0==VertID2Node.end()) throw new Fatal("FEM::Domain::Domain:: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   Could not find Node # %d corresponding to a Pin.",p->first->ID);
#ifdef HAS_MPI
            if (PARA)
            {
                // skip pins that aren't in this partition
                if (p->first->PartIDs.Find(MPI::COMM_WORLD.Get_rank())<0) continue;
            }
#endif
            Array<Node*> con_nods; // connected nodes to this pin
            for (size_t i=0; i<p->second.Size(); ++i)
            {
                std::map<int,Node*>::const_iterator it1 = VertID2Node.find (p->second[i]->ID);
                if (it1==VertID2Node.end()) throw new Fatal("FEM::Domain::Domain:: There is a problem with the Msh structure (probably the domain decomposition has failed).\n   Could not find Node # %d connected to Pin # %d. Both pins must be in the same partition.",p->second[i]->ID,p->first->ID);
                con_nods.Push (it1->second);
            }
            it0->second->SetPin (con_nods);
            NodsPins.Push (it0->second);
        }
    }

    // set subset of active nodes (must be after the creation of elements)
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        if (Nods[i]->NShares>0) ActNods.Push (Nods[i]);
    }

    // print header of output files
    OutResults (/*vtufname*/NULL, /*headeronly*/true);

    // write file with IDs and tags of output nodes
    if (FNKey!=NULL && OutV!=NULL)
    {
        String buf(FNKey);  buf.append("_out_nods.res");
        std::ofstream of(buf.CStr(), std::ios::out);
        buf.Printf("%6s %6s %16s %16s %16s\n","id","tag","x","y","z");
        of << buf;
        for (size_t i=0; i<OutNods.Size(); ++i)
        {
            buf.Printf("%6d %6d %16.8e %16.8e %16.8e\n",OutNods[i]->Vert.ID, OutNods[i]->Vert.Tag, OutNods[i]->Vert.C(0), OutNods[i]->Vert.C(1), OutNods[i]->Vert.C(2));
            of << buf;
        }
        of.close();
    }
}

inline Domain::~Domain()
{
    for (Domain::Models_t::iterator p=Mdls.begin(); p!=Mdls.end(); ++p) delete p->second;
    for (Domain::Models_t::iterator p=XMdls.begin(); p!=XMdls.end(); ++p) delete p->second;
    for (size_t i=0; i<Nods   .Size(); ++i) if (Nods   [i]!=NULL) delete Nods   [i];
    for (size_t i=0; i<Eles   .Size(); ++i) if (Eles   [i]!=NULL) delete Eles   [i];
    for (size_t i=0; i<FilNods.Size(); ++i) if (FilNods[i]!=NULL) { FilNods[i]->close(); delete FilNods[i]; }
}

// Methods

inline void Domain::SetBCsNew (Array<BCfunc> const & Bfuns, Array<FaceBc> const & Fbcs, Array<NodeBc> const & Nbcs)
{
}

inline void Domain::SetBCs (Dict const & BCs)
{
    // info
#ifdef HAS_MPI
    if (PARA && MPI::COMM_WORLD.Get_rank()!=0) WithInfo = false;
#endif
    Util::Stopwatch stopwatch(/*activated*/WithInfo);
    if (WithInfo) printf("\n%s--- Domain --- setting BCs and activating/deactivating elements --------------------%s\n",TERM_CLR1,TERM_RST);

    // clear previous BCs
    for (size_t i=0; i<NodsWithPU.Size(); ++i) NodsWithPU[i]->DelPUs    ();
    for (size_t i=0; i<NodsWithPF.Size(); ++i) NodsWithPF[i]->DelPFs    ();
    for (size_t i=0; i<NodsIncSup.Size(); ++i) NodsIncSup[i]->DelIncSup ();
    for (size_t i=0; i<ElesWithBC.Size(); ++i) ElesWithBC[i]->ClrBCs    ();

    // map used to track whether BC was already set or not
    std::map<int,BryTagType_t> bc_set; // maps: btag => type of BC: none, element, border, node
    for (Dict_t::const_iterator it=BCs.begin(); it!=BCs.end(); ++it) bc_set[it->first] = None_t;

    // 'BCs' at the element itself (ex: 'activate', 's':source, 'qn' at beams, etc.).
    // This must come before the setting up of 'real' boundary conditions, since, for example,
    // some nodes may become activated or deactivated
    ElesWithBC.Resize (0);
    for (Dict_t::const_iterator dict_it=BCs.begin(); dict_it!=BCs.end(); ++dict_it)
    {
        // auxiliar variables
        int            btag = dict_it->first;  // the element/boundary tag
        SDPair const & bcs  = dict_it->second; // the boundary conditions

        // check if tag corresponds to element tag (and not boundary tag)
        Tag2Eles_t::const_iterator it = Tag2Eles.find(btag);
        if (it!=Tag2Eles.end()) // found elements with this btag
        {
            bc_set[btag] = Element_t; // flag BC set
            for (size_t k=0; k<bcs.Keys.Size(); ++k)
            {
                // callback function
                BCFuncs * bcfun = NULL;
                if (AllExtraKeysBCF.Has(bcs.Keys[k])) // callback specified
                {
                    String mfkey;
                    mfkey.Printf ("%d_%s", btag, bcs.Keys[k].ToLowerCpy().CStr());
                    MDatabase_t::const_iterator q = MFuncs.find(mfkey);
                    if (q==MFuncs.end()) throw new Fatal("FEM::Domain::SetBCs: Multiplier function with Element Tag=%d was not found in MFuncs database ('%s' key)",btag,mfkey.CStr());
                    else bcfun = q->second;
                }

                // check if bc key is not U or F key
                if (AllUKeys.Has(bcs.Keys[k])) throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' with tag==%d cannot be applied to the Element", bcs.Keys[k].CStr(), btag);
                if (AllFKeys.Has(bcs.Keys[k])) throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' with tag==%d cannot be applied to the Element", bcs.Keys[k].CStr(), btag);


                printf("\n\n");
                for (size_t ib=0; ib<bcs.Keys.Size(); ++ib)
                {
                    printf("[1;35m%s  =>  %g[0m\n", bcs.Keys[ib].CStr(), bcs(bcs.Keys[ib]));
                }
                printf("\n\n");
















                // TODO: 'deactivate' and 'gravity' must be together
                //       currently this code is splitting them => major rewrite needed

















                // set BCs: ex: 'activate', 's', 'qn', etc.
                SDPair bc;
                bc.Set (bcs.Keys[k].ToLowerCpy().CStr(), bcs(bcs.Keys[k]));

                // TODO: extra arguments need to be passed onto some elements => fix this
                //if (bcs.HasKey("activate"))
                //{
                    //bc.Set("activate", bcs("activate"));
                    //bc.Set("gravity",  bcs("gravity"));
                //}
                //if (bcs.HasKey("deactivate")) // environment temperature
                //{
                    //bc.Set("deactivate", bcs("deactivate"));
                    //bc.Set("gravity",  bcs("gravity"));
                //}

                for (size_t i=0; i<it->second.Size(); ++i) // for each element
                {
                    Element * ele = it->second[i];
                    ele->SetBCs (/*idx_side(ignored)*/0, bc, bcfun);
                    ElesWithBC.Push (ele);
                }
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
                Mesh::BryTag_t const & eftags = TgdEles[j]->Cell.BryTags;
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
                            if (bcs.Keys[k]=="incsup") // is inclined support
                            {
                                eleside2tag_to_ubc[es].first  = btag;
                                eleside2tag_to_ubc[es].second = bcs; // has to be copied because of the extra parameters such as 'alpha' and so on
                                break;
                            }
                            else
                            {
                                if (AllUKeys.Has(bcs.Keys[k]) || AllUKeysBCF.Has(bcs.Keys[k])) // is U key or U callback function
                                {
                                    eleside2tag_to_ubc[es].first = btag;
                                    eleside2tag_to_ubc[es].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k])); // cannot copy bcs here because we don't want F values
                                }
                                else // is F key or 'qn', 'qt', etc.
                                {
                                    eleside2tag_to_fbc[es].first = btag;
                                    eleside2tag_to_fbc[es].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k])); // cannot copy bcs here because we don't want U values
                                }
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
                        if (bcs.Keys[k]=="incsup") // inclined support
                        {
                            nod2tag_to_ubc[TgdNods[j]].first  = btag;
                            nod2tag_to_ubc[TgdNods[j]].second = bcs; // has to be copied because of the extra parameters such as 'alpha' and so on
                            break;
                        }
                        else
                        {
                            if (AllUKeys.Has(bcs.Keys[k]) || AllUKeysBCF.Has(bcs.Keys[k])) // is U key or U callback function
                            {
                                nod2tag_to_ubc[TgdNods[j]].first = btag;
                                nod2tag_to_ubc[TgdNods[j]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k])); // cannot copy bcs here because we don't want F values
                            }
                            else if (AllFKeys.Has(bcs.Keys[k]) || AllFKeysBCF.Has(bcs.Keys[k])) // is F key or F callback function
                            {
                                nod2tag_to_fbc[TgdNods[j]].first = btag;
                                nod2tag_to_fbc[TgdNods[j]].second.Set (bcs.Keys[k].CStr(), bcs(bcs.Keys[k])); // cannot copy bcs here because we don't want F values
                            }
                            else throw new Fatal("FEM::Domain::SetBCs: Boundary condition '%s' (tag==%d) cannot be specified to Nodes (node # %d)", bcs.Keys[k].CStr(), btag, TgdNods[j]->Vert.ID);
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
        for (size_t k=0; k<bcs.Keys.Size(); ++k)
        {
            String key = bcs.Keys[k].ToLowerCpy();
            SDPair bc;
            bc.Set (key.CStr(), bcs(bcs.Keys[k]));
            
            // TODO: extra arguments need to be passed onto some elements => fix this
            if (bcs.HasKey("h")) // convection transmissivity coefficient
            {
                bc.Set("h", bcs("h"));
            }
            if (bcs.HasKey("Tinf")) // environment temperature
            {
                bc.Set("Tinf", bcs("Tinf"));
            }

            if (AllFKeysBCF.Has(bcs.Keys[k]) || AllExtraKeysBCF.Has(bcs.Keys[k])) // callback specified
            {
                String mfkey;
                mfkey.Printf ("%d_%s", bc_tag, key.CStr());
                MDatabase_t::const_iterator q = MFuncs.find(mfkey);
                if (q==MFuncs.end()) throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary (edge/face) Tag=%d was not found in MFuncs database ('%s' key)",bc_tag,mfkey.CStr());
                ele->SetBCs (idx_side, bc, q->second);
            }
            else ele->SetBCs (idx_side, bc, NULL);
        }
    }

    // set F bcs at nodes
    for (std::map<Node*,tagbcs_t>::iterator p=nod2tag_to_fbc.begin(); p!=nod2tag_to_fbc.end(); ++p)
    {
        Node         * nod    = p->first;
        int            bc_tag = p->second.first;
        SDPair const & bcs    = p->second.second;
        for (size_t k=0; k<bcs.Keys.Size(); ++k)
        {
            if (AllFKeysBCF.Has(bcs.Keys[k])) // callback specified
            {
                String key = bcs.Keys[k].ToLowerCpy();
                String mfkey;
                mfkey.Printf ("%d_%s", bc_tag, key.CStr());
                MDatabase_t::const_iterator q = MFuncs.find(mfkey);
                if (q==MFuncs.end()) throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary (nodes) Tag=%d was not found in MFuncs database ('%s' key)",bc_tag,mfkey.CStr());
                nod->AddToPF (key, bcs(bcs.Keys[k]), q->second);
            }
            else nod->AddToPF (bcs.Keys[k], bcs(bcs.Keys[k]), NULL);
        }
    }

    // set U bcs at sides (edges/faces) of elements
    for (std::map<eleside_t,tagbcs_t>::iterator p=eleside2tag_to_ubc.begin(); p!=eleside2tag_to_ubc.end(); ++p)
    {
        Element      * ele      = p->first.first;
        int            idx_side = p->first.second;
        int            bc_tag   = p->second.first;
        SDPair const & bcs      = p->second.second;
        for (size_t k=0; k<bcs.Keys.Size(); ++k)
        {
            String key = bcs.Keys[k].ToLowerCpy();
            SDPair bc;
            bc.Set (key.CStr(), bcs(bcs.Keys[k]));
            if (AllUKeysBCF.Has(bcs.Keys[k])) // callback specified
            {
                String mfkey;
                mfkey.Printf ("%d_%s", bc_tag, key.CStr());
                MDatabase_t::const_iterator q = MFuncs.find(mfkey);
                if (q==MFuncs.end()) throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary (edge/face) Tag=%d was not found in MFuncs database ('%s' key)",bc_tag,mfkey.CStr());
                ele->SetBCs (idx_side, bc, q->second);
            }
            else ele->SetBCs (idx_side, bc, NULL);
        }
    }

    // set U bcs at nodes
    for (std::map<Node*,tagbcs_t>::iterator p=nod2tag_to_ubc.begin(); p!=nod2tag_to_ubc.end(); ++p)
    {
        Node         * nod    = p->first;
        int            bc_tag = p->second.first;
        SDPair const & bcs    = p->second.second;
        for (size_t k=0; k<bcs.Keys.Size(); ++k)
        {
            if (bcs.HasKey("incsup"))
            {
                if (NDim!=2) throw new Fatal("FEM::Domain::SetBCs: Inclined support is not implemented for 3D problems yet");
                nod->SetIncSup (bcs("alpha"));
                break;
            }
            else
            {
                if (AllUKeysBCF.Has(bcs.Keys[k])) // callback specified
                {
                    String key = bcs.Keys[k].ToLowerCpy();
                    String mfkey;
                    mfkey.Printf ("%d_%s", bc_tag, key.CStr());
                    MDatabase_t::const_iterator q = MFuncs.find(mfkey);
                    if (q==MFuncs.end()) throw new Fatal("FEM::Domain::SetBCs: Multiplier function with boundary (nodes) Tag=%d was not found in MFuncs database ('%s' key)",bc_tag,mfkey.CStr());
                    nod->SetPU (key, bcs(bcs.Keys[k]), q->second);
                }
                else nod->SetPU (bcs.Keys[k], bcs(bcs.Keys[k]), NULL);
            }
        }
    }

    // set subsets of active nodes, nodes with prescribed U and/or F, and nodes with inclined supports
    Array<Node*> prev_act_nods(ActNods); // previous active nodes
    StgNewNods.Resize (0);               // new nodes activated here
    ActNods   .Resize (0);
    NodsWithPU.Resize (0);
    NodsWithPF.Resize (0);
    NodsIncSup.Resize (0);
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        if (Nods[i]->NShares>0)
        {
            ActNods.Push (Nods[i]); 
            if (!prev_act_nods.Has(Nods[i])) StgNewNods.Push (Nods[i]); // node was just activated
        }
        if (Nods[i]->NPU()>0)     NodsWithPU.Push (Nods[i]);
        if (Nods[i]->HasIncSup()) NodsIncSup.Push (Nods[i]);
        if (Nods[i]->NPF()>0)
        {
            NodsWithPF.Push  (Nods[i]);
            Nods[i]->AccumPF (); // accumulate PF inside Node (to calculate Reactions later)
        }
    }

    //std::cout << "prev_act_nods = "; for (size_t i=0; i<prev_act_nods.Size(); ++i) std::cout << prev_act_nods[i]->Vert.ID << " "; std::cout << std::endl;
    //std::cout << "StgNewNods    = "; for (size_t i=0; i<StgNewNods   .Size(); ++i) std::cout << StgNewNods   [i]->Vert.ID << " "; std::cout << std::endl;

    // set subsets of active elements
    ActEles.Resize (0);
    for (size_t i=0; i<Eles.Size(); ++i) { if (Eles[i]->Active) ActEles.Push (Eles[i]); }
}

inline void Domain::NewNodsClearU ()
{
    for (size_t i=0; i<StgNewNods.Size(); ++i) StgNewNods[i]->ClearU ();
}

inline void Domain::PrintResults (char const * NF, bool OnlySummary, bool WithElems, double Tol) const
{
    printf("\n%s--- Results ------------------------------------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    // number format for text
    String fmt;
    fmt.TextFmt (NF);

    if (!OnlySummary)
    {
        // nodes: header
        String buf;
        std::cout << TERM_CLR2 << Util::_6 << "Node";
        for (size_t i=0; i<AllUKeys.Size(); ++i) { buf.Printf(fmt, AllUKeys[i].CStr());  std::cout<<buf; }
        for (size_t i=0; i<AllFKeys.Size(); ++i) { buf.Printf(fmt, AllFKeys[i].CStr());  std::cout<<buf; }
        printf("%s\n",TERM_RST);

        // nodes: data
        for (size_t i=0; i<ActNods.Size(); ++i)
        {
            std::cout << Util::_6 << ActNods[i]->Vert.ID;
            for (size_t j=0; j<AllUKeys.Size(); ++j)
            {
                buf.Printf(NF, ActNods[i]->UOrZero(AllUKeys[j]));
                std::cout << buf;
            }
            for (size_t j=0; j<AllFKeys.Size(); ++j)
            {
                buf.Printf(NF, ActNods[i]->FOrZero(AllFKeys[j]));
                std::cout << buf;
            }
            printf("\n");
        }
        printf("\n");
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
        buf.Printf (fmt, tmp.CStr());
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

inline void Domain::WriteMPY (char const * FNKey, FEM::MPyPrms const & Prms) const
{
    // pointer to array of elements
    Array<Element*> const * eles = (Prms.OnlyBeams ? &Beams : &ActEles);

    // find min and max M and N
    if (Prms.FindMLimits && Beams.Size()>0)
    {
        FEM::Beam const * e0 = static_cast<FEM::Beam const *>(Beams[0]);
        double N,V,M;
        e0->CalcRes (0.0, N,V,M);
        Prms.EleMmin = e0;
        Prms.EleMmax = e0;
        Prms.EleNmin = e0;
        Prms.EleNmax = e0;
        Prms.EleVmin = e0;
        Prms.EleVmax = e0;
        Prms.Mmin    = M;
        Prms.Mmax    = M;
        Prms.Nmin    = N;
        Prms.Nmax    = N;
        Prms.Vmin    = V;
        Prms.Vmax    = V;
        Prms.rMmin   = 0.0;
        Prms.rMmax   = 0.0;
        for (size_t i=0; i<Beams.Size(); ++i)
        {
            FEM::Beam const * e = static_cast<FEM::Beam const *>(Beams[i]);
            for (size_t j=0; j<=Prms.NDiv; ++j)
            {
                double r = static_cast<double>(j)/static_cast<double>(Prms.NDiv);
                e->CalcRes (r, N,V,M);
                if (M<Prms.Mmin) { Prms.Mmin=M;  Prms.EleMmin=e;  Prms.rMmin=r; }
                if (M>Prms.Mmax) { Prms.Mmax=M;  Prms.EleMmax=e;  Prms.rMmax=r; }
                if (N<Prms.Nmin) { Prms.Nmin=N;  Prms.EleNmin=e; }
                if (N>Prms.Nmax) { Prms.Nmax=N;  Prms.EleNmax=e; }
                if (V<Prms.Vmin) { Prms.Vmin=V;  Prms.EleVmin=e; }
                if (V>Prms.Vmax) { Prms.Vmax=V;  Prms.EleVmax=e; }
            }
        }
    }

    // calculate scale factor for diagrams
    if (Prms.AutoLimits)
    {
        // bounding box
        Vec3_t min, max, del;
        min = (*eles)[0]->Con[0]->Vert.C(0), (*eles)[0]->Con[0]->Vert.C(1), (*eles)[0]->Con[0]->Vert.C(2);
        max = (*eles)[0]->Con[1]->Vert.C(0), (*eles)[0]->Con[1]->Vert.C(1), (*eles)[0]->Con[1]->Vert.C(2);
        for (size_t i=1; i<eles->Size(); ++i)
        {
            for (size_t j=0; j<(*eles)[i]->Con.Size(); ++j)
            {
                if ((*eles)[i]->Con[j]->Vert.C(0)<min(0)) min(0) = (*eles)[i]->Con[j]->Vert.C(0);
                if ((*eles)[i]->Con[j]->Vert.C(1)<min(1)) min(1) = (*eles)[i]->Con[j]->Vert.C(1);
                if ((*eles)[i]->Con[j]->Vert.C(2)<min(2)) min(2) = (*eles)[i]->Con[j]->Vert.C(2);
                if ((*eles)[i]->Con[j]->Vert.C(0)>max(0)) max(0) = (*eles)[i]->Con[j]->Vert.C(0);
                if ((*eles)[i]->Con[j]->Vert.C(1)>max(1)) max(1) = (*eles)[i]->Con[j]->Vert.C(1);
                if ((*eles)[i]->Con[j]->Vert.C(2)>max(2)) max(2) = (*eles)[i]->Con[j]->Vert.C(2);
            }
        }

        // max distance in mesh
        del = max - min;
        Prms.MaxDist = Norm(del);

        // scale factor
        double max_abs = 0;
        double N,V,M;
        if (Prms.DrawN)
        {
            if (Prms.EleNmin!=NULL)
            {
                static_cast<FEM::Beam const *>(Prms.EleNmin)->CalcRes (0.0, N,V,M);
                if (fabs(N)>max_abs) max_abs = fabs(N);
            }
            if (Prms.EleNmax!=NULL)
            {
                static_cast<FEM::Beam const *>(Prms.EleNmax)->CalcRes (0.0, N,V,M);
                if (fabs(N)>max_abs) max_abs = fabs(N);
            }
        }
        else if (Prms.DrawV)
        {
            if (Prms.EleVmin!=NULL)
            {
                static_cast<FEM::Beam const *>(Prms.EleVmin)->CalcRes (0.0, N,V,M);
                if (fabs(V)>max_abs) max_abs = fabs(V);
            }
            if (Prms.EleVmax!=NULL)
            {
                static_cast<FEM::Beam const *>(Prms.EleVmax)->CalcRes (0.0, N,V,M);
                if (fabs(V)>max_abs) max_abs = fabs(V);
            }
        }
        else
        {
            if (Prms.EleMmin!=NULL)
            {
                static_cast<FEM::Beam const *>(Prms.EleMmin)->CalcRes (Prms.rMmin, N,V,M);
                if (fabs(M)>max_abs) max_abs = fabs(M);
            }
            if (Prms.EleMmax!=NULL)
            {
                static_cast<FEM::Beam const *>(Prms.EleMmax)->CalcRes (Prms.rMmax, N,V,M);
                if (fabs(M)>max_abs) max_abs = fabs(M);
            }
        }
        Prms.SF = (Prms.PctMaxDist * Prms.MaxDist) / (max_abs>0.0 ? max_abs : 1.0);
    }

    // elements and diagrams
    std::ostringstream oss;
    MPL::Header (oss);
    for (size_t i=0; i<eles->Size(); ++i) (*eles)[i]->Draw (oss, Prms);
    MPL::AddPatch (oss);

    // reactions
    if (Prms.WithReac)
    {
        Table  reac;
        SDPair sum;
        CalcReactions (reac, sum);
        for (size_t i=0; i<reac.NRows; ++i)
        {
            int id = static_cast<int>(reac("Node",i));
            if (Prms.ReacNodes.Has(id))
            {
                std::map<int,Node*>::const_iterator it = VertID2Node.find(id);
                if (it!=VertID2Node.end())
                {
                    FEM::Node const * nod = const_cast<FEM::Node const *>(it->second);
                    String buf;
                    buf.Printf ("text(%g,%g, 'Rx=%g, Ry=%g', fontsize=9, color='blue', ha='center', va='top', backgroundcolor='white')\n",
                                nod->Vert.C(0), nod->Vert.C(1), (reac.Keys.Has("ux") ? reac("ux",i) : 0),
                                                                (reac.Keys.Has("uy") ? reac("uy",i) : 0));
                    oss << buf;
                }
            }
        }
    }

    // output string
    if (Prms.Extra!=NULL) oss << Prms.Extra;
    if (Prms.PNG)         MPL::SaveFig (FNKey, oss);
    else                  MPL::Show    (oss);

    // write file
    String fn(FNKey);  fn.append(".mpy");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close ();
}

inline void Domain::WriteVTU (char const * FNKey, bool DoExtrapolation) const
{
    // extrapolate results
    if (DoExtrapolation) NodalResults ();

    // data
    size_t nn = ActNods.Size(); // number of active nodes
    size_t ne = ActEles.Size(); // number of active elements

    // header
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\"?>\n";
    oss << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    oss << "  <UnstructuredGrid>\n";
    oss << "    <Piece NumberOfPoints=\"" << nn << "\" NumberOfCells=\"" << ne << "\">\n";

    // nodes: coordinates
    oss << "      <Points>\n";
    oss << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    size_t k = 0; oss << "        ";
    std::map<int,int> VertID2LocID;
    for (size_t i=0; i<nn; ++i)
    {
        VertID2LocID[ActNods[i]->Vert.ID] = i;
        oss << "  " << Util::_8s <<          ActNods[i]->Vert.C(0) << " ";
        oss <<         Util::_8s <<          ActNods[i]->Vert.C(1) << " ";
        oss <<         Util::_8s << (NDim==3?ActNods[i]->Vert.C(2):0.0);
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
        size_t nne = (NDim==2 ? NVertsToVTKNVerts2D[ActEles[i]->Con.Size()] : NVertsToVTKNVerts3D[ActEles[i]->Con.Size()]);
        oss << "  ";
        for (size_t j=0; j<nne; ++j) oss << VertID2LocID[ActEles[i]->Con[j]->Vert.ID] << " ";
        k++;
        VTU_NEWLINE (i,k,ne,40/ActEles[i]->Con.Size(),oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    size_t offset = 0;
    for (size_t i=0; i<ne; ++i)
    {
        size_t nne = (NDim==2 ? NVertsToVTKNVerts2D[ActEles[i]->Con.Size()] : NVertsToVTKNVerts3D[ActEles[i]->Con.Size()]);
        offset += nne;
        oss << (k==0?"  ":" ") << offset;
        k++;
        VTU_NEWLINE (i,k,ne,40,oss);
    }
    oss << "        </DataArray>\n";
    oss << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t i=0; i<ne; ++i)
    {
        if (NDim==2) oss << (k==0?"  ":" ") << NVertsToVTKCell2D[ActEles[i]->Con.Size()];
        else         oss << (k==0?"  ":" ") << NVertsToVTKCell3D[ActEles[i]->Con.Size()];
        k++;
        VTU_NEWLINE (i,k,ne,40,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </Cells>\n";

    // data -- nodes
    oss << "      <PointData Scalars=\"point_scalars\">\n";
    for (size_t i=0; i<AllUKeys.Size(); ++i)
    {
        if (!(AllUKeys[i]=="ux" || AllUKeys[i]=="uy" || AllUKeys[i]=="uz"))
        {
            oss << "        <DataArray type=\"Float64\" Name=\"" << AllUKeys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                oss << (k==0?"  ":" ") << ActNods[j]->UOrZero(AllUKeys[i]);
                k++;
                VTU_NEWLINE (j,k,nn,6-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }

    // data -- nodes -- IDs
    oss << "        <DataArray type=\"Float64\" Name=\"" << "tag,id" << "\" NumberOfComponents=\"2\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t j=0; j<nn; ++j)
    {
        oss << "  " << ActNods[j]->Vert.Tag << " " << ActNods[j]->Vert.ID << " ";
        k++;
        VTU_NEWLINE (j,k,nn,40/2-1,oss);
    }
    oss << "        </DataArray>\n";

    // data -- nodes -- displacements
    if (HasDisps)
    {
        oss << "        <DataArray type=\"Float64\" Name=\"" << "u" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t i=0; i<nn; ++i)
        {
            oss << "  " << Util::_8s <<            ActNods[i]->U("ux") << " ";
            oss <<         Util::_8s <<            ActNods[i]->U("uy") << " ";
            oss <<         Util::_8s << (NDim==3 ? ActNods[i]->U("uz") : 0.0);
            k++;
            VTU_NEWLINE (i,k,nn,6/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }

    // data -- nodes -- extrapolated values
    for (size_t i=0; i<AllEKeys.Size(); ++i)
    {
        if (!(AllEKeys[i]=="vx" || AllEKeys[i]=="vy" || AllEKeys[i]=="vz"))
        {
            oss << "        <DataArray type=\"Float64\" Name=\"" << AllEKeys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                double cnt = NodResCount[ActNods[j]][i];
                oss << (k==0?"  ":" ") << NodResults[ActNods[j]][i] / (cnt>0.0 ? cnt : 1.0);
                k++;
                VTU_NEWLINE (j,k,nn,6-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }

    // data -- nodes -- velocities
    if (HasVeloc)
    {
        size_t idx_vx=-1, idx_vy=-1, idx_vz=-1;
        for (size_t i=0; i<AllEKeys.Size(); ++i)
        {
            if (AllEKeys[i]=="vx") idx_vx = i;
            if (AllEKeys[i]=="vy") idx_vy = i;
            if (AllEKeys[i]=="vz") idx_vz = i;
        }
        oss << "        <DataArray type=\"Float64\" Name=\"" << "v" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t j=0; j<nn; ++j)
        {
            double cnt_x =            NodResCount[ActNods[j]][idx_vx];
            double cnt_y =            NodResCount[ActNods[j]][idx_vy];
            double cnt_z = (NDim==3 ? NodResCount[ActNods[j]][idx_vz] : 0.0);
            oss << "  " << Util::_8s <<            NodResults[ActNods[j]][idx_vx] / (cnt_x>0.0 ? cnt_x : 1.0) << " ";
            oss <<         Util::_8s <<            NodResults[ActNods[j]][idx_vy] / (cnt_y>0.0 ? cnt_y : 1.0) << " ";
            oss <<         Util::_8s << (NDim==3 ? NodResults[ActNods[j]][idx_vz] / (cnt_z>0.0 ? cnt_z : 1.0) : 0.0);
            k++;
            VTU_NEWLINE (j,k,nn,6/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }

    // data -- nodes -- specific water discharge
    if (HasWDisch)
    {
        size_t idx_vx=-1, idx_vy=-1, idx_vz=-1;
        for (size_t i=0; i<AllEKeys.Size(); ++i)
        {
            if (AllEKeys[i]=="qwx") idx_vx = i;
            if (AllEKeys[i]=="qwy") idx_vy = i;
            if (AllEKeys[i]=="qwz") idx_vz = i;
        }
        oss << "        <DataArray type=\"Float64\" Name=\"" << "qw" << "\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        k = 0; oss << "        ";
        for (size_t j=0; j<nn; ++j)
        {
            double cnt_x =            NodResCount[ActNods[j]][idx_vx];
            double cnt_y =            NodResCount[ActNods[j]][idx_vy];
            double cnt_z = (NDim==3 ? NodResCount[ActNods[j]][idx_vz] : 0.0);
            oss << "  " << Util::_8s <<            NodResults[ActNods[j]][idx_vx] / (cnt_x>0.0 ? cnt_x : 1.0) << " ";
            oss <<         Util::_8s <<            NodResults[ActNods[j]][idx_vy] / (cnt_y>0.0 ? cnt_y : 1.0) << " ";
            oss <<         Util::_8s << (NDim==3 ? NodResults[ActNods[j]][idx_vz] / (cnt_z>0.0 ? cnt_z : 1.0) : 0.0);
            k++;
            VTU_NEWLINE (j,k,nn,6/3-1,oss);
        }
        oss << "        </DataArray>\n";
    }

    // extra values from callback (solution?)
    if (VTUfunc!=NULL)
    {
        SDPair pairs;
        (*VTUfunc) (Vec3_t(0,0,0),0, pairs);
        for (size_t i=0; i<pairs.Keys.Size(); ++i)
        {
            oss << "        <DataArray type=\"Float64\" Name=\"" << pairs.Keys[i] << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
            k = 0; oss << "        ";
            for (size_t j=0; j<nn; ++j)
            {
                SDPair p;
                (*VTUfunc) (Nods[j]->Vert.C, Time, p);
                oss << (k==0?"  ":" ") << p(pairs.Keys[i]);
                k++;
                VTU_NEWLINE (j,k,nn,6-1,oss);
            }
            oss << "        </DataArray>\n";
        }
    }

    // data -- nodes -- end
    oss << "      </PointData>\n";

    // data -- elements
    oss << "      <CellData Scalars=\"cell_scalars\">\n";
    oss << "        <DataArray type=\"Float64\" Name=\"" << "Tag,ID" << "\" NumberOfComponents=\"2\" format=\"ascii\">\n";
    k = 0; oss << "        ";
    for (size_t j=0; j<ne; ++j)
    {
        oss << "  " << ActEles[j]->Cell.Tag << " " << ActEles[j]->Cell.ID << " ";
        k++;
        VTU_NEWLINE (j,k,ne,40/2-1,oss);
    }
    oss << "        </DataArray>\n";
    oss << "      </CellData>\n";

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
    bool reactions_calculated  = false;
    bool nodresults_calculated = false;
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
                    if      (AllUKeys.Has(key)) err[j] = fabs(Nods[j]->U(key.CStr()) - NodSol(key,j));
                    else if (AllFKeys.Has(key)) err[j] = fabs(Nods[j]->F(key.CStr()) - NodSol(key,j));
                    else
                    {
                        if (!nodresults_calculated)
                        {
                            NodalResults ();
                            nodresults_calculated = true;
                        }
                        long pos = AllEKeys.Find(key);
                        if (pos<0) throw new Fatal("Domain::CheckErrorNods: Could not find key=%s in AllEKeys (extrapolated variables) array",key.CStr());
                        double cnt = NodResCount[Nods[j]][pos];
                        err[j] = fabs(NodResults[Nods[j]][pos]/(cnt>0.0?cnt:1.0) - NodSol(key,j));
                    }
                }
                else err[j] = 0.0;
            }
        }

        // summary
        if (!NodTol.HasKey(key)) throw new Fatal("Domain::CheckErrorNods: NodTol structure with the tolerances must have key=%s",key.CStr());
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
            Array<SDPair> loc_res; // local results: size==number of nodes in element
            ActEles[j]->StateAtNodes (loc_res);
            SDPair ave; // average
            ave = loc_res[0];
            for (size_t k=1; k<loc_res.Size(); ++k)
            {
                ave.AddVal (key.CStr(), loc_res[k](key));
            }
            err [j] = fabs(ave(key)/loc_res.Size() - EleSol(key,j));
        }

        // summary
        if (!EleTol.HasKey(key)) throw new Fatal("Domain::CheckErrorEles: EleTol structure with the tolerances must have key=%s",key.CStr());
        double max_err = err.TheMax();
        double tol     = EleTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }

    return error;
}

inline bool Domain::CheckErrorIPs (Table const & EleSol, SDPair const & EleTol) const
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
            if (Eles[j]->GE==NULL) throw new Fatal("Domain::CheckErrorIPs: This method works only when GE (geometry element) is not NULL");
            Array<SDPair> res;
            Eles[j]->StateAtIPs (res);
            for (size_t k=0; k<Eles[j]->GE->NIP; ++k)
            {
                size_t row = k + j*Eles[j]->GE->NIP;
                err [row] = fabs(res[k](key) - EleSol(key,row));
            }
        }

        // summary
        if (!EleTol.HasKey(key)) throw new Fatal("Domain::CheckErrorIPs: EleTol structure with the tolerances must have key=%s",key.CStr());
        double max_err = err.TheMax();
        double tol     = EleTol(key);
        std::cout << Util::_4<< key << Util::_8s<<err.TheMin() << Util::_8s<<err.Mean();
        std::cout << (max_err>tol ? TERM_RED : TERM_GREEN) << Util::_8s<<max_err << TERM_RST << Util::_8s<<err.Norm() << "\n";
        if (max_err>tol) error = true;
    }
    std::cout << "\n";

    return error;
}

inline void Domain::SaveState (char const * FileKey) const
{
#ifdef USE_HDF5
    // open HDF5 file
    String fn(FileKey);
    fn.append(".hdf5");
    hid_t hdf = H5Fcreate (fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // number of nodes and elements
    int     nods_sz[1] = { static_cast<int>(Nods.Size()) };
    int     eles_sz[1] = { static_cast<int>(Eles.Size()) };
    double  time   [1] = { Time };
    hsize_t unit_sz[1] = { 1 };
    H5LTmake_dataset_int    (hdf, "/NNods", 1, unit_sz, nods_sz);
    H5LTmake_dataset_int    (hdf, "/NEles", 1, unit_sz, eles_sz);
    H5LTmake_dataset_double (hdf, "/Time",  1, unit_sz, time);

    // define HDF5 table with node results
    size_t ncols    = AllUKeys.Size() + AllFKeys.Size();
    size_t row_size = sizeof(double) * ncols;
    size_t col_offsets[ncols];
    //size_t col_sizes  [ncols];
    hid_t  col_types  [ncols];
    for (size_t i=0; i<ncols; ++i)
    {
        col_offsets[i] = i * sizeof(double);
        //col_sizes  [i] =     sizeof(double);
        col_types  [i] = H5T_NATIVE_DOUBLE;
    }
    char const * col_labels[ncols];
    for (size_t i=0; i<AllUKeys.Size(); ++i) col_labels[                i] = AllUKeys[i].CStr();
    for (size_t i=0; i<AllFKeys.Size(); ++i) col_labels[AllUKeys.Size()+i] = AllFKeys[i].CStr();

    // create HDF5 table with node results
    H5TBmake_table ("U and F values at nodes", hdf, "/Nods", /*nfields*/ncols, /*nrecords*/Nods.Size(), row_size, col_labels, col_offsets, col_types, /*chunk_size*/10, /*fill_data*/NULL, /*compress*/false, /*data*/NULL);

    // fill HDF5 table with node results
    int    col_idx[1];
    double dbl_dat[1];
    size_t dbl_siz[1] = { sizeof(double) };
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        for (size_t j=0; j<AllUKeys.Size(); ++j)
        {
            col_idx[0] = j;
            dbl_dat[0] = Nods[i]->UOrZero(AllUKeys[j]);
            H5TBwrite_fields_index (hdf, "/Nods", /*nfields*/1, col_idx, /*row_idx*/i, /*nrecords*/1, sizeof(double), /*offset*/0, dbl_siz, dbl_dat);
        }
        for (size_t j=0; j<AllFKeys.Size(); ++j)
        {
            col_idx[0] = AllFKeys.Size() + j;
            dbl_dat[0] = Nods[i]->FOrZero(AllFKeys[j]);
            H5TBwrite_fields_index (hdf, "/Nods", /*nfields*/1, col_idx, /*row_idx*/i, /*nrecords*/1, sizeof(double), /*offset*/0, dbl_siz, dbl_dat);
        }
    }

    // save elements' states
    String  buf;
    hsize_t sta_size[1];
    hid_t eles = H5Gcreate (hdf, "/Eles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        buf.Printf ("ele_%d",i);
        hid_t ele = H5Gcreate (eles, buf.CStr(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        for (size_t j=0; j<Eles[i]->Sta.Size(); ++j)
        {
            buf.Printf ("sta_%d",j);
            Array<double> pck;
            Eles[i]->Sta[j]->Pack (pck);
            sta_size[0] = pck.Size();
            H5LTmake_dataset_double (ele, buf.CStr(), /*rank*/1, sta_size, pck.GetPtr());
        }
    }

    // close file
    H5Fclose (hdf);
#else
    throw new Fatal("Domain::SaveState: This method needs USE_HDF5 defined");
#endif
}

inline void Domain::LoadState (char const * FileKey)
{
#ifdef USE_HDF5
    // open HDF5 file
    String fn(FileKey);
    fn.append(".hdf5");
    hid_t hdf = H5Fopen (fn.CStr(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // number of nodes and elements
    int     nods_sz[1];
    int     eles_sz[1];
    double  time   [1];
    H5LTread_dataset_int    (hdf, "/NNods", nods_sz);
    H5LTread_dataset_int    (hdf, "/NEles", eles_sz);
    H5LTread_dataset_double (hdf, "/Time",  time);
    if (nods_sz[0]!=(int)Nods.Size()) throw new Fatal("Domain::LoadState(%s) Number of nodes in file (%d) is different from number of nodes in Domain (%zd)",FileKey,nods_sz[0],Nods.Size());
    if (eles_sz[0]!=(int)Eles.Size()) throw new Fatal("Domain::LoadState(%s) Number of elemnts in file (%d) is different from number of elements in Domain (%zd)",FileKey,eles_sz[0],Eles.Size());
    Time = time[0];

    // read table with node results
    int    col_idx[1];
    double dbl_dat[1];
    size_t dbl_siz[1] = { sizeof(double) };
    for (size_t i=0; i<Nods.Size(); ++i)
    {
        for (size_t j=0; j<AllUKeys.Size(); ++j)
        {
            col_idx[0] = j;
            H5TBread_fields_index (hdf, "/Nods", /*nfields*/1, col_idx, /*row_idx*/i, /*nrecords*/1, sizeof(double), /*offset*/0, dbl_siz, dbl_dat);
            Nods[i]->TrySetU (AllUKeys[j], dbl_dat[0]);
        }
        for (size_t j=0; j<AllFKeys.Size(); ++j)
        {
            col_idx[0] = AllFKeys.Size() + j;
            H5TBread_fields_index (hdf, "/Nods", /*nfields*/1, col_idx, /*row_idx*/i, /*nrecords*/1, sizeof(double), /*offset*/0, dbl_siz, dbl_dat);
            Nods[i]->TrySetF (AllFKeys[j], dbl_dat[0]);
        }
    }

    // read elements' states
    String  buf;
    //hsize_t sta_size[1];
    hid_t eles = H5Gopen (hdf, "/Eles", H5P_DEFAULT);
    for (size_t i=0; i<Eles.Size(); ++i)
    {
        buf.Printf ("ele_%d",i);
        hid_t ele = H5Gopen (eles, buf.CStr(), H5P_DEFAULT);
        for (size_t j=0; j<Eles[i]->Sta.Size(); ++j)
        {
            buf.Printf ("sta_%d",j);
            Array<double> pck(Eles[i]->Sta[j]->PckSize());
            //sta_size[0] = pck.Size();
            H5LTread_dataset_double (ele, buf.CStr(), pck.GetPtr());
            Eles[i]->Sta[j]->Unpack (pck);
        }
    }

    // close file
    H5Fclose (hdf);

#else
    throw new Fatal("Domain::LoadState: This method needs USE_HDF5 defined");
#endif
}

inline void Domain::SetGeostatic (SDPair const & Data)
{
    // Prp must have: geosta, surf, K0, gamNat, gamSat
    // Prp may  have: wtable, gamW, pospw (to force only positive values of pw)

    /*
        double z_surf    = Prp("surf");
        double K0        = Prp("K0");
        bool   pos_pw    = (Prp.HasKey("pospw") ? Prp("pospw")>0 : false);
        bool   has_water = Prp.HasKey("water");
        double z_water   = 0;
        if (GTy==pse_t)          throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, geometry cannot be of 'plane-stress' (pse) type");
        if (!Prp.HasKey("K0"))   throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, 'K0' must be provided in 'Prp' dictionary");
        if (!Prp.HasKey("surf")) throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, 'surf' must be provided in 'Prp' dictionary");
        if (has_water)
        {
            z_water = Prp("water");
            if (z_water>z_surf) throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, 'water' must be smaller than or equal to 'surf'");
            // TODO: this last condition can be removed, but sv calculation in the next lines must be corrected
        }
        if (Mdl->GamW  <0.0) throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, 'gamW' must be positive");
        if (Mdl->GamNat<0.0) throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, 'gamNat' must be positive");
        if (Mdl->GamSat<0.0) throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, 'gamSat' must be positive");
        for (size_t i=0; i<GE->NIP; ++i)
        {
            // elevation of point
            Vec_t X;
            CoordsOfIP (i, X);
            double z = (NDim==2 ? X(1) : X(2));
            if (z>z_surf) throw new Fatal("EquilibElem::EquilibElem: 'surf' must be greater than any point in the domain.\n\tThere is a point [%g,%g,%g] above surf=%g.",X(0),X(1),(NDim==3?X(2):0.0),z_surf);

            // pore-water pressure and total vertical stress
            double pw = 0.0;
            double sv;
            if (has_water)
            {
                double hw = z_water-z; // column of water
                pw = (hw>0.0 ? Mdl->GamW*hw : (pos_pw ? 0.0 : Mdl->GamW*hw));
                if (z>z_water) sv = (z_surf-z)*Mdl->GamNat;
                else sv = (z_surf-z_water)*Mdl->GamNat + (z_water-z)*Mdl->GamSat;
            }
            else sv = (z_surf-z)*Mdl->GamNat;   
            sv *= (-1.0); // convert soil-mech. convention to classical mech. convention

            // vertical and horizontal effective stresses
            double sv_ = sv + pw;  // effective vertical stresss
            double sh_ = K0*sv_;   // effective horizontal stress
            double sh  = sh_ - pw; // total horizontal stress

            // set stress tensor with __total__ stresses
            Vec_t & sig = static_cast<EquilibState *>(Sta[i])->Sig;
            if (NDim==2) sig = sh, sv, sh, 0.0;
            else         sig = sh, sh, sv, 0.0,0.0,0.0;
        }


        // pw at IP
        if (geosta)
        {
            // elevation of point
            Vec_t X;
            CoordsOfIP (i, X);
            double z = (NDim==2 ? X(1) : X(2));

            // pore-water pressure
            double hw = Prp("water")-z; // column of water
            double pw = (hw>0.0 ? gamW*hw : (pos_pw ? 0.0 : gamW*hw));
            ini.Set ("pw", pw);
        }
    */
}

// Internal methods

inline void Domain::NodalResults (bool OnlyOutNods) const
{
    // extrapolate from elements' IPs to nodes
    Array<Node*>    const * nods;
    Array<Element*> const * eles;
    if (OnlyOutNods)
    {
        nods = &OutNods;
        eles = &OutEles;
    }
    else
    {
        nods = &ActNods;
        eles = &ActEles;
    }
    for (size_t i=0; i<nods->Size(); ++i)
    {
        NodResults [(*nods)[i]].Resize    (AllEKeys.Size()); // results at nodes
        NodResCount[(*nods)[i]].Resize    (AllEKeys.Size()); // count how many times a varable was added to a node
        NodResults [(*nods)[i]].SetValues (0.0);
        NodResCount[(*nods)[i]].SetValues (0.0);
    }
    for (size_t i=0; i<eles->Size(); ++i)
    {
        Array<SDPair> loc_res; // local results: size==number of nodes in element
        (*eles)[i]->StateAtNodes (loc_res);
        for (size_t j=0; j<AllEKeys.Size(); ++j)
        {
            if (loc_res[0].HasKey(AllEKeys[j]))
            {
                for (size_t k=0; k<(*eles)[i]->Con.Size(); ++k)
                {
                    if (OnlyOutNods)
                    { 
                        // skip nodes that are not in array OutNods
                        if (!nods->Has((*eles)[i]->Con[k])) continue;
                    }
                    NodResults [(*eles)[i]->Con[k]][j] += loc_res[k](AllEKeys[j]);
                    NodResCount[(*eles)[i]->Con[k]][j] += 1.0;
                }
            }
        }
    }
}

inline void Domain::OutResults (char const * VTUFName, bool HeaderOnly)
{
    // VTU
    if (VTUFName!=NULL)
    {
        // extrapolate values to (all) nodes
        NodalResults ();

        // write file
        String fkey;
        fkey.Printf ("%s_%08d", VTUFName, IdxOut);
        WriteVTU    (fkey.CStr(), /*do_extrapolation*/false);

        // next index of files
        IdxOut++;
    }
    else if (OutNods.Size()>0)
    {
        // extrapolate values to (some) nodes
        if (!HeaderOnly) NodalResults (/*only_outnods*/true);
    }

    // nodes
    for (size_t i=0; i<OutNods.Size(); ++i)
    {
        if (HeaderOnly)
        {
            (*FilNods[i]) << Util::_8s << "Time";
            for (size_t j=0; j<AllUKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << AllUKeys[j];
            for (size_t j=0; j<AllUKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << ("v"+AllUKeys[j].substr(1,AllUKeys[j].size()));
            for (size_t j=0; j<AllFKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << AllFKeys[j];
            for (size_t j=0; j<AllEKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << AllEKeys[j];
            (*FilNods[i]) << "\n";
        }
        else
        {
            (*FilNods[i]) << Util::_8s << Time;
            for (size_t j=0; j<AllUKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->UOrZero(AllUKeys[j]);
            for (size_t j=0; j<AllUKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->VOrZero(AllUKeys[j]);
            for (size_t j=0; j<AllFKeys.Size(); ++j) (*FilNods[i]) << Util::_8s << OutNods[i]->FOrZero(AllFKeys[j]);
            for (size_t j=0; j<AllEKeys.Size(); ++j)
            {
                double cnt = NodResCount[OutNods[i]][j];
                (*FilNods[i]) << Util::_8s << NodResults[OutNods[i]][j] / (cnt>0.0 ? cnt : 1.0);
            }
            (*FilNods[i]) << "\n";
        }
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
    NodesReactions.SetZero (buf.CStr(), NodsWithPU.Size());

    // find reactions
    size_t k = 0;
    for (size_t i=0; i<NodsWithPU.Size(); ++i)
    {
        Node const * nod = NodsWithPU[i];
        NodesReactions("Node", k) = nod->Vert.ID;
        std::map<String,double> R;
        nod->Reactions (R);
        for (std::map<String,double>::const_iterator it=R.begin(); it!=R.end(); ++it)
        {
            NodesReactions(it->first, k) = it->second;
            SumReactions[it->first] += it->second;
        }
        k++;
    }
}


}; // namespace FEM

#endif // MECHSYS_DOMAIN_H
