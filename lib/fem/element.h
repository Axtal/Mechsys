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

#ifndef MECHSYS_FEM_ELEMENT_H
#define MECHSYS_FEM_ELEMENT_H

// MechSys
#include <mechsys/geomtype.h>
#include <mechsys/mesh/mesh.h>
#include <mechsys/fem/node.h>
#include <mechsys/fem/geomelem.h>
#include <mechsys/models/model.h>
#include <mechsys/util/maps.h>
#include <mechsys/util/fatal.h>
#include <mechsys/linalg/matvec.h>

namespace FEM
{

typedef double (*pCalcM) (double t);              ///< M(ultiplier) callback
typedef std::pair<IntDbl_t,pCalcM>      BCData_t; ///< Bry conds data (idxDOF,valueBCs,multiplier)
typedef std::map<Node const *,BCData_t> NodBCs_t; ///< Map Node => bry conds. (idxDOF,valueBCS)


class Element
{
public:
    // Constructor
    Element (int                  NDim,   ///< Space dimension
             Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
             Model        const * Mdl,    ///< Model
             SDPair       const & Prp,    ///< Properties
             SDPair       const & Ini,    ///< Initial values
             Array<Node*> const & Nodes); ///< Array with all nodes (used to set the connectivity)

    // Destructor
    ~Element ();

    // Methods
    virtual void IncNLocDOF   (size_t & NEq)                             const {} ///< Increment the number of local DOFs
    virtual void GetLoc       (Array<size_t> & Loc)                      const;   ///< Get location vector for mounting K/M matrices
    virtual void BackupState  ()                                         const;   ///< Backup element state
    virtual void RestoreState ()                                         const;   ///< Restore element state
    virtual void SetBCs       (size_t IdxEdgeOrFace, SDPair const & BCs,
                               NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM) {}     ///< Set boundary conditions
    virtual void ClrBCs       ()                                           {}     ///< Clear boundary conditions
    virtual void Gravity      (NodBCs_t & pF, pCalcM CalcM, double gAccel) {}     ///< Apply gravity (body forces)
    virtual void Deactivate   (NodBCs_t & pF, pCalcM CalcM, double gAccel,
                               NodBCs_t & pU)                              {}     ///< Deactivate element
    virtual void CalcFint     (Vec_t * F_int=NULL)                       const { if (F_int!=NULL) for (size_t i=0; i<size((*F_int)); ++i) (*F_int)(i) = 0.0; }
    virtual void CalcK        (Mat_t & K)                                const { throw new Fatal("Element::CalcK: Method not implement for this element"); }
    virtual void CalcM        (Mat_t & M)                                const { throw new Fatal("Element::CalcM: Method not implement for this element"); }
    virtual void CalcC        (Mat_t & C)                                const { throw new Fatal("Element::CalcC: Method not implement for this element"); }
    virtual void CalcK        (Vec_t const & U, double Alpha, double dt, Mat_t & KK, Vec_t & dF) const { throw new Fatal("Element::CalcK (with U, Alpha, Dt => HydroMech): Method not implement for this element"); }
    virtual void CalcKCM      (Mat_t & KK, Mat_t & CC, Mat_t & MM)       const { throw new Fatal("Element::CalcKCM: Method not implement for this element"); }
    virtual void UpdateState  (Vec_t const & dU, Vec_t * F_int=NULL)     const {}
    virtual void StateKeys    (Array<String> & Keys)                     const {} ///< Get state keys, ex: sx, sy, sxy, ex, ey, exy
    virtual void StateAtIP    (SDPair & KeysVals, int IdxIP)             const {} ///< Get state at IP
    virtual void StateAtIPs   (Array<SDPair> & Results)                  const;   ///< Get state (internal values: sig, eps) at all integration points
    virtual void StateAtCt    (SDPair & Results)                         const;   ///< Get state at centroid
    virtual void StateAtNodes (Array<SDPair> & Results)                  const;   ///< Get state (internal values: sig, eps) at all nodes (applies extrapolation)
    virtual void Centroid     (Vec_t & X)                                const;   ///< Centroid of element
    virtual void Draw         (std::ostream & os, double MaxDist)        const;   ///< Draw element with MatPlotLib

    // Methods that depend on GE
    void CoordMatrix   (Mat_t & C)                 const; ///< Matrix with coordinates of nodes
    void FCoordMatrix  (size_t IdxFace, Mat_t & C) const; ///< Matrix with coordinates of face nodes
    void ShapeMatrix   (Mat_t & M)                 const; ///< Matrix with shape functions evaluated at IPs
    void CalcShape     (Mat_t const & C,  IntegPoint const & IP,  double & detJ, double & Coef) const; ///< Results are in GE->N
    void CalcFaceShape (Mat_t const & FC, IntegPoint const & FIP, double & detJ, double & Coef) const; ///< Results are in GE->FN
    void CoordsOfIP    (size_t IdxIP, Vec_t & X)   const; ///< x-y-z coordinates of integration point (IP)

    // Data
    int                NDim;   ///< Space dimension
    Mesh::Cell const & Cell;   ///< Geometric information: ID, Tag, connectivity
    Model      const * Mdl;    ///< The model
    GeomElem         * GE;     ///< Element for geometry definition
    bool               Active; ///< Active element ?
    GeomType           GTy;    ///< Geometry type
    Array<Node*>       Con;    ///< Connectivity
    Array<State*>      Sta;    ///< State at the centre of the element or at each IP
    Array<String>      UKeys;  ///< DOF keys such as 'ux', 'uy', 'uz'
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Element::Element (int TheNDim, Mesh::Cell const & TheCell, Model const * TheMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : NDim(TheNDim), Cell(TheCell), Mdl(TheMdl), GE(NULL), Active(Prp.HasKey("active") ? Prp("active") : true),
      GTy(SDPairToGType(Prp,(NDim==3?"d3d":"d2d")))
{
    // connectivity
    Con.Resize (Cell.V.Size());
    for (size_t i=0; i<Con.Size(); ++i)
    {
        Con[i] = Nodes[Cell.V[i]->ID];
        if (Active) Con[i]->NShares++;
    }

    // geometry element
    if (Prp.HasKey("geom"))
    {
        String geom_name;
        GEOM.Val2Key (Prp("geom"), geom_name);
        GE = AllocGeomElem (geom_name, NDim);

        // set number of integration points
        if (Prp.HasKey("nip")) GE->SetIPs (static_cast<int>(Prp("nip")));

        // check number of nodes
        if (Con.Size()!=GE->NN) throw new Fatal("Element::Element: Number of vertices (%d) of an element (%s) in mesh is incorrect. It must be equal to %d for %s", Con.Size(),geom_name.CStr(),GE->NN,geom_name.CStr());
    }

    // check model
    if (Mdl!=NULL) { if (GTy!=Mdl->GTy) throw new Fatal("Element::Element: Element geometry type (%s) must be equal to Model geometry type (%s)", GTypeToStr(GTy).CStr(), GTypeToStr(Mdl->GTy).CStr()); }
}

inline Element::~Element ()
{
    if (GE!=NULL) delete GE;
    for (size_t i=0; i<Sta.Size(); ++i) { if (Sta[i]!=NULL) delete Sta[i]; }
}

inline void Element::GetLoc (Array<size_t> & Loc) const
{
    size_t nnods = Con.Size();
    size_t ndofs = UKeys.Size();
    Loc.Resize (nnods*ndofs);
    for (size_t i=0; i<nnods; ++i)
    {
        for (size_t j=0; j<ndofs; ++j)
        {
            size_t idx = Con[i]->UMap(UKeys[j]); // index in Node corresponding to each DOF
            Loc[i*ndofs+j] = Con[i]->EQ[idx];
        }
    }
}

inline void Element::BackupState () const
{
    for (size_t i=0; i<Sta.Size(); ++i) Sta[i]->Backup();
}

inline void Element::RestoreState () const 
{
    for (size_t i=0; i<Sta.Size(); ++i) Sta[i]->Restore();
}

inline void Element::CoordMatrix (Mat_t & C) const
{
    if (GE==NULL) throw new Fatal("Element::CoordMatrix: This method works only when GE (geometry element) is not NULL");

    C.change_dim (GE->NN, NDim);
    for (size_t i=0; i<GE->NN; ++i)
    for (int    j=0; j<NDim;   ++j)
        C(i,j) = Con[i]->Vert.C[j];
}

inline void Element::FCoordMatrix (size_t IdxFace, Mat_t & C) const
{
    if (GE==NULL) throw new Fatal("Element::FCoordMatrix: This method works only when GE (geometry element) is not NULL");

    C.change_dim (GE->NFN, NDim);
    for (size_t i=0; i<GE->NFN; ++i)
    for (int    j=0; j<NDim;    ++j)
        C(i,j) = Con[GE->FNode(IdxFace,i)]->Vert.C[j];
}

inline void Element::ShapeMatrix (Mat_t & M) const
{
    if (GE==NULL) throw new Fatal("Element::ShapeMatrix: This method works only when GE (geometry element) is not NULL");

    /* Shape functions matrix:
              _                                              _
             |   N11      N12      N13      ...  N1(NN)       |
             |   N21      N22      N23      ...  N2(NN)       |
         M = |   N31      N32      N33      ...  N3(NN)       |
             |          ...............................       |
             |_  N(NIP)1  N(NIP)2  N(NIP)3  ...  N(NIP)(NN)  _| (NIP x NN)
    */
    M.change_dim (GE->NIP, GE->NN);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        GE->Shape (GE->IPs[i].r, GE->IPs[i].s, GE->IPs[i].t);
        for (size_t j=0; j<GE->NN; ++j) M(i,j) = GE->N(j);
    }
}

inline void Element::CalcShape (Mat_t const & C, IntegPoint const & IP, double & detJ, double & Coef) const
{
    if (GE==NULL) throw new Fatal("Element::ShapeMatrix: This method works only when GE (geometry element) is not NULL");

    // geometric data
    GE->Shape  (IP.r, IP.s, IP.t);
    GE->Derivs (IP.r, IP.s, IP.t);

    // Jacobian and its determinant
    Mat_t J(GE->dNdR * C); // J = dNdR * C
    detJ = Det(J);

    // coefficient used during integration
    Coef = detJ*IP.w;
}

inline void Element::CalcFaceShape (Mat_t const & FC, IntegPoint const & FIP, double & detJ, double & Coef) const
{
    if (GE==NULL) throw new Fatal("Element::ShapeMatrix: This method works only when GE (geometry element) is not NULL");

    // geometric data
    GE->FaceShape  (FIP.r, FIP.s);
    GE->FaceDerivs (FIP.r, FIP.s);

    // face/edge Jacobian and its determinant
    Mat_t J(GE->FdNdR * FC);
    detJ = Det(J);

    // coefficient used during integration
    Coef = detJ*FIP.w;
}

inline void Element::StateAtIPs (Array<SDPair> & Results) const
{
    if (GE==NULL) throw new Fatal("Element::StateAtIPs: This method works only when GE (geometry element) is not NULL");

    // one set of results per IP
    Results.Resize (GE->NIP);
    for (size_t i=0; i<GE->NIP; ++i) StateAtIP (Results[i], i);
}

inline void Element::StateAtCt (SDPair & Results) const
{
    if (GE==NULL) throw new Fatal("Element::StateAtCt: This method works only when GE (geometry element) is not NULL");

    // keys
    Array<String> keys;
    StateKeys (keys);

    // state at IPs
    Array<SDPair> state_at_IPs;
    StateAtIPs (state_at_IPs);

    // average
    Results.clear();
    for (size_t k=0; k<keys.Size(); ++k)
    {
        // average values at IPs
        double val = 0.0;
        for (size_t i=0; i<GE->NIP; ++i) val += state_at_IPs[i](keys[k]);
        val /= GE->NIP;

        // set results
        Results.Set (keys[k].CStr(), val);
    }
}

inline void Element::StateAtNodes (Array<SDPair> & Results) const
{
    if (GE==NULL) throw new Fatal("Element::StateAtNodes: This method works only when GE (geometry element) is not NULL");

    // resize results array (one set of results per node)
    Results.Resize (GE->NN);

    // shape func matrix
    Mat_t M, Mi;
    ShapeMatrix (M);
    Inv (M, Mi);

    // keys
    Array<String> keys;
    StateKeys (keys);

    // state at IPs
    Array<SDPair> state_at_IPs;
    StateAtIPs (state_at_IPs);

    // extrapolate
    Vec_t val_at_IPs (GE->NIP); // values at IPs
    Vec_t val_at_Nods(GE->NN);  // values at nodes
    for (size_t k=0; k<keys.Size(); ++k)
    {
        // gather values from IPs
        for (size_t i=0; i<GE->NIP; ++i) val_at_IPs(i) = state_at_IPs[i](keys[k]);

        // extrapolate to nodes
        val_at_Nods = Mi * val_at_IPs;

        // scatter values to nodes
        for (size_t i=0; i<GE->NN; ++i) Results[i].Set (keys[k].CStr(), val_at_Nods(i));
    }
}

inline void Element::Centroid (Vec_t & X) const
{
    if (GE==NULL) throw new Fatal("Element::Centroid: This method works only when GE (geometry element) is not NULL");

    GE->Shape (GE->Rct.r, GE->Rct.s, GE->Rct.t);
    X.change_dim (NDim);
    set_to_zero  (X);
    for (size_t i=0; i<GE->NN; ++i)
    for (int    j=0; j<NDim;   ++j)
        X(j) += GE->N(i)*Con[i]->Vert.C[j];
}

inline void Element::CoordsOfIP (size_t IdxIP, Vec_t & X) const
{
    if (GE==NULL) throw new Fatal("Element::CoordsOfIP: This method works only when GE (geometry element) is not NULL");

    //std::cout << "\n" << GE->IPs[IdxIP].r << "  " << GE->IPs[IdxIP].s << "\n\n";
    GE->Shape (GE->IPs[IdxIP].r, GE->IPs[IdxIP].s, GE->IPs[IdxIP].t);
    X.change_dim (NDim);
    set_to_zero  (X);
    for (size_t i=0; i<GE->NN; ++i)
    for (int    j=0; j<NDim;   ++j)
        X(j) += GE->N(i)*Con[i]->Vert.C[j];
}

inline void Element::Draw (std::ostream & os, double MaxDist) const
{
    if (GE==NULL) throw new Fatal("Element::CoordsIP: This method works only when GE (geometry element) is not NULL");

    // draw shape
    os << "dat.append((PH.MOVETO, (" << Con[0]->Vert.C[0] << "," << Con[0]->Vert.C[1] << ")))\n";
    if (GE->NN<=4)
    {
        for (size_t i=1; i<GE->NN; ++i) os << "dat.append((PH.LINETO,    (" << Con[i]->Vert.C[0] << "," << Con[i]->Vert.C[1] << ")))\n";
        os << "dat.append((PH.CLOSEPOLY, (" << Con[0]->Vert.C[0] << "," << Con[0]->Vert.C[1] << ")))\n";
    }
    else if (GE->NN==6)
    {
        os << "dat.append((PH.LINETO,    (" << Con[3]->Vert.C(0) << "," << Con[3]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[1]->Vert.C(0) << "," << Con[1]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[4]->Vert.C(0) << "," << Con[4]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[2]->Vert.C(0) << "," << Con[2]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[5]->Vert.C(0) << "," << Con[5]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.CLOSEPOLY, (" << Con[0]->Vert.C(0) << "," << Con[0]->Vert.C(1) << ")))\n";
    }
    else if (GE->NN==8)
    {
        os << "dat.append((PH.LINETO,    (" << Con[4]->Vert.C(0) << "," << Con[4]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[1]->Vert.C(0) << "," << Con[1]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[5]->Vert.C(0) << "," << Con[5]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[2]->Vert.C(0) << "," << Con[2]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[6]->Vert.C(0) << "," << Con[6]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[3]->Vert.C(0) << "," << Con[3]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.LINETO,    (" << Con[7]->Vert.C(0) << "," << Con[7]->Vert.C(1) << ")))\n";
        os << "dat.append((PH.CLOSEPOLY, (" << Con[0]->Vert.C(0) << "," << Con[0]->Vert.C(1) << ")))\n";
    }

    // model has deviatoric plastic strain ?
    bool has_edp = false;
    if (Mdl->IvNames.Find("edp")>=0) has_edp = true;

    // draw IPs
    os << "x_ips, y_ips = [], []\n";
    os << "x_edp, y_edp = [], []\n";
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Vec_t X;
        CoordsOfIP (i,X);
        os << "x_ips.append(" << X(0) << ")\n";
        os << "y_ips.append(" << X(1) << ")\n";
        if (has_edp)
        {
            SDPair res;
            StateAtIP (res, i);
            if (res("edp")>0.0)
            {
                os << "x_edp.append(" << X(0) << ")\n";
                os << "y_edp.append(" << X(1) << ")\n";
            }
        }
    }
    os << "plot(x_ips,y_ips,'r*',zorder=1)\n";
    if (has_edp) os << "plot(x_edp,y_edp,'ko',ms=9,zorder=2)\n";
}

std::ostream & operator<< (std::ostream & os, Element const & E)
{
    os << Util::_4 << E.Cell.ID << " ";
    os << "[";
    for (size_t i=0; i<E.UKeys.Size(); ++i)
    {
        os << E.UKeys[i];
        if (i==E.UKeys.Size()-1) os << "]";
        else                     os << ",";
    }
    os << (E.Active?" active":" inactive") << " ";
    os << GTypeToStr(E.GTy) << " ";
    if (E.GE!=NULL)  os << E.GE->Name  << " " << "NIP=" << E.GE->NIP << " ";
    else             os << "GE=NULL"   << " ";
    if (E.Mdl!=NULL) os << E.Mdl->Name << " ";
    else             os << "Mdl=NULL";
    os << "(";
    for (size_t i=0; i<E.Con.Size(); ++i)
    {
        os << E.Con[i]->Vert.ID;
        if (i==E.Con.Size()-1) os << ") ";
        else                   os << ",";
    }
    os << E.Cell.Tag;
    return os;
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


SDPair PROB;

typedef Element * (*ElementMakerPtr)(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes);

typedef std::map<String, ElementMakerPtr> ElementFactory_t;

ElementFactory_t ElementFactory;

Element * AllocElement(String const & Name, int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
{
    ElementFactory_t::iterator it = ElementFactory.find(Name);
    if (it==ElementFactory.end()) throw new Fatal("AllocElement: '%s' is not available", Name.CStr());

    Element * ptr = (*it->second)(NDim,Cell,Mdl,Prp,Ini,Nodes);

    return ptr;
}

}; // namespace FEM

#ifdef USE_BOOST_PYTHON
double PyPROB (BPy::str const & Key) { return FEM::PROB(BPy::extract<char const *>(Key)()); }
#endif

#endif // MECHSYS_FEM_ELEMENT
