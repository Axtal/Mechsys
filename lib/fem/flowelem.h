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

#ifndef MECHSYS_FEM_FLOWELEM_H
#define MECHSYS_FEM_FLOWELEM_H

// Std Lib
#include <map>

// MechSys
#include <mechsys/fem/element.h>
#include <mechsys/models/flowstate.h>
#include <mechsys/models/flowupdate.h>

namespace FEM
{

class FlowElem : public Element
{
public:
    // Constructor
    FlowElem (int                  NDim,   ///< Space dimension
              Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
              Model        const * Mdl,    ///< Model
              SDPair       const & Prp,    ///< Properties
              SDPair       const & Ini,    ///< Initial values
              Array<Node*> const & Nodes); ///< Connectivity

    // Methods
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF); ///< If setting body forces, IdxEdgeOrFace is ignored
    void ClrBCs      ();                                                        ///< Clear BCs
    void GetLoc      (Array<size_t> & Loc)                               const; ///< Get location vector for mounting K/M matrices
    void CalcK       (Mat_t & K)                                         const; ///< Stiffness matrix
    void CalcM       (Mat_t & M)                                         const; ///< Mass matrix
	void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL)              const; ///< Update state at IPs
    void StateKeys   (Array<String> & Keys)                              const; ///< Get state keys, ex: vx, vy, gx, gy
    void StateAtIP   (SDPair & KeysVals, int IdxIP)                      const; ///< State at IP

    // Internal methods
    void CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const; ///< Strain-displacement matrix. Coef: coefficient used during integration

    // Data
    bool                 HasConv; // Has convection ?
    std::map<int,double> Bry2h;   // Map: boundary (edge/face) ID ==> to convection coefficient (h)
    double               m;       // Coefficient for mass matrix
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline FlowElem::FlowElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes), HasConv(false)
{
    // check
    if (GE==NULL)  throw new Fatal("FlowElem::FlowElem: GE (geometry element) must be defined");
    if (Mdl==NULL) throw new Fatal("FlowElem::FlowElem: Model must be defined");

    // check GTy
         if (NDim==2 && GTy==d2d_t) { /*OK*/ }
    else if (NDim==3 && GTy==d3d_t) { /*OK*/ }
    else throw new Fatal("FlowElem::FlowElem: For NDim=%d, geometry type (GTy) must be equal to 'd%dd'. GTy=%s is invalid",NDim,NDim,GTypeToStr(GTy).CStr());

    // properties
    m = (Prp.HasKey("m") ? Prp("m") : 0.0);

    // allocate and initialize state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Sta.Push (new FlowState(NDim));
        Mdl->InitIvs (Ini, Sta[i]);
    }
}

inline void FlowElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF)
{
    bool has_s    = BCs.HasKey("s");    // source term
    bool has_H    = BCs.HasKey("H");    // potential
    bool has_flux = BCs.HasKey("flux"); // flux
    bool has_conv = BCs.HasKey("conv"); // convection

    if (has_s || has_flux || has_conv)
    {
        // source term
        if (has_s)
        {
            double detJ, coef;
            double s = BCs("s");
            Mat_t C;
            CoordMatrix (C);
            for (size_t i=0; i<GE->NIP; ++i)
            {
                CalcShape (C, GE->IPs[i], detJ, coef);
                for (size_t j=0; j<GE->NN; ++j)
                {
                    Con[j]->AddToPF("Q", s*coef*GE->N(j), BCF);
                }
            }
        }

        // flux
        if (has_flux)
        {
            double qn = BCs("flux");
            double detJ, coef;
            Mat_t FC;
            FCoordMatrix (IdxEdgeOrFace, FC);
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (FC, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    size_t k = GE->FNode(IdxEdgeOrFace,j);
                    Con[k]->AddToPF("Q", coef*GE->FN(j)*qn, BCF);
                }
            }
        }

        // convection
        else if (has_conv)
        {
            // data
            double h    = BCs("h");    // convection coefficient
            double Tinf = BCs("Tinf"); // temperature of surrounding environment
            HasConv     = true;        // has convection
            Bry2h[IdxEdgeOrFace] = h;  // Map: bry => h

            // add to pF
            double detJ, coef;
            Mat_t FC;
            FCoordMatrix (IdxEdgeOrFace, FC);
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (FC, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    size_t k = GE->FNode(IdxEdgeOrFace,j);
                    Con[k]->AddToPF("Q", coef*GE->FN(j)*h*Tinf, BCF);
                }
            }
        }
    }

    // potential
    else if (has_H)
    {
        double H = BCs("H");
        for (size_t j=0; j<GE->NFN; ++j)
        {
            size_t k = GE->FNode(IdxEdgeOrFace,j);
            Con[k]->SetPU("H", H, BCF);
        }
    }
}

inline void FlowElem::ClrBCs ()
{
    HasConv = false;
    Bry2h.clear();
}

inline void FlowElem::GetLoc (Array<size_t> & Loc) const
{
    Loc.Resize (GE->NN);
    for (size_t i=0; i<GE->NN; ++i) Loc[i] = Con[i]->Eq("H");
}

inline void FlowElem::CalcK (Mat_t & K) const
{
    double detJ, coef;
    Mat_t C, D, B;
    int nrows = Con.Size(); // number of rows in local K matrix
    K.change_dim (nrows,nrows);
    set_to_zero  (K);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Mdl->Stiffness (Sta[i], D);
        CalcB (C, GE->IPs[i], B, detJ, coef);
        Mat_t BtDB(trans(B)*D*B);
        K += coef * BtDB;
    }
    if (HasConv)
    {
        for (std::map<int,double>::const_iterator p=Bry2h.begin(); p!=Bry2h.end(); ++p)
        {
            int    idx_bry = p->first;
            double h       = p->second;

            // add to K
            Mat_t FC;
            FCoordMatrix (idx_bry, FC);
            //Mat_t Kh(nrows,nrows);  set_to_zero(Kh);
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (FC, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    int row = GE->FNode(idx_bry,j);
                    for (size_t k=0; k<GE->NFN; ++k)
                    {
                        int col = GE->FNode(idx_bry,k);
                        K(row,col) += h*coef * GE->FN(j)*GE->FN(k);
                    }
                }
            }
        }
        
    }
}

inline void FlowElem::CalcM (Mat_t & M) const
{
    double detJ, coef;
    Mat_t C;
    int nrows = Con.Size(); // number of rows in local K matrix
    M.change_dim (nrows,nrows);
    set_to_zero  (M);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        CalcShape (C, GE->IPs[i], detJ, coef);
        for (size_t j=0; j<GE->NN; ++j)
        {
            for (size_t k=0; k<GE->NN; ++k)
            {
                M(j,k) += m*coef*GE->N(j)*GE->N(k);
            }
        }
    }
}

inline void FlowElem::CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const
{
    // deriv of shape func w.r.t natural coordinates
    GE->Derivs (IP.r, IP.s, IP.t);

    // Jacobian and its determinant
    Mat_t J(GE->dNdR * C); // J = dNdR * C
    detJ = Det(J);

    // deriv of shape func w.r.t real coordinates
    Mat_t Ji;
    Inv (J, Ji);

    // coefficient used during integration
    Coef = detJ*IP.w;

    // B matrix
    int nrows = Con.Size(); // number of rows in local K matrix
    B.change_dim (NDim,nrows);
    B = Ji * GE->dNdR; // B = dNdX = Inv(J) * dNdR
}

inline void FlowElem::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal potential
    int nrows = Con.Size(); // number of rows in local K matrix
    Vec_t dUe(nrows);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // update state at each IP
    FlowUpdate fu(Mdl);
    double detJ, coef;
    Mat_t  C, B;
    Vec_t  dFe(nrows), dvel(NDim), dgra(NDim);
    set_to_zero (dFe);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // B matrix
        CalcB (C, GE->IPs[i], B, detJ, coef);

        // velocity and gradient increments
        dgra = B * dUe;
        fu.Update (dgra, Sta[i], dvel);

        // element nodal forces
        Vec_t Btdvel(trans(B)*dvel);
        dFe -= coef * (Btdvel); // '-=' because dvel is -k*dgra
    }

    // add contribution to dFe due to convection term
    if (HasConv)
    {
        for (std::map<int,double>::const_iterator p=Bry2h.begin(); p!=Bry2h.end(); ++p)
        {
            int    idx_bry = p->first;
            double h       = p->second;

            // add to dFe
            Mat_t FC;
            FCoordMatrix (idx_bry, FC);
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (FC, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    int row = GE->FNode(idx_bry,j);
                    for (size_t k=0; k<GE->NFN; ++k)
                    {
                        int col = GE->FNode(idx_bry,k);
                        dFe(row) += h*coef * GE->FN(j)*GE->FN(k) * dUe(col);
                    }
                }
            }
        }
    }

    // add results to Fint (internal forces)
    if (F_int!=NULL) for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
}

inline void FlowElem::StateKeys (Array<String> & Keys) const
{
    Keys.Resize (2*NDim + Mdl->NIvs);
    size_t k=0;
    Keys[k++] = "vx";
    Keys[k++] = "vy";
    if (NDim==3)
    {
        Keys[k++] = "vz";
    }
    Keys[k++] = "gx";
    Keys[k++] = "gy";
    if (NDim==3)
    {
        Keys[k++] = "gz";
    }
    for (size_t i=0; i<Mdl->NIvs; ++i) Keys[k++] = Mdl->IvNames[i];
}

inline void FlowElem::StateAtIP (SDPair & KeysVals, int IdxIP) const
{
    Vec_t const & vel = static_cast<FlowState const *>(Sta[IdxIP])->Vel;
    Vec_t const & gra = static_cast<FlowState const *>(Sta[IdxIP])->Gra;

    if (NDim==2)
    {
        KeysVals.Set("vx vy  gx gy",
                     vel(0), vel(1),
                     gra(0), gra(1));
    }
    else
    {
        KeysVals.Set("vx vy vz  gx gy gz",
                     vel(0), vel(1), vel(2),
                     gra(0), gra(1), gra(2));
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * FlowElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new FlowElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int FlowElemRegister()
{
    ElementFactory["Flow"]   = FlowElemMaker;
    ElementVarKeys["Flow2D"] = std::make_pair ("H", "Q");
    ElementVarKeys["Flow3D"] = std::make_pair ("H", "Q");
    PROB.Set ("Flow", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __FlowElem_dummy_int  = FlowElemRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_FLOWELEM
