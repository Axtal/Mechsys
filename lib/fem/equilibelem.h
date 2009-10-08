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

#ifndef MECHSYS_FEM_EQUILIBELEM_H
#define MECHSYS_FEM_EQUILIBELEM_H

// MechSys
#include "fem/element.h"
#include "models/equilibstate.h"
#include "models/stressupdate.h"

namespace FEM
{

class EquilibElem : public Element
{
public:
    // Constructor
    EquilibElem (int                  NDim,   ///< Space dimension
                 Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                 Model        const * Mdl,    ///< Model
                 SDPair       const & Prp,    ///< Properties
                 SDPair       const & Ini,    ///< Initial values
                 Array<Node*> const & Nodes); ///< Array with all nodes (used to set the connectivity)

    // Methods
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs);   ///< If setting body forces, IdxEdgeOrFace is ignored
    void CalcFint    (Vec_t * F_int=NULL)                   const; ///< Calculate or set Fint. Set nodes if F_int==NULL
    void CalcK       (Mat_t & K)                            const; ///< Stiffness matrix
    void CalcM       (Mat_t & M)                            const; ///< Mass matrix
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL) const; ///< Update state at IPs
    void GetState    (SDPair & KeysVals, int IdxIP=-1)      const; ///< IdxIP<0 => At centroid
    void GetState    (Array<SDPair> & Results)              const; ///< At each integration point (IP)

    // Internal methods
    void CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const; ///< Strain-displacement matrix. Coef: coefficient used during integration
    void CalcN (Mat_t const & C, IntegPoint const & IP, Mat_t & N, double & detJ, double & Coef) const; ///< Shape functions matrix

    // Constants
    double h;   ///< Thickness of the element
    double rho; ///< Density
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline EquilibElem::EquilibElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes)
{
    // check GE
    if (GE==NULL) throw new Fatal("EquilibElem::EquilibElem: GE (geometry element) must be defined");

    // parameters/properties
    h   = (Prp.HasKey("h")   ? Prp("h")   : 1.0);
    rho = (Prp.HasKey("rho") ? Prp("rho") : 1.0);

    // allocate and initialize state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Sta.Push (new EquilibState(NDim));
        Mdl->InitIvs (Ini, Sta[i]);
    }

    // allocate state at centroid
    Sta.Push (new EquilibState(NDim));
    Mdl->InitIvs (Ini, Sta[Sta.Size()-1]);

    if (NDim==2)
    {
        // set UKeys in parent element
        UKeys.Resize (NDim);
        UKeys = "ux", "uy";

        // initialize DOFs
        for (size_t i=0; i<GE->NN; ++i) Con[i]->AddDOF("ux uy", "fx fy");
    }
    else // 3D
    {
        // set UKeys in parent element
        UKeys.Resize (NDim);
        UKeys = "ux", "uy", "uz";

        // initialize DOFs
        for (size_t i=0; i<GE->NN; ++i) Con[i]->AddDOF("ux uy uz", "fx fy fz");
    }

    // set F in nodes due to Fint
    CalcFint ();
}

inline void EquilibElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs)
{
    bool has_bx  = BCs.HasKey("bx");
    bool has_by  = BCs.HasKey("by");
    bool has_bz  = BCs.HasKey("bz");
    bool has_cbx = BCs.HasKey("cbx"); // cbx: centrifugal body force along x (in axisymmetric problems)
    if (has_bx || has_by || has_bz || has_cbx) // prescribed body forces
    {
        double bx = (has_bx  ? BCs("bx")  : 0.0);
        double by = (has_by  ? BCs("by")  : 0.0);
        double bz = (has_bz  ? BCs("bz")  : 0.0);
               bx = (has_cbx ? BCs("cbx") : bx );
        Mat_t C;
        CoordMatrix (C);
        //Vec_t Fb(GE->NN*NDim); set_to_zero(Fb);
        for (size_t i=0; i<GE->NIP; ++i)
        {
            // geometric data
            GE->Shape  (GE->IPs[i].r, GE->IPs[i].s, GE->IPs[i].t);
            GE->Derivs (GE->IPs[i].r, GE->IPs[i].s, GE->IPs[i].t);

            // Jacobian and its determinant
            Mat_t J(GE->dNdR * C); // J = dNdR * C
            double detJ = Det(J);

            // coefficient used during integration
            double coef = h*detJ*GE->IPs[i].w;
            if (GTy==axs_t)
            {
                // calculate radius=x at this IP
                double radius = 0.0;
                for (size_t j=0; j<GE->NN; ++j) radius += GE->N(j)*Con[j]->Vert.C[0];

                // correct coef
                if (has_cbx) coef *= radius*radius;
                else         coef *= radius;
            }

            // add to dF
            for (size_t j=0; j<GE->NN; ++j)
            {
                Con[j]->DF[Con[j]->FMap("fx")] += coef*GE->N(j)*bx;
                Con[j]->DF[Con[j]->FMap("fy")] += coef*GE->N(j)*by;  if (NDim==3)
                Con[j]->DF[Con[j]->FMap("fz")] += coef*GE->N(j)*bz;
                //Fb(0+j*NDim) += coef*GE->N(j)*bx;
                //Fb(1+j*NDim) += coef*GE->N(j)*by;  if (NDim==3)
                //Fb(2+j*NDim) += coef*GE->N(j)*bz;
            }
        }
        //std::cout << bx << "   Fb(" << Cell.ID << ") = " << PrintVector(Fb);
    }
    else
    {
        bool has_ux = BCs.HasKey("ux");
        bool has_uy = BCs.HasKey("uy");
        bool has_uz = BCs.HasKey("uz");
        if (has_ux || has_uy || has_uz) // prescribed displacements
        {
            for (size_t i=0; i<GE->NFN; ++i)
            {
                Node & nod = (*Con[GE->FNode(IdxEdgeOrFace,i)]);
                nod.SetBCs (BCs);
            }
        }
        else // surface loading
        {
            // matrix of coordinates of edge/face
            Mat_t C;
            FCoordMatrix (IdxEdgeOrFace, C);

            // loading
            double qx = (BCs.HasKey("qx") ? BCs("qx") : 0.0);
            double qy = (BCs.HasKey("qy") ? BCs("qy") : 0.0);
            double qz = (BCs.HasKey("qz") ? BCs("qz") : 0.0);
            double qn = (BCs.HasKey("qn") ? BCs("qn") : 0.0);
            double qt = (BCs.HasKey("qt") ? BCs("qt") : 0.0);

            for (size_t i=0; i<GE->NFIP; ++i)
            {
                // geometric data
                GE->FaceShape  (GE->FIPs[i].r, GE->FIPs[i].s);
                GE->FaceDerivs (GE->FIPs[i].r, GE->FIPs[i].s);

                // face/edge Jacobian and its determinant
                Mat_t J(GE->FdNdR * C);

                // normal to edge/face
                Vec_t n(NDim); // normal multiplied by detJ
                if (NDim==2) n = J(0,1), -J(0,0);
                else
                {
                    // vectorial product
                    Vec_t a(3);  a = J(0,0), J(0,1), J(0,2);
                    Vec_t b(3);  b = J(1,0), J(1,1), J(1,2);
                    n = a(1)*b(2) - a(2)*b(1),
                        a(2)*b(0) - a(0)*b(2),
                        a(0)*b(1) - a(1)*b(0);
                }

                // loading
                if (NDim==2)
                {
                    qx = n(0)*qn - n(1)*qt;
                    qy = n(1)*qn + n(0)*qt;
                }
                else
                {
                    qx = n(0)*qn;
                    qy = n(1)*qn;
                    qz = n(2)*qn;
                }

                // coefficient used during integration
                double coef = h*GE->FIPs[i].w; // *detJ is not neccessary since qx,qy,qz are already multiplied by detJ (due to normal)
                if (GTy==axs_t)
                {
                    // calculate radius=x at this FIP
                    double radius = 0.0;
                    for (size_t j=0; j<GE->NFN; ++j) radius += GE->FN(j)*Con[GE->FNode(IdxEdgeOrFace,j)]->Vert.C[0];
                    coef *= radius; // correct coef
                }

                // add to dF
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    Node & nod = (*Con[GE->FNode(IdxEdgeOrFace,j)]);
                    nod.DF[nod.FMap("fx")] += coef*GE->FN(j)*qx;
                    nod.DF[nod.FMap("fy")] += coef*GE->FN(j)*qy;  if (NDim==3)
                    nod.DF[nod.FMap("fz")] += coef*GE->FN(j)*qz;
                }
            }
        }
    }
}

inline void EquilibElem::CalcFint (Vec_t * F_int) const
{
    Array<size_t> loc;
    GetLoc (loc);
    double detJ, coef;
    Mat_t  C, B;
    Vec_t  Fe(GE->NN*NDim);
    set_to_zero (Fe);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        CalcB (C, GE->IPs[i], B, detJ, coef);
        Vec_t Btsig(trans(B)*static_cast<EquilibState const *>(Sta[i])->Sig);
        Fe += coef * (Btsig);
    }
    if (F_int==NULL) // set nodes
    {
        // add to F
        for (size_t i=0; i<GE->NN; ++i)
        {
            Con[i]->F[Con[i]->FMap("fx")] += Fe(0+i*NDim);
            Con[i]->F[Con[i]->FMap("fy")] += Fe(1+i*NDim);  if (NDim==3)
            Con[i]->F[Con[i]->FMap("fz")] += Fe(2+i*NDim);
        }
    }
    else for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += Fe(i);
}

inline void EquilibElem::CalcK (Mat_t & K) const
{
    double detJ, coef;
    Mat_t C, D, B;
    int nrows = GE->NN*NDim; // number of rows in local K matrix
    K.change_dim (nrows,nrows);
    set_to_zero  (K);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Mdl->Stiffness (Sta[i], D);
        CalcB (C, GE->IPs[i], B, detJ, coef);
        Mat_t BtDB(trans(B)*D*B);
        K += coef * (BtDB);
    }
    //std::cout << "D = \n" << PrintMatrix(D);
    //std::cout << "K = \n" << PrintMatrix(K);
}

inline void EquilibElem::CalcM (Mat_t & M) const
{
    double detJ, coef;
    Mat_t  C, N;
    int nrows = GE->NN*NDim; // number of rows in local M matrix
    M.change_dim (nrows,nrows);
    set_to_zero  (M);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        CalcN (C, GE->IPs[i], N, detJ, coef);
        Mat_t NtN(trans(N)*N);
        M += (rho*coef) * (NtN);
    }
    //std::cout << "M = \n" << PrintMatrix(M);
}

inline void EquilibElem::CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const
{
    // deriv of shape func w.r.t natural coordinates
    GE->Derivs (IP.r, IP.s, IP.t);

    // Jacobian and its determinant
    Mat_t J(GE->dNdR * C); // J = dNdR * C
    detJ = Det(J);

    //std::cout << "J = \n" << PrintMatrix(J);
    //std::cout << "detJ = " << detJ << "\n";

    // deriv of shape func w.r.t real coordinates
    Mat_t Ji;
    Inv (J, Ji);
    Mat_t dNdX(Ji * GE->dNdR); // dNdX = Inv(J) * dNdR

    // coefficient used during integration
    Coef = h*detJ*IP.w;

    // B matrix
    int nrows = GE->NN*NDim; // number of rows in local K matrix
    int ncomp = 2*NDim;          // number of stress/strain components
    B.change_dim (ncomp,nrows);
    set_to_zero  (B);
    if (NDim==2)
    {
        if (GTy==axs_t)
        {
            // shape functions
            GE->Shape (IP.r, IP.s, IP.t);

            // correct Coef
            double radius = 0.0; // radius=x at this IP
            for (size_t j=0; j<GE->NN; ++j) radius += GE->N(j)*Con[j]->Vert.C[0];
            Coef *= radius;

            // B matrix
            for (size_t i=0; i<GE->NN; ++i)
            {
                B(0,0+i*NDim) = dNdX(0,i);
                B(1,1+i*NDim) = dNdX(1,i);
                B(2,0+i*NDim) = GE->N(i)/radius;
                B(3,0+i*NDim) = dNdX(1,i)/sqrt(2.0);
                B(3,1+i*NDim) = dNdX(0,i)/sqrt(2.0);
            }
            //std::cout << "B = \n" << PrintMatrix(B);
        }
        else // pse_t, psa_t
        {
            for (size_t i=0; i<GE->NN; ++i)
            {
                B(0,0+i*NDim) = dNdX(0,i);
                B(1,1+i*NDim) = dNdX(1,i);
                B(3,0+i*NDim) = dNdX(1,i)/sqrt(2.0);
                B(3,1+i*NDim) = dNdX(0,i)/sqrt(2.0);
            }
        }
    }
    else // 3D
    {
        for (size_t i=0; i<GE->NN; ++i)
        {
            B(0,0+i*NDim) = dNdX(0,i);
            B(1,1+i*NDim) = dNdX(1,i);
            B(2,2+i*NDim) = dNdX(2,i);
            B(3,0+i*NDim) = dNdX(1,i)/sqrt(2.0);   B(3,1+i*NDim) = dNdX(0,i)/sqrt(2.0);
            B(4,1+i*NDim) = dNdX(2,i)/sqrt(2.0);   B(4,2+i*NDim) = dNdX(1,i)/sqrt(2.0);
            B(5,2+i*NDim) = dNdX(0,i)/sqrt(2.0);   B(5,0+i*NDim) = dNdX(2,i)/sqrt(2.0);
        }
    }
}

inline void EquilibElem::CalcN (Mat_t const & C, IntegPoint const & IP, Mat_t & N, double & detJ, double & Coef) const
{
    // deriv of shape func w.r.t natural coordinates
    GE->Shape  (IP.r, IP.s, IP.t);
    GE->Derivs (IP.r, IP.s, IP.t);

    // Jacobian and its determinant
    Mat_t J(GE->dNdR * C); // J = dNdR * C
    detJ = Det(J);

    // coefficient used during integration
    Coef = h*detJ*IP.w;

    // N matrix
    int nrows = GE->NN*NDim; // number of rows in local K matrix
    N.change_dim (NDim,nrows);
    set_to_zero  (N);
    if (GTy==axs_t)
    {
        double radius = 0.0; // radius=x at this IP
        for (size_t j=0; j<GE->NN; ++j) radius += GE->N(j)*Con[j]->Vert.C[0];
        Coef *= radius; // correct coef for axisymmetric problems
    }
    for (int    i=0; i<NDim;   ++i)
    for (size_t j=0; j<GE->NN; ++j)
        N(i,i+j*NDim) = GE->N(j);
}

inline void EquilibElem::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal displacements
    int nrows = GE->NN*NDim; // number of rows in local K matrix
    Vec_t dUe(nrows);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // update state at each IP
    StressUpdate su(Mdl);
    int ncomp = 2*NDim; // number of stress/strain components
    double detJ, coef;
    Mat_t  C, B;
    Vec_t  dFe(nrows), dsig(ncomp), deps(ncomp);
    set_to_zero (dFe);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // B matrix
        CalcB (C, GE->IPs[i], B, detJ, coef);

        // strain and stress increments
        deps = B * dUe;
        su.Update (deps, Sta[i], dsig);

        // element nodal forces
        Vec_t Btdsig(trans(B)*dsig);
        dFe += coef * (Btdsig);
    }

    // update state at centroid
    CalcB (C, GE->Rct, B, detJ, coef);
    deps = B * dUe;
    su.Update (deps, Sta[Sta.Size()-1], dsig);

    // add results to Fint (internal forces)
    if (F_int!=NULL) for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
}

inline void EquilibElem::GetState (SDPair & KeysVals, int IdxIP) const
{
    double sq2 = sqrt(2.0);
    if (false) // extrapolate, then interpolate at centroid
    {
        // shape func matrix
        Mat_t M, Mi;
        ShapeMatrix (M);
        Inv (M, Mi);

        // shape at centroid
        GE->Shape (GE->Rct.r, GE->Rct.s, GE->Rct.t);

        // values at centroid
        int ncomp = 2*NDim; // number of stress/strain components
        Vec_t sig_ct(ncomp), eps_ct(ncomp);
        set_to_zero (sig_ct);
        set_to_zero (eps_ct);
        for (int i=0; i<ncomp; ++i)
        {
            // extrapolate to nodes
            Vec_t sig_at_IPs (GE->NIP);
            Vec_t sig_at_Nods(GE->NN);
            Vec_t eps_at_IPs (GE->NIP);
            Vec_t eps_at_Nods(GE->NN);
            for (size_t j=0; j<GE->NIP; ++j)
            {
                sig_at_IPs(j) = static_cast<EquilibState const *>(Sta[j])->Sig(i);
                eps_at_IPs(j) = static_cast<EquilibState const *>(Sta[j])->Eps(i);
            }
            sig_at_Nods = Mi * sig_at_IPs;
            eps_at_Nods = Mi * eps_at_IPs;

            // interpolate to centroid
            for (size_t j=0; j<GE->NN; ++j)
            {
                sig_ct(i) += GE->N(j)*sig_at_Nods(j);
                eps_ct(i) += GE->N(j)*eps_at_Nods(j);
            }
        }
        KeysVals.Set("sx sy sz sxy  ex ey ez exy",//  p q ev ed",
                     sig_ct(0),sig_ct(1),sig_ct(2),sig_ct(3)/sq2,
                     eps_ct(0),eps_ct(1),eps_ct(2),eps_ct(3)/sq2);
    }
    else if (false) // average
    {
        // average values
        int ncomp = 2*NDim; // number of stress/strain components
        Vec_t sig_ave(ncomp), eps_ave(ncomp);
        set_to_zero (sig_ave);
        set_to_zero (eps_ave);
        for (size_t i=0; i<Sta.Size(); ++i)
        {
            EquilibState const * sta = static_cast<EquilibState const *>(Sta[i]);
            sig_ave += sta->Sig;
            eps_ave += sta->Eps;
        }
        sig_ave /= Sta.Size();
        eps_ave /= Sta.Size();
        KeysVals.Set("sx sy sz sxy  ex ey ez exy",//  p q ev ed",
                     sig_ave(0),sig_ave(1),sig_ave(2),sig_ave(3)/sq2,
                     eps_ave(0),eps_ave(1),eps_ave(2),eps_ave(3)/sq2);
    }
    else
    {
        Vec_t const * sig;
        Vec_t const * eps;
        if (IdxIP<0) // centroid
        {
            sig = &(static_cast<EquilibState const *>(Sta[Sta.Size()-1])->Sig);
            eps = &(static_cast<EquilibState const *>(Sta[Sta.Size()-1])->Eps);
        }
        else
        {
            sig = &(static_cast<EquilibState const *>(Sta[IdxIP])->Sig);
            eps = &(static_cast<EquilibState const *>(Sta[IdxIP])->Eps);
        }
        if (NDim==2)
        {
            KeysVals.Set("sx sy sz sxy  ex ey ez exy",
                         (*sig)(0),(*sig)(1),(*sig)(2),(*sig)(3)/sq2,
                         (*eps)(0),(*eps)(1),(*eps)(2),(*eps)(3)/sq2);
        }
        else
        {
            KeysVals.Set("sx sy sz sxy syz szx  ex ey ez exy eyz ezx",
                         (*sig)(0),(*sig)(1),(*sig)(2),(*sig)(3)/sq2,(*sig)(4)/sq2,(*sig)(5)/sq2,
                         (*eps)(0),(*eps)(1),(*eps)(2),(*eps)(3)/sq2,(*eps)(4)/sq2,(*eps)(5)/sq2);
        }
    }
}

inline void EquilibElem::GetState (Array<SDPair> & Results) const
{
    double sq2 = sqrt(2.0);
    Results.Resize (GE->NIP);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Vec_t const & sig = static_cast<EquilibState const *>(Sta[i])->Sig;
        Vec_t const & eps = static_cast<EquilibState const *>(Sta[i])->Eps;
        if (NDim==2)
        {
            Results[i].Set("sx sy sz sxy  ex ey ez exy",
                           sig(0),sig(1),sig(2),sig(3)/sq2,
                           eps(0),eps(1),eps(2),eps(3)/sq2);
        }
        else
        {
            Results[i].Set("sx sy sz sxy syz szx  ex ey ez exy eyz ezx",
                           sig(0),sig(1),sig(2),sig(3)/sq2,sig(4)/sq2,sig(5)/sq2,
                           eps(0),eps(1),eps(2),eps(3)/sq2,eps(4)/sq2,eps(5)/sq2);
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * EquilibElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new EquilibElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int EquilibElemRegister()
{
    ElementFactory["Equilib"] = EquilibElemMaker;
    PROB.Set ("Equilib", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __EquilibElem_dummy_int  = EquilibElemRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIBELEM
