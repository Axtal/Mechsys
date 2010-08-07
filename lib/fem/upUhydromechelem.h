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

#ifndef MECHSYS_FEM_HYDROMECHELEM_H
#define MECHSYS_FEM_HYDROMECHELEM_H

// MechSys
#include <mechsys/fem/element.h>
#include <mechsys/models/equilibstate.h>
#include <mechsys/models/stressupdate.h>

namespace FEM
{

class HydroMechElem : public Element
{
public:
    // static
    static size_t NDn; ///< Number of DOFs per node = NDim*2+1
    static size_t NDt; ///< Total number of DOFs = GE->NN*NDn
    static size_t NDs; ///< Total number of DOFs of solid matrix = GE->NN*NDim
    static size_t NDp; ///< Total number of DOFs of pressure = GE->NN
    static size_t NCo; ///< Number of stress/strain components = 2*NDim
    static Mat_t  Im;  ///< Identity tensor (NCo,1)
    static Vec_t  Iv;  ///< Identity tensor (NCo)

    // Constructor
    HydroMechElem (int                  NDim,   ///< Space dimension
                   Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                   Model        const * Mdl,    ///< Model
                   SDPair       const & Prp,    ///< Properties
                   SDPair       const & Ini,    ///< Initial values
                   Array<Node*> const & Nodes); ///< Connectivity

    // Methods
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs,
                      NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM); ///< If setting body forces, IdxEdgeOrFace is ignored
    void Matrices    (Mat_t & M, Mat_t & C, Mat_t & K)      const;
    void CalcFint    (Vec_t * F_int=NULL)                   const; ///< Calculate or set Fint. Set nodes if F_int==NULL
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL) const; ///< Update state at IPs
    void GetState    (SDPair & KeysVals, int IdxIP=-1)      const; ///< IdxIP<0 => At centroid
    void GetState    (Array<SDPair> & Results)              const; ///< At each integration point (IP)

    // Internal methods
    void Interp (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & Bp, Mat_t & N, Mat_t & Np, double & detJ, double & Coef) const; ///< Interpolation matrices

    // Constants
    double h;    ///< Thickness of the element
    double alp;  ///< Compressibility coefficient
    double nn;   ///< Porosity
    double rhoS; ///< Density of solids
    double rhoF; ///< Density of fluid
    double kk;   ///< Permeability
    double Qs;   ///< Storage due to compressibility
    double rho;  ///< Density of mixture
};

size_t HydroMechElem::NDn = 0;
size_t HydroMechElem::NDt = 0;
size_t HydroMechElem::NDs = 0;
size_t HydroMechElem::NDp = 0;
size_t HydroMechElem::NCo = 0;
Mat_t  HydroMechElem::Im;
Vec_t  HydroMechElem::Iv;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline HydroMechElem::HydroMechElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes)
{
    // check GE
    if (GE==NULL) throw new Fatal("HydroMechElem::HydroMechElem: GE (geometry element) must be defined");

    // parameters/properties
    h    = (Prp.HasKey("h") ? Prp("h") : 1.0);
    alp  = Prp("alp");
    nn   = Prp("n");
    rhoS = Prp("rhoS");
    rhoF = Prp("rhoF");
    kk   = Prp("k");
    Qs   = Prp("Q");
    rho  = nn*rhoF+(1.0-nn)*rhoS;

    // allocate and initialize state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Sta.Push (new EquilibState(NDim));
        Mdl->InitIvs (Ini, Sta[i]);
    }

    // set constants of this class (once)
    if (NDn==0)
    {
        NDn = NDim*2+1;
        NDt = GE->NN*NDn;
        NDs = GE->NN*NDim;
        NDp = GE->NN;
        NCo = 2*NDim;
        Im.change_dim (NCo,1);
        Iv.change_dim (NCo);
        set_to_zero (Im);
        set_to_zero (Iv);
        Im(0,0)=1.0;   Im(1,0)=1.0;   Im(2,0)=1.0;
        Iv(0)  =1.0;   Iv(1)  =1.0;   Iv(2)  =1.0;
    }

    // set UKeys in parent element and initialize DOFs
    UKeys.Resize (NDn);
    if (NDim==2)
    {
        UKeys = "ux", "uy",  "pw",  "Ux", "Uy";
        for (size_t i=0; i<GE->NN; ++i) Con[i]->AddDOF("ux uy  pw  Ux Uy", "fx fy  Qw  Fx Fy");
    }
    else // 3D
    {
        UKeys = "ux", "uy", "uz",  "pw",  "Ux", "Uy", "Uz";
        for (size_t i=0; i<GE->NN; ++i) Con[i]->AddDOF("ux uy uz  pw  Ux Uy Uz", "fx fy fz  Qw  Fx Fy Fz");
    }

    // set F in nodes due to Fint
    CalcFint ();
}

inline void HydroMechElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM)
{
    bool has_bx   = BCs.HasKey("bx");   // x component of body force
    bool has_by   = BCs.HasKey("by");   // y component of body force
    bool has_bz   = BCs.HasKey("bz");   // z component of body force
    bool has_cbx  = BCs.HasKey("cbx");  // centrifugal body force along x (in axisymmetric problems)
    bool has_qx   = BCs.HasKey("qx");   // x component of distributed loading
    bool has_qy   = BCs.HasKey("qy");   // y component of distributed loading
    bool has_qz   = BCs.HasKey("qz");   // z component of distributed loading
    bool has_qn   = BCs.HasKey("qn");   // normal distributed loading
    bool has_qt   = BCs.HasKey("qt");   // tangential distributed loading (2D only)
    bool has_ux   = BCs.HasKey("ux");   // x displacement
    bool has_uy   = BCs.HasKey("uy");   // y displacement
    bool has_uz   = BCs.HasKey("uz");   // z displacement
    bool has_pw   = BCs.HasKey("pw");   // water pressure
    bool has_Ux   = BCs.HasKey("Ux");   // x displacement of fluid
    bool has_Uy   = BCs.HasKey("Uy");   // y displacement of fluid
    bool has_Uz   = BCs.HasKey("Uz");   // z displacement of fluid
    bool has_flux = BCs.HasKey("flux"); // normal flux to boundary

    // body forces
    if (has_bx || has_by || has_bz || has_cbx) // prescribed body forces
    {
        // matrix of coordinates of nodes
        Mat_t C;
        CoordMatrix (C);
        
        // loading
        double bx = (has_bx  ? BCs("bx")  : 0.0);
        double by = (has_by  ? BCs("by")  : 0.0);
        double bz = (has_bz  ? BCs("bz")  : 0.0);
               bx = (has_cbx ? BCs("cbx") : bx );

        // set
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
                if (has_bx || has_cbx) pF[Con[j]].first[Con[j]->FMap("fx")] += coef*GE->N(j)*bx;
                if (has_by           ) pF[Con[j]].first[Con[j]->FMap("fy")] += coef*GE->N(j)*by;
                if (has_bz           ) pF[Con[j]].first[Con[j]->FMap("fz")] += coef*GE->N(j)*bz;
            }
        }

        // set CalcM
        for (size_t j=0; j<GE->NN; ++j) pF[Con[j]].second = CalcM;
    }

    // surface loading
    if (has_qx || has_qy || has_qz || has_qn  || has_qt)
    {
        // matrix of coordinates of edge/face
        Mat_t Cf;
        FCoordMatrix (IdxEdgeOrFace, Cf);

        // loading
        double qx = (has_qx ? BCs("qx") : 0.0);
        double qy = (has_qy ? BCs("qy") : 0.0);
        double qz = (has_qz ? BCs("qz") : 0.0);
        double qn = (has_qn ? BCs("qn") : 0.0);
        double qt = (has_qt ? BCs("qt") : 0.0);

        // set
        for (size_t i=0; i<GE->NFIP; ++i)
        {
            // geometric data
            GE->FaceShape  (GE->FIPs[i].r, GE->FIPs[i].s);
            GE->FaceDerivs (GE->FIPs[i].r, GE->FIPs[i].s);

            // face/edge Jacobian and its determinant
            Mat_t J(GE->FdNdR * Cf);

            // coefficient used during integration
            double coef = h*GE->FIPs[i].w; // *detJ is not neccessary since qx,qy,qz are already multiplied by detJ (due to normal)
            if (GTy==axs_t)
            {
                // calculate radius=x at this FIP
                double radius = 0.0;
                for (size_t j=0; j<GE->NFN; ++j) radius += GE->FN(j)*Con[GE->FNode(IdxEdgeOrFace,j)]->Vert.C[0];
                coef *= radius; // correct coef
            }

            // calculate qx, qy and qz from qn and qt
            if (has_qn || has_qt)
            {
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
            }

            // add to dF
            for (size_t j=0; j<GE->NFN; ++j)
            {
                size_t k = GE->FNode(IdxEdgeOrFace,j);
                pF[Con[k]].first[Con[k]->FMap("fx")] += coef*GE->FN(j)*qx;
                pF[Con[k]].first[Con[k]->FMap("fy")] += coef*GE->FN(j)*qy;  if (NDim==3)
                pF[Con[k]].first[Con[k]->FMap("fz")] += coef*GE->FN(j)*qz;
            }
        }

        // set CalcM
        for (size_t j=0; j<GE->NFN; ++j) pF[Con[GE->FNode(IdxEdgeOrFace,j)]].second = CalcM;
    }

    // flux
    if (has_flux)
    {
        double Qn = BCs("flux");
        Mat_t Cf;
        FCoordMatrix (IdxEdgeOrFace, Cf);
        double detJ, coef;
        for (size_t i=0; i<GE->NFIP; ++i)
        {
            CalcFaceShape (Cf, GE->FIPs[i], detJ, coef);
            if (GTy==axs_t) // correct Coef for axisymmetric problems
            {
                double radius = 0.0; // calculate radius=x at this FIP
                for (size_t j=0; j<GE->NFN; ++j) radius += GE->FN(j)*Con[GE->FNode(IdxEdgeOrFace,j)]->Vert.C[0];
                coef *= radius; // correct coef
            }
            for (size_t j=0; j<GE->NFN; ++j)
            {
                size_t k = GE->FNode(IdxEdgeOrFace,j);
                pF[Con[k]].first[Con[k]->FMap("Qw")] += coef*GE->FN(j)*Qn;
            }
        }

        // set CalcM
        for (size_t j=0; j<GE->NFN; ++j) pF[Con[GE->FNode(IdxEdgeOrFace,j)]].second = CalcM;
    }

    // prescribed displacements
    if (has_ux || has_uy || has_uz || has_pw ||
        has_Ux || has_Uy || has_Uz)
    {
        double ux = (has_ux ? BCs("ux") : 0.0);
        double uy = (has_uy ? BCs("uy") : 0.0);
        double uz = (has_uz ? BCs("uz") : 0.0);
        double pw = (has_pw ? BCs("pw") : 0.0);
        double Ux = (has_Ux ? BCs("Ux") : 0.0);
        double Uy = (has_Uy ? BCs("Uy") : 0.0);
        double Uz = (has_Uz ? BCs("Uz") : 0.0);
        for (size_t j=0; j<GE->NFN; ++j)
        {
            size_t k = GE->FNode(IdxEdgeOrFace,j);
            if (has_ux) pU[Con[k]].first[Con[k]->UMap("ux")] = ux;
            if (has_uy) pU[Con[k]].first[Con[k]->UMap("uy")] = uy;
            if (has_uz) pU[Con[k]].first[Con[k]->UMap("uz")] = uz;
            if (has_pw) pU[Con[k]].first[Con[k]->UMap("pw")] = pw;
            if (has_Ux) pU[Con[k]].first[Con[k]->UMap("Ux")] = Ux;
            if (has_Uy) pU[Con[k]].first[Con[k]->UMap("Uy")] = Uy;
            if (has_Uz) pU[Con[k]].first[Con[k]->UMap("Uz")] = Uz;
        }
    }
}

inline void HydroMechElem::Matrices (Mat_t & MM, Mat_t & CC, Mat_t & KK) const
{
    // permeability tensor
    /*
    Mat_t km(NDim,NDim);
    set_to_zero (km);
    km(0,0)=kk;  km(1,1)=kk;
    */

    // submatrices
    Mat_t Ms(NDs,NDs); // mass matrix (solids)
    Mat_t Mf(NDs,NDs); // mass matrix (fluid)
    Mat_t C1(NDs,NDs); // permeability matrix ?
    Mat_t K (NDs,NDs); // stiffness matrix
    Mat_t G1(NDs,NDp); // coupling matrix (solids)
    Mat_t G2(NDs,NDp); // coupling matrix (fluid)
    Mat_t S (NDp,NDp); // compressibility matrix
    set_to_zero (Ms);
    set_to_zero (Mf);
    set_to_zero (C1);
    set_to_zero (K);
    set_to_zero (G1);
    set_to_zero (G2);
    set_to_zero (S);

    // auxiliar matrices
    double detJ, coef;
    Mat_t C, D, B, Bp, N, Np;
    Mat_t NtN   (NDs,NDs);
    Mat_t BtDB  (NDs,NDs);
    Mat_t BtINp (NDs,NDp);
    Mat_t NptNp (NDp,NDp);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Mdl->Stiffness (Sta[i], D);
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);
        NtN    = trans(N)*N;
        BtDB   = trans(B)*D*B;
        BtINp  = trans(B)*Im*Np;
        NptNp  = trans(Np)*Np;
        Ms    += (coef*(1.0-nn)*rhoS) * NtN;
        Mf    += (coef*nn*rhoF)       * NtN;
        C1    += (coef*nn*nn/kk)      * NtN;
        K     += (coef)               * BtDB;
        G1    += (coef*(alp-nn))      * BtINp;
        G2    += (coef*nn)            * BtINp;
        S     += (coef/Qs)            * NptNp;
    }

    // assemble
    MM.change_dim (NDt,NDt);
    CC.change_dim (NDt,NDt);
    KK.change_dim (NDt,NDt);
    set_to_zero (MM);
    set_to_zero (CC);
    set_to_zero (KK);
    for (size_t i=0; i<NDs; ++i)
    {
        for (size_t j=0; j<NDs; ++j)
        {
            MM(        i,         j) =  Ms(i,j);
            MM(NDs+NDp+i, NDs+NDp+j) =  Mf(i,j);
            CC(        i,         j) =  C1(i,j);
            CC(        i, NDs+NDp+j) = -C1(i,j);
            CC(NDs+NDp+i,         j) = -C1(j,i);
            CC(NDs+NDp+i, NDs+NDp+j) =  C1(i,j);
            KK(        i,         j) =  K (i,j);
        }
        for (size_t j=0; j<NDp; ++j)
        {
            KK(        i, NDs+j) = -G1(i,j);
            KK(NDs+NDp+i, NDs+j) = -G2(i,j);
        }
    }
    for (size_t i=0; i<NDp; ++i)
    {
        for (size_t j=0; j<NDs; ++j)
        {
            KK(NDs+i,         j) = -G1(j,i);
            KK(NDs+i, NDs+NDp+j) = -G2(j,i);
        }
        for (size_t j=0; j<NDp; ++j)
        {
            KK(NDs+i, NDs+j) = S(i,j);
        }
    }
}

inline void HydroMechElem::Interp (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & Bp, Mat_t & N, Mat_t & Np, double & detJ, double & Coef) const
{
    // deriv of shape func w.r.t natural coordinates
    GE->Shape  (IP.r, IP.s, IP.t);
    GE->Derivs (IP.r, IP.s, IP.t);

    // Jacobian and its determinant
    Mat_t J(GE->dNdR * C); // J = dNdR * C
    detJ = Det(J);

    // deriv of shape func w.r.t real coordinates
    Mat_t Ji;
    Inv (J, Ji);
    Mat_t dNdX(Ji * GE->dNdR); // dNdX = Inv(J) * dNdR

    // coefficient used during integration
    Coef = h*detJ*IP.w;

    // B matrix
    B.change_dim (NCo,NDs);
    set_to_zero  (B);
    if (NDim==2)
    {
        if (GTy==axs_t)
        {
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

    // Bp matrix
    Bp.change_dim (NDim,NDp);
    Bp = Ji * GE->dNdR; // Bp = dNdX = Inv(J) * dNdR

    // N matrix
    N.change_dim (NDim,NDs);
    set_to_zero  (N);
    for (int    i=0; i<NDim;   ++i)
    for (size_t j=0; j<GE->NN; ++j)
        N(i,i+j*NDim) = GE->N(j);

    // Np matrix
    Np.change_dim (1,NDp);
    for (size_t j=0; j<GE->NN; ++j) Np(0,j) = GE->N(j);
}

inline void HydroMechElem::CalcFint (Vec_t * F_int) const
{
    return;
    // calc element forces
    double detJ, coef;
    Mat_t C, B, Bp, N, Np;
    Vec_t BtSig (NDs);
    Vec_t fs    (NDs);
    Vec_t fp    (NDp);
    set_to_zero (fs);
    set_to_zero (fp);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // interpolation matrices
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);

        // force due to stress in solids
        BtSig = trans(B)*static_cast<EquilibState const *>(Sta[i])->Sig;
        fs   += coef * BtSig;
    }

    // assemble Fe vector (element force)
    Vec_t Fe(NDt);
    for (size_t i=0; i<NDs; ++i) Fe(i)     = fs(i);
    for (size_t i=0; i<NDp; ++i) Fe(NDs+i) = fp(i);

    // set nodes directy
    if (F_int==NULL)
    {
        // add to F
        for (size_t i=0; i<GE->NN; ++i)
        {
            size_t k = 0;
            Con[i]->F[Con[i]->FMap("fx")] += Fe(k+i*NDn);  k++;
            Con[i]->F[Con[i]->FMap("fy")] += Fe(k+i*NDn);  k++;  if (NDim==3) {
            Con[i]->F[Con[i]->FMap("fz")] += Fe(k+i*NDn);  k++; }
            Con[i]->F[Con[i]->FMap("Qw")] += Fe(k+i*NDn);  k++;
        }
    }

    // add to F_int
    else
    {
        Array<size_t> loc;
        GetLoc (loc);
        for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += Fe(i);
    }
}

inline void HydroMechElem::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    /*
    // permeability tensor
    Mat_t km(NDim,NDim);
    set_to_zero (km);
    km(0,0)=kk;  km(1,1)=kk;
    */

    // element nodal displacements
    Array<size_t> loc;
    GetLoc (loc);
    Vec_t dUe(NDt);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // split dUe
    Vec_t dus(NDs); // solid
    Vec_t dup(NDp); // pressure
    Vec_t duf(NDs); // fluid
    for (size_t i=0; i<NDs; ++i)
    {
        dus(i) = dUe(i);
        duf(i) = dUe(NDs+NDp+i);
    }
    for (size_t i=0; i<NDp; ++i) dup(i) = dUe(NDs+i);

    // update state at each IP
    StressUpdate su(Mdl);
    double detJ, coef;
    Mat_t C, B, Bp, N, Np;
    Vec_t dfs(NDs), dsig(NCo),  deps(NCo);
    Vec_t dfp(NDp);
    Vec_t dff(NDs);
    Vec_t BtdSig(NDs);
    set_to_zero (dfs);
    set_to_zero (dfp);
    set_to_zero (dff);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // interpolation matrices
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);

        // strain and (effective) stress increments
        deps = B * dus;
        su.Update (deps, Sta[i], dsig);

        // pressure increment
        Vec_t tmp1(Np*dup);
        double dp = tmp1(0);

        // element nodal forces (solid)
        dsig  -= ((alp-nn)*dp)*Iv;
        BtdSig = trans(B)*dsig;
        dfs   += coef * BtdSig;

        // element nodal forces (pressure)
        Vec_t dE(B * duf);
        for (size_t i=0; i<NDp; ++i) dfp(i) += coef*Np(0,i)*(dp/Qs-(alp-nn)*dot(deps,Iv)-nn*dot(dE,Iv));

        // element nodal forces (fluid)
        Vec_t tmp(trans(B)*Iv);
        dff -= (coef*nn) * tmp;
    }

    //std::cout << "dfp = " << dfp << std::endl;
    //throw new Fatal("stop");

    if (F_int!=NULL)
    {
        // assemble dFe vector (element force)
        Vec_t dFe(NDt);
        for (size_t i=0; i<NDs; ++i)
        {
            dFe(i)         = dfs(i);
            dFe(NDs+NDp+i) = dff(i);
        }
        for (size_t i=0; i<NDp; ++i) dFe(NDs+i) = dfp(i);

        // add results to Fint (internal forces)
        for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
    }
}

inline void HydroMechElem::GetState (SDPair & KeysVals, int IdxIP) const
{
    // average values
    Vec_t sig_ave(NCo), eps_ave(NCo);
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
    if (NDim==2)
    {
        KeysVals.Set("sx sy sz sxy  ex ey ez exy",
                     sig_ave(0),sig_ave(1),sig_ave(2),sig_ave(3)/Util::SQ2,
                     eps_ave(0),eps_ave(1),eps_ave(2),eps_ave(3)/Util::SQ2);
    }
    else
    {
        KeysVals.Set("sx sy sz sxy syz szx  ex ey ez exy eyz ezx",
                     sig_ave(0),sig_ave(1),sig_ave(2),sig_ave(3)/Util::SQ2,sig_ave(4)/Util::SQ2,sig_ave(5)/Util::SQ2,
                     eps_ave(0),eps_ave(1),eps_ave(2),eps_ave(3)/Util::SQ2,eps_ave(4)/Util::SQ2,eps_ave(5)/Util::SQ2);
    }
}

inline void HydroMechElem::GetState (Array<SDPair> & Results) const
{
    Results.Resize (GE->NIP);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Vec_t const & sig = static_cast<EquilibState const *>(Sta[i])->Sig;
        Vec_t const & eps = static_cast<EquilibState const *>(Sta[i])->Eps;
        if (NDim==2)
        {
            Results[i].Set("sx sy sz sxy  ex ey ez exy",
                           sig(0),sig(1),sig(2),sig(3)/Util::SQ2,
                           eps(0),eps(1),eps(2),eps(3)/Util::SQ2);
        }
        else
        {
            Results[i].Set("sx sy sz sxy syz szx  ex ey ez exy eyz ezx",
                           sig(0),sig(1),sig(2),sig(3)/Util::SQ2,sig(4)/Util::SQ2,sig(5)/Util::SQ2,
                           eps(0),eps(1),eps(2),eps(3)/Util::SQ2,eps(4)/Util::SQ2,eps(5)/Util::SQ2);
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * HydroMechElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new HydroMechElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int HydroMechElemRegister()
{
    ElementFactory["HMEquilib"] = HydroMechElemMaker;
    PROB.Set ("HMEquilib", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __HydroMechElem_dummy_int  = HydroMechElemRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_HYDROMECHELEM
