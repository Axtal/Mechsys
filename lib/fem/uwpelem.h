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

#ifndef MECHSYS_UWPELEM_H
#define MECHSYS_UWPELEM_H

// MechSys
#include <mechsys/fem/equilibelem.h>
#include <mechsys/models/unsatflow.h>

namespace FEM
{

class UWPElem : public EquilibElem
{
public:
    // typedefs
    typedef std::map<size_t,SDPair  const *> IdxToBcs_t;
    typedef std::map<size_t,BCFuncs       *> IdxToBcf_t;

    // Static
    static size_t NDp;  ///< Number of DOFs of pressure = GEp->NN
    static size_t NDt;  ///< Total number of DOFs = NDu + NDp
    static Mat_t  M;    ///< NDu,NDu
    static Mat_t  Mw;   ///< NDu,NDu
    static Mat_t  Hw;   ///< NDp,NDp
    static Vec_t  f;    ///< NDu
    static Vec_t  fw;   ///< NDu
    static Vec_t  hw;   ///< NDu
    static Vec_t  cw;   ///< NDu
    static Vec_t  ew;   ///< NDp
    static Mat_t  dNdX; ///< NDim,GE->NN
    static Mat_t  B;    ///< NCo,NDu
    static Mat_t  Bp;   ///< NDim,GE->NN  == dNdX
    static Mat_t  N;    ///< NDim,NDu
    static Mat_t  Np;   ///< 1,NDp
    static Vec_t  Npv;  ///< NDp
    static Mat_t  J;    ///< NDim,NDim     Jacobian
    static Mat_t  Ji;   ///< NDim,NDim     Inverse of jacobian
    static Mat_t  Jf;   ///< NDim-1,NDim   Jacobian of face
    static double DetJ; ///< Determinant of Jacobian
    static double Coef; ///< Coef for integration
    static Vec_t  V;    ///< NDu
    static Vec_t  Ww;   ///< NDu
    static Vec_t  Pw;   ///< NDp
    static Mat_t  Co;   ///< GE->NN,NDim   Coordinates matrix
    static Mat_t  Cf;   ///< GE->NFN,NDim  Face coordinates matrix
    static Mat_t  lw;   ///< NDim,NDim     O2 tensor: velocity gradient == Sum_n (vsn+wwn) dy Gn
    static Vec_t  ww;   ///< NDim          ww = N * Ww
    static Vec_t  d;    ///< NCo           O2 tensor: rate of deformation of solids == B*v
    static Vec_t  fwd;  ///< NDim          drag force vector

    // Constructor
    UWPElem (int                  NDim,   ///< Space dimension
             Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
             Model        const * Mdl,    ///< Model
             Model        const * XMdl,   ///< Extra Model
             SDPair       const & Prp,    ///< Properties
             SDPair       const & Ini,    ///< Initial values
             Array<Node*> const & Nodes); ///< Connectivity

    // Destructor
    ~UWPElem () { for (size_t i=0; i<GE->NIP; ++i) delete FSta[i]; }

    // Methods
    void SetBCs        (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF);
    void CalcSurfLoads (double Time, Vec_t & Sl, Vec_t & Sq)  const;
    void GetLoc        (Array<size_t> & LocV, Array<size_t> & LocWw, Array<size_t> & LocPw) const;
    void ElemEqs       (double Time, Vec_t const & V_g, Vec_t const & Ww_g, Vec_t const & Pw_g) const;
    void StateKeys     (Array<String> & Keys)                 const;
    void StateAtIP     (SDPair & KeysVals, int IdxIP)         const;

    // Internal Methods
    void Interp (Mat_t const & C, IntegPoint const & IP) const; ///< Interpolation matrices

    // Data
    Array<UnsatFlowState*>         FSta;    ///< Flow state
    UnsatFlow              const * FMdl;    ///< Flow model
    IdxToBcs_t                     LBCs;    ///< Loads boundary conditions
    IdxToBcf_t                     LBCf;    ///< Loads callbacks
    bool                           HasGrav; ///< Has gravity?
    Vec_t                          g;       ///< Gravity vector: [0,-g]^T  or  [0,0,-g]^T
};

size_t UWPElem::NDp = 0;
size_t UWPElem::NDt = 0;
Mat_t  UWPElem::M;
Mat_t  UWPElem::Mw;
Mat_t  UWPElem::Hw;
Vec_t  UWPElem::f;
Vec_t  UWPElem::fw;
Vec_t  UWPElem::hw;
Vec_t  UWPElem::cw;
Vec_t  UWPElem::ew;
Mat_t  UWPElem::dNdX;
Mat_t  UWPElem::B;
Mat_t  UWPElem::Bp;
Mat_t  UWPElem::N;
Mat_t  UWPElem::Np;
Vec_t  UWPElem::Npv;
Mat_t  UWPElem::J;
Mat_t  UWPElem::Ji;
Mat_t  UWPElem::Jf;
double UWPElem::DetJ;
double UWPElem::Coef;
Vec_t  UWPElem::V;
Vec_t  UWPElem::Ww;
Vec_t  UWPElem::Pw;
Mat_t  UWPElem::Co;
Mat_t  UWPElem::Cf;
Mat_t  UWPElem::lw;
Vec_t  UWPElem::ww;
Vec_t  UWPElem::d;
Vec_t  UWPElem::fwd;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline UWPElem::UWPElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : EquilibElem(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes),
      FMdl(static_cast<UnsatFlow const*>(XMdl)), HasGrav(false)
{
    // check
    if (XMdl==NULL) throw new Fatal("UWPElem::UWPElem: E(x)tra Model (flow model) must be defined by means of 'xname' key");

    // gravity vector
    g.change_dim(NDim);
    set_to_zero(g);

    // set constants of this class (just once)
    if (NDp==0)
    {
        NDp = GE->NN;
        NDt = NDu + NDp;

        M .change_dim (NDu,NDu);
        Mw.change_dim (NDu,NDu);
        Hw.change_dim (NDp,NDp);
        f .change_dim (NDu);
        fw.change_dim (NDu);
        hw.change_dim (NDu);
        cw.change_dim (NDu);
        ew.change_dim (NDp);

        dNdX.change_dim (NDim, GE->NN);
        B   .change_dim (NCo,  NDu);
        Bp  .change_dim (NDim, GE->NN);
        N   .change_dim (NDim, NDu);
        Np  .change_dim (1,    NDp);
        Npv .change_dim (NDp);

        J .change_dim (NDim,   NDim);
        Ji.change_dim (NDim,   NDim);
        Jf.change_dim (NDim-1, NDim);

        V .change_dim (NDu);
        Ww.change_dim (NDu);
        Pw.change_dim (NDp);

        Co.change_dim (GE->NN,  NDim);
        Cf.change_dim (GE->NFN, NDim);

        lw .change_dim (NDim, NDim);
        ww .change_dim (NDim);
        d  .change_dim (NCo);
        fwd.change_dim (NDim);
    }

    // initial data
    bool   geosta = (Prp.HasKey("geosta") ? Prp("geosta")>0 : false);
    bool   pos_pw = (Prp.HasKey("pospw")  ? Prp("pospw")>0  : false);
    double gamW   = XMdl->Prms("gamW");

    // allocate and initialize flow state at each IP
    SDPair ini(Ini);
    for (size_t i=0; i<GE->NIP; ++i)
    {
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

        // init flow state
        FSta.Push     (new UnsatFlowState(NDim));
        FMdl->InitIvs (ini, FSta[i]);
    }

    // pore-water pressure at nodes (extrapolated)
    Vec_t pwe(NDp);
    Array<SDPair> res;
    StateAtNodes (res);
    for (size_t i=0; i<NDp; ++i)
    {
        double pw = -res[i]("pc");
        Con[i]->U("pw") = pw;
        pwe(i) = pw;
    }
}

inline void UWPElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF)
{
    // has gravity?
    HasGrav = BCs.HasKey("gravity");
    if (HasGrav)
    {
        if (NDim==2) g = 0.0,      -BCs("gravity");
        else         g = 0.0, 0.0, -BCs("gravity");
    }

    // set essential BCs (do not depend on geometry shape)
    if (BCs.HasKey("ux"))     { for (size_t j=0; j<GE->NFN; ++j) Con[GE->FNode(IdxEdgeOrFace,j)]->SetPU("ux", BCs("ux"), BCF); }
    if (BCs.HasKey("uy"))     { for (size_t j=0; j<GE->NFN; ++j) Con[GE->FNode(IdxEdgeOrFace,j)]->SetPU("uy", BCs("uy"), BCF); }
    if (BCs.HasKey("uz"))     { for (size_t j=0; j<GE->NFN; ++j) Con[GE->FNode(IdxEdgeOrFace,j)]->SetPU("uz", BCs("uz"), BCF); }
    if (BCs.HasKey("pw"))     { for (size_t j=0; j<GE->NFN; ++j) Con[GE->FNode(IdxEdgeOrFace,j)]->SetPU("pw", BCs("pw"), BCF); }
    if (BCs.HasKey("incsup")) { for (size_t j=0; j<GE->NFN; ++j) Con[GE->FNode(IdxEdgeOrFace,j)]->SetIncSup(BCs("alpha"));     }

    // set maps for prescribed loads
    bool has_qx = BCs.HasKey("qx");  // x component of distributed loading
    bool has_qy = BCs.HasKey("qy");  // y component of distributed loading
    bool has_qz = BCs.HasKey("qz");  // z component of distributed loading
    bool has_qn = BCs.HasKey("qn");  // normal distributed loading
    bool has_qt = BCs.HasKey("qt");  // tangential distributed loading (2D only)
    bool has_qw = BCs.HasKey("qw");  // water flux
    if (has_qx || has_qy || has_qz || has_qn || has_qt || has_qw)
    {
        LBCs[IdxEdgeOrFace] = &BCs;
        LBCf[IdxEdgeOrFace] = BCF;
    }
}

inline void UWPElem::CalcSurfLoads (double Time, Vec_t & Sl, Vec_t & Sq) const
{
    // surface load
    //Sl.change_dim(NDu);
    set_to_zero(Sl);

    // surface flux
    //Sq.change_dim(NDp);
    set_to_zero(Sq);

    for (IdxToBcs_t::const_iterator it=LBCs.begin(); it!=LBCs.end(); ++it)
    {
        // local index of edge of face
        size_t ief = it->first; // IdxEdgeOrFace

        // load multiplier
        IdxToBcf_t::const_iterator itt = LBCf.find(ief);
        double lm = (itt->second==NULL ? 1.0 : itt->second->fm(Time));

        // prescribed quantities
        bool has_qx = it->second->HasKey("qx");  // x component of distributed loading
        bool has_qy = it->second->HasKey("qy");  // y component of distributed loading
        bool has_qz = it->second->HasKey("qz");  // z component of distributed loading
        bool has_qn = it->second->HasKey("qn");  // normal distributed loading
        bool has_qt = it->second->HasKey("qt");  // tangential distributed loading (2D only)
        bool has_qw = it->second->HasKey("qw");  // water flux

        // matrix of coordinates of edge/face
        for (size_t i=0; i<GE->NFN; ++i)
        for (int    j=0; j<NDim;    ++j)
            Cf(i,j) = Con[GE->FNode(ief,i)]->Vert.C[j];

        // surface loading
        if (has_qx || has_qy || has_qz || has_qn || has_qt)
        {
            // loading
            double qx = (has_qx ? (*it->second)("qx") : 0.0);
            double qy = (has_qy ? (*it->second)("qy") : 0.0);
            double qz = (has_qz ? (*it->second)("qz") : 0.0);
            double qn = (has_qn ? (*it->second)("qn") : 0.0);
            double qt = (has_qt ? (*it->second)("qt") : 0.0);

            // set
            double detJ, coef;
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (Cf, GE->FIPs[i], Jf, detJ, coef);
                coef /= detJ; // *detJ is not necessary since qx,qy,qz are already multiplied by detJ (due to normal)

                // calculate qx, qy and qz from qn and qt
                if (has_qn || has_qt)
                {
                    // normal to edge/face
                    Vec_t n(NDim); // normal multiplied by detJ
                    if (NDim==2) n = Jf(0,1), -Jf(0,0);
                    else
                    {
                        // vectorial product
                        Vec_t a(3);  a = Jf(0,0), Jf(0,1), Jf(0,2);
                        Vec_t b(3);  b = Jf(1,0), Jf(1,1), Jf(1,2);
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

                // set boundary conditions
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    size_t k = GE->FNode(ief,j);
                    Sl(k*NDim+0) += coef*GE->FN(j)*qx*lm;
                    Sl(k*NDim+1) += coef*GE->FN(j)*qy*lm;  if (NDim==3)
                    Sl(k*NDim+2) += coef*GE->FN(j)*qz*lm;
                }
            }
        }

        // water flux
        if (has_qw)
        {
            double qw = (*it->second)("qw");
            double detJ, coef;
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (Cf, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    size_t k = GE->FNode(ief,j);
                    Sq(k) += coef*GE->FN(j)*(-qw)*lm;
                }
            }
        }
    }
}

inline void UWPElem::GetLoc (Array<size_t> & LocV, Array<size_t> & LocWw, Array<size_t> & LocPw) const
{
    LocV .Resize (NDu);
    LocWw.Resize (NDu);
    LocPw.Resize (NDp);
    for (size_t i=0; i<GE->NN; ++i)
    {
        LocV [i*NDim+0] = Con[i]->Eq("ux");
        LocV [i*NDim+1] = Con[i]->Eq("uy");  if (NDim==3)
        LocV [i*NDim+2] = Con[i]->Eq("uz");
        LocWw[i*NDim+0] = Con[i]->Eq("wwx");
        LocWw[i*NDim+1] = Con[i]->Eq("wwy");  if (NDim==3)
        LocWw[i*NDim+2] = Con[i]->Eq("wwz");
        LocPw[i]        = Con[i]->Eq("pw");
    }
}

inline void UWPElem::ElemEqs (double Time, Vec_t const & V_g, Vec_t const & Ww_g, Vec_t const & Pw_g) const
{
    // element vectors
    for (size_t i=0; i<GE->NN; ++i)
    {
        V (i*NDim+0) = V_g (Con[i]->Eq("ux"));
        V (i*NDim+1) = V_g (Con[i]->Eq("uy"));  if (NDim==3)
        V (i*NDim+2) = V_g (Con[i]->Eq("uz"));
        Ww(i*NDim+0) = Ww_g(Con[i]->Eq("wwx"));
        Ww(i*NDim+1) = Ww_g(Con[i]->Eq("wwy"));  if (NDim==3)
        Ww(i*NDim+2) = Ww_g(Con[i]->Eq("wwz"));
        Pw(i)        = Pw_g(Con[i]->Eq("pw"));
    }

    // clear output matrices and vectors
    set_to_zero (M);
    set_to_zero (Mw);
    set_to_zero (Hw);
    set_to_zero (fw);
    set_to_zero (hw);
    set_to_zero (cw);

    // surface loads
    CalcSurfLoads (Time, f, ew); // f=Sl  ew=Sq

    // coordinates matrix
    for (size_t i=0; i<GE->NN; ++i)
    for (int    j=0; j<NDim;   ++j)
        Co(i,j) = Con[i]->Vert.C[j];

    // integrate
    double Cpw, Cvs, trd;
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Interp (Co, GE->IPs[i]);

        // water relative velocity at IP
        ww = N * Ww;

        // calc Cpw, Cvs and fwd
        FMdl->TgVars (FSta[i], ww,  Cpw, Cvs, fwd);

        // values at IP
        double n    = FSta[i]->n;          // porosity
        double Sw   = FSta[i]->Sw;         // water saturation
        double nw   = n*Sw;                // water volume fraction
        double RhoW = FSta[i]->RhoW;       // real water density
        double RhoS = FSta[i]->RhoS;       // real solids density
        double rhow = nw*RhoW;             // water partial density
        double rho  = (1.0-n)*RhoS + rhow; // mixture density

        // calc lw
        set_to_zero (lw);
        for (size_t j=0; j<GE->NN; ++j)
        {
            for (int r=0; r<NDim; ++r)
            for (int s=0; s<NDim; ++s)
                lw(r,s) += (V(j*NDim+r)+Ww(j*NDim+r)) * dNdX(s,j);
        }

        // solids rate of deformation
        d   = B * V;
        trd = Tra(d);

        // vectors
        f  += (Coef * rho)     * trans(N)  * g;
        f  -= (Coef)           * trans(B)  * static_cast<EquilibState*>(Sta[i])->Sig;
        cw += (Coef * rhow)    * trans(N)  * lw * ww;
        hw += (Coef * nw)      * trans(N)  * Bp * Pw;
        fw += (Coef)           * trans(N)  * fwd;
        fw += (Coef * rhow)    * trans(N)  * g;
        ew += (Coef * rhow)    * trans(Bp) * ww;
        ew -= (Coef * Cvs*trd) * Npv;

        // matrices
        M  += (Coef * rho)  * trans(N)  * N;
        Mw += (Coef * rhow) * trans(N)  * N;
        Hw += (Coef * Cpw)  * trans(Np) * Np;
    }
}

inline void UWPElem::StateKeys (Array<String> & Keys) const
{
    EquilibElem::StateKeys (Keys);
    Keys.Push ("n");
    Keys.Push ("pc");
    Keys.Push ("Sw");
    Keys.Push ("rw");
    Keys.Push ("qwx");
    Keys.Push ("qwy"); if (NDim==3)
    Keys.Push ("qwz");
    Keys.Push ("H");
}

inline void UWPElem::StateAtIP (SDPair & KeysVals, int IdxIP) const
{
    // check
    if (FSta[IdxIP]->Sw<0.0) throw new Fatal("UWPElem::StateAtIP: Sw<0");

    // output
    double rw = FMdl->rw(FSta[IdxIP]->Sw);
    EquilibElem::StateAtIP (KeysVals, IdxIP);
    KeysVals.Set ("n",  FSta[IdxIP]->n);
    KeysVals.Set ("pc", FSta[IdxIP]->pc);
    KeysVals.Set ("Sw", FSta[IdxIP]->Sw);
    KeysVals.Set ("rw", rw);

    // elevation of point
    Vec_t X;
    CoordsOfIP (IdxIP, X);
    double z = (NDim==2 ? X(1) : X(2)); // the Datum will be z=0 always

    // total water head
    double pw = -FSta[IdxIP]->pc;
    KeysVals.Set ("H", z + pw/FMdl->GamW);

    // vector of current pw at nodes of element
    Vec_t pwe(NDp);
    for (size_t i=0; i<GE->NN; ++i) pwe(i) = Con[i]->U("pw");

    // pore-water pressure gradient at IP
    //double detJ, coef;
    //Mat_t  C, dNdX, B, Bp, N, Np;
    //CoordMatrix (C);
    //Interp (C, GE->IPs[IdxIP], dNdX, B, Bp, N, Np, detJ, coef);
    //Vec_t grad_pw(Bp * pwe);

    // relative specific discharge
    //grad_pw += FMdl->GamW*zv;          // grad_pw = grad_pw + gamW*z
    //Vec_t mqw(FMdl->kwsatb * grad_pw); // -qw
    //KeysVals.Set ("qwx", -rw*mqw(0));
    //KeysVals.Set ("qwy", -rw*mqw(1)); if (NDim==3)
    //KeysVals.Set ("qwz", -rw*mqw(2));
}

inline void UWPElem::Interp (Mat_t const & C, IntegPoint const & IP) const
{
    // deriv of shape func w.r.t natural coordinates
    GE->Shape  (IP.r, IP.s, IP.t);
    GE->Derivs (IP.r, IP.s, IP.t);

    // Jacobian and its determinant
    J    = GE->dNdR * C;
    DetJ = Det(J);

    // deriv of shape func w.r.t real coordinates
    Inv (J, Ji);
    dNdX = Ji * GE->dNdR; // dNdX = Inv(J) * dNdR

    // coefficient used during integration
    Coef = h*DetJ*IP.w;

    // B matrix
    set_to_zero (B);
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
    Bp = dNdX;

    // N matrix
    set_to_zero (N);
    for (int    i=0; i<NDim;   ++i)
    for (size_t j=0; j<GE->NN; ++j)
        N(i,i+j*NDim) = GE->N(j);

    // Np matrix
    for (size_t j=0; j<GE->NN; ++j)
    {
        Np (0,j) = GE->N(j);
        Npv(j)   = GE->N(j);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * UWPElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new UWPElem(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes); }

// Register element
int UWPElemRegister()
{
    ElementFactory["UWP"]   = UWPElemMaker;
    ElementVarKeys["UWP2D"] = std::make_pair ("ux uy wwx wwy pw",        "fx fy fwx fwy qw");
    ElementVarKeys["UWP3D"] = std::make_pair ("ux uy uz wwx wwy wwz pw", "fx fy fz fwx fwy fwz qw");
    PROB.Set ("UWP", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __UWPElem_dummy_int  = UWPElemRegister();

}; // namespace FEM

#endif // MECHSYS_UWPELEM_H
