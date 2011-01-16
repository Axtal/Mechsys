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

#ifndef MECHSYS_HYDROMECH_H
#define MECHSYS_HYDROMECH_H

// MechSys
#include <mechsys/fem/equilibelem.h>
#include <mechsys/models/unsatflow.h>

namespace FEM
{

class HydroMechElem : public EquilibElem
{
public:
    // Static
    static size_t NDp; ///< Number of DOFs of pressure = GEp->NN
    static size_t NDt; ///< Total number of DOFs = NDu + NDp
    static Mat_t  Im;  ///< Identity column matrix (NCo,1)
    static Vec_t  Iv;  ///< Identity vector (NCo)
    static Vec_t  zv;  ///< z-vector: pointing up along y (2) or z (3D)

    // Constructor
    HydroMechElem (int                  NDim,   ///< Space dimension
                   Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                   Model        const * Mdl,    ///< Model
                   SDPair       const & Prp,    ///< Properties
                   SDPair       const & Ini,    ///< Initial values
                   Array<Node*> const & Nodes); ///< Connectivity

    // Destructor
    ~HydroMechElem () { if (FMdl!=NULL) delete FMdl; for (size_t i=0; i<GE->NIP; ++i) delete FSta[i]; }

    // Methods
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF);
    void ClrBCs      ();
    void AddToF      (double Time, Vec_t & F)               const;
    void CalcKCM     (Mat_t & KK, Mat_t & CC, Mat_t & MM)   const;
    void GetLoc      (Array<size_t> & Loc)                  const;
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL) const;
    void StateKeys   (Array<String> & Keys)                 const;
    void StateAtIP   (SDPair & KeysVals, int IdxIP)         const;

    // Internal Methods
    void Interp (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & Bp, Mat_t & N, Mat_t & Np, double & detJ, double & Coef) const; ///< Interpolation matrices

    // Methods for the Runge-Kutta method
    size_t NIVs       ()                                                            const;
    double GetIV      (size_t i)                                                    const;
    void   SetIV      (size_t i, double Val);
    void   CalcIVRate (double Time, Vec_t const & U, Vec_t const & V, Vec_t & Rate) const;

    // Data
    UnsatFlow              * FMdl;     ///< Flow model
    Array<UnsatFlowState*>   FSta;     ///< Flow state
    bool                     HasGrav;  ///< Has gravity ?
    BCFuncs                * GravMult; ///< gravity multiplier
};

size_t HydroMechElem::NDp = 0;
size_t HydroMechElem::NDt = 0;
Mat_t  HydroMechElem::Im;
Vec_t  HydroMechElem::Iv;
Vec_t  HydroMechElem::zv;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline HydroMechElem::HydroMechElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : EquilibElem(NDim,Cell,Mdl,Prp,Ini,Nodes), HasGrav(false), GravMult(NULL)
{
    // set constants of this class (just once)
    if (NDp==0)
    {
        NDp = GE->NN;
        NDt = NDu + NDp;

        Im.change_dim (NCo,1);
        Iv.change_dim (NCo);
        set_to_zero (Im);
        set_to_zero (Iv);
        Im(0,0)=1.0;   Im(1,0)=1.0;   Im(2,0)=1.0;
        Iv(0)  =1.0;   Iv(1)  =1.0;   Iv(2)  =1.0;

        zv.change_dim (NDim);
        if (NDim==2) zv = 0.0, 1.0;
        else         zv = 0.0, 0.0, 1.0;
    }

    // initial data
    double pos_pw  = (Ini.HasKey("only_positive_pw") ? true : false);
    double z_water = 0.0;
    double z_surf  = 0.0;
    double gamW    = Mdl->Prms("gamW");
    bool   has_pw  = Ini.HasKey("pw");
    bool   has_Sw  = Ini.HasKey("Sw");
    bool   has_geo = Ini.HasKey("geostatic");

    // check
    if ((has_pw || has_Sw) && has_geo) throw new Fatal("HydroMechElem::HydroMechElem: 'geostatic' cannot be specified together with pw or Sw (in 'Inis')");
    if (has_geo)
    {
        z_water = Ini("water");
        z_surf  = Ini("surf");
    }

    // allocate flow model TODO: move this to the outside, like Mdl, otherwise we're going to allocate one Mdl per element
    FMdl = new UnsatFlow (NDim, Mdl->Prms);

    // initial porosity and density of mixture
    double n    = Ini("n");
    double rhoW = Mdl->Prms("rhoW");
    double rhoS = Mdl->Prms("rhoS");
    rho = n*rhoW + (1.0-n)*rhoS;
    SDPair ini(Ini);

    // allocate and initialize flow state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // pw at IP
        if (has_geo)
        {
            // elevation of point
            Vec_t X;
            CoordsOfIP (i, X);
            double z = (NDim==2 ? X(1) : X(2));

            // pore-water pressure
            double hw = z_water-z; // column of water
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

    // set F in nodes due to initial (hydraulic) values
    double detJ, coef;
    Mat_t  C, B, Bp, N, Np;
    CoordMatrix (C);
    Vec_t  fe(NDp);
    set_to_zero (fe);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);
        FMdl->TgVars (FSta[i]); // set c, C, chi, and kwb
        Mat_t kBp    (FMdl->kwb * Bp);
        Vec_t kBppwe (kBp * pwe);
        fe += (coef) * trans(Bp)*kBppwe;
    }
    for (size_t i=0; i<GE->NN; ++i) Con[i]->F("qw") += fe(i);
}

inline void HydroMechElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF)
{
    bool has_flux = BCs.HasKey("flux");    // normal flux to boundary
    bool has_srcw = BCs.HasKey("srcw");    // water source term
    bool has_pw   = BCs.HasKey("pw");      // water pressure

    if (has_flux || has_srcw || has_pw)
    {
        // flux
        if (has_flux)
        {
            double qn = BCs("flux");
            double detJ, coef;
            Mat_t  Cf;
            FCoordMatrix (IdxEdgeOrFace, Cf);
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
                    Con[k]->AddToPF("qw", -coef*GE->FN(j)*qn, BCF);
                }
            }
        }

        // source
        else if (has_srcw)
        {
            double srcw = BCs("srcw");
            double detJ, coef;
            Mat_t  C, B, Bp, N, Np;
            CoordMatrix (C);
            for (size_t i=0; i<GE->NIP; ++i)
            {
                Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);
                for (size_t j=0; j<GE->NN; ++j)
                {
                    Con[j]->AddToPF("qw", srcw*coef*Np(0,j), BCF);
                }
            }
        }

        // prescribed pw
        else if (has_pw)
        {
            double pw = BCs("pw");
            for (size_t j=0; j<GE->NFN; ++j)
            {
                size_t k = GE->FNode(IdxEdgeOrFace,j);
                Con[k]->SetPU("pw", pw, BCF);
            }
        }
    }

    // other BCs => set by EquilibElem
    EquilibElem::SetBCs (IdxEdgeOrFace, BCs, BCF);

    // gravity
    if (BCs.HasKey("fgravity"))
    {
        HasGrav  = true;
        GravMult = BCF;
    }
}

inline void HydroMechElem::ClrBCs ()
{
    HasGrav  = false;
    GravMult = NULL;
}

inline void HydroMechElem::AddToF (double Time, Vec_t & F) const
{
    if (HasGrav)
    {
        // get location array
        Array<size_t> loc;
        GetLoc (loc);

        // force vector
        Vec_t  fe(NDp);
        set_to_zero (fe);
        double detJ, coef;
        Mat_t  C, B, Bp, N, Np;
        CoordMatrix (C);
        for (size_t i=0; i<GE->NIP; ++i)
        {
            Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);
            FMdl->TgVars (FSta[i]); // set c, C, chi, and kwb
            Mat_t kw     (FMdl->kwb * FMdl->gamW);
            Vec_t kwz    (kw * zv);
            fe += (-coef) * trans(Bp) * kwz;
        }

        // add results to F (external forces)
        for (size_t i=0; i<NDp; ++i)
        {
            double fm = (GravMult==NULL ? 1.0 : GravMult->fm(Time));
            F(loc[NDu+i]) += fe(i) * fm;
        }
    }
}

inline void HydroMechElem::GetLoc (Array<size_t> & Loc) const
{
    Loc.Resize (NDt);
    for (size_t i=0; i<GE->NN; ++i)
    {
        Loc[i*NDim+0] = Con[i]->Eq("ux");
        Loc[i*NDim+1] = Con[i]->Eq("uy");  if (NDim==3)
        Loc[i*NDim+2] = Con[i]->Eq("uz");
        Loc[NDu+i]    = Con[i]->Eq("pw");
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
    B.change_dim (NCo,NDu);
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
    Bp = dNdX;

    // N matrix
    N.change_dim (NDim,NDu);
    set_to_zero  (N);
    for (int    i=0; i<NDim;   ++i)
    for (size_t j=0; j<GE->NN; ++j)
        N(i,i+j*NDim) = GE->N(j);

    // Np matrix
    Np.change_dim (1,NDp);
    for (size_t j=0; j<GE->NN; ++j) Np(0,j) = GE->N(j);
}

inline void HydroMechElem::CalcKCM (Mat_t & KK, Mat_t & CC, Mat_t & MM) const
{
    // mechanical matrices
    Mat_t M, K, D;
    M.change_dim (NDu,NDu); // mass matrix
    K.change_dim (NDu,NDu); // stiffness matrix
    set_to_zero (M);
    set_to_zero (K);

    // hydraulic matrices
    Mat_t Q, Qb, H, S;
    Q .change_dim (NDu,NDp); // coupling matrix
    Qb.change_dim (NDu,NDp); // coupling matrix
    H .change_dim (NDp,NDp); // permeability matrix
    S .change_dim (NDp,NDp); // compressibility matrix
    set_to_zero (Q);
    set_to_zero (Qb);
    set_to_zero (H);
    set_to_zero (S);

    // auxiliar matrices
    Mat_t C, B, Bp, N, Np;
    double detJ, coef;
    Mat_t NtN    (NDu,NDu);
    Mat_t BtDB   (NDu,NDu);
    Mat_t BtmNp  (NDu,NDp);
    Mat_t NptNp  (NDp,NDp);
    Mat_t BptkBp (NDp,NDp);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        FMdl -> TgVars    (FSta[i]); // set c, C, chi, and kwb
        Mdl  -> Stiffness (Sta[i], D);
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);
        NtN    = trans(N)*N;
        BtDB   = trans(B)*D*B;
        BtmNp  = trans(B)*Im*Np;
        NptNp  = trans(Np)*Np;
        BptkBp = trans(Bp)*(FMdl->kwb)*Bp;
        M     += (coef*rho)       * NtN;
        K     += (coef)           * BtDB;
        Q     += (coef*FMdl->c)   * BtmNp;
        Qb    += (coef*FMdl->chi) * BtmNp;
        S     += (coef*FMdl->C)   * NptNp;
        H     += (coef)           * BptkBp;
    }

    //std::cout << "S = \n" << PrintMatrix(S) << std::endl;

    // assemble
    MM.change_dim (NDt,NDt);
    CC.change_dim (NDt,NDt);
    KK.change_dim (NDt,NDt);
    set_to_zero (MM);
    set_to_zero (CC);
    set_to_zero (KK);
    bool schemeA = true;
    for (size_t i=0; i<NDu; ++i)
    {
        for (size_t j=0; j<NDu; ++j)
        {
            MM(i,j) = M(i,j);
            CC(i,j) = 0.0; // Rayleigh: Am*M + Ak*K
            KK(i,j) = K(i,j);
        }
        for (size_t j=0; j<NDp; ++j) KK(i,NDu+j) = -Qb(i,j);
    }
    for (size_t i=0; i<NDp; ++i)
    {
        for (size_t j=0; j<NDu; ++j) CC(NDu+i,j) = Q(j,i);
        for (size_t j=0; j<NDp; ++j)
        {
            KK(NDu+i,NDu+j) = H(i,j);
            if (schemeA) CC(NDu+i,NDu+j) = S(i,j);
            else         MM(NDu+i,NDu+j) = S(i,j);
        }
    }
}

inline void HydroMechElem::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal displacements
    Vec_t dUe(NDu);
    for (size_t i=0; i<NDu; ++i) dUe(i) = dU(loc[i]);

    // element nodal pore-water pressure
    Vec_t dpwe(NDp);
    for (size_t i=0; i<NDp; ++i) dpwe(i) = dU(loc[NDu+i]);

    // auxiliary matrices and vectors
    Mat_t kBp     (NDim, NDp);
    Vec_t kBpdpwe (NDim);

    // update state at each IP
    double detJ, coef;
    Mat_t  C, B, Bp, N, Np;
    Vec_t  dFe(NDu), dsig(NCo), deps(NCo), dfe(NDp);
    set_to_zero (dFe);
    set_to_zero (dfe);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // interpolation functions
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);

        // calculate dpw at IP
        Vec_t Npdpwe(Np * dpwe);
        double dpw = Npdpwe(0);

        // strain and effective stress increments
        deps = B * dUe;
        // TODO: The model expects __effective__ stresses => Sta[i] has to be corrected here
        Mdl->SUp.Update (deps, Sta[i], dsig);

#ifdef DO_DEBUG
        double normdsig = Norm(dsig);
        if (Util::IsNan(normdsig)) throw new Fatal("HydroMechElem::UpdateState: normdsig is NaN");
#endif

        // update flow state at IP
        double dev = deps(0)+deps(1)+deps(2);
        FMdl->Update (dpw, dev, FSta[i]);

#ifdef DO_DEBUG
        if (Util::IsNan(FSta[i]->Sw))       throw new Fatal("HydroMechElem::UpdateState: Sw is NaN");
        if (Util::IsNan(FSta[i]->pc))       throw new Fatal("HydroMechElem::UpdateState: pc is NaN");
        if (Util::IsNan(FSta[i]->kwb(0,0))) throw new Fatal("HydroMechElem::UpdateState: kwb(0,0) is NaN");
#endif

        // updated tg variables
        FMdl->TgVars (FSta[i]); // set c, C, chi, and kwb

        // total stress increments
        dsig -= (FMdl->chi)*dpw*Iv;

        // element nodal forces
        kBp     = (FMdl->kwb) * Bp;
        kBpdpwe = kBp * dpwe;
        dFe += (coef) * trans(B) *dsig;
        dfe += (coef) * trans(Bp)*kBpdpwe;
    }

    // add results to Fint (internal forces)
    if (F_int!=NULL)
    {
        for (size_t i=0; i<NDu; ++i) (*F_int)(loc[i])     += dFe(i);
        for (size_t i=0; i<NDp; ++i) (*F_int)(loc[NDu+i]) += dfe(i);
    }
}

inline void HydroMechElem::StateKeys (Array<String> & Keys) const
{
    EquilibElem::StateKeys (Keys);
    Keys.Push ("n");
    Keys.Push ("pc");
    Keys.Push ("Sw");
    Keys.Push ("kw");
    Keys.Push ("qwx");
    Keys.Push ("qwy"); if (NDim==3)
    Keys.Push ("qwz");
    Keys.Push ("H");
}

inline void HydroMechElem::StateAtIP (SDPair & KeysVals, int IdxIP) const
{
    // check
    if (FSta[IdxIP]->Sw<0.0) throw new Fatal("HydroMechElem::StateAtIP: Sw<0");

    // output
    EquilibElem::StateAtIP (KeysVals, IdxIP);
    KeysVals.Set ("n",  FSta[IdxIP]->n);
    KeysVals.Set ("pc", FSta[IdxIP]->pc);
    KeysVals.Set ("Sw", FSta[IdxIP]->Sw);
    KeysVals.Set ("kw", FSta[IdxIP]->kwb(0,0)*FMdl->gamW);

    // elevation of point
    Vec_t X;
    CoordsOfIP (IdxIP, X);
    double z = (NDim==2 ? X(1) : X(2)); // the Datum will be z=0 always

    // total water head
    double pw = -FSta[IdxIP]->pc;
    KeysVals.Set ("H", z + pw/FMdl->gamW);

    // vector of current pw at nodes of element
    Vec_t pwe(NDp);
    for (size_t i=0; i<GE->NN; ++i) pwe(i) = Con[i]->U("pw");

    // pore-water pressure gradient at IP
    double detJ, coef;
    Mat_t  C, B, Bp, N, Np;
    CoordMatrix (C);
    Interp (C, GE->IPs[IdxIP], B, Bp, N, Np, detJ, coef);
    Vec_t grad_pw(Bp * pwe);

    // relative specific discharge
    grad_pw += FMdl->gamW*zv;              // grad_pw = grad_pw + gamW*z
    Vec_t mqw(FSta[IdxIP]->kwb * grad_pw); // -qw
    KeysVals.Set ("qwx", -mqw(0));
    KeysVals.Set ("qwy", -mqw(1)); if (NDim==3)
    KeysVals.Set ("qwz", -mqw(2));
}

inline size_t HydroMechElem::NIVs () const 
{
    return NCo*GE->NIP;
}

inline double HydroMechElem::GetIV (size_t i) const 
{
    int iip = static_cast<int>(i) / static_cast<int>(NCo); // index of IP
    int ico = static_cast<int>(i) % static_cast<int>(NCo); // index of component
    return static_cast<EquilibState const*>(Sta[iip])->Sig(ico);
}

inline void HydroMechElem::SetIV (size_t i, double Val) 
{
    int iip = static_cast<int>(i) / static_cast<int>(NCo); // index of IP
    int ico = static_cast<int>(i) % static_cast<int>(NCo); // index of component
    static_cast<EquilibState*>(Sta[iip])->Sig(ico) = Val;
}

inline void HydroMechElem::CalcIVRate (double Time, Vec_t const & U, Vec_t const & V, Vec_t & Rate) const
{
    // resize rate
    Rate.change_dim (NCo*GE->NIP);

    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal velocities
    Vec_t Ve(NDu);
    for (size_t i=0; i<NDu; ++i) Ve(i) = V(loc[i]);

    // element dpwdt
    Vec_t dpwedt(NDp);
    for (size_t i=0; i<NDp; ++i) dpwedt(i) = V(loc[NDu+i]);

    // calc dsigdt
    double detJ, coef;
    Mat_t  C, B, D, Bp, N, Np;
    Vec_t  depsdt(NCo), dsigdt(NCo);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // interpolation functions
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);

        // calculate dpwdt at IP
        Vec_t Npdpwedt(Np * dpwedt);
        double dpwdt = Npdpwedt(0);

        // effective stress rate
        Mdl->Stiffness (Sta[i], D);
        depsdt = B * Ve;
        dsigdt = D * depsdt;

        // total stress rate
        dsigdt -= (FMdl->chi)*dpwdt*Iv;

        // set rate vector
        for (size_t j=0; j<NCo; ++j) Rate(j+i*NCo) = dsigdt(j);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * HydroMechElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new HydroMechElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int HydroMechElemRegister()
{
    ElementFactory["HydroMech"]   = HydroMechElemMaker;
    ElementVarKeys["HydroMech2D"] = std::make_pair ("ux uy pw",    "fx fy qw");
    ElementVarKeys["HydroMech3D"] = std::make_pair ("ux uy uz pw", "fx fy fz qw");
    PROB.Set ("HydroMech", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __HydroMechElem_dummy_int  = HydroMechElemRegister();

}; // namespace FEM

#endif // MECHSYS_HYDROMECH_H
