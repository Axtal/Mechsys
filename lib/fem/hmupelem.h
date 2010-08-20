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

#ifndef MECHSYS_FEM_HMUPELEM_H
#define MECHSYS_FEM_HMUPELEM_H

// MechSys
#include <mechsys/fem/element.h>
#include <mechsys/models/hmstate.h>
#include <mechsys/models/hmupdate.h>

using std::cout;
using std::endl;

namespace FEM
{

class HMupElem : public Element
{
public:
    // static
    static size_t NDn; ///< Number of DOFs per node = NDim+1
    static size_t NDt; ///< Total number of DOFs = GE->NN*NDn
    static size_t NDs; ///< Total number of DOFs of solid matrix = GE->NN*NDim
    static size_t NDp; ///< Total number of DOFs of pressure = GE->NN
    static size_t NCo; ///< Number of stress/strain components = 2*NDim
    static Mat_t  Im;  ///< Identity tensor (NCo,1)
    static Vec_t  Iv;  ///< Identity tensor (NCo)

    // Constructor
    HMupElem (int                  NDim,   ///< Space dimension
              Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
              Model        const * Mdl,    ///< Model
              SDPair       const & Prp,    ///< Properties
              SDPair       const & Ini,    ///< Initial values
              Array<Node*> const & Nodes); ///< Connectivity

    // Methods
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs,
                      NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM); ///< If setting body forces, IdxEdgeOrFace is ignored
    void CalcK       (Vec_t const & U, double Alpha, double dt, Mat_t & KK, Vec_t & dF) const;
    void CalcKC      (Mat_t & KK, Mat_t & CC)               const; ///< Augmented stiffness(K) and damping(C) matrices
    void CalcKCM     (Mat_t & KK, Mat_t & CC, Mat_t & MM)   const; ///< Augmented stiffness(K), damping(C), and mass(M) matrices
    void CalcFint    (Vec_t * F_int=NULL)                   const; ///< Calculate or set Fint. Set nodes if F_int==NULL
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL) const; ///< Update state at IPs
    void StateKeys   (Array<String> & Keys)                 const; ///< Get state keys, ex: sx, sy, sxy, ex, ey, exy
    void StateAtIP   (SDPair & KeysVals, int IdxIP)         const; ///< Get state at IP

    // Internal methods
    void Matrices (Mat_t & M, Mat_t & K, Mat_t & Q, Mat_t & Qb, Mat_t & H, Mat_t & S) const;
    void Interp   (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & Bp, Mat_t & N, Mat_t & Np, double & detJ, double & Coef) const; ///< Interpolation matrices

    // Constants
    double alp;  ///< Compressibility coefficient
    double nn;   ///< Porosity
    double rhoS; ///< Density of solids
    double rhoF; ///< Density of fluid
    double rho;  ///< Density of mixture
};

size_t HMupElem::NDn = 0;
size_t HMupElem::NDt = 0;
size_t HMupElem::NDs = 0;
size_t HMupElem::NDp = 0;
size_t HMupElem::NCo = 0;
Mat_t  HMupElem::Im;
Vec_t  HMupElem::Iv;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline HMupElem::HMupElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes)
{
    // check GE
    if (GE==NULL)   throw new Fatal("HMupElem::HMupElem: GE (geometry element) must be defined");
    if (GTy==pse_t) throw new Fatal("HMupElem::HMupElem: This element does not work for plane-stress (pse)");
    if (Mdl==NULL)  throw new Fatal("HMupElem::HMupElem: Model must be defined");

    // parameters/properties
    alp  = (Prp.HasKey("alp") ? Prp("alp") : 1.0);
    nn   = Prp("n");
    rhoS = Prp("rhoS");
    rhoF = Prp("rhoF");
    rho  = nn*rhoF+(1.0-nn)*rhoS;

    // allocate and initialize state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Sta.Push (new HMState(NDim));
        Mdl->InitIvs (Ini, Sta[i]); // initialize with effective stresses
    }

    // set constants of this class (once)
    if (NDn==0)
    {
        NDn = NDim+1;
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
        UKeys = "ux", "uy", "pw";
        for (size_t i=0; i<GE->NN; ++i) Con[i]->AddDOF("ux uy pw", "fx fy Qw");
    }
    else // 3D
    {
        UKeys = "ux", "uy", "uz", "pw";
        for (size_t i=0; i<GE->NN; ++i) Con[i]->AddDOF("ux uy uz pw", "fx fy fz Qw");
    }

    // set initial values
    if (Ini.HasKey("geostatic"))
    {
        if (!Ini.HasKey("K0"))                 throw new Fatal("HMupElem::HMupElem: For geostatic stresses, 'K0' must be provided in 'Ini' dictionary");
        if (!Ini.HasKey("gamT"))               throw new Fatal("HMupElem::HMupElem: For geostatic stresses, 'gamT' must be provided in 'Ini' dictionary");
        if (!Ini.HasKey("gamW"))               throw new Fatal("HMupElem::HMupElem: For geostatic stresses, 'gamW' must be provided in 'Ini' dictionary");
        if (NDim==2 && !Ini.HasKey("y_surf"))  throw new Fatal("HMupElem::HMupElem: For geostatic stresses in 2D, 'y_surf' (domain surface) must be provided in 'Ini' dictionary");
        if (NDim==3 && !Ini.HasKey("z_surf"))  throw new Fatal("HMupElem::HMupElem: For geostatic stresses in 3D, 'z_surf' (domain surface) must be provided in 'Ini' dictionary");
        if (NDim==2 && !Ini.HasKey("y_water")) throw new Fatal("HMupElem::HMupElem: For geostatic stresses in 2D, 'y_water' (water table) must be provided in 'Ini' dictionary");
        if (NDim==3 && !Ini.HasKey("z_water")) throw new Fatal("HMupElem::HMupElem: For geostatic stresses in 3D, 'z_water' (water table) must be provided in 'Ini' dictionary");
        double K0    = Ini("K0");
        double gamT  = Ini("gamT");
        double gamW  = Ini("gamW");
        double surf  = (NDim==2 ? Ini("y_surf")  : Ini("z_surf"));
        double water = (NDim==2 ? Ini("y_water") : Ini("z_water"));
        Vec_t X;                  // coords of IP
        Vec_t pw_at_IPs(GE->NIP); // pw at each IP
        for (size_t i=0; i<GE->NIP; ++i)
        {
            // get coordinates of IP
            CoordsOfIP (i, X);
            double z = (NDim==2 ? X(1) : X(2)); // position of point
            if (z>surf) throw new Fatal("HMupElem::HMupElem: y_surf(2D) or z_surf(3D) must be higher than every point in the domain.\n\tThere is a point [%g,%g,%g] above surf=%g.",X(0),X(1),(NDim==3?X(2):0.0),surf);

            // calculate compressive stress and pressure
            double hs   = fabs(surf-z); // column of (wet) soil
            double hw   = water-z;      // column of water
            double sv   = gamT*hs;      // compressive total vertical stress
            double pw   = gamW*hw;      // pore-water pressure
            if (Ini.HasKey("force_zero_pw")) { if (pw<0.0) pw = 0.0; }
            double sv_  = sv-pw;        // compressive effective vertical stress
            double sh_  = K0*sv_;       // compressive effective horizontal stress

            // set effective stress and water pore-pressure
            Vec_t & sig = static_cast<HMState*>(Sta[i])->Sig;
            if (NDim==2) sig = -sh_, -sv_, -sh_, 0.0;
            else         sig = -sh_, -sh_, -sv_, 0.0, 0.0, 0.0;
            static_cast<HMState*>(Sta[i])->pw = pw;

            // pw at IP
            pw_at_IPs(i) = pw;
        }

        // extrapolate pw to nodes
        Mat_t M, Mi;
        ShapeMatrix (M);
        Inv (M, Mi);
        Vec_t pw_at_Nods(Mi * pw_at_IPs);
        for (size_t i=0; i<GE->NN; ++i)
        {
            Con[i]->U    [Con[i]->UMap("pw")] += pw_at_Nods(i);
            Con[i]->NaddU[Con[i]->UMap("pw")] += 1.0;
        }
    }

    // set F in nodes due to Fint
    CalcFint ();
}

inline void HMupElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM)
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
            double coef = detJ*GE->IPs[i].w;
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
            double coef = GE->FIPs[i].w; // *detJ is not neccessary since qx,qy,qz are already multiplied by detJ (due to normal)
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
    if (has_ux || has_uy || has_uz || has_pw)
    {
        double ux = (has_ux ? BCs("ux") : 0.0);
        double uy = (has_uy ? BCs("uy") : 0.0);
        double uz = (has_uz ? BCs("uz") : 0.0);
        double pw = (has_pw ? BCs("pw") : 0.0);
        for (size_t j=0; j<GE->NFN; ++j)
        {
            size_t k = GE->FNode(IdxEdgeOrFace,j);
            if (has_ux) pU[Con[k]].first[Con[k]->UMap("ux")] = ux;
            if (has_uy) pU[Con[k]].first[Con[k]->UMap("uy")] = uy;
            if (has_uz) pU[Con[k]].first[Con[k]->UMap("uz")] = uz;
            if (has_pw) pU[Con[k]].first[Con[k]->UMap("pw")] = pw;
        }
    }
}

inline void HMupElem::Matrices (Mat_t & M, Mat_t & K, Mat_t & Q, Mat_t & Qb, Mat_t & H, Mat_t & S) const
{
    // submatrices
    M .change_dim (NDs,NDs); // mass matrix
    K .change_dim (NDs,NDs); // stiffness matrix
    Q .change_dim (NDs,NDp); // coupling matrix
    Qb.change_dim (NDp,NDs); // coupling matrix (bar)
    H .change_dim (NDp,NDp); // permeability matrix
    S .change_dim (NDp,NDp); // compressibility matrix
    set_to_zero (M);
    set_to_zero (K);
    set_to_zero (Q);
    set_to_zero (Qb);
    set_to_zero (H);
    set_to_zero (S);

    // auxiliar matrices
    double Sw=1, chiw, Cs;
    double detJ, coef;
    Vec_t Dwv;
    Mat_t C, D, Dw, B, Bp, N, Np, Kw;
    Mat_t X      (Im.num_rows(),1);
    Mat_t NtN    (NDs,NDs);
    Mat_t BtDB   (NDs,NDs);
    Mat_t BtXNp  (NDs,NDp);
    Mat_t NptItB (NDp,NDs);
    Mat_t BptkBp (NDp,NDp);
    Mat_t NptNp  (NDp,NDp);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        //Sw = static_cast<HMState const *>(Sta[i])->Sw;
        Interp         (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);
        Mdl->Hydraulic (Sta[i], Kw, chiw, Cs);
        Mdl->Stiffness (Sta[i], D, Dwv);
        Vec2ColMat     (Dwv, Dw);
        X      = Dw - chiw*Im;
        NtN    = trans(N)*N;
        BtDB   = trans(B)*D*B;
        BtXNp  = trans(B)*X*Np;
        NptItB = trans(Np)*trans(Im)*B;
        BptkBp = trans(Bp)*Kw*Bp;
        NptNp  = trans(Np)*Np;
        M     += (coef*rho) * NtN;
        K     += (coef)     * BtDB;
        Q     += (coef)     * BtXNp;
        Qb    += (coef*Sw)  * NptItB;
        H     += (coef)     * BptkBp;
        S     += (coef*Cs)  * NptNp;
    }

    //cout << "M  = \n" << PrintMatrix(M, "%18.8e");
    //cout << "K  = \n" << PrintMatrix(K, "%14.4e");
    //cout << "Q  = \n" << PrintMatrix(Q, "%14.4e");
    //cout << "Qb = \n" << PrintMatrix(Mat_t(-1.0*Qb),"%14.4e");
    //cout << "H  = \n" << PrintMatrix(Mat_t(-1.0*H), "%14.4e");
    //cout << "S  = \n" << PrintMatrix(S, "%14.4e");
}

inline void HMupElem::CalcK (Vec_t const & U, double Alpha, double dt, Mat_t & KK, Vec_t & dF) const
{
    // submatrices
    Mat_t     M, K, Q, Qb, H, S;
    Matrices (M, K, Q, Qb, H, S);

    // pressure vector
    Vec_t Up(NDp);
    for (size_t i=0; i<GE->NN; ++i) Up(i) = U(Con[i]->EQ[Con[i]->FMap("Qw")]);

    // auxiliar vector
    Vec_t HUp(H*Up);

    // assemble
    KK.change_dim (NDt,NDt);
    set_to_zero   (KK);
    for (size_t i=0; i<NDs; ++i)
    {
        for (size_t j=0; j<NDs; ++j) KK(i,j)     = K(i,j);
        for (size_t j=0; j<NDp; ++j) KK(i,NDs+j) = Q(i,j);
    }
    for (size_t i=0; i<NDp; ++i)
    {
        for (size_t j=0; j<NDs; ++j) KK(NDs+i,j)     = Qb(i,j);
        for (size_t j=0; j<NDp; ++j) KK(NDs+i,NDs+j) = S(i,j) + Alpha*dt*H(i,j);
        dF(Con[i]->EQ[Con[i]->FMap("Qw")]) += -dt*HUp(i);
    }

    //cout << "KK = \n" << PrintMatrix(KK);
    //cout << "CC = \n" << PrintMatrix(CC);
    //cout << "MM = \n" << PrintMatrix(MM);
}

inline void HMupElem::CalcKCM (Mat_t & KK, Mat_t & CC, Mat_t & MM) const
{
    // submatrices
    Mat_t     M, K, Q, Qb, H, S;
    Matrices (M, K, Q, Qb, H, S);

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
            MM(i,j) = M(i,j);
            CC(i,j) = 0.0; // Rayleigh: Am*M + Ak*K
            KK(i,j) = K(i,j);
        }
        for (size_t j=0; j<NDp; ++j) KK(i,NDs+j) = Q(i,j);
    }
    for (size_t i=0; i<NDp; ++i)
    {
        for (size_t j=0; j<NDs; ++j) CC(NDs+i,j) = Qb(i,j);
        for (size_t j=0; j<NDp; ++j)
        {
            CC(NDs+i,NDs+j) = S(i,j);
            KK(NDs+i,NDs+j) = H(i,j);
        }
    }

    //cout << "KK = \n" << PrintMatrix(KK);
    //cout << "CC = \n" << PrintMatrix(CC);
    //cout << "MM = \n" << PrintMatrix(MM);
}

inline void HMupElem::Interp (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & Bp, Mat_t & N, Mat_t & Np, double & detJ, double & Coef) const
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

    //cout << "dNdR =\n" << PrintMatrix(GE->dNdR,"%14.4e");
    //cout << "C =\n" << PrintMatrix(C,"%14.4e");
    //cout << "J =\n" << PrintMatrix(J,"%14.4e");
    //cout << "dNdR =\n" << PrintMatrix(GE->dNdR,"%18.8e");
    //cout << "C =\n" << PrintMatrix(C,"%18.8e");
    //cout << "J =\n" << PrintMatrix(J,"%18.8e");
    //cout << "Ji =\n" << PrintMatrix(Ji,"%18.8e");
    //cout << "IP = " << IP.r << " " << IP.s << " " << IP.t << endl;
    //cout << "dNdX =\n" << PrintMatrix(dNdX,"%18.8e");

    // coefficient used during integration
    Coef = detJ*IP.w;

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

inline void HMupElem::CalcFint (Vec_t * F_int) const
{
    // calc element forces
    double chiw, invQs;
    double detJ, coef;
    Mat_t C, B, Bp, N, Np, Kw;
    Vec_t BtSig (NDs);
    Vec_t fu    (NDs);
    Vec_t fp    (NDp);
    set_to_zero (fu);
    set_to_zero (fp);
    CoordMatrix (C);
    Vec_t sig(NDim*2); // total stress
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // interpolation matrices
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);

        // hydraulic variables
        Mdl->Hydraulic (Sta[i], Kw, chiw, invQs);

        // force due to displacements
        HMState const * sta = static_cast<HMState const *>(Sta[i]);
        sig   = sta->Sig - (alp*chiw*sta->pw)*Iv;
        BtSig = trans(B)*sig;
        fu   += coef * BtSig;
    }

    // set nodes directly
    if (F_int==NULL)
    {
        // add to F
        for (size_t i=0; i<GE->NN; ++i)
        {
            Con[i]->F[Con[i]->FMap("fx")] += fu(0+i*NDim);
            Con[i]->F[Con[i]->FMap("fy")] += fu(1+i*NDim); if (NDim==3) {
            Con[i]->F[Con[i]->FMap("fz")] += fu(2+i*NDim); }
            Con[i]->F[Con[i]->FMap("Qw")] += fp(i);
        }
    }

    // add to F_int
    else
    {
        for (size_t i=0; i<GE->NN; ++i)
        {
            (*F_int)(Con[i]->EQ[Con[i]->FMap("fx")]) += fu(0+i*NDim);
            (*F_int)(Con[i]->EQ[Con[i]->FMap("fy")]) += fu(1+i*NDim); if (NDim==3) {
            (*F_int)(Con[i]->EQ[Con[i]->FMap("fz")]) += fu(2+i*NDim); }
            (*F_int)(Con[i]->EQ[Con[i]->FMap("Qw")]) += fp(i);
        }
    }
}

inline void HMupElem::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // element nodal displacements
    Array<size_t> loc;
    GetLoc (loc);
    Vec_t dUe(NDt);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // split dUe
    Vec_t dus(NDs);
    Vec_t dup(NDp);
    for (size_t i=0; i<NDs; ++i) dus(i) = dUe(i);
    for (size_t i=0; i<NDp; ++i) dup(i) = dUe(NDs+i);

    // update state at each IP
    HMUpdate su(Mdl);
    double chiw, invQs, dpw;
    double detJ, coef;
    Mat_t C, B, Bp, N, Np, Kw;
    Vec_t dfs(NDs), dsig(NCo),  deps(NCo);
    Vec_t dfp(NDp), dvel(NDim), dgra(NDim);
    Vec_t Bptdvel(NDp);
    Vec_t BtdSig(NDs);
    set_to_zero (dfs);
    set_to_zero (dfp);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // interpolation matrices
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);

        // strain and (effective) stress increments
        Vec_t tmp(Np*dup);
        deps = B * dus;
        dpw  = tmp(0);
        su.Update (dpw, deps, Sta[i], dsig);
        dsig -= (alp*chiw*dpw)*Iv; // dsig is now total stress increment

        // element nodal forces (solid)
        BtdSig = trans(B)*dsig;
        dfs   += coef * BtdSig;

        // hydraulic variables
        Mdl->Hydraulic (Sta[i], Kw, chiw, invQs);

        // element nodal forces (fluid)
        dgra    = Bp * dup;
        dvel    = Kw * dgra;
        Bptdvel = trans(Bp)*dvel;
        dfp    += coef * Bptdvel; // add here S
    }

    if (F_int!=NULL)
    {
        //Mat_t M, K, Q, Qb, H, S;
        //Matrices (M, K, Q, Qb, H, S);
        //dfp += (1.0/0.1)*Qb*dus + (1.0/0.1)*S*dup;

        // assemble dFe vector (element force)
        Vec_t dFe(NDt);
        for (size_t i=0; i<NDs; ++i) dFe(i)     = dfs(i);
        for (size_t i=0; i<NDp; ++i) dFe(NDs+i) = dfp(i);

        // add results to Fint (internal forces)
        for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
    }
}

/*
inline void HMupElem::GetLoc (Array<size_t> & uLoc, Array<size_t> & pLoc) const
{
    size_t nnods = Con.Size();
    uLoc.Resize (nnods*NDim);
    pLoc.Resize (nnods);
    for (size_t i=0; i<nnods; ++i)
    {
        for (int j=0; j<NDim; ++j)
        {
            size_t idx = Con[i]->UMap(UKeys[j]); // index in Node corresponding to each DOF
            uLoc[i*NDim+j] = Con[i]->EQ[idx];
        }
        size_t idx = Con[i]->UMap(UKeys[NDim]); // index in Node corresponding to each DOF
        pLoc[i] = Con[i]->EQ[idx];
    }
}
*/

inline void HMupElem::StateKeys (Array<String> & Keys) const
{
    //           sig,eps       p,q,ev,ed       pw
    Keys.Resize (2*Mdl->NCps + 4 + Mdl->NIvs + 1);
    size_t k=0;
    Keys[k++] = "sx'";
    Keys[k++] = "sy'";
    Keys[k++] = "sz'";
    Keys[k++] = "sxy'";
    if (NDim==3)
    {
        Keys[k++] = "syz'";
        Keys[k++] = "szx'";
    }
    Keys[k++] = "ex";
    Keys[k++] = "ey";
    Keys[k++] = "ez";
    Keys[k++] = "exy";
    if (NDim==3)
    {
        Keys[k++] = "eyz";
        Keys[k++] = "ezx";
    }
    Keys[k++] = "pcam'";
    Keys[k++] = "qcam";
    Keys[k++] = "ev";
    Keys[k++] = "ed";
    for (size_t i=0; i<Mdl->NIvs; ++i) Keys[k++] = Mdl->IvNames[i];
    Keys[k++] = "pwIP"; // cannot be pw, since ParaView SegFaults in case there are point and cell data with the same name
}

inline void HMupElem::StateAtIP (SDPair & KeysVals, int IdxIP) const
{
    Vec_t  const & sig = static_cast<HMState const *>(Sta[IdxIP])->Sig;
    Vec_t  const & eps = static_cast<HMState const *>(Sta[IdxIP])->Eps;
    Vec_t  const & ivs = static_cast<HMState const *>(Sta[IdxIP])->Ivs;
    double const & pw  = static_cast<HMState const *>(Sta[IdxIP])->pw;

    if (NDim==2)
    {
        KeysVals.Set("sx' sy' sz' sxy'  ex ey ez exy  pcam' qcam  ev ed  pwIP",
                     sig(0), sig(1), sig(2), sig(3)/Util::SQ2,
                     eps(0), eps(1), eps(2), eps(3)/Util::SQ2,
                     Calc_pcam(sig), Calc_qcam(sig), Calc_ev(eps), Calc_ed(eps), pw);
    }
    else
    {
        KeysVals.Set("sx' sy' sz' sxy' syz' szx'  ex ey ez exy eyz ezx  pcam' qcam  ev ed  pwIP",
                     sig(0), sig(1), sig(2), sig(3)/Util::SQ2, sig(4)/Util::SQ2, sig(5)/Util::SQ2,
                     eps(0), eps(1), eps(2), eps(3)/Util::SQ2, eps(4)/Util::SQ2, eps(5)/Util::SQ2,
                     Calc_pcam(sig), Calc_qcam(sig), Calc_ev(eps), Calc_ed(eps), pw);
    }
    for (size_t k=0; k<Mdl->NIvs; ++k) KeysVals.Set(Mdl->IvNames[k].CStr(), ivs(k));
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * HMupElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new HMupElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int HMupElemRegister()
{
    ElementFactory["HMup"] = HMupElemMaker;
    PROB.Set ("HMup", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __HMupElem_dummy_int  = HMupElemRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_HMUPELEM
