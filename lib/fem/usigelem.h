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

#ifndef MECHSYS_FEM_USIGELEM_H
#define MECHSYS_FEM_USIGELEM_H

// MechSys
#include <mechsys/fem/element.h>
#include <mechsys/models/equilibstate.h>
#include <mechsys/models/strainupdate.h>
#include <mechsys/models/stressupdate.h>

namespace FEM
{

class USigElem : public Element
{
public:
    // Static
    static size_t NCo; ///< Number of stress/strain components == 2*NDim
    static size_t NDs; ///< Number of DOFs of stress (sigma) == NCo*NumNodesSig
    static size_t NDu; ///< Number of DOFs (displacements) == NN*NDim

    // Constructor
    USigElem (int                  NDim,   ///< Space dimension
              Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
              Model        const * Mdl,    ///< Model
              SDPair       const & Prp,    ///< Properties
              SDPair       const & Ini,    ///< Initial values
              Array<Node*> const & Nodes); ///< Connectivity

    // Destructor
    ~USigElem () { if (GEs!=NULL) delete GEs; }

    // Methods
    void IncNLocDOF  (size_t & NEq)                         const { FirstEQ = NEq;   NEq += NDs; }
    void GetLoc      (Array<size_t> & Loc)                  const; ///< Get location vector for mounting K/M matrices
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs,
                      NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM); ///< If setting body forces, IdxEdgeOrFace is ignored
    void CalcK       (Mat_t & K)                            const; ///< Stiffness
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL) const; ///< Update state at IPs
    void StateKeys   (Array<String> & Keys)                 const; ///< Get state keys, ex: sx, sy, sxy, ex, ey, exy
    void StateAtIP   (SDPair & KeysVals, int IdxIP)         const; ///< Get state at IP

    // Internal methods
    void Matrices (Mat_t & A, Mat_t & Q) const;
    void Interp   (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & N, Mat_t & Ns, double & detJ, double & Coef) const; ///< Interpolation matrices

    // Constants
    bool         StrainResid; ///< Strain-residual formulation ?
    GeomElem   * GEs;         ///< Local nodes
    mutable long FirstEQ;     ///< First equation of sig DOF
};

size_t USigElem::NCo = 0;
size_t USigElem::NDs = 0;
size_t USigElem::NDu = 0;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline USigElem::USigElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes), StrainResid(true)
{
    // check GE
    if (GE==NULL)   throw new Fatal("USigElem::USigElem: GE (geometry element) must be defined");
    if (GTy==pse_t) throw new Fatal("USigElem::USigElem: This element does not work for plane-stress (pse)");
    if (Mdl==NULL)  throw new Fatal("USigElem::USigElem: Model must be defined");

    // properties
    if (Prp.HasKey("strain_resid")) StrainResid = true;

    // local nodes
    if (strcmp(GE->Name.CStr(),"Quad8")==0) GEs = AllocGeomElem ("Quad4", NDim);
    else throw new Fatal("USigElem::USigElem: GE must be Quad8 for the time being");

    // set constants of this class (just once)
    if (NDs==0)
    {
        NCo = 2*NDim;
        NDs = NCo*GEs->NN;
        NDu = NDim*GE->NN;
    }

    // allocate and initialize state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Sta.Push (new EquilibState(NDim));
        Mdl->InitIvs (Ini, Sta[i]);
    }

    // set F in nodes due to initial stresses
    double detJ, coef;
    Mat_t C, B, N, Ns;
    Vec_t Fe(NDu);
    CoordMatrix (C);
    set_to_zero (Fe);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Interp (C, GE->IPs[i], B, N, Ns, detJ, coef);
        Vec_t const & sig = static_cast<EquilibState const *>(Sta[i])->Sig;
        Fe += (coef) * trans(B)*sig;
    }
    for (size_t i=0; i<GE->NN; ++i)
    {
        Con[i]->F[Con[i]->FMap("fx")] += Fe(0+i*NDim);
        Con[i]->F[Con[i]->FMap("fy")] += Fe(1+i*NDim);  if (NDim==3)
        Con[i]->F[Con[i]->FMap("fz")] += Fe(2+i*NDim);
    }
}

inline void USigElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM)
{
    // check
    if (!Active) throw new Fatal("USigElem::SetBCs: Element %d is inactive",Cell.ID);

    bool has_bx  = BCs.HasKey("bx");  // x component of body force
    bool has_by  = BCs.HasKey("by");  // y component of body force
    bool has_bz  = BCs.HasKey("bz");  // z component of body force
    bool has_cbx = BCs.HasKey("cbx"); // centrifugal body force along x (in axisymmetric problems)
    bool has_qx  = BCs.HasKey("qx");  // x component of distributed loading
    bool has_qy  = BCs.HasKey("qy");  // y component of distributed loading
    bool has_qz  = BCs.HasKey("qz");  // z component of distributed loading
    bool has_qn  = BCs.HasKey("qn");  // normal distributed loading
    bool has_qt  = BCs.HasKey("qt");  // tangential distributed loading (2D only)
    bool has_ux  = BCs.HasKey("ux");  // x displacement
    bool has_uy  = BCs.HasKey("uy");  // y displacement
    bool has_uz  = BCs.HasKey("uz");  // z displacement

    // force components specified
    if (has_bx || has_by || has_bz || has_cbx ||
        has_qx || has_qy || has_qz || has_qn  || has_qt)
    {
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
                double coef = GE->FIPs[i].w; // *detJ is not necessary since qx,qy,qz are already multiplied by detJ (due to normal)

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

                    //std::cout << "n = " << PrintVector(n);

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

            /*
            for (size_t j=0; j<GE->NFN; ++j)
            {
                size_t k = GE->FNode(IdxEdgeOrFace,j);
                std::cout << pF[Con[GE->FNode(IdxEdgeOrFace,j)]].first[Con[k]->FMap("fx")] << "  ";
                std::cout << pF[Con[GE->FNode(IdxEdgeOrFace,j)]].first[Con[k]->FMap("fy")] << std::endl;
            }
            */
        }
    }

    // prescribed displacements
    else if (has_ux || has_uy || has_uz)
    {
        double ux = (has_ux ? BCs("ux") : 0.0);
        double uy = (has_uy ? BCs("uy") : 0.0);
        double uz = (has_uz ? BCs("uz") : 0.0);
        for (size_t j=0; j<GE->NFN; ++j)
        {
            size_t k = GE->FNode(IdxEdgeOrFace,j);
            if (has_ux) pU[Con[k]].first[Con[k]->UMap("ux")] = ux;
            if (has_uy) pU[Con[k]].first[Con[k]->UMap("uy")] = uy;
            if (has_uz) pU[Con[k]].first[Con[k]->UMap("uz")] = uz;
        }
    }
}

inline void USigElem::GetLoc (Array<size_t> & Loc) const
{
    // Sig DOFs
    Loc.Resize (NDs + NDu);
    for (size_t i=0; i<NDs; ++i) Loc[i] = FirstEQ + i;

    // U DOFs
    //for (size_t i=0; i<GE->NN; ++i)
    //{
        //for (size_t j=0; j<UKeys.Size(); ++j)
        //{
            //size_t idx = Con[i]->UMap(UKeys[j]); // index in Node corresponding to each DOF
            //Loc[NDs + i*NDim+j] = Con[i]->EQ[idx];
        //}
    //}
}

inline void USigElem::Matrices (Mat_t & A, Mat_t & Q) const
{
    // submatrices
    A.change_dim (NDs,NDs);
    Q.change_dim (NDs,NDu);
    set_to_zero  (A);
    set_to_zero  (Q);

    // auxiliar matrices
    double detJ, coef;
    Mat_t N, Ns, B, C, D, Di;
    CoordMatrix (C);

    for (size_t i=0; i<GE->NIP; ++i)
    {
        // interpolation
        Interp (C, GE->IPs[i], B, N, Ns, detJ, coef);

        // calc D
        Mdl->Stiffness (Sta[i], D);

        // matrices
        Q += (coef) * trans(Ns)*B;
        if (StrainResid)
        {
            Inv (D, Di);
            A += (-coef) * trans(Ns)*Di*Ns;
        }
        else
        {
            A += (-coef) * trans(B)*D*B;
        }
    }

    //cout << "A    = \n" << PrintMatrix(A, "%14.3e");
    //cout << "detA = "   << Det(A) << endl;
    //cout << "Q    = \n" << PrintMatrix(Q, "%14.3e");
    //cout << "detQ = "   << Det(Q) << endl;
}

inline void USigElem::CalcK (Mat_t & K) const
{
    Mat_t A, Q;
    Matrices (A, Q);
    K.change_dim (NDs+NDu,NDs+NDu);
    set_to_zero  (K);
    if (StrainResid)
    {
        for (size_t i=0; i<NDs; ++i)
        {
            for (size_t j=0; j<NDs; ++j) K(i,j) = A(i,j);
            for (size_t j=0; j<NDu; ++j)
            {
                K(i,NDs+j) = Q(i,j);
                K(NDs+j,i) = Q(i,j);
            }
        }
    }
    else
    {
        for (size_t i=0; i<NDu; ++i)
        {
            for (size_t j=0; j<NDs; ++j)
            {
                K(i,    j) = Q(j,i);
                K(NDu+i,j) = Q(j,i);
            }
            for (size_t j=0; j<NDu; ++j) K(i,NDs+j) = A(i,j);
        }
        /*
        for (size_t i=0; i<NDu; ++i)
        {
            for (size_t j=0; j<NDu; ++j) K(NDs+i,NDs+j) = A(i,j);
            for (size_t j=0; j<NDs; ++j)
            {
                K(j,NDs+i) = Q(j,i);
                K(NDs+i,j) = Q(j,i);
            }
        }
        */
    }
}

inline void USigElem::Interp (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & N, Mat_t & Ns, double & detJ, double & Coef) const
{
    // deriv of shape func w.r.t natural coordinates
    GE ->Shape  (IP.r, IP.s, IP.t);
    GE ->Derivs (IP.r, IP.s, IP.t);
    GEs->Shape  (IP.r, IP.s, IP.t);

    // Jacobian and its determinant
    Mat_t J(GE->dNdR * C); // J = dNdR * C
    detJ = Det(J);

    // deriv of shape func w.r.t real coordinates
    Mat_t Ji;
    Inv (J, Ji);
    Mat_t dNdX(Ji * GE->dNdR); // dNdX = Inv(J) * dNdR

    // coefficient used during integration
    Coef = detJ*IP.w;

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

    // N matrix
    N.change_dim (NDim,NDu);
    set_to_zero  (N);
    for (int    i=0; i<NDim;   ++i)
    for (size_t j=0; j<GE->NN; ++j)
        N(i,i+j*NDim) = GE->N(j);

    // Ns matrix
    Ns.change_dim (NCo,NDs);
    set_to_zero   (Ns);
    for (size_t i=0; i<NCo;     ++i)
    for (size_t j=0; j<GEs->NN; ++j)
        Ns(i,i+j*NCo) = GEs->N(j);
}

inline void USigElem::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // stress increment
    Vec_t dSe(NDs);
    for (size_t i=0; i<NDs; ++i) dSe(i) = dU(FirstEQ+i);

    // displacement increment
    Vec_t dUe(NDu);
    //for (size_t i=0; i<GE->NN; ++i)
    //for (int    j=0; j<NDim;   ++j)
    //{
        //size_t idx = Con[i]->UMap(UKeys[j]);
        //dUe(i*NDim+j) = dU(Con[i]->EQ[idx]);
    //}

    // update F_int
    if (F_int!=NULL)
    {
        double detJ, coef;
        Mat_t C, B, N, Ns;
        Vec_t dsig(NCo), deps(NCo), dF1(NDs), dF2(NDu);
        CoordMatrix (C);
        set_to_zero (dF1);
        set_to_zero (dF2);
        for (size_t i=0; i<GE->NIP; ++i)
        {
            // interpolation
            Interp (C, GE->IPs[i], B, N, Ns, detJ, coef);

            // (external) stress and strain increment
            dsig = Ns * dSe;
            deps = B  * dUe;

            if (StrainResid)
            {
                // (internal) strain increment
                Vec_t deps_i(NCo);
                StrainUpdate staup(Mdl);
                staup.Update (dsig, Sta[i], deps_i);

                // (internal) forces
                dF1 += (-coef) * trans(Ns)*deps_i;
                dF1 +=  (coef) * trans(Ns)*deps;
                dF2 +=  (coef) * trans(B)*dsig;
            }
            else
            {
                // (internal) stress increment
                Vec_t dsig_i(NCo);
                StressUpdate steup(Mdl);
                steup.Update (deps, Sta[i], dsig_i);

                // (internal) forces
                dF1 +=  (coef) * trans(Ns)*deps;
                dF2 +=  (coef) * trans(B)*dsig;
                dF2 += (-coef) * trans(B)*dsig_i;
            }
        }

        // add to F_int
        for (size_t i=0; i<NDs; ++i) (*F_int)(FirstEQ+i) += dF1(i);
        for (size_t i=0; i<GE->NN; ++i)
        for (int    j=0; j<NDim;   ++j)
        {
            //size_t idx = Con[i]->UMap(UKeys[j]);
            //(*F_int)(Con[i]->EQ[idx]) += dF2(i*NDim+j);
        }
    }
}

inline void USigElem::StateKeys (Array<String> & Keys) const
{
    Keys = EquilibState::Keys;
    for (size_t i=0; i<Mdl->NIvs; ++i) Keys.Push (Mdl->IvNames[i]);
}

inline void USigElem::StateAtIP (SDPair & KeysVals, int IdxIP) const
{
    Vec_t const & sig = static_cast<EquilibState const *>(Sta[IdxIP])->Sig;
    Vec_t const & eps = static_cast<EquilibState const *>(Sta[IdxIP])->Eps;
    Vec_t const & ivs = static_cast<EquilibState const *>(Sta[IdxIP])->Ivs;

    if (NDim==2)
    {
        KeysVals.Set("sx sy sz sxy  ex ey ez exy  pcam qcam  ev ed",
                     sig(0), sig(1), sig(2), sig(3)/Util::SQ2,
                     eps(0), eps(1), eps(2), eps(3)/Util::SQ2,
                     Calc_pcam(sig), Calc_qcam(sig), Calc_ev(eps), Calc_ed(eps));
    }
    else
    {
        KeysVals.Set("sx sy sz sxy syz szx  ex ey ez exy eyz ezx  pcam qcam  ev ed",
                     sig(0), sig(1), sig(2), sig(3)/Util::SQ2, sig(4)/Util::SQ2, sig(5)/Util::SQ2,
                     eps(0), eps(1), eps(2), eps(3)/Util::SQ2, eps(4)/Util::SQ2, eps(5)/Util::SQ2,
                     Calc_pcam(sig), Calc_qcam(sig), Calc_ev(eps), Calc_ed(eps));
    }
    for (size_t k=0; k<Mdl->NIvs; ++k) KeysVals.Set(Mdl->IvNames[k].CStr(), ivs(k));
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * USigElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new USigElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int USigElemRegister()
{
    ElementFactory["USig"]   = USigElemMaker;
    ElementVarKeys["USig2D"] = std::make_pair ("ux uy",    "fx fy");
    ElementVarKeys["USig3D"] = std::make_pair ("ux uy uz", "fx fy fz");
    PROB.Set ("USig", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __USigElem_dummy_int  = USigElemRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_USIGELEM_H
