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

// This is the mixed u-sig element (usigelem) with the condensation of the local DOFs

#ifndef MECHSYS_FEM_USIGCONDELEM_H
#define MECHSYS_FEM_USIGCONDELEM_H

// MechSys
#include <mechsys/fem/element.h>
#include <mechsys/models/equilibstate.h>
#include <mechsys/models/stressupdate.h>

namespace FEM
{

class USigCondElem : public Element
{
public:
    // static
    static size_t NCo; ///< Number of stress/strain components == 2*NDim
    static size_t NDs; ///< Number of DOFs of stress (sigma) == NCo*NumNodesSig
    static size_t NDu; ///< Number of DOFs (displacements) == NN*NDim

    // Constructor
    USigCondElem (int                  NDim,   ///< Space dimension
                  Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                  Model        const * Mdl,    ///< Model
                  Model        const * XMdl,   ///< Extra Model
                  SDPair       const & Prp,    ///< Properties
                  SDPair       const & Ini,    ///< Initial values
                  Array<Node*> const & Nodes); ///< Array with all nodes (used to set the connectivity)

    // Destructor
    ~USigCondElem () { if (GEs!=NULL) delete GEs; }

    // Methods
    void GetLoc      (Array<size_t> & Loc)                               const; ///< Get location vector for mounting K/M matrices
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF); ///< If setting body forces, IdxEdgeOrFace is ignored
    void CalcK       (Mat_t & K)                                         const; ///< Stiffness
    void CalcFint    (Vec_t * F_int=NULL)                                const; ///< Calculate or set Fint. Set nodes if F_int==NULL
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL)              const; ///< Update state at IPs
    void StateKeys   (Array<String> & Keys)                              const; ///< Get state keys, ex: sx, sy, sxy, ex, ey, exy
    void StateAtIP   (SDPair & KeysVals, int IdxIP)                      const; ///< Get state at IP

    // Internal methods
    void Matrices (Mat_t & Ai, Mat_t & Q) const;
    void Interp   (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & N, Mat_t & Nsig, double & detJ, double & Coef) const; ///< Interpolation matrices

    // Constants
    GeomElem   * GEs;        ///< Local nodes
    mutable long FirstEQsig; ///< First equation of sig DOF
};

size_t USigCondElem::NCo = 0;
size_t USigCondElem::NDs = 0;
size_t USigCondElem::NDu = 0;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline USigCondElem::USigCondElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes)
{
    // check GE
    if (GE==NULL)   throw new Fatal("USigCondElem::USigCondElem: GE (geometry element) must be defined");
    if (GTy==pse_t) throw new Fatal("USigCondElem::USigCondElem: This element does not work for plane-stress (pse)");
    if (Mdl==NULL)  throw new Fatal("USigCondElem::USigCondElem: Model must be defined");

    // local nodes
    if (strcmp(GE->Name.CStr(),"Quad8")==0) GEs = AllocGeomElem ("Quad4", NDim);
    else throw new Fatal("USigCondElem::USigCondElem: GE must be Quad8 for the time being");

    // set constants of this class (once)
    if (NDs==0)
    {
        NCo = 2*NDim;
        NDs = NCo*GEs->NN;
        NDu = NDim*GE->NN;
    }

    // allocate and initialize state at each local node
    for (size_t i=0; i<GEs->NN; ++i)
    {
        Sta.Push (new EquilibState(NDim));
        Mdl->InitIvs (Ini, Sta[i]); // initialize with effective stresses
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
        Con[i]->F("fx") += Fe(0+i*NDim);
        Con[i]->F("fy") += Fe(1+i*NDim);  if (NDim==3)
        Con[i]->F("fz") += Fe(2+i*NDim);
    }
}

inline void USigCondElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF)
{
    // deactivate/activate commands
    if (BCs.HasKey("deactivate")) throw new Fatal("USigCondElem::SetBCs: 'deactivate' command failed: method not implemented yet");
    if (BCs.HasKey("activate"))   throw new Fatal("USigCondElem::SetBCs: 'activate' command failed: method not implemented yet");

    // check
    if (!Active) throw new Fatal("USigCondElem::SetBCs: Element %d is inactive",Cell.ID);

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
    bool has_sup = BCs.HasKey("incsup");  // inclined support

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

                // set boundary conditions
                for (size_t j=0; j<GE->NN; ++j)
                {
                    Con[j]->AddToPF("fx", coef*GE->N(j)*bx, BCF);
                    Con[j]->AddToPF("fy", coef*GE->N(j)*by, BCF);  if (NDim==3)
                    Con[j]->AddToPF("fz", coef*GE->N(j)*bz, BCF);
                }
            }
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
                    size_t k = GE->FNode(IdxEdgeOrFace,j);
                    Con[k]->AddToPF("fx", coef*GE->FN(j)*qx, BCF);
                    Con[k]->AddToPF("fy", coef*GE->FN(j)*qy, BCF);  if (NDim==3)
                    Con[k]->AddToPF("fz", coef*GE->FN(j)*qz, BCF);
                }
            }
        }
    }

    // prescribed displacements
    else if (has_ux || has_uy || has_uz || has_sup)
    {
        double ux = (has_ux ? BCs("ux") : 0.0);
        double uy = (has_uy ? BCs("uy") : 0.0);
        double uz = (has_uz ? BCs("uz") : 0.0);
        double alpha = (has_sup ? BCs("alpha") : 0.0);
        for (size_t j=0; j<GE->NFN; ++j)
        {
            size_t k = GE->FNode(IdxEdgeOrFace,j);
            if (has_ux) Con[k]->SetPU("ux", ux, BCF);
            if (has_uy) Con[k]->SetPU("uy", uy, BCF);
            if (has_uz) Con[k]->SetPU("uz", uz, BCF);
            if (has_sup) Con[k]->SetIncSup (alpha);
        }
    }
}

inline void USigCondElem::GetLoc (Array<size_t> & Loc) const
{
    Loc.Resize (NDu);
    for (size_t i=0; i<GE->NN; ++i)
    {
        Loc[i*NDim+0] = Con[i]->Eq("ux");
        Loc[i*NDim+1] = Con[i]->Eq("uy");  if (NDim==3)
        Loc[i*NDim+2] = Con[i]->Eq("uz");
    }
}

inline void USigCondElem::Matrices (Mat_t & Ai, Mat_t & Q) const
{
    // submatrices
    Q.change_dim (NDs,NDu);
    set_to_zero  (Q);

    // auxiliar matrices
    double detJ, coef;
    Mat_t NtDiN  (NDs,NDs);
    Mat_t NtB    (NDs,NDu);
    Mat_t N, Nsig, A(NDs,NDs), B, C, D(NCo,NCo), Di;
    CoordMatrix  (C);
    set_to_zero  (A);

    double E  = 1000.0;
    double nu = 0.2;
    double c  = E/((1.0+nu)*(1.0-2.0*nu));
    D = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,
             c*nu ,  c*(1.0-nu),      c*nu ,            0.0,
             c*nu ,       c*nu , c*(1.0-nu),            0.0,
              0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu);

    for (size_t i=0; i<GE->NIP; ++i)
    {
        Interp (C, GE->IPs[i], B, N, Nsig, detJ, coef);

        Inv (D, Di);

        NtDiN = trans(Nsig)*Di*Nsig;
        NtB   = trans(Nsig)*B;
        A    -= (coef)     * NtDiN;
        Q    += (coef)     * NtB;
    }
    
    Inv (A, Ai);
}

inline void USigCondElem::CalcK (Mat_t & K) const
{
  Mat_t Ai, Q;
  Matrices (Ai, Q);
  K.change_dim (NDu,NDu);
  
  K = (-1.)*trans(Q)*Ai*Q;
}

inline void USigCondElem::Interp (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & N, Mat_t & Nsig, double & detJ, double & Coef) const
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

    // Nsig matrix
    Nsig.change_dim (NCo,NDs);
    set_to_zero     (Nsig);
    for (size_t i=0; i<NCo;     ++i)
    for (size_t j=0; j<GEs->NN; ++j)
        Nsig(i,i+j*NCo) = GEs->N(j);
}

inline void USigCondElem::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
}

inline void USigCondElem::StateKeys (Array<String> & Keys) const
{
}

inline void USigCondElem::StateAtIP (SDPair & KeysVals, int IdxIP) const
{
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * USigCondMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new USigCondElem(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes); }

// Register element
int USigCondRegister()
{
    ElementFactory["USigCond"]   = USigCondMaker;
    ElementVarKeys["USigCond2D"] = std::make_pair ("ux uy",    "fx fy");
    ElementVarKeys["USigCond3D"] = std::make_pair ("ux uy uz", "fx fy fz");
    PROB.Set ("USigCond", (double)PROB.Keys.Size());
    return 0;
}    
 
// Call register
int __USigCond_dummy_int  = USigCondRegister();
 
}; // namespace FEM

#endif // MECHSYS_FEM_USIGCONDELEM_H
