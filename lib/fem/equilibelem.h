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
#include <mechsys/fem/element.h>
#include <mechsys/models/equilibstate.h>
#include <mechsys/models/stressupdate.h>

namespace FEM
{

class EquilibElem : public Element
{
public:
    // Static
    static size_t NCo; ///< Number of stress/strain components == 2*NDim
    static size_t NDu; ///< Number of DOFs (displacements) == NN*NDim

    // Constructor & Destructor
    EquilibElem (int                  NDim,   ///< Space dimension
                 Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                 Model        const * Mdl,    ///< Model
                 Model        const * XMdl,   ///< Extra Model
                 SDPair       const & Prp,    ///< Properties
                 SDPair       const & Ini,    ///< Initial values
                 Array<Node*> const & Nodes); ///< Connectivity
    virtual ~EquilibElem () {}

    // Methods
    virtual void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF); ///< If setting body forces, IdxEdgeOrFace is ignored
    virtual void GetLoc      (Array<size_t> & Loc)                               const; ///< Get location vector for mounting K/M matrices
    void         CalcK       (Mat_t & K)                                         const; ///< Stiffness matrix
    void         CalcM       (Mat_t & M)                                         const; ///< Mass matrix
    virtual void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL)              const; ///< Update state at IPs
    virtual void SetFint     (Vec_t * Fint=NULL)                                 const; ///< Set Fint=K*Ue or set nodes if Fint==NULL
    virtual void StateKeys   (Array<String> & Keys)                              const; ///< Get state keys, ex: sx, sy, sxy, ex, ey, exy
    virtual void StateAtIP   (SDPair & KeysVals, int IdxIP)                      const; ///< Get state at IP

    // Internal methods
    void CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const; ///< Strain-displacement matrix. Coef: coefficient used during integration
    void CalcN (Mat_t const & C, IntegPoint const & IP, Mat_t & N, double & detJ, double & Coef) const; ///< Shape functions matrix

    // Methods for the Runge-Kutta method
    virtual size_t NIVs       ()                                                            const;
    virtual double GetIV      (size_t i)                                                    const;
    virtual void   SetIV      (size_t i, double Val);
    virtual void   CalcIVRate (double Time, Vec_t const & U, Vec_t const & V, Vec_t & Rate) const;
    virtual void   CorrectIVs ();

    // Constants
    double h;   ///< Thickness of the element
    double rho; ///< Density
};

size_t EquilibElem::NCo = 0;
size_t EquilibElem::NDu = 0;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline EquilibElem::EquilibElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes)
{
    // check
    if (GE==NULL)  throw new Fatal("EquilibElem::EquilibElem: GE (geometry element) must be defined");
    if (Mdl==NULL) throw new Fatal("EquilibElem::EquilibElem: Model must be defined");

    // set constants of this class (just once)
    if (NDu==0)
    {
        NCo = 2*NDim;
        NDu = NDim*GE->NN;
    }

    // parameters/properties
    h   = (Prp.HasKey("h")   ? Prp("h")   : 1.0);
    rho = Mdl->Rho;//(Prp.HasKey("rho") ? Prp("rho") : 1.0);
    if (h<1.0e-8)   throw new Fatal("EquilibElem::EquilibElem: The thickness of the element must be greater than 1.0e-8. h=%g is invalid",h);
    if (rho<1.0e-8) throw new Fatal("EquilibElem::EquilibElem: 'rho' must be greater than zero (rho=%g is invalid)",rho);

    // allocate and initialize state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Sta.Push (new EquilibState(NDim));
        Mdl->InitIvs (Ini, Sta[i]);
    }

    // set initial values
    bool geosta = (Prp.HasKey("geosta") ? Prp("geosta")>0 : false);
    if (geosta)
    {
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
    }

    // set F in nodes due to initial stresses
    SetFint ();
}

inline void EquilibElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, BCFuncs * BCF)
{
    // deactivate/activate commands
    if (BCs.HasKey("deactivate"))
    {
        // check
        if (!Active) throw new Fatal("EquilibElem::SetBCs: 'deactivate' command failed: Element # %d is already inactive",Cell.ID);

        // calc force to be applied after removal of element
        double gra = (BCs.HasKey("gravity") ? BCs("gravity") : 0.0);
        double detJ, coef;
        Mat_t  C, B;
        Vec_t  Fi(NDu), Fb(NDu);
        set_to_zero (Fi);
        set_to_zero (Fb);
        CoordMatrix (C);
        size_t idx_grav = (NDim==2 ? 1 : 2);
        for (size_t i=0; i<GE->NIP; ++i)
        {
            // internal forces
            GE->Shape (GE->IPs[i].r, GE->IPs[i].s, GE->IPs[i].t);
            CalcB     (C, GE->IPs[i], B, detJ, coef);
            Vec_t const & sig = static_cast<EquilibState const *>(Sta[i])->Sig;
            Fi += (coef) * trans(B)*sig;

            // body forces
            for (size_t j=0; j<GE->NN; ++j) Fb(idx_grav+j*NDim) += -coef*GE->N(j)*rho*gra;
        }

        // set nodes
        for (size_t i=0; i<GE->NN; ++i)
        {
            // remove sharing information
            Con[i]->NShares--;
            if (Con[i]->NShares<0) throw new Fatal("EquilibElem::SetBCs: __internal_error__: 'deactivate' command failed: NShares==%d must be positive",Con[i]->NShares<0);

            // clear all values in node in case it becomes inactive
            if (Con[i]->NShares==0) Con[i]->Clear();
            else
            {
                // set boundary conditions
                Con[i]->AddToPF("fx", Fi(0+i*NDim) - Fb(0+i*NDim), BCF);
                Con[i]->AddToPF("fy", Fi(1+i*NDim) - Fb(1+i*NDim), BCF);  if (NDim==3)
                Con[i]->AddToPF("fz", Fi(2+i*NDim) - Fb(2+i*NDim), BCF);
            }
        }

        // clear state at IPs
        SDPair inis;
        inis.Set("zero",1.);
        for (size_t i=0; i<GE->NIP; ++i) Sta[i]->Init (inis, size(static_cast<EquilibState const *>(Sta[i])->Ivs));

        // deactivate element
        Active = false;
        return; // must return otherwise 'gravity' may be set in the following code
    }
    else if (BCs.HasKey("activate"))
    {
        // check
        if (Active) throw new Fatal("EquilibElem::SetBCs: 'activate' command failed: Element # %d is already active",Cell.ID);

        // add information to shares array in nodes
        for (size_t i=0; i<GE->NN; ++i) Con[i]->NShares++;

        // clear state at IPs
        SDPair inis;
        inis.Set("zero",1.);
        for (size_t i=0; i<GE->NIP; ++i) Sta[i]->Init (inis, size(static_cast<EquilibState const *>(Sta[i])->Ivs));

        // activate
        Active = true;
        // continue to set gravity
    }

    // check
    if (!Active) throw new Fatal("EquilibElem::SetBCs: Element %d is inactive",Cell.ID);

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
    bool has_gra = BCs.HasKey("gravity"); // gravity

    // force components specified
    if (has_bx || has_by || has_bz || has_cbx || has_gra ||
        has_qx || has_qy || has_qz || has_qn  || has_qt)
    {
        // body forces
        if (has_bx || has_by || has_bz || has_cbx || has_gra) // prescribed body forces
        {
            // matrix of coordinates of nodes
            Mat_t C;
            CoordMatrix (C);
            
            // loading
            double bx = (has_bx  ? BCs("bx")  : 0.0);
            double by = (has_by  ? BCs("by")  : 0.0);
            double bz = (has_bz  ? BCs("bz")  : 0.0);
                   bx = (has_cbx ? BCs("cbx") : bx );
            if (has_gra)
            {
                if (NDim==2) by = -rho*BCs("gravity");
                else         bz = -rho*BCs("gravity");
            }

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
                double coef = h*GE->FIPs[i].w; // *detJ is not necessary since qx,qy,qz are already multiplied by detJ (due to normal)

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

inline void EquilibElem::GetLoc (Array<size_t> & Loc) const
{
    Loc.Resize (NDu);
    for (size_t i=0; i<GE->NN; ++i)
    {
        Loc[i*NDim+0] = Con[i]->Eq("ux");
        Loc[i*NDim+1] = Con[i]->Eq("uy");  if (NDim==3)
        Loc[i*NDim+2] = Con[i]->Eq("uz");
    }
}

inline void EquilibElem::CalcK (Mat_t & K) const
{
    double detJ, coef;
    Mat_t C, D, B;
    K.change_dim (NDu, NDu);
    set_to_zero  (K);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Mdl->Stiffness (Sta[i], D);
        CalcB (C, GE->IPs[i], B, detJ, coef);
        K += (coef) * trans(B)*D*B;
    }
}

inline void EquilibElem::CalcM (Mat_t & M) const
{
    double detJ, coef;
    Mat_t  C, N;
    M.change_dim (NDu, NDu);
    set_to_zero  (M);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        CalcN (C, GE->IPs[i], N, detJ, coef);
        M += (rho*coef) * trans(N)*N;
    }
}

inline void EquilibElem::CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const
{
    // deriv of shape func w.r.t natural coordinates
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



    //std::cout << "J = \n" << PrintMatrix(J);
    //std::cout << "dNdR = \n" << PrintMatrix(GE->dNdR);
    //std::cout << "C = \n" << PrintMatrix(C);
    //std::cout << "dNdX = \n" << PrintMatrix(dNdX);
    //printf("detJ = %g\n", detJ);


    // B matrix
    B.change_dim (NCo, NDu);
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
    N.change_dim (NDim, NDu);
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
    Vec_t dUe(NDu);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // update state at each IP
    double detJ, coef;
    Mat_t  C, B;
    Vec_t  dFe(NDu), dsig(NCo), deps(NCo);
    set_to_zero (dFe);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // B matrix
        CalcB (C, GE->IPs[i], B, detJ, coef);

        // strain and stress increments
        deps = B * dUe;
        std::cout << "deps = " << PrintVector(deps);
        Mdl->SUp.Update (deps, Sta[i], dsig);

        // element nodal forces
        dFe += (coef) * trans(B)*dsig;
    }

    // add results to Fint (internal forces)
    if (F_int!=NULL) for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
}

inline void EquilibElem::SetFint (Vec_t * Fint) const
{
    // element force
    double detJ, coef;
    Mat_t  C, B;
    Vec_t  Fe(NDu);
    CoordMatrix (C);
    set_to_zero (Fe);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        CalcB (C, GE->IPs[i], B, detJ, coef);
        Vec_t const & sig = static_cast<EquilibState const *>(Sta[i])->Sig;
        Fe += (coef) * trans(B)*sig;
    }

    // set nodes
    if (Fint==NULL)
    {
        for (size_t i=0; i<GE->NN; ++i)
        {
            Con[i]->F("fx") += Fe(0+i*NDim);
            Con[i]->F("fy") += Fe(1+i*NDim);  if (NDim==3)
            Con[i]->F("fz") += Fe(2+i*NDim);
        }
    }

    // return Fint
    else
    {
        Array<size_t> loc;
        GetLoc (loc);
        for (size_t i=0; i<loc.Size(); ++i) (*Fint)(loc[i]) += Fe(i);
    }
}

inline void EquilibElem::StateKeys (Array<String> & Keys) const
{
    Keys = EquilibState::Keys;
    for (size_t i=0; i<Mdl->NIvs; ++i) Keys.Push (Mdl->IvNames[i]);
}

inline void EquilibElem::StateAtIP (SDPair & KeysVals, int IdxIP) const
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

inline size_t EquilibElem::NIVs () const 
{
    size_t niv = size(static_cast<EquilibState const*>(Sta[0])->Ivs);
    return (NCo+niv)*GE->NIP;
}

inline double EquilibElem::GetIV (size_t i) const 
{
    size_t niv = size(static_cast<EquilibState const*>(Sta[0])->Ivs);
    size_t iip = static_cast<int>(i) / static_cast<int>(NCo+niv); // index of IP
    size_t ico = static_cast<int>(i) % static_cast<int>(NCo+niv); // index of component/iv
    if (ico<NCo) return static_cast<EquilibState const*>(Sta[iip])->Sig(ico);
    else         return static_cast<EquilibState const*>(Sta[iip])->Ivs(ico-NCo);
}

inline void EquilibElem::SetIV (size_t i, double Val) 
{
    size_t niv = size(static_cast<EquilibState const*>(Sta[0])->Ivs);
    size_t iip = static_cast<int>(i) / static_cast<int>(NCo+niv); // index of IP
    size_t ico = static_cast<int>(i) % static_cast<int>(NCo+niv); // index of component
    if (ico<NCo) static_cast<EquilibState*>(Sta[iip])->Sig(ico)     = Val;
    else         static_cast<EquilibState*>(Sta[iip])->Ivs(ico-NCo) = Val;
}

inline void EquilibElem::CalcIVRate (double Time, Vec_t const & U, Vec_t const & V, Vec_t & Rate) const
{
    // resize rate
    size_t niv = size(static_cast<EquilibState const*>(Sta[0])->Ivs);
    Rate.change_dim ((NCo+niv)*GE->NIP);

    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal velocities
    Vec_t Ve(NDu);
    for (size_t i=0; i<loc.Size(); ++i) Ve(i) = V(loc[i]);

    // calc dsigdt
    double detJ, coef;
    Mat_t  C, B, D;
    Vec_t  depsdt(NCo), dsigdt(NCo), divsdt(niv);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        CalcB (C, GE->IPs[i], B, detJ, coef);
        depsdt = B * Ve;
        Mdl->Rate (Sta[i], depsdt,  dsigdt, divsdt);
        for (size_t j=0; j<NCo; ++j) Rate(    j+i*(NCo+niv)) = dsigdt(j);
        for (size_t j=0; j<niv; ++j) Rate(NCo+j+i*(NCo+niv)) = divsdt(j);
    }
}

inline void EquilibElem::CorrectIVs ()
{
    if (Mdl->SUp.CDrift)
    {
        for (size_t i=0; i<GE->NIP; ++i) Mdl->CorrectDrift (Sta[i]);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * EquilibElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new EquilibElem(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes); }

// Register element
int EquilibElemRegister()
{
    ElementFactory  ["Equilib"]   = EquilibElemMaker;
    ElementVarKeys  ["Equilib2D"] = std::make_pair ("ux uy",    "fx fy");
    ElementVarKeys  ["Equilib3D"] = std::make_pair ("ux uy uz", "fx fy fz");
    ElementExtraKeys["Equilib2D"] = Array<String>  ("activate", "deactivate", "grav", "qn", "qt");
    ElementExtraKeys["Equilib3D"] = Array<String>  ("activate", "deactivate", "grav", "qn");
    PROB.Set ("Equilib", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __EquilibElem_dummy_int  = EquilibElemRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIBELEM
