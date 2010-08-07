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

    // Constructor
    EquilibElem (int                  NDim,   ///< Space dimension
                 Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                 Model        const * Mdl,    ///< Model
                 SDPair       const & Prp,    ///< Properties
                 SDPair       const & Ini,    ///< Initial values
                 Array<Node*> const & Nodes); ///< Connectivity

    // Methods
    void SetBCs       (size_t IdxEdgeOrFace, SDPair const & BCs, 
                       NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM); ///< If setting body forces, IdxEdgeOrFace is ignored
    void Gravity      (NodBCs_t & pF, pCalcM CalcM, double gAccel); ///< Apply gravity
    void Deactivate   (NodBCs_t & pF, pCalcM CalcM, double gAccel,
                       NodBCs_t & pU);                              ///< Deactivate element
    void CalcK        (Mat_t & K)                            const; ///< Stiffness matrix
    void CalcM        (Mat_t & M)                            const; ///< Mass matrix
    void UpdateState  (Vec_t const & dU, Vec_t * F_int=NULL) const; ///< Update state at IPs
    void StateKeys    (Array<String> & Keys)                 const; ///< Get state keys, ex: sx, sy, sxy, ex, ey, exy
    void StateAtIP    (SDPair & KeysVals, int IdxIP)         const; ///< Get state at IP

    // Internal methods
    void CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const; ///< Strain-displacement matrix. Coef: coefficient used during integration
    void CalcN (Mat_t const & C, IntegPoint const & IP, Mat_t & N, double & detJ, double & Coef) const; ///< Shape functions matrix

    // Constants
    double h;   ///< Thickness of the element
    double rho; ///< Density
};

size_t EquilibElem::NCo = 0;
size_t EquilibElem::NDu = 0;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline EquilibElem::EquilibElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes)
{
    // check GE
    if (GE==NULL) throw new Fatal("EquilibElem::EquilibElem: GE (geometry element) must be defined");

    // set constants of this class (just once)
    if (NDu==0)
    {
        NCo = 2*NDim;
        NDu = NDim*GE->NN;
    }

    // parameters/properties
    h   = (Prp.HasKey("h")   ? Prp("h")   : 1.0);
    rho = (Prp.HasKey("rho") ? Prp("rho") : 1.0);
    if (h<1.0e-8) throw new Fatal("EquilibElem::EquilibElem: The thickness of the element must be greater than 1.0e-8. h=%g is invalid",h);

    // allocate and initialize state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Sta.Push (new EquilibState(NDim));
        Mdl->InitIvs (Ini, Sta[i]);
    }

    // set initial values
    if (Ini.HasKey("geostatic"))
    {
        if (GTy==pse_t)                       throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, geometry cannot be of 'plane-stress' (pse) type");
        if (!Ini.HasKey("K0"))                throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, 'K0' must be provided in 'Ini' dictionary");
        if (!Ini.HasKey("gam"))               throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses, 'gam' must be provided in 'Ini' dictionary");
        if (NDim==2 && !Ini.HasKey("y_surf")) throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses in 2D, 'y_surf' must be provided in 'Ini' dictionary");
        if (NDim==3 && !Ini.HasKey("z_surf")) throw new Fatal("EquilibElem::EquilibElem: For geostatic stresses in 3D, 'z_surf' must be provided in 'Ini' dictionary");
        double K0   = Ini("K0");
        double gam  = Ini("gam");
        double surf = (NDim==2 ? Ini("y_surf") : Ini("z_surf"));
        Vec_t X; // coords of IP
        for (size_t i=0; i<GE->NIP; ++i)
        {
            CoordsOfIP (i, X);
            Vec_t & sig = static_cast<EquilibState *>(Sta[i])->Sig;
            if (NDim==2)
            {
                double h = fabs(surf-X(1));
                sig(1) = -gam*h;    // sy
                sig(0) = K0*sig(1); // sx
                sig(2) = K0*sig(1); // sz
                sig(3) = 0.0;       // sxy*sq2
            }
            else // 3D
            {
                double h = fabs(surf-X(2));
                sig(2) = -gam*h;    // sz
                sig(0) = K0*sig(2); // sx
                sig(1) = K0*sig(2); // sy
                sig(3) = 0.0;       // sxy*sq2
                sig(4) = 0.0;       // syz*sq2
                sig(5) = 0.0;       // szx*sq2
            }
        }
    }

    // set UKeys in parent element and initialize DOFs
    if (NDim==2)
    {
        UKeys.Resize (NDim);
        UKeys = "ux", "uy";
        for (size_t i=0; i<GE->NN; ++i) Con[i]->AddDOF("ux uy", "fx fy");
    }
    else // 3D
    {
        UKeys.Resize (NDim);
        UKeys = "ux", "uy", "uz";
        for (size_t i=0; i<GE->NN; ++i) Con[i]->AddDOF("ux uy uz", "fx fy fz");
    }

    // set F in nodes due to initial stresses
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
    for (size_t i=0; i<GE->NN; ++i)
    {
        Con[i]->F[Con[i]->FMap("fx")] += Fe(0+i*NDim);
        Con[i]->F[Con[i]->FMap("fy")] += Fe(1+i*NDim);  if (NDim==3)
        Con[i]->F[Con[i]->FMap("fz")] += Fe(2+i*NDim);
    }
}

inline void EquilibElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, NodBCs_t & pF, NodBCs_t & pU, pCalcM CalcM)
{
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

inline void EquilibElem::Gravity (NodBCs_t & pF, pCalcM CalcM, double gAccel)
{
    // check
    if (!Active) throw new Fatal("EquilibElem::Gravity: Element %d is inactive",Cell.ID);

    // matrix of coordinates of nodes
    Mat_t C;
    CoordMatrix (C);
    
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
            coef *= radius;
        }

        // add to dF
        for (size_t j=0; j<GE->NN; ++j)
        {
            if (NDim==2) pF[Con[j]].first[Con[j]->FMap("fy")] += -coef*GE->N(j)*rho*gAccel;
            else         pF[Con[j]].first[Con[j]->FMap("fz")] += -coef*GE->N(j)*rho*gAccel;
        }
    }

    // set CalcM
    for (size_t j=0; j<GE->NN; ++j) pF[Con[j]].second = CalcM;
}

inline void EquilibElem::Deactivate (NodBCs_t & pF, pCalcM CalcM, double gAccel, NodBCs_t & pU)
{
    // check
    if (!Active) throw new Fatal("EquilibElem::Deactivate: Element %d is already inactive",Cell.ID);

    // calc force to be applied after removal of element
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
        for (size_t j=0; j<GE->NN; ++j)
            Fb(idx_grav+j*NDim) += coef*GE->N(j)*rho*gAccel;
    }

    /*
    for (size_t j=0; j<GE->NN; ++j) // debug
    {
        std::cout << Util::_4 << Con[j]->Vert.ID << ": ";
        for (size_t k=0; k<NDim; ++k) std::cout << Util::_6_3 << Fi(k+j*NDim);
        std::cout << "  :   ";
        for (size_t k=0; k<NDim; ++k) std::cout << Util::_6_3 << Fb(k+j*NDim);
        std::cout << std::endl;
    }
    */

    // set nodes
    for (size_t i=0; i<GE->NN; ++i)
    {
        // remove internal force contribution from external forces vector
        Con[i]->F[Con[i]->FMap("fx")] -= Fi(0+i*NDim);
        Con[i]->F[Con[i]->FMap("fy")] -= Fi(1+i*NDim);  if (NDim==3)
        Con[i]->F[Con[i]->FMap("fz")] -= Fi(2+i*NDim);

        // remove sharing information
        Con[i]->NShares--;
        if (Con[i]->NShares<0) throw new Fatal("EquilibElem::Deactivate: __internal_error__: NShares==%d must be positive",Con[i]->NShares<0);

        // delete node from pF and pU in case node becomes inactive
        if (Con[i]->NShares==0)
        {
            pF.erase (Con[i]);
            pU.erase (Con[i]);
        }
        else // set boundary conditions (prescribed F: pF)
        {
            // add to F
            pF[Con[i]].first[Con[i]->FMap("fx")] += Fi(0+i*NDim) + Fb(0+i*NDim);
            pF[Con[i]].first[Con[i]->FMap("fy")] += Fi(1+i*NDim) + Fb(1+i*NDim);  if (NDim==3)
            pF[Con[i]].first[Con[i]->FMap("fz")] += Fi(2+i*NDim) + Fb(2+i*NDim);

            // set CalcM (multiplier callback)
            pF[Con[i]].second = CalcM;
        }
    }

    // deactivate element
    Active = false;
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
    StressUpdate su(Mdl);
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
        su.Update (deps, Sta[i], dsig);

        // element nodal forces
        dFe += (coef) * trans(B)*dsig;
    }

    // add results to Fint (internal forces)
    if (F_int!=NULL) for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
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


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * EquilibElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new EquilibElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int EquilibElemRegister()
{
    ElementFactory["Equilib"]   = EquilibElemMaker;
    ElementVarKeys["Equilib2D"] = std::make_pair ("ux uy",    "fx fy");
    ElementVarKeys["Equilib3D"] = std::make_pair ("ux uy uz", "fx fy fz");
    PROB.Set ("Equilib", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __EquilibElem_dummy_int  = EquilibElemRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_EQUILIBELEM
