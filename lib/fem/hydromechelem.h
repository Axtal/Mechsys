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

    // Constructor
    HydroMechElem (int                  NDim,   ///< Space dimension
                   Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                   Model        const * Mdl,    ///< Model
                   SDPair       const & Prp,    ///< Properties
                   SDPair       const & Ini,    ///< Initial values
                   Array<Node*> const & Nodes); ///< Connectivity

    // Methods
    void CalcKCM     (Mat_t & KK, Mat_t & CC, Mat_t & MM)   const;
    void GetLoc      (Array<size_t> & Loc)                  const;
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL) const;

    // Internal Methods
    void Interp (Mat_t const & C, IntegPoint const & IP, Mat_t & B, Mat_t & Bp, Mat_t & N, Mat_t & Np, double & detJ, double & Coef) const; ///< Interpolation matrices

    // Data
    double cval;
    double Cval;
    double chi;
    Mat_t  kb; // == k/gammaW
};

size_t HydroMechElem::NDp = 0;
size_t HydroMechElem::NDt = 0;
Mat_t  HydroMechElem::Im;
Vec_t  HydroMechElem::Iv;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline HydroMechElem::HydroMechElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : EquilibElem(NDim,Cell,Mdl,Prp,Ini,Nodes)
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
    }

    // set initial pw values at Nodes
    if (Ini.HasKey("geostatic"))
    {
        if (!Ini.HasKey("gamW"))               throw new Fatal("HydroMechElem::HydroMechElem: For geostatic stresses, 'gamW' must be provided in 'Ini' dictionary");
        if (NDim==2 && !Ini.HasKey("y_water")) throw new Fatal("HydroMechElem::HydroMechElem: For geostatic stresses in 2D, 'y_water' must be provided in 'Ini' dictionary");
        if (NDim==3 && !Ini.HasKey("z_water")) throw new Fatal("HydroMechElem::HydroMechElem: For geostatic stresses in 3D, 'z_water' must be provided in 'Ini' dictionary");
        bool   pos_pw  = (Ini.HasKey("only_positive_pw") ? true : false);
        double gamW    = Ini("gamW");
        double z_water = (NDim==2 ? Ini("y_water") : Ini("z_water"));
        for (size_t i=0; i<GE->NN; ++i)
        {
            // elevation of node
            double z = (NDim==2 ? Con[i]->Vert.C(1) : Con[i]->Vert.C(2));

            // pore-water pressure
            double hw = z_water-z; // column of water
            double pw = (hw>0.0 ? gamW*hw : (pos_pw ? 0.0 : gamW*hw));

            // set U in node
            Con[i]->U("pw") = pw;
        }
    }

    // hydraulic conductivity tensor
    kb.change_dim(NDim,NDim);
    set_to_zero (kb);
    kb(0,0)=1;  kb(1,1)=1;  if (NDim==3) kb(2,2)=1;
    double gammaW = 1.;
    kb  /= gammaW;
    cval = 1.0;
    Cval = 0.0;
    chi  = 1.0;
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
        Mdl->Stiffness (Sta[i], D);
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);
        NtN    = trans(N)*N;
        BtDB   = trans(B)*D*B;
        BtmNp  = trans(B)*Im*Np;
        NptNp  = trans(Np)*Np;
        BptkBp = trans(Bp)*kb*Bp;
        M     += (coef*rho)  * NtN;
        K     += (coef)      * BtDB;
        Q     += (coef*cval) * BtmNp;
        Qb    += (coef*chi)  * BtmNp;
        S     += (coef*Cval) * NptNp;
        H     += (coef)      * BptkBp;
    }

    // assemble
    MM.change_dim (NDt,NDt);
    CC.change_dim (NDt,NDt);
    KK.change_dim (NDt,NDt);
    set_to_zero (MM);
    set_to_zero (CC);
    set_to_zero (KK);
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
            CC(NDu+i,NDu+j) = S(i,j);
            KK(NDu+i,NDu+j) = H(i,j);
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
    StressUpdate su(Mdl);
    double detJ, coef;
    Mat_t  C, B, Bp, N, Np;
    Vec_t  dFe(NDu), dsig(NCo), deps(NCo), dfe(NDp);
    set_to_zero (dFe);
    set_to_zero (dfe);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // B matrix
        Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);

        // strain and effective stress increments
        deps = B * dUe;
        su.Update (deps, Sta[i], dsig);

        // total stress increments
        Vec_t Npdpwe(Np * dpwe);
        double dpw = Npdpwe(0);
        dsig -= chi*dpw*Iv;

        // element nodal forces
        kBp     = kb * Bp;
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
