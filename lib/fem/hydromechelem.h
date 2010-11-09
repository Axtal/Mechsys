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

    // Constructor
    HydroMechElem (int                  NDim,   ///< Space dimension
                   Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
                   Model        const * Mdl,    ///< Model
                   SDPair       const & Prp,    ///< Properties
                   SDPair       const & Ini,    ///< Initial values
                   Array<Node*> const & Nodes); ///< Connectivity

    // Destructor
    ~HydroMechElem () { if (GEp!=NULL) delete GEp; }

    // Methods
    void CalcKCM (Mat_t & KK, Mat_t & CC, Mat_t & MM) const;
    void GetLoc  (Array<size_t> & Loc)                const; ///< Get location vector for mounting KCM matrices

    // Data
    GeomElem * GEp; ///< Element for the pressure DOFs
};

size_t HydroMechElem::NDp = 0;
size_t HydroMechElem::NDt = 0;


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline HydroMechElem::HydroMechElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : EquilibElem(NDim,Cell,Mdl,Prp,Ini,Nodes)
{
    // allocate GEp
    if (GE->Name=="Hex20") GEp = AllocGeomElem ("Hex8", NDim);
    else throw new Fatal("HydroMechElem::HydroMechElem: Geom element==%s is not implemented yet",GE->Name.CStr());

    // set constants of this class (just once)
    if (NDp==0)
    {
        NDp = GEp->NN;
        NDt = NDu + NDp;
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
        for (size_t i=0; i<GEp->NN; ++i)
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
}

inline void HydroMechElem::GetLoc (Array<size_t> & Loc) const
{
    // u DOFs
    Loc.Resize (NDt);
    for (size_t i=0; i<GE->NN; ++i)
    {
        Loc[i*NDim+0] = Con[i]->Eq("ux");
        Loc[i*NDim+1] = Con[i]->Eq("uy");  if (NDim==3)
        Loc[i*NDim+2] = Con[i]->Eq("uz");
    }

    // pw DOFs
    for (size_t i=0; i<GEp->NN; ++i)
    {
        Loc[NDu+i] = Con[i]->Eq("pw");
    }
}

inline void HydroMechElem::CalcKCM (Mat_t & KK, Mat_t & CC, Mat_t & MM) const
{
    // hydraulic conductivity tensor
    Mat_t kb(NDim,NDim);
    set_to_zero (kb);
    kb(0,0)=1;  kb(1,1)=1;  if (NDim==3) kb(2,2)=1;
    double gammaW = 1.;
    kb /= gammaW;

    // mechanical matrices
    Mat_t M, K;
    EquilibElem::CalcM (M);
    EquilibElem::CalcK (K);

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
    Mat_t C;
    double detJ, coef;
    //Mat_t NtN    (NDu,NDu);
    //Mat_t BtDB   (NDu,NDu);
    //Mat_t BtINp  (NDu,NDp);
    //Mat_t BptkBp (NDp,NDp);
    //Mat_t NptNp  (NDp,NDp);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        //Mdl->Stiffness (Sta[i], D);
        //Interp (C, GE->IPs[i], B, Bp, N, Np, detJ, coef);
        //NtN    = trans(N)*N;
        //BtDB   = trans(B)*D*B;
        //BtINp  = trans(B)*Im*Np;
        //BptkBp = trans(Bp)*km*Bp;
        //NptNp  = trans(Np)*Np;
        //M     += (coef*rho) * NtN;
        //K     += (coef)     * BtDB;
        //L     += (coef*alp) * BtINp;
        //H     += (coef)     * BptkBp;
        //S     += (coef/Qs)  * NptNp;
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


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * HydroMechElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new HydroMechElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int HydroMechElemRegister()
{
    ElementFactory["HydroMech"]   = HydroMechElemMaker;
    ElementVarKeys["HydroMech2D"] = std::make_pair ("ux uy pw",    "fx fy qw");
    //ElementVarKeys["HydroMech3D"] = std::make_pair ("ux uy uz", "fx fy fz");
    ElementVarKeys["HydroMech3D"] = std::make_pair ("ux uy uz pw", "fx fy fz qw");
    PROB.Set ("HydroMech", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __HydroMechElem_dummy_int  = HydroMechElemRegister();

}; // namespace FEM

#endif // MECHSYS_HYDROMECH_H
