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

#ifndef MECHSYS_FEM_FLOWELEM_H
#define MECHSYS_FEM_FLOWELEM_H

// Std Lib
#include <map>

// MechSys
#include "fem/element.h"
#include "models/flowstate.h"
#include "models/flowupdate.h"

namespace FEM
{

class FlowElem : public Element
{
public:
    // Constructor
    FlowElem (int                  NDim,   ///< Space dimension
              Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
              Model        const * Mdl,    ///< Model
              SDPair       const & Prp,    ///< Properties
              SDPair       const & Ini,    ///< Initial values
              Array<Node*> const & Nodes); ///< Array with all nodes (used to set the connectivity)

    // Methods
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs);       ///< If setting body forces, IdxEdgeOrFace is ignored
    void ClrBCs      ();                                               ///< Clear BCs
    void CalcK       (Mat_t & K)                                const; ///< Stiffness matrix
    void CalcM       (Mat_t & M)                                const; ///< Mass matrix
	void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL)     const; ///< Update state at IPs
    void GetState    (SDPair & KeysVals, int IdxIP=-1)          const; ///< IdxIP<0 => At centroid

    // Internal methods
    void CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const; ///< Strain-displacement matrix. Coef: coefficient used during integration

    // Data
    bool                 HasConv; // Has convection ?
    std::map<int,double> Bry2h;   // Map: boundary (edge/face) ID ==> to convection coefficient (h)
    double               m;       // Coefficient for mass matrix
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline FlowElem::FlowElem (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes), HasConv(false)
{
    // check GE
    if (GE==NULL) throw new Fatal("FlowElem::FlowElem: GE (geometry element) must be defined");

    // check GTy
         if (NDim==2 && GTy==d2d_t) { /*OK*/ }
    else if (NDim==3 && GTy==d3d_t) { /*OK*/ }
    else throw new Fatal("FlowElem::FlowElem: For NDim=%d, geometry type (GTy) must be equal to 'd%dd'. GTy=%s is invalid",NDim,NDim,GTypeToStr(GTy).CStr());

    // properties
    m = (Prp.HasKey("m") ? Prp("m") : 0.0);

    // allocate and initialize state at each IP
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Sta.Push (new FlowState(NDim));
        Mdl->InitIvs (Ini, Sta[i]);
    }

    // allocate state at centroid
    Sta.Push (new FlowState(NDim));
    Mdl->InitIvs (Ini, Sta[Sta.Size()-1]);

    // set UKeys in parent element
    UKeys.Resize (1);
    UKeys = "H"; // total head, temperature, voltage, etc.

    // initialize DOFs
    for (size_t i=0; i<Con.Size(); ++i) Con[i]->AddDOF("H", "Q");
}

inline void FlowElem::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs)
{
    if (BCs.HasKey("s")) // prescribed source term
    {
        double detJ, coef;
        double s = BCs("s");
        Mat_t C;
        CoordMatrix (C);
        //Vec_t Fs(Con.Size()); set_to_zero(Fs);
        for (size_t i=0; i<GE->NIP; ++i)
        {
            CalcShape (C, GE->IPs[i], detJ, coef);
            for (size_t j=0; j<GE->NN; ++j)
            {
                Con[j]->DF[Con[j]->FMap("Q")] += s*coef*GE->N(j);
                //Fs(j) += coef*GE->N(j)*s;
            }
        }
        //std::cout << s << "   Fs(" << Cell.ID << ") = " << PrintVector(Fs);
    }
    else
    {
        if (BCs.HasKey("H")) // prescribed potential (total head, temperature, ...)
        {
            for (size_t i=0; i<GE->NFN; ++i)
            {
                Node & nod = (*Con[GE->FNode(IdxEdgeOrFace,i)]);
                nod.SetBCs (BCs);
            }
        }
        else if (BCs.HasKey("flux"))
        {
            // flux value
            double qn = BCs("flux");

            // add to dF
            double detJ, coef;
            Mat_t FC;
            FCoordMatrix (IdxEdgeOrFace, FC);
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (FC, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    Node & nod = (*Con[GE->FNode(IdxEdgeOrFace,j)]);
                    nod.DF[nod.FMap("Q")] += coef*GE->FN(j)*qn;
                }
            }
        }
        else if (BCs.HasKey("conv"))
        {
            // data
            double h    = BCs("h");    // convection coefficient
            double Tinf = BCs("Tinf"); // temperature of surrounding environment
            HasConv     = true;        // has convection
            Bry2h[IdxEdgeOrFace] = h;  // Map: bry => h

            // add to dF
            double detJ, coef;
            Mat_t FC;
            FCoordMatrix (IdxEdgeOrFace, FC);
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (FC, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    Node & nod = (*Con[GE->FNode(IdxEdgeOrFace,j)]);
                    nod.DF[nod.FMap("Q")] += coef*GE->FN(j)*h*Tinf;
                }
            }
        }
    }
}

inline void FlowElem::ClrBCs ()
{
    HasConv = false;
    Bry2h.clear();
}

inline void FlowElem::CalcK (Mat_t & K) const
{
    double detJ, coef;
    Mat_t C, D, B;
    int nrows = Con.Size(); // number of rows in local K matrix
    K.change_dim (nrows,nrows);
    set_to_zero  (K);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        Mdl->Stiffness (Sta[i], D);
        CalcB (C, GE->IPs[i], B, detJ, coef);
        Mat_t BtDB(trans(B)*D*B);
        K += coef * BtDB;
    }
    if (HasConv)
    {
        for (std::map<int,double>::const_iterator p=Bry2h.begin(); p!=Bry2h.end(); ++p)
        {
            int    idx_bry = p->first;
            double h       = p->second;

            // add to K
            Mat_t FC;
            FCoordMatrix (idx_bry, FC);
            //Mat_t Kh(nrows,nrows);  set_to_zero(Kh);
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (FC, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    int row = GE->FNode(idx_bry,j);
                    for (size_t k=0; k<GE->NFN; ++k)
                    {
                        int col = GE->FNode(idx_bry,k);
                        K(row,col) += h*coef * GE->FN(j)*GE->FN(k);
                        //Kh(row,col) += h*coef * GE->FN(j)*GE->FN(k);
                    }
                }
            }
            //std::cout << "Kh = \n" << PrintMatrix(Kh);
        }
        
    }
    //std::cout << "D = \n" << PrintMatrix(D);
    //std::cout << "K = \n" << PrintMatrix(K);
}

inline void FlowElem::CalcM (Mat_t & M) const
{
    double detJ, coef;
    Mat_t C;
    int nrows = Con.Size(); // number of rows in local K matrix
    M.change_dim (nrows,nrows);
    set_to_zero  (M);
    CoordMatrix  (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        CalcShape (C, GE->IPs[i], detJ, coef);
        for (size_t j=0; j<GE->NN; ++j)
        {
            for (size_t k=0; k<GE->NN; ++k)
            {
                M(j,k) += m*coef*GE->N(j)*GE->N(k);
            }
        }
    }
}

inline void FlowElem::CalcB (Mat_t const & C, IntegPoint const & IP, Mat_t & B, double & detJ, double & Coef) const
{
    // deriv of shape func w.r.t natural coordinates
    GE->Derivs (IP.r, IP.s, IP.t);

    // Jacobian and its determinant
    Mat_t J(GE->dNdR * C); // J = dNdR * C
    detJ = Det(J);

    //std::cout << "J = \n" << PrintMatrix(J);
    //std::cout << "detJ = " << detJ << "\n";

    // deriv of shape func w.r.t real coordinates
    Mat_t Ji;
    Inv (J, Ji);

    // coefficient used during integration
    Coef = detJ*IP.w;

    // B matrix
    int nrows = Con.Size(); // number of rows in local K matrix
    B.change_dim (NDim,nrows);
    B = Ji * GE->dNdR; // B = dNdX = Inv(J) * dNdR
}

inline void FlowElem::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal potential
    int nrows = Con.Size(); // number of rows in local K matrix
    Vec_t dUe(nrows);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // update state at each IP
    FlowUpdate fu(Mdl);
    double detJ, coef;
    Mat_t  C, B;
    Vec_t  dFe(nrows), dvel(NDim), dgra(NDim);
    set_to_zero (dFe);
    CoordMatrix (C);
    for (size_t i=0; i<GE->NIP; ++i)
    {
        // B matrix
        CalcB (C, GE->IPs[i], B, detJ, coef);

        // velocity and gradient increments
        dgra = B * dUe;
        fu.Update (dgra, Sta[i], dvel);

        // element nodal forces
        Vec_t Btdvel(trans(B)*dvel);
        dFe -= coef * (Btdvel); // '-=' because dvel is -k*dgra
    }

    // add contribution to dFe due to convection term
    if (HasConv)
    {
        for (std::map<int,double>::const_iterator p=Bry2h.begin(); p!=Bry2h.end(); ++p)
        {
            int    idx_bry = p->first;
            double h       = p->second;

            // add to dFe
            Mat_t FC;
            FCoordMatrix (idx_bry, FC);
            for (size_t i=0; i<GE->NFIP; ++i)
            {
                CalcFaceShape (FC, GE->FIPs[i], detJ, coef);
                for (size_t j=0; j<GE->NFN; ++j)
                {
                    int row = GE->FNode(idx_bry,j);
                    for (size_t k=0; k<GE->NFN; ++k)
                    {
                        int col = GE->FNode(idx_bry,k);
                        dFe(row) += h*coef * GE->FN(j)*GE->FN(k) * dUe(col);
                    }
                }
            }
        }
    }

    // update state at centroid
    CalcB (C, GE->Rct, B, detJ, coef);
    dgra = B * dUe;
    fu.Update (dgra, Sta[Sta.Size()-1], dvel);

    // add results to Fint (internal forces)
    if (F_int!=NULL) for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
}

inline void FlowElem::GetState (SDPair & KeysVals, int IdxIP) const
{
    Vec_t const * vel;
    Vec_t const * gra;
    if (IdxIP<0) // centroid
    {
        vel = &(static_cast<FlowState const *>(Sta[Sta.Size()-1])->Vel);
        gra = &(static_cast<FlowState const *>(Sta[Sta.Size()-1])->Gra);
    }
    else
    {
        vel = &(static_cast<FlowState const *>(Sta[IdxIP])->Vel);
        gra = &(static_cast<FlowState const *>(Sta[IdxIP])->Gra);
    }
    if (NDim==2)
    {
        KeysVals.Set("vx vy  gx gy",
                     (*vel)(0),(*vel)(1),
                     (*gra)(0),(*gra)(1));
    }
    else
    {
        KeysVals.Set("vx vy vz  gx gy gz",
                     (*vel)(0),(*vel)(1),(*vel)(2),
                     (*gra)(0),(*gra)(1),(*gra)(2));
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * FlowElemMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new FlowElem(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int FlowElemRegister()
{
    ElementFactory["Flow"] = FlowElemMaker;
    PROB.Set ("Flow", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __FlowElem_dummy_int  = FlowElemRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_FLOWELEM
