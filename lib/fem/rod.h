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

#ifndef MECHSYS_FEM_ROD_H
#define MECHSYS_FEM_ROD_H

// MechSys
#include "fem/element.h"
#include "models/model.h"

using std::cout;
using std::endl;

namespace FEM
{

class RodState : public State
{
public:
    RodState (int NDim) : State(NDim) {}
    void Init    (SDPair const & Ini) { ea = (Ini.HasKey("ea") ? Ini("ea") : 0.0); }
    void Backup  () { ea_bkp = ea; }
    void Restore () { ea = ea_bkp; }
    double ea, ea_bkp; ///< Axial strain
};

class Rod : public Element
{
public:
    // Constructor
    Rod (int                  NDim,   ///< Space dimension
         Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
         Model        const * Mdl,    ///< Model
         SDPair       const & Prp,    ///< Properties
         SDPair       const & Ini,    ///< Initial values
         Array<Node*> const & Nodes); ///< Array with all nodes (used to set the connectivity)

    // Methods
    void CalcK       (Mat_t & K)                            const; ///< Stiffness matrix
    void CalcT       (Mat_t & T, double & l)                const; ///< Transformation matrix
	void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL) const;
    void GetState    (SDPair & KeysVals, int none=-1)       const;
    void Centroid    (Vec_t & X)                            const; ///< Centroid of element

    // Constants
    double E; ///< Young modulus
    double A; ///< Cross-sectional area
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Rod::Rod (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes)
{
    // check GTy
    if (GTy!=fra_t) throw new Fatal("Rod::Rod: Geometry type (GTy) must be equal to 'fra' (Frame). GTy=%s is invalid",GTypeToStr(GTy).CStr());

    // parameters/properties
    E = Prp("E");
    A = Prp("A");

    // allocate and initialize state (constant along element)
    Sta.Push (new RodState(NDim));
    Sta[0]->Init (Ini);
    
    // set UKeys in parent element
    UKeys.Resize (NDim);
    if (NDim==2) UKeys = "ux", "uy";
    else         UKeys = "ux", "uy", "uz";

    // initialize DOFs
    if (NDim==2) for (size_t i=0; i<Con.Size(); ++i) Con[i]->AddDOF("ux uy",    "fx fy");
    else         for (size_t i=0; i<Con.Size(); ++i) Con[i]->AddDOF("ux uy uz", "fx fy fz");
}

inline void Rod::CalcK (Mat_t & K) const
{
    // rod length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // local K
    double m = E*A/l;
    Mat_t Kl(2,2);
    Kl = m, -m,
        -m,  m;

    // K matrix
    int nrows = 2*NDim; // number of rows in local K matrix
    K.change_dim (nrows,nrows);
    K = trans(T)*Kl*T;
}

inline void Rod::CalcT (Mat_t & T, double & l) const
{
    // coordinates
    double x0 = Con[0]->Vert.C[0];
    double y0 = Con[0]->Vert.C[1];
    double x1 = Con[1]->Vert.C[0];
    double y1 = Con[1]->Vert.C[1];

    if (NDim==2)
    {
        // derived variables
        l = sqrt(pow(x1-x0,2.0)+pow(y1-y0,2.0)); // rod length
        double c = (x1-x0)/l;                    // cosine
        double s = (y1-y0)/l;                    // sine

        // transformation matrix
        T.change_dim (2,4);
        T =   c,   s, 0.0, 0.0,
            0.0, 0.0,   c,   s;
    }
    else
    {
        // derived variables
        double z0 = Con[0]->Vert.C[2];
        double z1 = Con[1]->Vert.C[2];
        l = sqrt(pow(x1-x0,2.0)+pow(y1-y0,2.0)+pow(z1-z0,2.0)); // rod length
        double c = (x1-x0)/l;
        double s = (y1-y0)/l;
        double t = (z1-z0)/l;

        // transformation matrix
        T.change_dim (2,6);
        T =   c,   s,   t, 0.0, 0.0, 0.0,
            0.0, 0.0, 0.0,   c,   s,   t;
    }
}

inline void Rod::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal displacements
    int nrows = Con.Size()*NDim; // number of rows in local K matrix
    Vec_t dUe(nrows);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // K matrix
    Mat_t K;
    CalcK (K);

    // element nodal forces
    Vec_t dFe(nrows);
    dFe = K * dUe;

    // add results to Fint (internal forces)
    if (F_int!=NULL) for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);

    // rod length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // axial strain
    Vec_t dUa(T * dUe); // dUa = T * dUe
    double dea = (dUa(1)-dUa(0))/l;
    static_cast<RodState*>(Sta[0])->ea += dea;
}

inline void Rod::GetState (SDPair & KeysVals, int none) const
{
    double ea = static_cast<RodState*>(Sta[0])->ea;
    double sa = E*ea;
    double fa = A*sa;
    KeysVals.Set("fa sa ea",fa,sa,ea);
}

inline void Rod::Centroid (Vec_t & X) const
{
    X.change_dim (NDim);
    X(0) = (Con[0]->Vert.C[0] + Con[1]->Vert.C[0])/2.0;
    X(1) = (Con[0]->Vert.C[1] + Con[1]->Vert.C[1])/2.0;  if (NDim==3)
    X(2) = (Con[0]->Vert.C[2] + Con[1]->Vert.C[2])/2.0;
}

////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * RodMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new Rod(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int RodRegister()
{
    ElementFactory["Rod"] = RodMaker;
    PROB.Set ("Rod", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __Rod_dummy_int = RodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_ROD
