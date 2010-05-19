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
#include <mechsys/fem/element.h>

namespace FEM
{

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
    void CalcK        (Mat_t & K)                            const; ///< Stiffness matrix
    void CalcT        (Mat_t & T, double & l)                const; ///< Transformation matrix
	void UpdateState  (Vec_t const & dU, Vec_t * F_int=NULL) const;
    void StateKeys    (Array<String> & Keys)                 const; ///< Get state keys
    void StateAtCt    (SDPair & KeysVals)                    const; ///< State at centroid
    void StateAtNodes (Array<SDPair> & Results)              const; ///< State at nodes
    void Centroid     (Vec_t & X)                            const; ///< Centroid of element
    void Draw         (std::ostream & os, double SF)         const;

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
    if (F_int!=NULL)
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
        for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
    }
}

inline void Rod::StateKeys (Array<String> & Keys) const
{
    Keys.Resize (1);
    Keys[0] = "N";
}

inline void Rod::StateAtCt (SDPair & KeysVals) const
{
    // rod length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // displacements in global coordinates
    Vec_t U(2*NDim);
    for (size_t j=0; j<2; ++j)
    {
        U(0+j*NDim) = Con[j]->U[Con[j]->UMap("ux")];
        U(1+j*NDim) = Con[j]->U[Con[j]->UMap("uy")];  if (NDim==3)
        U(2+j*NDim) = Con[j]->U[Con[j]->UMap("uz")];
    }

    // displacements in local coordinates
    Vec_t Ul(T * U);

    // axial force
    double N = E*A*(Ul(1)-Ul(0))/l;
    KeysVals.Set ("N", N);
}

inline void Rod::StateAtNodes (Array<SDPair> & Results) const
{
    SDPair res;
    StateAtCt (res);
    Results.Resize (2);
    Results[0] = res;
    Results[1] = res;
}

inline void Rod::Centroid (Vec_t & X) const
{
    X.change_dim (NDim);
    X(0) = (Con[0]->Vert.C[0] + Con[1]->Vert.C[0])/2.0;
    X(1) = (Con[0]->Vert.C[1] + Con[1]->Vert.C[1])/2.0;  if (NDim==3)
    X(2) = (Con[0]->Vert.C[2] + Con[1]->Vert.C[2])/2.0;
}

inline void Rod::Draw (std::ostream & os, double SF) const
{
    // coordinates
    double x0 = Con[0]->Vert.C[0];
    double y0 = Con[0]->Vert.C[1];
    double x1 = Con[1]->Vert.C[0];
    double y1 = Con[1]->Vert.C[1];

    if (NDim==2)
    {
        // draw shape
        os << "XY = array([["<<x0<<","<<y0<<"],["<<x1<<","<<y1<<"]])\n";
        os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=orange, lw=4))\n";
    }
    else throw new Fatal("Rod::Draw: Method not available for 3D yet");
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
