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
         Model        const * XMdl,   ///< Extra Model
         SDPair       const & Prp,    ///< Properties
         SDPair       const & Ini,    ///< Initial values
         Array<Node*> const & Nodes); ///< Connectivity

    // Methods
    void GetLoc       (Array<size_t> & Loc)                     const; ///< Get location vector for mounting K/M matrices
    void CalcK        (Mat_t & K)                               const; ///< Stiffness matrix
    void CalcM        (Mat_t & M)                               const; ///< Mass matrix
    void CalcT        (Mat_t & T, double & l)                   const; ///< Transformation matrix
    void UpdateState  (Vec_t const & dU, Vec_t * F_int=NULL)    const;
    void StateKeys    (Array<String> & Keys)                    const; ///< Get state keys
    void StateAtNodes (Array<SDPair> & Results)                 const; ///< State at nodes
    void Draw         (std::ostream & os, MPyPrms const & Prms) const;

    // Constants
    double E; ///< Young modulus
    double A; ///< Cross-sectional area
    double rho; ///< density
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Rod::Rod (int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes)
{
    // check GTy
    if (GTy!=fra_t) throw new Fatal("Rod::Rod: Geometry type (GTy) must be equal to 'fra' (Frame). GTy=%s is invalid",GTypeToStr(GTy).CStr());

    // parameters/properties
    E = Prp("E");
    A = Prp("A");
    rho = (Prp.HasKey("rho") ? Prp("rho") : 1.0);
}

inline void Rod::GetLoc (Array<size_t> & Loc) const
{
    Loc.Resize (2*NDim);
    for (size_t i=0; i<2; ++i)
    {
        Loc[i*NDim+0] = Con[i]->Eq("ux");
        Loc[i*NDim+1] = Con[i]->Eq("uy");  if (NDim==3)
        Loc[i*NDim+2] = Con[i]->Eq("uz");
    }
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

inline void Rod::CalcM (Mat_t & M) const
{
    // length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // local M
    Mat_t Ml(2,2);
    Ml = 2., 1.,
         1., 2.;

    // M matrix
    double m = rho*A*l/6.;
    M.change_dim (4,4);
    M = m*trans(T)*Ml*T;
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

inline void Rod::StateAtNodes (Array<SDPair> & Results) const
{
    // rod length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // displacements in global coordinates
    Vec_t U(2*NDim);
    for (size_t j=0; j<2; ++j)
    {
        U(0+j*NDim) = Con[j]->U("ux");
        U(1+j*NDim) = Con[j]->U("uy");  if (NDim==3)
        U(2+j*NDim) = Con[j]->U("uz");
    }

    // displacements in local coordinates
    Vec_t Ul(T * U);

    // axial force
    double N = E*A*(Ul(1)-Ul(0))/l;
    Results.Resize (2);
    Results[0].Set ("N", N);
    Results[1].Set ("N", N);
}

inline void Rod::Draw (std::ostream & os, MPyPrms const & Prms) const
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
Element * RodMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new Rod(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes); }

// Register element
int RodRegister()
{
    ElementFactory  ["Rod"]   = RodMaker;
    ElementVarKeys  ["Rod2D"] = std::make_pair ("ux uy",    "fx fy");
    ElementVarKeys  ["Rod3D"] = std::make_pair ("ux uy uz", "fx fy fz");
    ElementExtraKeys["Rod2D"] = Array<String>  ();
    ElementExtraKeys["Rod3D"] = Array<String>  ();
    PROB.Set ("Rod", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __Rod_dummy_int = RodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_ROD
