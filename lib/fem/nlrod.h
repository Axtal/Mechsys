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

#ifndef MECHSYS_FEM_NLROD_H
#define MECHSYS_FEM_NLROD_H

// MechSys
#include <mechsys/fem/element.h>
#include <mechsys/models/model.h>

namespace FEM
{

class NLRodState : public State
{
public:
    NLRodState (int NDim) : State(NDim) {}
    void   Init    (SDPair const & Ini, size_t NIvs=0) { sa = (Ini.HasKey("sa") ? Ini("sa") : 0.0);  ea = (Ini.HasKey("ea") ? Ini("ea") : 0.0); }
    void   Backup  () { sa_bkp = sa;  ea_bkp = ea; }
    void   Restore () { sa = sa_bkp;  ea = ea_bkp; }
    size_t PckSize ()                        const { return 2; }
    void   Pack    (Array<double>       & V) const { V.Resize(2);  V = sa, ea; }
    void   Unpack  (Array<double> const & V)       { sa=V[0]; ea=V[1]; }
    double sa, sa_bkp; ///< Axial stress
    double ea, ea_bkp; ///< Axial strain
};

class NLRod : public Element
{
public:
    // Constructor
    NLRod (int                  NDim,   ///< Space dimension
           Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
           Model        const * Mdl,    ///< Model
           Model        const * XMdl,   ///< Extra Model
           SDPair       const & Prp,    ///< Properties
           SDPair       const & Ini,    ///< Initial values
           Array<Node*> const & Nodes); ///< Connectivity

    // Methods
    void GetLoc      (Array<size_t> & Loc)                  const; ///< Get location vector for mounting K/M matrices
    void CalcK       (Mat_t & K)                            const; ///< Stiffness matrix
    void CalcT       (Mat_t & T, double & l)                const; ///< Transformation matrix
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL) const;
    void GetState    (SDPair & KeysVals, int none=-1)       const;

    // Constants
    double E0;  ///< Initial Young modulus
    double alp; ///< Nonlinear parameter
    double A;   ///< Cross-sectional area
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline NLRod::NLRod (int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes)
{
    // check GTy
    if (GTy!=fra_t) throw new Fatal("NLRod::NLRod: Geometry type (GTy) must be equal to 'fra' (Frame). GTy=%s is invalid",GTypeToStr(GTy).CStr());

    // parameters/properties
    E0  = Prp("E0");
    alp = Prp("alp");
    A   = Prp("A");

    // allocate and initialize state (constant along element)
    Sta.Push (new NLRodState(NDim));
    Sta[0]->Init (Ini);
}

inline void NLRod::GetLoc (Array<size_t> & Loc) const
{
    Loc.Resize (2*NDim);
    for (size_t i=0; i<2; ++i)
    {
        Loc[i*NDim+0] = Con[i]->Eq("ux");
        Loc[i*NDim+1] = Con[i]->Eq("uy");  if (NDim==3)
        Loc[i*NDim+2] = Con[i]->Eq("uz");
    }
}

inline void NLRod::CalcK (Mat_t & K) const
{
    // rod length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // local K
    double dgde = (1.0-2.0*alp*static_cast<NLRodState*>(Sta[0])->ea);
    double m    = dgde*E0*A/l;
    Mat_t Kl(2,2);
    Kl = m, -m,
        -m,  m;

    // K matrix
    int nrows = 2*NDim; // number of rows in local K matrix
    K.change_dim (nrows,nrows);
    K = trans(T)*Kl*T;
}

inline void NLRod::CalcT (Mat_t & T, double & l) const
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

inline void NLRod::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal displacements
    int nrows = Con.Size()*NDim; // number of rows in local K matrix
    Vec_t dUe(nrows);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // rod length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // axial strain increment
    Vec_t dUa(T * dUe); // dUa = T * dUe
    double dea = (dUa(1)-dUa(0))/l;

    // update strain and stress
    NLRodState * sta = static_cast<NLRodState*>(Sta[0]);
    double dsa = sta->sa;
    sta->ea += dea;
    sta->sa  = E0*(sta->ea - alp*pow(sta->ea,2.0));
    dsa      = sta->sa - dsa;

    // element nodal forces
    Vec_t dFl(2); // local
    dFl = -A*(dsa), A*(dsa);
    Vec_t dFe(nrows);
    dFe = trans(T)*dFl;

    // add results to Fint (internal forces)
    if (F_int!=NULL) for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
}

inline void NLRod::GetState (SDPair & KeysVals, int none) const
{
    double ea = static_cast<NLRodState*>(Sta[0])->ea;
    double sa = static_cast<NLRodState*>(Sta[0])->sa;
    double fa = A*sa;
    KeysVals.Set("fa sa ea",fa,sa,ea);
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * NLRodMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, Model const * XMdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new NLRod(NDim,Cell,Mdl,XMdl,Prp,Ini,Nodes); }

// Register element
int NLRodRegister()
{
    ElementFactory["NLRod"]   = NLRodMaker;
    ElementVarKeys["NLRod2D"] = std::make_pair ("ux uy",    "fx fy");
    ElementVarKeys["NLRod3D"] = std::make_pair ("ux uy uz", "fx fy fz");
    PROB.Set ("NLRod", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __NLRod_dummy_int = NLRodRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_NLROD
