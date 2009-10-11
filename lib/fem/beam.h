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

#ifndef MECHSYS_FEM_BEAM_H
#define MECHSYS_FEM_BEAM_H

// Std Lib
#include <ostream>

// MechSys
#include "fem/element.h"
#include "models/model.h"
#include "util/util.h"

namespace FEM
{

class Beam : public Element
{
public:
    // Constructor
    Beam (int                  NDim,   ///< Space dimension
          Mesh::Cell   const & Cell,   ///< Geometric information: ID, Tag, connectivity
          Model        const * Mdl,    ///< Model
          SDPair       const & Prp,    ///< Properties
          SDPair       const & Ini,    ///< Initial values
          Array<Node*> const & Nodes); ///< Array with all nodes (used to set the connectivity)

    // Methods
    void SetBCs      (size_t IdxEdgeOrFace, SDPair const & BCs);           ///< IdxEdgeOrFace is ignored
    void CalcK       (Mat_t & K)                                    const; ///< Stiffness matrix
    void CalcM       (Mat_t & M)                                    const; ///< Mass matrix
    void CalcT       (Mat_t & T, double & l)                        const; ///< Transformation matrix
    void UpdateState (Vec_t const & dU, Vec_t * F_int=NULL)         const;
    void GetState    (SDPair & KeysVals, int none=-1)               const;
    void Centroid    (Vec_t & X)                                    const; ///< Centroid of element
    void CalcRes     (double r, double & P, double & V, double & M) const; ///< Resultants: Axial force P, Shear force V, Bending moment M
    void Draw        (std::ostream & os, double SF)                 const;

    // Constants
    double E;     ///< Young modulus
    double A;     ///< Cross-sectional area
    double Izz;   ///< Inertia
    double rho;   ///< Density
    double qnl;   ///< Normal load (left)
    double qnr;   ///< Normal load (right)
    bool   HasQn; ///< Has normal load ?
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Beam::Beam (int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes)
    : Element(NDim,Cell,Mdl,Prp,Ini,Nodes), qnl(0.0), qnr(0.0), HasQn(false)
{
    // check GTy
    if (GTy!=fra_t) throw new Fatal("Beam::Beam: Geometry type (GTy) must be equal to 'fra' (Frame). GTy=%s is invalid",GTypeToStr(GTy).CStr());

    // parameters/properties
    E   = Prp("E");
    A   = Prp("A");
    Izz = Prp("Izz");
    rho = (Prp.HasKey("rho") ? Prp("rho") : 1.0);

    // set UKeys in parent element
    if (NDim==2) { UKeys.Resize(3);  UKeys = "ux", "uy", "wz"; }
    else throw new Fatal("Beam::Beam: 3D Beam is not available yet");

    // initialize DOFs
    if (NDim==2) for (size_t i=0; i<Con.Size(); ++i) Con[i]->AddDOF("ux uy wz", "fx fy mz");
    else {}
}

inline void Beam::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs)
{
    // length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    if (BCs.HasKey("qn") || BCs.HasKey("qnl") || BCs.HasKey("qnr"))
    {
        if (NDim==3) throw new Fatal("Beam::SetBCs: Method not available for 3D yet (with 'qn' or 'qnl' or 'qnr')");

        // loading
        qnl = (BCs.HasKey("qn") ? BCs("qn") : BCs("qnl")); // qn at left
        qnr = (BCs.HasKey("qn") ? BCs("qn") : BCs("qnr")); // qn at right
        HasQn = true;

        // is node 0 leftmost ?
        bool n0_is_left = true;
        if (fabs(Con[1]->Vert.C[0]-Con[0]->Vert.C[0])<1.0e-7) { // vertical segment
             if (Con[1]->Vert.C[1]<Con[0]->Vert.C[1]) n0_is_left = false; }
        else if (Con[1]->Vert.C[0]<Con[0]->Vert.C[0]) n0_is_left = false;
        if (!n0_is_left) { Util::Swap(qnl,qnr);  qnl *= -1.;  qnr *= -1.; }

        // local and global forces
        Vec_t Fe(6);
        Fe = 0.0, l*(7.0*qnl+3.0*qnr)/20.0,  l*l*(3.0*qnl+2.0*qnr)/60.0,
             0.0, l*(3.0*qnl+7.0*qnr)/20.0, -l*l*(2.0*qnl+3.0*qnr)/60.0;
        Vec_t F(trans(T)*Fe);

        // add to nodes
        for (size_t j=0; j<2; ++j)
        {
            Con[j]->DF[Con[j]->FMap("fx")] += F(0+j*3);
            Con[j]->DF[Con[j]->FMap("fy")] += F(1+j*3);
            Con[j]->DF[Con[j]->FMap("mz")] += F(2+j*3);
        }
    }
    else
    {
        std::ostringstream oss;
        oss << BCs;
        throw new Fatal("Beam::SetBCs: This method does not work yet with BCs=%s",oss.str().c_str());
    }
}

inline void Beam::CalcK (Mat_t & K) const
{
    // length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // local K
    double ll = l*l;
    double m  = E*A/l;
    double n  = E*Izz/(ll*l);
    Mat_t Kl(6,6);
    Kl =  m,        0.,       0.,   -m,        0.,       0.,
         0.,   12.  *n,  6.*l *n,   0.,  -12.  *n,  6.*l *n,
         0.,    6.*l*n,  4.*ll*n,   0.,   -6.*l*n,  2.*ll*n,
         -m,        0.,       0.,    m,        0.,       0.,
         0.,  -12.  *n, -6.*l *n,   0.,   12.  *n, -6.*l *n,
         0.,    6.*l*n,  2.*ll*n,   0.,   -6.*l*n,  4.*ll*n;

    // K matrix
    K.change_dim (6,6);
    K = trans(T)*Kl*T;
}

inline void Beam::CalcM (Mat_t & M) const
{
    // length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    // local M
    double ll = l*l;
    double m  = rho*A*l/420.;
    Mat_t Ml(6,6);
    Ml = 140.*m ,    0.     ,   0.      ,   70.*m ,    0.     ,    0.      ,
           0.   ,  156.*m   ,  22.*l*m  ,    0.   ,   54.*m   ,  -13.*l*m  ,
           0.   ,   22.*l*m ,   4.*ll*m ,    0.   ,   13.*l*m ,   -3.*ll*m ,
          70.*m ,    0.     ,   0.      ,  140.*m ,    0.     ,    0.      ,
           0.   ,   54.*m   ,  13.*l*m  ,    0.   ,  156.*m   ,  -22.*l*m  ,
           0.   ,  -13.*l*m ,  -3.*ll*m ,    0.   ,  -22.*l*m ,    4.*ll*m ;

    // M matrix
    M.change_dim (6,6);
    M = trans(T)*Ml*T;
}

inline void Beam::CalcT (Mat_t & T, double & l) const
{
    // coordinates
    double x0 = Con[0]->Vert.C[0];
    double y0 = Con[0]->Vert.C[1];
    double x1 = Con[1]->Vert.C[0];
    double y1 = Con[1]->Vert.C[1];

    if (NDim==2)
    {
        // derived variables
        l = sqrt(pow(x1-x0,2.0)+pow(y1-y0,2.0)); // length
        double c = (x1-x0)/l;                    // cosine
        double s = (y1-y0)/l;                    // sine

        // transformation matrix
        T.change_dim (6,6);
        T =   c,  s,  0.,  0.,  0.,  0.,
             -s,  c,  0.,  0.,  0.,  0.,
             0., 0.,  1.,  0.,  0.,  0.,
             0., 0.,  0.,   c,   s,  0.,
             0., 0.,  0.,  -s,   c,  0.,
             0., 0.,  0.,  0.,  0.,  1.;
    }
    else throw new Fatal("Beam::CalcT: 3D Beam is not available yet");
}

inline void Beam::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    // get location array
    Array<size_t> loc;
    GetLoc (loc);

    // element nodal displacements
    Vec_t dUe(6);
    for (size_t i=0; i<loc.Size(); ++i) dUe(i) = dU(loc[i]);

    // K matrix
    Mat_t K;
    CalcK (K);

    // element nodal forces
    Vec_t dFe(6);
    dFe = K * dUe;

    // add results to Fint (internal forces)
    if (F_int!=NULL) for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
}

inline void Beam::GetState (SDPair & KeysVals, int none) const
{
    double P, V, M;
    CalcRes (/*rct*/0.5, P, V, M);
    KeysVals.Set ("P V M",P,V,M);
}

inline void Beam::Centroid (Vec_t & X) const
{
    X.change_dim (NDim);
    X(0) = (Con[0]->Vert.C[0] + Con[1]->Vert.C[0])/2.0;
    X(1) = (Con[0]->Vert.C[1] + Con[1]->Vert.C[1])/2.0;  if (NDim==3)
    X(2) = (Con[0]->Vert.C[2] + Con[1]->Vert.C[2])/2.0;
}

inline void Beam::CalcRes (double r, double & P, double & V, double & M) const
{
    // rod length and T matrix
    double l;
    Mat_t  T;
    CalcT (T, l);

    if (NDim==2)
    {
        // displacements in global coordinates
        Vec_t U(6);
        for (size_t j=0; j<2; ++j)
        {
            U(0+j*3) = Con[j]->U[Con[j]->UMap("ux")];
            U(1+j*3) = Con[j]->U[Con[j]->UMap("uy")];
            U(2+j*3) = Con[j]->U[Con[j]->UMap("wz")];
        }

        // displacements in local coordinates
        Vec_t Ul(T * U);

        // local (natural) coordinate
        double s   = r*l;
        double ll  = l*l;
        double lll = ll*l;

        // axial force
        P = E*A*(Ul(3)-Ul(0))/l;

        // shear force
        V = E*Izz*((12.*Ul(1))/lll + (6.*Ul(2))/ll - (12.*Ul(4))/lll + (6.*Ul(5))/ll);

        // bending moment
        M = E*Izz*(Ul(1)*((12.*s)/lll-6./ll) + Ul(2)*((6.*s)/ll-4./l) + Ul(4)*(6./ll-(12.*s)/lll) + Ul(5)*((6.*s)/ll-2./l));

        if (HasQn)
        {
            double ss  = s*s;
            double sss = ss*s;
            V += -(3.*qnr*ll +7.*qnl*ll-20.*qnl*s*l -10.*qnr*ss  +10.*qnl*ss)/(20.*l);
            M +=  (2.*qnr*lll+3.*qnl*lll-9.*qnr*s*ll-21.*qnl*s*ll+30.*qnl*ss*l+10.*qnr*sss-10.*qnl*sss)/(60.*l);
        }
    }
}

inline void Beam::Draw (std::ostream & os, double SF) const
{
    // constants
    size_t ndiv = 10;

    // coordinates
    double x0 = Con[0]->Vert.C[0];
    double y0 = Con[0]->Vert.C[1];
    double x1 = Con[1]->Vert.C[0];
    double y1 = Con[1]->Vert.C[1];

    if (NDim==2)
    {
        // draw shape
        os << "XY = array([["<<x0<<","<<y0<<"],["<<x1<<","<<y1<<"]])\n";
        os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=dred, lw=8))\n";

        // derived variables
        double l = sqrt(pow(x1-x0,2.0)+pow(y1-y0,2.0)); // length
        double c = (x1-x0)/l;                           // cosine
        double s = (y1-y0)/l;                           // sine

        // normal
        double xn = -s;
        double yn =  c;

        // results
        double P, V, M, Mmax;
        double r, x, y, rMmax;
        double sf, xf, yf;
        rMmax = 0.0;
        CalcRes (rMmax, P, V, Mmax);
        os << "dat_beam = []\n";
        for (size_t i=0; i<ndiv+1; ++i)
        {
            r = static_cast<double>(i)/static_cast<double>(ndiv);
            CalcRes (r, P, V, M);
            x  = x0 + r*(x1-x0);
            y  = y0 + r*(y1-y0);
            sf = SF*M;
            xf = x - sf*xn;
            yf = y - sf*yn;
            os << "XY = array([["<<x<<","<<y<<"],["<<xf<<","<<yf<<"]])\n";
            os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor='blue', lw=1))\n";
            if (i>0) os << "dat_beam.append((PH.LINETO, (" << xf << "," << yf << ")))\n";
            else     os << "dat_beam.append((PH.MOVETO, (" << xf << "," << yf << ")))\n";
            if (fabs(M)>fabs(Mmax)) { rMmax = r;  Mmax = M; }
        }
        os << "cmd_beam,vert_beam = zip(*dat_beam)\n";
        os << "ph_beam       = PH (vert_beam, cmd_beam)\n";
        os << "pc_beam       = PC (ph_beam, facecolor='none', edgecolor=dblue, linewidth=1)\n";
        os << "ax.add_patch  (pc_beam)\n\n";

        // max M
        x  = x0 + rMmax*(x1-x0);
        y  = y0 + rMmax*(y1-y0);
        sf = SF*Mmax;
        xf = x - sf*xn;
        yf = y - sf*yn;
        String buf;
        buf.Printf ("%g",Mmax);
        os << "XY = array([["<<x<<","<<y<<"],["<<xf<<","<<yf<<"]])\n";
        os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=dblue, lw=4))\n";
        os << "ax.text ("<<(x+xf)/2.<<","<<(y+yf)/2.<<", " << buf << ", backgroundcolor=pink, va='top', ha='center', fontsize=12)\n";
        //std::cout << "r(MaxM) = " << rMmax << "    MaxM = " << Mmax << std::endl;
    }
    else throw new Fatal("Beam::GetState: 3D Beam is not available yet");
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * BeamMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new Beam(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int BeamRegister()
{
    ElementFactory["Beam"] = BeamMaker;
    PROB.Set ("Beam", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __Beam_dummy_int = BeamRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_BEAM