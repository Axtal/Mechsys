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
#include <mechsys/fem/element.h>
#include <mechsys/util/util.h>

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
          Array<Node*> const & Nodes); ///< Connectivity

    // Methods
    void SetBCs       (size_t IdxEdgeOrFace, SDPair const & BCs, PtBCMult MFun); ///< IdxEdgeOrFace is ignored
    void ClrBCs       ();                                                        ///< Clear boundary conditions
    void GetLoc       (Array<size_t> & Loc)                               const; ///< Get location vector for mounting K/M matrices
    void CalcK        (Mat_t & K)                                         const; ///< Stiffness matrix
    void CalcM        (Mat_t & M)                                         const; ///< Mass matrix
    void CalcT        (Mat_t & T, double & l, Vec3_t * normal=NULL)       const; ///< Transformation matrix
    void UpdateState  (Vec_t const & dU, Vec_t * F_int=NULL)              const;
    void StateKeys    (Array<String> & Keys)                              const; ///< Get state keys
    void StateAtNodes (Array<SDPair> & Results)                           const; ///< State at nodes
    void CalcRes      (double r, double & N, double & V, double & M)      const; ///< Resultants: Axial force N, Shear force V, Bending moment M
    void Draw         (std::ostream & os, MPyPrms const & Prms)           const; ///< Draw beam and diagram of M, V, or N

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
}

inline void Beam::SetBCs (size_t IdxEdgeOrFace, SDPair const & BCs, PtBCMult MFun)
{
    if (NDim==3) throw new Fatal("Beam::SetBCs: Method not available for 3D yet");

    // length and T matrix
    double l;
    Mat_t  T;
    Vec3_t normal;
    CalcT (T, l, &normal);

    // boundary conditions
    bool has_qx  = BCs.HasKey("qx");  // x component of distributed loading
    bool has_qy  = BCs.HasKey("qy");  // y component of distributed loading
    bool has_qn  = BCs.HasKey("qn");  // normal distributed loading
    bool has_qnl = BCs.HasKey("qnl"); // normal distributed loading (left)
    bool has_qnr = BCs.HasKey("qnr"); // normal distributed loading (right)
    bool has_qt  = BCs.HasKey("qt");  // tangential distributed loading (2D only)

    // loads
    HasQn = false;
    qnl   = 0.0;
    qnr   = 0.0;
    double qt = 0.0;
    if (BCs.HasKey("gravity"))
    {
        double g   = BCs("gravity");
        double gam = rho*g;
        double qy  = -gam*A;
        double qn  =  normal(1)*qy;
        qt    = -normal(0)*qy;
        qnl   = qn;
        qnr   = qn;
        HasQn = true;
    }
    else if (has_qx || has_qy)
    {
        throw new Fatal("Beam::SetBCs: Setting up of boundary conditions with 'qx' and 'qy' is not available");
        // TODO: this need check (seems OK, though)
        /*
        double qx = (has_qx ? BCs("qx") : 0.0);
        double qy = (has_qy ? BCs("qy") : 0.0);
        double qn = normal(0)*qx + normal(1)*qy;
        qt    = normal(1)*qx - normal(0)*qy;
        qnl   = qn;
        qnr   = qn;
        HasQn = true;
        */
    }
    else if (has_qn || has_qnl || has_qnr || has_qt)
    {
        if (has_qnl || has_qnr)
        {
            qnl = (has_qnl ? BCs("qnl") : 0.0);
            qnr = (has_qnr ? BCs("qnr") : 0.0);
            if (qnl*qnr<0.0) throw new Fatal("Beam::SetBCs: qnl=%g and qnr=%g must have the same sign",qnl,qnr);
        }
        else
        {
            double qn = (has_qn ? BCs("qn") : 0.0);
            qnl = qn;
            qnr = qn;
        }
        qt    = (has_qt ? BCs("qt") : 0.0);
        HasQn = true;

        // is node 0 leftmost ?
        bool n0_is_left = true;
        if (fabs(Con[1]->Vert.C[0]-Con[0]->Vert.C[0])<1.0e-7) { // vertical segment
             if (Con[1]->Vert.C[1]<Con[0]->Vert.C[1]) n0_is_left = false; }
        else if (Con[1]->Vert.C[0]<Con[0]->Vert.C[0]) n0_is_left = false;
        if (!n0_is_left)
        {
            Util::Swap(qnl,qnr);
            qnl *= -1.;
            qnr *= -1.; 
            qt  *= -1.;
        }
    }

    // set equivalent forces at nodes
    if (HasQn)
    {
        // local and global forces
        Vec_t Fe(6);
        Fe = qt*l/2.0, l*(7.0*qnl+3.0*qnr)/20.0,  l*l*(3.0*qnl+2.0*qnr)/60.0,
             qt*l/2.0, l*(3.0*qnl+7.0*qnr)/20.0, -l*l*(2.0*qnl+3.0*qnr)/60.0;
        Vec_t F(trans(T)*Fe);

        // add to nodes
        for (size_t j=0; j<2; ++j)
        {
            Con[j]->AddToPF("fx", F(0+j*3), MFun);
            Con[j]->AddToPF("fy", F(1+j*3), MFun);
            Con[j]->AddToPF("mz", F(2+j*3), MFun);
        }
    }
}

inline void Beam::ClrBCs ()
{
    qnl   = 0.0;
    qnr   = 0.0;
    HasQn = false;
}

inline void Beam::GetLoc (Array<size_t> & Loc) const
{
    Loc.Resize (6);
    for (size_t i=0; i<2; ++i)
    {
        Loc[i*3+0] = Con[i]->Eq("ux");
        Loc[i*3+1] = Con[i]->Eq("uy");
        Loc[i*3+2] = Con[i]->Eq("wz");
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

inline void Beam::CalcT (Mat_t & T, double & l, Vec3_t * normal) const
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

        // normal
        if (normal!=NULL) (*normal) = -s, c, 0.0;
    }
    else throw new Fatal("Beam::CalcT: 3D Beam is not available yet");
}

inline void Beam::UpdateState (Vec_t const & dU, Vec_t * F_int) const
{
    if (F_int!=NULL)
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

        /*
        // length and T matrix
        double l;
        Mat_t  T;
        CalcT (T, l);

        // increment of displacements in local coordinates
        Vec_t dUl(T * dUe);

        // axial and shear forces
        double ll  = l*l;
        double lll = ll*l;
        double dN  = E*A*(dUl(3)-dUl(0))/l;
        double dV  = E*Izz*((12.*dUl(1))/lll + (6.*dUl(2))/ll - (12.*dUl(4))/lll + (6.*dUl(5))/ll);
        double dM0 = E*Izz*(dUl(1)*(           -6./ll) + dUl(2)*(         -4./l) + dUl(4)*(6./ll            ) + dUl(5)*(         -2./l));
        double dM1 = E*Izz*(dUl(1)*((12.*l)/lll-6./ll) + dUl(2)*((6.*l)/ll-4./l) + dUl(4)*(6./ll-(12.*l)/lll) + dUl(5)*((6.*l)/ll-2./l));

        // element nodal forces in global coordinates
        Vec_t dFl(6);
        dFl = -dN, dV, -dM0,  dN, -dV, dM1;
        Vec_t dFe(trans(T) * dFl);
        */

        // add results to Fint (internal forces)
        for (size_t i=0; i<loc.Size(); ++i) (*F_int)(loc[i]) += dFe(i);
    }
}

inline void Beam::StateKeys (Array<String> & Keys) const
{
    Keys.Resize (3);
    Keys = "N", "V", "M";
}

inline void Beam::StateAtNodes (Array<SDPair> & Results) const
{
    Results.Resize (2);
    double N, V, M;
    CalcRes (0.0, N, V, M);
    Results[0].Set ("N V M",N,V,M);
    CalcRes (1.0, N, V, M);
    Results[1].Set ("N V M",N,V,M);
}

inline void Beam::CalcRes (double r, double & N, double & V, double & M) const
{
    // rod length and T matrix
    double l;
    Mat_t  T;
    Vec3_t normal;
    CalcT (T, l, &normal);

    if (NDim==2)
    {
        // displacements in global coordinates
        Vec_t U(6);
        for (size_t j=0; j<2; ++j)
        {
            U(0+j*3) = Con[j]->U("ux");
            U(1+j*3) = Con[j]->U("uy");
            U(2+j*3) = Con[j]->U("wz");
        }

        // displacements in local coordinates
        Vec_t Ul(T * U);

        // local (natural) coordinate
        double s   = r*l;
        double ll  = l*l;
        double lll = ll*l;

        // axial force
        N = E*A*(Ul(3)-Ul(0))/l;

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

        // correct the sign of M
        if (qnl>0.0) M = -M;
    }
}

inline void Beam::Draw (std::ostream & os, MPyPrms const & Prms) const
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
        os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor=dred, lw=8))\n";

        // derived variables
        double l = sqrt(pow(x1-x0,2.0)+pow(y1-y0,2.0)); // length
        double c = (x1-x0)/l;                           // cosine
        double s = (y1-y0)/l;                           // sine

        // normal
        double xn = -s;
        double yn =  c;

        // results
        double N, V, M,  Mmax,  Mmin,  Vmax,  Vmin;
        double r, x, y, rMmax, rMmin, rVmax, rVmin;
        double sf, xf, yf, val;
        rMmax = 0.0;
        rMmin = 0.0;
        rVmax = 0.0;
        rVmin = 0.0;
        CalcRes (rMmax, N, Vmax, Mmax);
        CalcRes (rMmin, N, Vmin, Mmin);
        os << "dat_beam = []\n";
        if (Prms.DrawN) os << "dat_beam_ax = []\n";
        for (size_t i=0; i<Prms.NDiv+1; ++i)
        {
            r = static_cast<double>(i)/static_cast<double>(Prms.NDiv);
            CalcRes (r, N, V, M);
            if (Prms.DrawN)
            {
                val = N;
                sf  = fabs(Prms.SF*N)/2.0;
                x   = x0 + r*(x1-x0);
                y   = y0 + r*(y1-y0);
                xf  = x - sf*xn;
                yf  = y - sf*yn;
                x   = x + sf*xn;
                y   = y + sf*yn;
            }
            else
            {
                if (Prms.DrawV)
                {
                    val = V;
                    sf  = Prms.SF*V;
                    if (V>Vmax) { rVmax = r;  Vmax = V; }
                    if (V<Vmin) { rVmin = r;  Vmin = V; }
                }
                else
                {
                    val = M;
                    if (qnl>0.0) sf = -Prms.SF*M;
                    else         sf =  Prms.SF*M;
                    if (M>Mmax) { rMmax = r;  Mmax = M; }
                    if (M<Mmin) { rMmin = r;  Mmin = M; }
                }
                x  = x0 + r*(x1-x0);
                y  = y0 + r*(y1-y0);
                xf = x - sf*xn;
                yf = y - sf*yn;
            }
            os << "XY = array([["<<x<<","<<y<<"],["<<xf<<","<<yf<<"]])\n";
            os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor="<<(val<0.0?"dpink":"dblue")<<", lw=1))\n";
            if (i>0) os << "dat_beam.append((PH.LINETO, (" << xf << "," << yf << ")))\n";
            else     os << "dat_beam.append((PH.MOVETO, (" << xf << "," << yf << ")))\n";
            if (Prms.DrawN)
            {
                if (i>0) os << "dat_beam_ax.append((PH.LINETO, (" << x << "," << y << ")))\n";
                else     os << "dat_beam_ax.append((PH.MOVETO, (" << x << "," << y << ")))\n";
            }
        }
        os << "cmd_beam,vert_beam = zip(*dat_beam)\n";
        os << "ph_beam       = PH (vert_beam, cmd_beam)\n";
        os << "pc_beam       = PC (ph_beam, facecolor='none', edgecolor=dblue, linewidth=1)\n";
        os << "ax.add_patch  (pc_beam)\n\n";
        if (Prms.DrawN)
        {
            os << "cmd_beam_ax,vert_beam_ax = zip(*dat_beam_ax)\n";
            os << "ph_beam_ax       = PH (vert_beam_ax, cmd_beam_ax)\n";
            os << "pc_beam_ax       = PC (ph_beam_ax, facecolor='none', edgecolor=dblue, linewidth=1)\n";
            os << "ax.add_patch    (pc_beam_ax)\n\n";
        }

        // Text
        if (Prms.WithTxt)
        {
            if (Prms.DrawN)
            {
                // max or min N
                bool is_max_or_min = (this==Prms.EleNmax || this==Prms.EleNmin);
                bool skip = (Prms.OnlyTxtLim && !is_max_or_min);
                CalcRes (0.5, N, V, M);
                if (fabs(N)>1.0e-13 && !skip)
                {
                    sf = fabs(Prms.SF*N)/2.0;
                    x  = x0 + 0.5*(x1-x0);
                    y  = y0 + 0.5*(y1-y0);
                    String buf;
                    buf.Printf ("%g",N);
                    os << "ax.text ("<<x<<","<<y<<", " << buf << ", backgroundcolor=pink, va='top', ha='center', fontsize="<<Prms.TxtSz<<")\n";
                }
            }
            else if (Prms.DrawV)
            {
                // max V
                bool skip = (Prms.OnlyTxtLim && Prms.EleVmax!=this);
                if (fabs(Vmax)>1.0e-13 && !skip)
                {
                    x  = x0 + rVmax*(x1-x0);
                    y  = y0 + rVmax*(y1-y0);
                    sf = Prms.SF*Vmax;
                    xf = x - sf*xn;
                    yf = y - sf*yn;
                    String buf;
                    buf.Printf ("%g",Vmax);
                    os << "XY = array([["<<x<<","<<y<<"],["<<xf<<","<<yf<<"]])\n";
                    os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor="<<(Vmax<0.0?"dpink":"dblue")<<", lw=4))\n";
                    os << "ax.text ("<<(x+xf)/2.<<","<<(y+yf)/2.<<", " << buf << ", backgroundcolor=pink, va='top', ha='center', fontsize="<<Prms.TxtSz<<")\n";
                }

                // min V
                skip = (Prms.OnlyTxtLim && Prms.EleVmin!=this);
                if (fabs(Vmin)>1.0e-13 && !skip)
                {
                    x  = x0 + rVmin*(x1-x0);
                    y  = y0 + rVmin*(y1-y0);
                    sf = Prms.SF*Vmin;
                    xf = x - sf*xn;
                    yf = y - sf*yn;
                    String buf;
                    buf.Printf ("%g",Vmin);
                    os << "XY = array([["<<x<<","<<y<<"],["<<xf<<","<<yf<<"]])\n";
                    os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor="<<(Vmin<0.0?"dpink":"dblue")<<", lw=4))\n";
                    os << "ax.text ("<<(x+xf)/2.<<","<<(y+yf)/2.<<", " << buf << ", backgroundcolor=pink, va='top', ha='center', fontsize="<<Prms.TxtSz<<")\n";
                }
            }
            else // DrawM
            {
                // max M
                bool skip = (Prms.OnlyTxtLim && Prms.EleMmax!=this);
                if (fabs(Mmax)>1.0e-13 && !skip)
                {
                    x  = x0 + rMmax*(x1-x0);
                    y  = y0 + rMmax*(y1-y0);
                    if (qnl>0.0) sf = -Prms.SF*Mmax;
                    else         sf =  Prms.SF*Mmax;
                    xf = x - sf*xn;
                    yf = y - sf*yn;
                    String buf;
                    buf.Printf ("%g",Mmax);
                    os << "XY = array([["<<x<<","<<y<<"],["<<xf<<","<<yf<<"]])\n";
                    os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor="<<(Mmax<0.0?"dpink":"dblue")<<", lw=4))\n";
                    os << "ax.text ("<<(x+xf)/2.<<","<<(y+yf)/2.<<", " << buf << ", backgroundcolor=pink, va='top', ha='center', fontsize="<<Prms.TxtSz<<")\n";
                }
                
                // min M
                skip = (Prms.OnlyTxtLim && Prms.EleMmin!=this);
                if (fabs(Mmin)>1.0e-13 && !skip)
                {
                    x  = x0 + rMmin*(x1-x0);
                    y  = y0 + rMmin*(y1-y0);
                    if (qnl>0.0) sf = -Prms.SF*Mmin;
                    else         sf =  Prms.SF*Mmin;
                    xf = x - sf*xn;
                    yf = y - sf*yn;
                    String buf;
                    buf.Printf ("%g",Mmin);
                    os << "XY = array([["<<x<<","<<y<<"],["<<xf<<","<<yf<<"]])\n";
                    os << "ax.add_patch (MPL.patches.Polygon(XY, closed=False, edgecolor="<<(Mmin<0.0?"dpink":"dblue")<<", lw=4))\n";
                    os << "ax.text ("<<(x+xf)/2.<<","<<(y+yf)/2.<<", " << buf << ", backgroundcolor=pink, va='top', ha='center', fontsize="<<Prms.TxtSz<<")\n";
                }
            }
        }
    }
    else throw new Fatal("Beam::GetState: 3D Beam is not available yet");
}


////////////////////////////////////////////////////////////////////////////////////////////////// Factory /////


// Allocate a new element
Element * BeamMaker(int NDim, Mesh::Cell const & Cell, Model const * Mdl, SDPair const & Prp, SDPair const & Ini, Array<Node*> const & Nodes) { return new Beam(NDim,Cell,Mdl,Prp,Ini,Nodes); }

// Register element
int BeamRegister()
{
    ElementFactory["Beam"]   = BeamMaker;
    ElementVarKeys["Beam2D"] = std::make_pair ("ux uy wz", "fx fy mz");
    PROB.Set ("Beam", (double)PROB.Keys.Size());
    return 0;
}

// Call register
int __Beam_dummy_int = BeamRegister();

}; // namespace FEM

#endif // MECHSYS_FEM_BEAM
