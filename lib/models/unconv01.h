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

#ifndef MECHSYS_UNCONV01_H
#define MECHSYS_UNCONV01_H

// MechSys
#include <mechsys/models/elastoplastic.h>

using std::cout;
using std::endl;

class Unconv01 : public Model
{
public:
    // Constructor
    Unconv01 (int NDim, SDPair const & Prms);

    // Derived methods
    void InitIvs   (SDPair const & Ini, State * Sta)                             const;
    void TgIncs    (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const;
    void Stiffness (State const * Sta, Mat_t & D)                                const;
    bool LoadCond  (State const * Sta, Vec_t const & DEps, double & alpInt)      const;

    // Data
    double         l0,l1,l3,betb,betbb; ///< isotropic compression parameters
    mutable double v0,xR10,xR30;        ///< initial values
    double         M;                   ///< CSL slope
    mutable double K,G,nu;              ///< bulk and shear moduli
    Vec_t          I;                   ///< 2nd order identity
    Mat_t          De;                  ///< Elastic stiffness

    // Scratchpad
    mutable Vec_t  V,Vb,VDe,DeW,devSig,depsEl;
    mutable double p,q,t,y0,lamb,lambb,hp;

private:
    // internal methods
    void _calc_invariants_and_gradients (EquilibState const * Sta) const ///< calculate: devSig,p,q,t,V,y0
    {
        // invariants
        OctInvs (Sta->Sig, p,q,t);
        Dev     (Sta->Sig, devSig);

        // gradients
        double z0 = Sta->Ivs(0);
        double pb = exp(z0);
        y0 = -M*M*p*exp(z0);
        V  = (M*M*(pb-2.0*p)/Util::SQ3)*I + 2.0*devSig;
    }
    void _calc_hardening (EquilibState const * Sta) const ///< calculate: lamb,lambb,hp,Vb (must be after _calc_invariants_and_gradients)
    {
        // internal variables
        double z0 = Sta->Ivs(0);
        double z1 = Sta->Ivs(1);
        double v  = Sta->Ivs(2);

        // isotropic coefficients
        double xR1    = xR10 + (log(v0/v))/l1;
        double distb  = z1 - xR1;
        double distbb = z1 - z0;
        double trW    = Tra(V);
        if (distb <0.0) distb  = 0.0;
        if (distbb<0.0) distbb = 0.0;
        lamb  = l3+(l1  -l3)*exp(-betb *distb);
        lambb = l0+(lamb-l0)*exp(-betbb*distbb);
        double H0 = -trW/lambb;
        hp = y0*H0;
        Vb = V + (y0*(-1.0/(3.0*lambb*K)))*I;
    }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Unconv01::Unconv01 (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"Unconv01")
{
    // parameters
    l0    = Prms("l0");
    l1    = Prms("l1");
    l3    = Prms("l3");
    betb  = Prms("betb");
    betbb = Prms("betbb");
    M     = Phi2M(Prms("phi"));
    K     = Prms("K");
    G     = Prms("G");

    // constants
    I.change_dim (NCps);
    if (NDim==2) I = 1.0, 1.0, 1.0, 0.0;
    else         I = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

    // internal values
    NIvs = 3;
    IvNames.Push ("z0");
    IvNames.Push ("z1");
    IvNames.Push ("v"); // specific volume = 1+e

    // elastic stiffness
    if (GTy==pse_t) throw new Fatal("Unconv01::Unconv01: This model does not work for plane-stress (pse)");
    Mat_t Psd,IdyI;
    Calc_Psd  (NCps,Psd);
    Calc_IdyI (NCps,IdyI);
    De = (2.0*G)*Psd + K*IdyI;

    /*
    double E  = 6000.0;
    double nu = 0.3;
    double c  = E/((1.0+nu)*(1.0-2.0*nu));
    De = c*(1.0-nu),       c*nu ,      c*nu ,            0.0,            0.0,            0.0,
              c*nu ,  c*(1.0-nu),      c*nu ,            0.0,            0.0,            0.0,
              c*nu ,       c*nu , c*(1.0-nu),            0.0,            0.0,            0.0,
               0.0 ,        0.0 ,       0.0 , c*(1.0-2.0*nu),            0.0,            0.0,
               0.0 ,        0.0 ,       0.0 ,            0.0, c*(1.0-2.0*nu),            0.0,
               0.0 ,        0.0 ,       0.0 ,            0.0,            0.0, c*(1.0-2.0*nu);
    cout <<"De(KG) =\n" << PrintMatrix(De);
    cout <<"De(Enu) =\n" << PrintMatrix(De);
    */

    // scratchpad
    V     .change_dim (NCps);
    Vb    .change_dim (NCps);
    VDe   .change_dim (NCps);
    DeW   .change_dim (NCps);
    devSig.change_dim (NCps);
    depsEl.change_dim (NCps);
}

inline void Unconv01::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

    // NCL position and specific void
    xR10 = Ini("xR10");
    xR30 = Ini("xR30");
    v0   = Ini("v0");

    // invariants
    double p,q,t;
    OctInvs (sta->Sig, p,q,t);

    // internal variables
    double pb = p+(q*q)/(p*M*M);
    sta->Ivs(0) = log(pb);
    sta->Ivs(1) = xR30;
    sta->Ivs(2) = v0;

    // check initial yield function
    double f = q*q + (p-pb)*p*M*M;
    if (f>1.0e-8) throw new Fatal("Unconv01:InitIvs: stress point (sig=(%g,%g,%g,%g], p=%g, q=%g) is outside yield surface (f=%g) with pb=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,p,q,f,pb);
}

inline void Unconv01::TgIncs (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // invariants and gradients
    _calc_invariants_and_gradients (sta); // calc: devSig,p,q,t,V,y0

    // volume strain increment and internal variables
    double dev = Calc_ev(DEps);
    //double z0  = sta->Ivs(0);
    double v   = sta->Ivs(2);
    //double pb  = exp(z0);
    //double f   = q*q + (p-pb)*p*M*M;
    //cout << "f(before) = " << f << endl;

    // increments
    if (sta->Ldg)
    {
        // hardening
        _calc_hardening (sta); // calc: lamb,lambb,hp,Vb

        // plastic multiplier
        Mult (Vb, De, VDe);
        double phi = dot(VDe,V) - hp;
        double gam = dot(VDe,DEps)/phi;

        // stress increment
        depsEl = DEps - gam*V;
        DSig   = De*depsEl;

        //double gam_ = -dot(Vb,DSig)/hp;
        //cout << "gam = " << gam << ", gam_ = " << gam_ << endl;

        // increment of internal values
        DIvs(0) = (-1.0/lambb)*dev;
        DIvs(1) = (-1.0/lamb )*dev;
        DIvs(2) = v*dev;

        //Vec_t hh0((-1.0/(3.0*lambb*K))*I);
        //double dz0 = dot(hh0,DSig) - Tra(V)*gam/lambb;
        //cout << "DIvs(0) = " << DIvs(0) << ", dz0 = " << dz0 << endl;
        //DIvs(0) = dz0;

        //Vec_t sigf(sta->Sig+DSig);
        //double z0f = z0+DIvs(0);
        //double pbf = exp(z0f);
        //double pf  = Calc_poct(sigf);
        //double qf  = Calc_qoct(sigf);
        //double ff  = qf*qf + (pf-pbf)*pf*M*M;
        //cout << "f(after)  = " << ff << endl; cout << endl;
    }
    else
    {
        // stress increment
        DSig = De*DEps;

        // increment of internal values
        Vec_t VDe;
        Mult (V, De, VDe);
        DIvs(0) = -dot(V,DSig)/y0;
        DIvs(1) = 0.0;
        DIvs(2) = v*dev;
    }
}

inline void Unconv01::Stiffness (State const * Sta, Mat_t & D) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // stiffness
    if (sta->Ldg)
    {
        // invariants and gradients
        _calc_invariants_and_gradients (sta); // calc: devSig,p,q,t,V,y0

        // hardening
        _calc_hardening (sta); // calc: lamb,lambb,hp,Vb

        // auxiliar vectors
        Mult (Vb, De, VDe);
        double phi = dot(VDe,V) - hp;
        DeW = De*V;

        // elastoplastic stiffness
        D.change_dim (NCps, NCps);
        for (size_t i=0; i<NCps; ++i)
        for (size_t j=0; j<NCps; ++j)
            D(i,j) = De(i,j) - DeW(i)*VDe(j)/phi;
    }
    else D = De;
}

inline bool Unconv01::LoadCond (State const * Sta, Vec_t const & DEps, double & alpInt) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // invariants and gradients
    _calc_invariants_and_gradients (sta); // calc: devSig,p,q,t,V,y0

    // check loading condition
    Mult (V, De, VDe);
    double num = dot(VDe,DEps);
    if (num<0.0) throw new Fatal("Unconv03::LoadCond: num(%g)<0 not ready yet",num);

    alpInt = -1.0; // no intersection (never in this model)
    return (num>0.0);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * Unconv01Maker(int NDim, SDPair const & Prms) { return new Unconv01(NDim,Prms); }

int Unconv01Register()
{
    ModelFactory["Unconv01"] = Unconv01Maker;
    MODEL.Set ("Unconv01", (double)MODEL.Keys.Size());
    return 0;
}

int __Unconv01_dummy_int = Unconv01Register();

#endif // MECHSYS_UNCONV01_H
