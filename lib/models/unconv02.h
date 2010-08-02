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

#ifndef MECHSYS_UNCONV02_H
#define MECHSYS_UNCONV02_H

// MechSys
#include <mechsys/models/model.h>

using std::cout;
using std::endl;

class Unconv02 : public Model
{
public:
    // Constructor
    Unconv02 (int NDim, SDPair const & Prms);

    // Derived methods
    void InitIvs    (SDPair const & Ini, State * Sta)                             const;
    void TgIncs     (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const;
    void Stiffness  (State const * Sta, Mat_t & D)                                const;
    bool LoadCond   (State const * Sta, Vec_t const & DEps, double & alpInt)      const;
    void UpdatePath (State const * Sta, Vec_t const & DEps, Vec_t const & DSig)   const;

    // Data
    double         l0,l1,l3,betb,betbb; ///< isotropic compression parameters
    double         k0,k1,betk,ev1;      ///< deviatoric parameters
    mutable double v0,xR10,xR30;        ///< initial values
    mutable double K,G,nu;              ///< bulk and shear moduli, and Poisson's coefficient
    double         M;                   ///< CSL slope
    Vec_t          I;                   ///< 2nd order identity
    mutable Mat_t  De;                  ///< Elastic stiffness
    mutable double alpha;               ///< stress-path variable == atan(dq/dp)

    // Scratchpad
    mutable Vec_t  V,Vb,VDe,DeW,devSig,depsEl; ///< vectors
    mutable double p,q,t,R,A,B,y0;             ///< invariants and gradients
    mutable double lamb,lambb,kapbb,hp;        ///< hardening coefficients

private:
    // internal methods
    void _calc_invariants_and_gradients (EquilibState const * Sta) const ///< calculate: devSig,p,q,t,V,yb,y0,y2
    {
        // invariants
        OctInvs (Sta->Sig, p,q,t);
        Dev     (Sta->Sig, devSig);

        // gradients
        R  = q/(M*p);
        V  = ((R*R-1.0)/Util::SQ3)*I + (2.0/(M*M*p))*devSig;
        y0 = -exp(Sta->Ivs(0));
        //if (R>1.0) throw new Fatal("Unconv02::_calc_invariants_and_gradients: R(%g) must be smaller than 1",R);
    }
    void _calc_hardening (EquilibState const * Sta) const ///< calculate: lamb,lambb,kapbb,Vb,hp (must be after _calc_invariants_and_gradients)
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

        // Vb
        //double distk = (1.0-R);
        double ev    = Calc_ev    (Sta->Eps);
        double ed    = Calc_edoct (Sta->Eps);
        double evr   = ev1 - k1*ed;
        double distk = ev - evr;
        if (distk<0.0) distk = 0.0;
        kapbb = k0+(k1-k0)*exp(-betk*distk);
        //A     = -(1.01-sin(alpha))*lambb;
        //B     = -kapbb*sin(alpha);
        A     = -lambb;
        B     = -kapbb;
        Vb    = (1.0/(3.0*K))*I + (A/y0)*V;
        if (q>1.0e-10) Vb -= (B/(2.0*G*q))*devSig;

        // hardening coefficient
        double ndevW = 2.0*R/M;
        hp = -B*ndevW + trW;

        //cout << "R = " << R << endl;
        //cout << "Vb = " << PrintVector(Vb);
        //cout << "ndevW = " << ndevW << ", hp = " << hp << endl;
    }
    double _calc_yield_function (EquilibState const * Sta) const
    {
        OctInvs (Sta->Sig, p,q,t);
        double pb = exp(Sta->Ivs(0));
        R = q/(M*p);
        return p*(1.0+R*R) - pb;
    }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline Unconv02::Unconv02 (int NDim, SDPair const & Prms)
    : Model (NDim,Prms,"Unconv02"), alpha(0.0)
{
    // parameters
    k0    = Prms("k0");
    k1    = Prms("k1");
    betk  = Prms("betk");
    ev1   = Prms("ev1");
    l0    = Prms("l0");
    l1    = Prms("l1");
    l3    = Prms("l3");
    betb  = Prms("betb");
    betbb = Prms("betbb");
    M     = Phi2M(Prms("phi"));
    nu    = Prms("nu");
    if (l3<l1) throw new Fatal("Unconv02::Unconv02: l3(%g) must be greater than or equal to l1(%g)",l3,l1);

    // constants
    I.change_dim (NCps);
    if (NDim==2) I = 1.0, 1.0, 1.0, 0.0;
    else         I = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

    // internal values
    NIvs = 3;
    IvNames.Push ("z0"); // 0
    IvNames.Push ("z1"); // 1
    IvNames.Push ("v");  // 2: specific volume = 1+e

    // scratchpad
    V     .change_dim (NCps);
    Vb    .change_dim (NCps);
    VDe   .change_dim (NCps);
    DeW   .change_dim (NCps);
    devSig.change_dim (NCps);
    depsEl.change_dim (NCps);
}

inline void Unconv02::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

    // NCL position and specific void
    xR10 = Ini("xR10");
    xR30 = Ini("xR30");
    v0   = Ini("v0");

    // invariants
    OctInvs (sta->Sig, p,q,t);
    R = q/(M*p);

    // calculate K and G
    K = v0*p/(l0*Util::SQ3);
    G = 3.0*K*(1.0-2.0*nu)/(2.0*(1.0+nu));
    //cout << "K = " << K << ", G = " << G << endl;

    // elastic stiffness
    if (GTy==pse_t) throw new Fatal("Unconv02::Unconv02: This model does not work for plane-stress (pse)");
    Mat_t Psd, IdyI;
    Calc_Psd  (NCps,Psd);
    Calc_IdyI (NCps,IdyI);
    De = (2.0*G)*Psd + K*IdyI;

    // internal variables
    double pb = p*(1.0+R*R);
    sta->Ivs(0) = log(pb);
    sta->Ivs(1) = xR30;
    sta->Ivs(2) = v0;

    // check initial yield function
    double f = _calc_yield_function(sta);
    if (f>1.0e-15) throw new Fatal("Unconv02:InitIvs: stress point (sig=(%g,%g,%g,%g], p=%g, q=%g) is outside yield surface (f=%g) with pb=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,p,q,f,pb);
}

inline void Unconv02::TgIncs (State const * Sta, Vec_t & DEps, Vec_t & DSig, Vec_t & DIvs) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    //cout << "DEps = " << PrintVector(DEps);
    //cout << "f    = " << _calc_yield_function(sta) << endl;

    // invariants and gradients
    _calc_invariants_and_gradients (sta); // calculate: devSig,p,q,t,V,yb,y0,y2

    // volume strain increment and internal variables
    double dev = Calc_ev(DEps);
    double v   = sta->Ivs(2);

    // increments
    if (sta->Ldg)
    {
        // hardening
        _calc_hardening (sta);  // calculate: lamb,lambb,kapbb,Vb,hp

        // plastic multiplier
        Mult (Vb, De, VDe);
        double phi = dot(VDe,V) - hp;
        double gam = dot(VDe,DEps)/phi;
        //if (gam<0.0) throw new Fatal("Unconv02::TgIncs: Error: gam(%g)<0 for loading!",gam);
        // TODO: remove Loading method, since we can use this method from now on to check Ldg

        // stress increment
        depsEl = DEps - gam*V;
        DSig   = De*depsEl;

        //double gam_ = -dot(Vb,DSig)/hp;
        //cout << "gam = " << gam << ", gam_ = " << gam_ << endl << endl;;

        // increment of internal values
        double ded = Calc_edoct(DEps);
        DIvs(0) = (dev - B*ded)/A;
        DIvs(1) = (-1.0/lamb)*dev;
        DIvs(2) = v*dev;
    }
    else
    {
        // TODO: check this for unloading
        
        // stress increment
        DSig = De*DEps;

        // increment of internal values
        Vec_t VDe;
        Mult (V, De, VDe);
        DIvs(0) = 0.0;
        DIvs(1) = 0.0;
        DIvs(2) = v*dev;
    }
}

inline void Unconv02::Stiffness (State const * Sta, Mat_t & D) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // stiffness
    if (sta->Ldg)
    {
        // invariants and gradients
        _calc_invariants_and_gradients (sta); // calculate: devSig,p,q,t,V,yb,y0,y2

        // hardening
        _calc_hardening (sta); // calculate: lamb,lambb,kapbb,Vb,hp

        // auxiliar vectors
        Mult (Vb, De, VDe);
        double phi = dot(VDe,V) - hp;
        DeW = De*V;

        //cout << "V     = " << PrintVector(V);
        //cout << "VDe   = " << PrintVector(VDe);
        //cout << "VDe.V = " << dot(VDe,V) << ", phi = " << phi << endl;

        // elastoplastic stiffness
        D.change_dim (NCps, NCps);
        for (size_t i=0; i<NCps; ++i)
        for (size_t j=0; j<NCps; ++j)
            D(i,j) = De(i,j) - DeW(i)*VDe(j)/phi;

        //Mat_t Cep(NCps,NCps);
        //cout << "det(Dep) = " << Det(D) << endl;
        //cout << "Dep = \n" << PrintMatrix(D);
        //Inv (D, Cep);
        //cout << "Cep = \n" << PrintMatrix(Cep);
    }
    else D = De;
}

inline bool Unconv02::LoadCond (State const * Sta, Vec_t const & DEps, double & alpInt) const
{
    // no intersection (never in this model)
    alpInt = -1.0;

    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // invariants and gradients
    _calc_invariants_and_gradients (sta); // calculate: devSig,p,q,t,V,yb,y0,y2

    // check loading condition
    Mult (V, De, VDe);
    double num = dot(VDe,DEps);
    //cout << "VDe  = " << PrintVector(VDe);
    //cout << "DEps = " << PrintVector(DEps);
    //cout << "num  = " << num << endl;
    if (fabs(num)<1.0e-12) return true; // neutral loading
    if (num<0.0) throw new Fatal("Unconv02::LoadCond: num(%g)<0 not ready yet",num);
    return (num>0.0);
}

inline void Unconv02::UpdatePath (State const * Sta, Vec_t const & DEps, Vec_t const & DSig) const
{
    double dq  = Calc_qoct (DSig);
    double dp  = Calc_poct (DSig);
    //double dqc = Calc_qcam (DSig);
    //double dpc = Calc_pcam (DSig);
    alpha = atan2(dq,dp);
    //cout << "alpha = " << 180.0*alpha/Util::PI << ", sin(alpha) = " << sin(alpha) << endl;
    /*
    cout << "dqc = " << dqc << ",  dpc = " << dpc << ", dqc/dpc = " << dqc/dpc << ",  alphac = atan(dqc/dpc) = " << 180.0*atan2(dqc,dpc)/Util::PI << endl;
    cout << "dq  = " << dq  << ",  dp  = " << dp  << ", dq /dp  = " << dq /dp  << ",  alpha  = atan(dq /dp)  = " << 180.0*alpha/Util::PI          << endl;
    cout << endl;
    */
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * Unconv02Maker(int NDim, SDPair const & Prms) { return new Unconv02(NDim,Prms); }

int Unconv02Register()
{
    ModelFactory["Unconv02"] = Unconv02Maker;
    MODEL.Set ("Unconv02", (double)MODEL.Keys.Size());
    return 0;
}

int __Unconv02_dummy_int = Unconv02Register();

#endif // MECHSYS_UNCONV02_H
