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

#ifndef MECHSYS_BBMX_H
#define MECHSYS_BBMX_H

// MechSys
#include <mechsys/models/model.h>

class BBMx : public Model
{
public:
    // Constructor
    BBMx (int NDim, SDPair const & Prms);

    // Derived methods
    void   InitIvs      (SDPair const & Ini, State * Sta)                                                    const;
    void   TgIncs       (State const * Sta, double pw, Vec_t & DEps, double Dpw, Vec_t & DSig, Vec_t & DIvs) const;
    void   Stiffness    (State const * Sta, double pw, Mat_t & D, Vec_t & d)                                 const;
    bool   LoadCond     (State const * Sta, double pw, Vec_t const & DEps, double Dpw)                       const;
    size_t CorrectDrift (State * Sta, double & pw)                                                           const;

    // Internal methods
    void Gradients  (Vec_t const & Sig, Vec_t const & Ivs, double pc) const;
    void Hardening  (Vec_t const & Sig, Vec_t const & Ivs, double pc) const;
    void Calc_De_de (Vec_t const & Sig, Vec_t const & Ivs, double pc) const;

    // Parameters and constants
    double lam0,kap,phi,nu,r,bet,kc,lams,kaps,B,c0,c1;
    double Mcs,wcs,pref;

    // Initial values
    mutable double v0;

    // Data
    double qTol; ///< Tolerance for minimum qoct

    // Scratchpad
    mutable double p, q, t, th;  ///< Invariants
    mutable Vec_t  dthdsig;      ///< Derivative of theta w.r.t sigma
    mutable Vec_t  s;            ///< Deviator of sigma
    mutable Vec_t  V;            ///< NCps: Gradient of the yield surface
    mutable double S;            ///< dfdpc
    mutable Vec_t  W;            ///< NCps: Plastic flow rule direction
    mutable Vec_t  Y;            ///< NIvs: Derivative of the yield surface w.r.t internal variables
    mutable Vec_t  H;            ///< NIvs: Hardening coefficients, one for each internal variable
    mutable Mat_t  De;           ///< Elastic stiffness
    mutable Vec_t  de;           ///< Elastic stiffness
    mutable Mat_t  Dep;          ///< Elastoplastic stiffness
    mutable Vec_t  VDe;          ///< V*De
    mutable Vec_t  DeW;          ///< De*W
    mutable Vec_t  DEpsEl;       ///< Elastic strain increment
    mutable Vec_t  DEpsPl;       ///< Plastic strain increment == Lam*W
    mutable Vec_t  DSigTr;       ///< Trial (elastic) increment (used in LoadCond)

    // Auxiliary methods
    double Calc_lam (double pc)            const { return (pc>0.0 ? lam0*((1.0-r)*exp(-bet*pc)+r) : lam0);                  }
    double Calc_psi (double pc)            const { return (pc>0.0 ? (lam0-kap)/(Calc_lam(pc)-kap) : 1.0);                   }
    double Calc_ps  (double pc)            const { return (pc>0.0 ? kc*pc : pc);                                            }
    double Calc_p0  (double z0, double pc) const { return (pc>0.0 ? pref*pow(z0/pref, Calc_psi(pc)) : z0-pc);               }
    double Calc_C   (double z1, double pc) const { return (pc>0.0 ? pref*pref*(exp(B*(pc-z1)/pref)-exp(-B*z1/pref)) : 0.0); } 
    double Calc_M   (double sin3th)        const { return Mcs*pow(2.0*wcs/(1.0+wcs-(1.0-wcs)*sin3th),1.0/4.0);              }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline BBMx::BBMx (int NDim, SDPair const & Prms)
    : Model (NDim, Prms, /*niv*/4, "BBMx"), pref(1.0), qTol(1.0e-7)
{
    // set hm model
    HMCoup = true;

    // parameters
    lam0  = Prms("lam0");
    kap   = Prms("kap");
    phi   = Prms("phi");
    nu    = Prms("nu");
    r     = Prms("r");
    bet   = Prms("bet");
    kc    = Prms("kc");
    lams  = Prms("lams");
    kaps  = Prms("kaps");
    B     = Prms("B");
    c0    = Prms("c0");
    c1    = Prms("c1");

    // constants
    double sin_phi = sin(phi*Util::PI/180.0);
    Mcs = Phi2M(phi,"oct");
    wcs = pow((3.0-sin_phi)/(3.0+sin_phi),4.0);

    // internal values
    IvNames = "z0", "z1", "z2", "z3";

    // check
    if (GTy==pse_t) throw new Fatal("BBMx::BBMx: This model does not work for plane-stress (pse)");

    // allocate memory
    s      .change_dim (NCps);
    dthdsig.change_dim (NCps);
    V      .change_dim (NCps);
    W      .change_dim (NCps);
    Y      .change_dim (NIvs);
    H      .change_dim (NIvs);
    De     .change_dim (NCps,NCps);
    de     .change_dim (NCps);
    Dep    .change_dim (NCps,NCps);
    VDe    .change_dim (NCps);
    DeW    .change_dim (NCps);
    DEpsEl .change_dim (NCps);
    DEpsPl .change_dim (NCps);
    DSigTr .change_dim (NCps);
}

inline void BBMx::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);
    v0 = Ini("v0");

    // invariants
    OctInvs (sta->Sig, p,q,t,th,s, qTol);
    double M = Calc_M (t);

    // internal variables
    double pc0 = -Ini("pw");
    double p0  = p+(q*q)/(p*M*M);
    double OCR = (Ini.HasKey("OCR") ? Ini("OCR") : 1.0);
    double OSI = (Ini.HasKey("OSI") ? Ini("OSI") : 1.0);
    sta->Ivs(0) = p0;
    sta->Ivs(1) = pc0;
    sta->Ivs(2) = OCR*p0;
    sta->Ivs(3) = sta->Ivs(1) + OSI;
}

inline void BBMx::TgIncs (State const * Sta, double pw, Vec_t & DEps, double Dpw, Vec_t & DSig, Vec_t & DIvs) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // zero internal values
    DIvs.change_dim (NIvs);
    set_to_zero     (DIvs);

    // De and de: elastic stiffness
    double pc  = -pw;
    double Dpc = -Dpw;
    Calc_De_de (sta->Sig, sta->Ivs, pc);

    // increments
    if (sta->Ldg)
    {
        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs, pc);
        Hardening (sta->Sig, sta->Ivs, pc);
        double hp = Y(0)*H(0) + Y(1)*H(1);

        // plastic multiplier
        Mult (V, De, VDe);
        double phi = dot(VDe,W) - hp;
        double Lam = (dot(VDe,DEps) + dot(V,de)*Dpc + S*Dpc)/phi;

        // increments
        DEpsPl = Lam*W;
        DEpsEl = DEps - DEpsPl;
        DSig   = De*DEpsEl;

        // increment of internal values
        DIvs = Lam*H;
    }
    else
    {
        // stress increment
        DSig = De*DEps;

        // new stress update
        DIvs(0) = -dot(V, DSig) / Y(0);
        DIvs(1) = Dpc;
    }
}

inline void BBMx::Stiffness (State const * Sta, double pw, Mat_t & D, Vec_t & d) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // De and de: elastic stiffness
    double pc = -pw;
    Calc_De_de (sta->Sig, sta->Ivs, pc);

    // stiffness
    if (sta->Ldg)
    {
        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs, pc);
        Hardening (sta->Sig, sta->Ivs, pc);
        double hp = Y(0)*H(0) + Y(1)*H(1);

        // auxiliary vectors
        Mult (V, De, VDe);
        double phi = dot(VDe,W) - hp;
        DeW = De*W;

        // elastoplastic stiffness
        D.change_dim (NCps, NCps);
        for (size_t i=0; i<NCps; ++i)
        for (size_t j=0; j<NCps; ++j)
            D(i,j) = De(i,j) - DeW(i)*VDe(j)/phi;

        // dep
        d = de - ((dot(V,de)+S)/phi)*DeW;
    }
    else
    {
        D = De;
        d = de;
    }
}

inline bool BBMx::LoadCond (State const * Sta, double pw, Vec_t const & DEps, double Dpw) const
{
    // state
    EquilibState const * sta = static_cast<EquilibState const *>(Sta);

    // numerator of Lagrange multiplier
    double pc  = -pw;
    double dpc = -Dpw;
    Calc_De_de (sta->Sig, sta->Ivs, pc);
    Gradients  (sta->Sig, sta->Ivs, pc);
    DSigTr = De * DEps;
    double numL = dot(V, DSigTr) + dot(V,de)*dpc + S*dpc;
    if (numL>0.0) return true;
    else          return false;
}

inline size_t BBMx::CorrectDrift (State * Sta, double & pw) const
{
    // state
    //EquilibState * sta = static_cast<EquilibState *>(Sta);

    // iterations
    //double fnew = YieldFunc (sta->Sig, sta->Ivs, pw);
    size_t it   = 0;
    /*
    while (fnew>DCFTol && it<DCMaxIt)
    {
        // gradients, flow rule, hardening, and hp
        Gradients (sta->Sig, sta->Ivs, -pw);
        Hardening (sta->Sig, sta->Ivs, -pw);
        double hp = Y(0)*H(0) + Y(1)*H(1);

        // elastic stiffness
        if (it==0) Calc_De_de (sta->Sig, sta->Ivs, -pw);

        // auxiliary vectors
        Mult (V, De, VDe);
        double dgam = fnew/(dot(VDe,W)-hp);
        DeW = De*W;

        // update stress and ivs (only)
        sta->Sig -= dgam*DeW;
        sta->Ivs += dgam*H;
        fnew = YieldFunc (sta->Sig, sta->Ivs);
        //if (NewSU) fnew = fabs(fnew); // TODO: should do this, but it does not converge

        // check convergence
        if (fabs(fnew)<DCFTol) break;
        it++;
    }

    // check number of iterations
    if (it>=DCMaxIt) throw new Fatal("BBMx::CorrectDrift: Yield surface drift correction did not converge after %d iterations (fnew=%g, DCFTol=%g)",it,fnew,DCFTol);
    */
    return it;
}

inline void BBMx::Gradients (Vec_t const & Sig, Vec_t const & Ivs, double pc) const
{
    // variables
    OctInvs (Sig, p,q,t,th,s, qTol, &dthdsig);
    double lam = Calc_lam (pc);
    double psi = Calc_psi (pc);
    double ps  = Calc_ps  (pc);
    double p0  = Calc_p0  (Ivs(0), pc);

    // auxiliary variables
    double M    = Calc_M (t);
    double MM   = M*M;
    double dfdp = MM*(2.0*p+ps-p0);
    double m    = -dfdp/Util::SQ3;

    // gradients: pc
    double dCdpc  =  0.0;
    double dp0dpc = -1.0;
    double dpsdpc =  1.0;
    if (pc>0.0)
    {
        double dlamdpc = -bet*lam0*(1.0-r)*exp(-bet*pc);
        double dpsidpc = (-psi/(lam-kap))*dlamdpc;
        dCdpc  = pref*B*exp(B*(pc-Ivs(1))/pref);
        dp0dpc = p0*log(Ivs(0)/pref)*dpsidpc;
        dpsdpc = kc;
    }
    S = dCdpc - MM*dpsdpc*(p0-p) - MM*dp0dpc*(ps+p);

    // gradients: sig
    if (q>qTol)
    {
        double dMdt  = 0.25*M*(1.0-wcs)/(1.0+wcs-(1.0-wcs)*t);
        double dtdth = 3.0*cos(3.0*th);
        double dMdth = dMdt*dtdth;
        double dfdth = -2.0*M*dMdth*(p0-p)*(ps+p);
        V = 2.0*s + m*I + dfdth*dthdsig;
    }
    else V = m*I;
    W = V;

    // internal variables
    double dCdz1  = 0.0;
    double dp0dz0 = 1.0;
    if (pc>0.0)
    {
        dCdz1  = -B*Calc_C(Ivs(1),pc)/pref;
        dp0dz0 = pref*psi*pow(Ivs(0)/pref,psi)/Ivs(0);
    }
    Y(0) = -MM*dp0dz0*(ps+p);
    Y(1) = dCdz1;
    Y(2) = 0.0;
    Y(3) = 0.0;
}

inline void BBMx::Hardening (Vec_t const & Sig, Vec_t const & Ivs, double pc) const
{
    double trW   = W(0)+W(1)+W(2);
    double chi0  = -(lam0 - kap )/v0;
    double chis  = -(lams - kaps)/v0;
    double L0    = c0*pow(log(    Ivs(2)/     Ivs(0)) ,2.0);
    double L1    = c1*pow(log(1.0+Ivs(3)/(1.0+Ivs(1))),2.0);
    H(2) = Ivs(2)*trW/chi0;
    H(3) = Ivs(3)*trW/chis;
    H(0) = Ivs(0)*(trW+L0)/chi0;
    H(1) = Ivs(1)*(trW+L1)/chis;
}

inline void BBMx::Calc_De_de (Vec_t const & Sig, Vec_t const & Ivs, double pc) const
{
    // bulk modulus
    double pcam = fabs(Sig(0)+Sig(1)+Sig(2))/3.0;
    double K    = pcam*v0/kap;
    double Kc   = pc  *v0/kaps;

    // elastic stiffness
    double G = Calc_G_ (K, nu);
    De = (2.0*G)*Psd + K*IdyI;
    if (pc>0.0) de = (-K/Kc)*I;
    else set_to_zero(de);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * BBMxMaker(int NDim, SDPair const & Prms, Model const * O) { return new BBMx(NDim,Prms); }

int BBMxRegister()
{
    ModelFactory   ["BBMx"] = BBMxMaker;
    MODEL.Set      ("BBMx", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["BBMx"].Resize(12);
    MODEL_PRM_NAMES["BBMx"] = "lam0", "kap", "phi", "nu", "r", "bet", "kc", "lams", "kaps", "B", "c0", "c1";
    MODEL_IVS_NAMES["BBMx"].Resize(3);
    MODEL_IVS_NAMES["BBMx"] = "v0", "OCR", "OSI";
    return 0;
}

int __BBMx_dummy_int = BBMxRegister();

#endif // MECHSYS_BBMX_H
