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
    void InitIvs   (SDPair const & Ini, State * Sta)                                                    const;
    void TgIncs    (State const * Sta, double pw, Vec_t & DEps, double Dpw, Vec_t & DSig, Vec_t & DIvs) const;
    void Stiffness (State const * Sta, double pw, Mat_t & D, Vec_t & d)                                 const;
    bool LoadCond  (State const * Sta, double pw, Vec_t const & DEps, double Dpw, double & alpInt)      const;

    // Internal methods
    void Gradients  (EquilibState const * Sta, double pc) const;
    void Hardening  (EquilibState const * Sta, double pc) const;
    void Calc_De_de (EquilibState const * Sta, double pc) const;

    // Data
    double  qTol;       ///< Tolerance for minimum qoct
    double lam0,kap,phi,nu,r,bet,kc,lams,kaps,B,c0,c1;
    double Mcs,wcs,pref;

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

    // Auxiliary methods
    double Calc_lam (double pc)            const { return lam0*((1.0-r)*exp(-bet*pc)+r);                                    }
    double Calc_psi (double pc)            const { return (lam0-kap)/(Calc_lam(pc)-kap);                                    }
    double Calc_ps  (double pc)            const { return (pc>0.0 ? -kc*pc : -pc);                                          }
    double Calc_py  (double z0, double pc) const { return (pc>0.0 ? pref*pow(z0/pref, Calc_psi(pc)) : z0-pc);               }
    double Calc_C   (double z1, double pc) const { return (pc>0.0 ? pref*pref*(exp(B*(pc-z1)/pref)-exp(-B*z1/pref)) : 0.0); } 
    double Calc_M   (double sin3th)        const { return Mcs*pow(2.0*wcs/(1.0+wcs-(1.0-wcs)*sin3th),1.0/4.0);              }
    void   Calc_pqt (Vec_t const & Sig)    const { OctInvs(Sig,p,q,t);  th=asin(t)/3.0; }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline BBMx::BBMx (int NDim, SDPair const & Prms)
    : Model (NDim, Prms, /*niv*/4, "BBMx"), pref(1.0)
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
    Mcs = 6.0*sin_phi/(3.0-sin_phi);
    wcs = pow((3.0-sin_phi)/(3.0+sin_phi),4.0);

    // internal values
    IvNames = "z0", "z1", "z2", "z3";

    // check
    if (GTy==pse_t) throw new Fatal("BBMx::BBMx: This model does not work for plane-stress (pse)");

    // allocate memory
    V  .change_dim (NCps);
    W  .change_dim (NCps);
    Y  .change_dim (NIvs);
    H  .change_dim (NIvs);
    De .change_dim (NCps,NCps);
    Dep.change_dim (NCps,NCps);
    VDe.change_dim (NCps);
    DeW.change_dim (NCps);
}

inline void BBMx::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

}

inline void BBMx::TgIncs (State const * Sta, double pw, Vec_t & DEps, double Dpw, Vec_t & DSig, Vec_t & DIvs) const
{
    // state
    //EquilibState const * sta = static_cast<EquilibState const *>(Sta);
}

inline void BBMx::Stiffness (State const * Sta, double pw, Mat_t & D, Vec_t & d) const
{
    // state
    //EquilibState const * sta = static_cast<EquilibState const *>(Sta);
}

inline bool BBMx::LoadCond (State const * Sta, double pw, Vec_t const & DEps, double Dpw, double & alpInt) const
{
    // state
    //EquilibState const * sta = static_cast<EquilibState const *>(Sta);
    return false;
}

inline void BBMx::Gradients (EquilibState const * Sta, double pc) const
{
    // variables
    OctInvs (Sta->Sig, p,q,t,th,s, qTol, &dthdsig);
    double lam = Calc_lam (pc);
    double psi = Calc_psi (pc);
    double ps  = Calc_ps  (pc);
    double py  = Calc_py  (Sta->Ivs(0), pc);

    // auxiliary variables
    double M  = Calc_M (t);
    double MM = M*M;
    double m  = MM*(2.0*p-ps-py)/3.0;

    // gradients
    double dlamdpc = -bet*lam0*(1.0-r)*exp(-bet*pc);
    double dpsidpc = (-psi/(lam-kap))*dlamdpc;
    double dpsdpc  = (pc>0.0 ? -kc : -1.0);
    double dpydpc  = (pc>0.0 ? py*log(Sta->Ivs(0)/pref)*dpsidpc : -1.0);
    double dCdpc   = (pc>0.0 ? B*pref*exp(B*(pc-Sta->Ivs(1))/pref) : 0.0);
    S = MM*((py-p)*dpsdpc+(ps-p)*dpydpc) + dCdpc;

    // Gradients
    if (q>qTol)
    {
        double dfdM  = -2.0*M*(p-ps)*(py-p);
        double dMdt  = 0.25*M*(1.0-wcs)/(1.0+wcs-(1.0-wcs)*t);
        double dtdth = 3.0*cos(3.0*th);
        V = 3.0*s + m*I + (dfdM*dMdt*dtdth)*dthdsig;
    }
    else V = m*I;

    // internal variables
    double dpydz0 = (pc>0.0 ? pref*psi*pow(Sta->Ivs(0)/pref,psi)/Sta->Ivs(0) : 1.0);
    Y(0) = -MM*(p-ps)*dpydz0;
    Y(1) = (pc>0.0 ? -B*Calc_C(Sta->Ivs(1),pc)/pref : 0.0);
    W    = V;
}

inline void BBMx::Hardening (EquilibState const * Sta, double pc) const
{
    double trW   = W(0)+W(1)+W(2);
    double chi0  = (lam0 - kap )/nu;
    double chis  = (lams - kaps)/nu;
    H(0) =  Sta->Ivs(0)      *trW/chi0;
    //H(1) = (Sta->Ivs(1)+patm)*trW/chis;
    H(1) = Sta->Ivs(1)*trW/chis;
}

inline void BBMx::Calc_De_de (EquilibState const * Sta, double pc) const
{
    // Invariants
    double p = (Sta->Sig(0)+Sta->Sig(1)+Sta->Sig(2))/3.0;

    // Calculate Bulk modulus
    double K  =  p       *nu/kap;
    //double Ks = (pc+patm)*v/kaps;
    double Ks = pc*nu/kaps;

    // Calc De: elastic stiffness
    //Tensors::AddScaled(2.0*_G,Tensors::Psd, K,Tensors::IdyI, De);

    // Calc de: elastic stiffness
    de = (-K/Ks) * I;
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
