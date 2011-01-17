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

#ifndef MECHSYS_CAMCLAY_H
#define MECHSYS_CAMCLAY_H

// MechSys
#include <mechsys/models/elastoplastic.h>

class CamClay : public ElastoPlastic
{
public:
    // Constructor
    CamClay (int NDim, SDPair const & Prms);

    // Derived methods
    void   InitIvs   (SDPair const & Ini, State * Sta) const;
    void   Gradients (EquilibState const * Sta)        const;
    void   Hardening (EquilibState const * Sta)        const;
    double YieldFunc (EquilibState const * Sta)        const;
    double FailCrit  (EquilibState const * Sta)        const;
    double CalcE     (EquilibState const * Sta) const { return fabs(Sta->Sig(0)+Sta->Sig(1)+Sta->Sig(2))*(1.0-2.0*nu)*v0/kap; }
    //double CalcE     (EquilibState const * Sta) const { return 6000.0; }

    // Data
    double         lam;
    double         kap;
    double         phi;
    Vec_t          I;
    mutable double v0;
    mutable double chi;
    double         Mcs;
    double         wcs;


    double CalcM (double const & sin3th) const;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline CamClay::CamClay (int NDim, SDPair const & Prms)
    : ElastoPlastic (NDim,Prms,/*derived*/true)
{
    Name = "CamClay";

    // parameters
    lam = Prms("lam");
    kap = Prms("kap");
    nu  = Prms("nu");
    phi = Prms("phi");
    Mcs = Phi2M(phi,"oct");
    wcs = 0;//pow((3.0-sin(phiRad))/(3.0+sin(phiRad)),4.0);

    // constants
    I.change_dim (NCps);
    if (NDim==2) I = 1.0, 1.0, 1.0, 0.0;
    else         I = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;

    // internal values
    NIvs = 1;
    if (NewSU) NIvs = 2;
    Y.change_dim (NIvs);
    H.change_dim (NIvs);
    IvNames.Push ("z0");
    if (NewSU) IvNames.Push ("z1");
}

inline void CamClay::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

    // specific void
    v0  = Ini("v0");
    chi = (kap-lam)/v0;

    // invariants
    double p,q,t;
    OctInvs (sta->Sig, p,q,t);

    // internal variables
    double M   = CalcM(t);
    double p0  = p+(q*q)/(p*M*M);
    double OCR = (Ini.HasKey("OCR") ? Ini("OCR") : 1.0);
    sta->Ivs(0) = OCR*p0;

    if (NewSU)
    {
        sta->Ivs(1) = sta->Ivs(0);
        sta->Ivs(0) = p0;
    }

    // check initial yield function
    double f = YieldFunc (sta);
    if (f>1.0e-8) throw new Fatal("CamClay:InitIvs: stress point (sig=(%g,%g,%g,%g], p=%g, q=%g) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,p,q,f,sta->Ivs(0));
}

inline void CamClay::Gradients (EquilibState const * Sta) const
{
    // invariants
    double p,q,t;
    Vec_t dev_sig;
    OctInvs (Sta->Sig, p,q,t);
    Dev     (Sta->Sig, dev_sig);

    // gradients
    double M = CalcM(t);
    V    = (M*M*(Sta->Ivs(0)-2.0*p)/(Util::SQ3))*I + 2.0*dev_sig;
    Y(0) = -M*M*p;

    if (NewSU) Y(1) = 0.0;

    if (false)//q>1.0e-10)
    {
        double ss3th = pow(t,2.0);
        double cos3th;
        if (ss3th>1.0) cos3th = 0.0;
        else           cos3th = sqrt(1.0-ss3th);
        if (cos3th>1.0e-10)
        {
            Vec_t SS,dev_SS,dth_dsig;
            double dM_dth = (0.75*M*(1.0-wcs)*cos3th) / (1.0+wcs+(wcs-1.0)*t);
            Pow2 (dev_sig, SS);
            dev_SS   = SS - (I * (SS(0)+SS(1)+SS(2))/3.0);
            dth_dsig = (1.5/(q*q*cos3th)) * (dev_SS*(3.0/q) - dev_sig*t);
            V       += dth_dsig*(2.0*M*p*(Sta->Ivs(0)-p)*dM_dth);
        }
    }

}

inline void CamClay::Hardening (EquilibState const * Sta) const
{
    H(0) = Sta->Ivs(0)*Tra(W)/chi;

    if (NewSU)
    {
        double D = 2.0*Sta->Ivs(1)/(Sta->Ivs(1)+Sta->Ivs(0))-1.0;
        if (D<0.0) D = 0.0;
        H(1) = exp(-BetSU*D)*H(0); //  (D>0.0 ? 0.0 : H(0));
        H(0) = H(0) + AlpSU*(1.0-exp(-BetSU*D));
        //printf("z0=%g, z1=%g, D=%g, H0=%g, H1=%g\n",Sta->Ivs(0),Sta->Ivs(1),D,H(0),H(1));
    }
}

inline double CamClay::YieldFunc (EquilibState const * Sta) const
{
    double p,q,t;
    OctInvs (Sta->Sig, p,q,t);
    double M = CalcM(t);
    return q*q + (p - Sta->Ivs(0))*p*M*M;
}

inline double CamClay::FailCrit (EquilibState const * Sta) const
{
    double p,q,t;
    OctInvs (Sta->Sig, p,q,t);
    double M = CalcM(t);
    return q - M*p;
}

inline double CamClay::CalcM (double const & sin3th) const
{
    return Mcs;//*pow(2.0*wcs/(1.0+wcs+(wcs-1.0)*sin3th),0.25);
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * CamClayMaker(int NDim, SDPair const & Prms) { return new CamClay(NDim,Prms); }

int CamClayRegister()
{
    ModelFactory   ["CamClay"] = CamClayMaker;
    MODEL.Set      ("CamClay", (double)MODEL.Keys.Size());
    MODEL_PRM_NAMES["CamClay"].Resize(7);
    MODEL_PRM_NAMES["CamClay"] = "lam", "kap", "nu", "phi", "newsu", "betsu", "alpsu";
    MODEL_IVS_NAMES["CamClay"].Resize(2);
    MODEL_IVS_NAMES["CamClay"] = "v0", "OCR";
    return 0;
}

int __CamClay_dummy_int = CamClayRegister();

#endif // MECHSYS_CAMCLAY_H
