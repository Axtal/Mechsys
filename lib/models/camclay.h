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
    void   InitIvs   (SDPair const & Ini, State * Sta)                            const;
    void   Gradients (Vec_t const & Sig, Vec_t const & Ivs, bool Potential=false) const;
    void   Hardening (Vec_t const & Sig, Vec_t const & Ivs)                       const;
    double YieldFunc (Vec_t const & Sig, Vec_t const & Ivs)                       const;
    double CalcE     (Vec_t const & Sig, Vec_t const & Ivs)                       const
    { return fabs(Sig(0)+Sig(1)+Sig(2))*(1.0-2.0*nu)*v0/kap; }
    /*
    {
        double pcam = fabs(Sig(0)+Sig(1)+Sig(2))/3.0;
        double K    = pcam*v0/kap;
        return Calc_E_ (K, nu);
    }
    */

    // Data
    double         lam;
    double         kap;
    double         phi;
    mutable double v0;
    mutable double chi;
    mutable double Mcs;
    mutable double wcs;

    // Internal methods
    double Calc_M (double sin3th) const { return Mcs*pow(2.0*wcs/(1.0+wcs-(1.0-wcs)*sin3th),1.0/4.0); }
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline CamClay::CamClay (int NDim, SDPair const & Prms)
    : ElastoPlastic (NDim, Prms, /*niv*/4, "CamClay", /*derived*/true)
{
    // parameters
    lam = Prms("lam");
    kap = Prms("kap");
    nu  = Prms("nu");
    phi = Prms("phi");
    Mcs = Phi2M(phi,"oct");
    double phi_rad = phi*Util::PI/180.0;
    wcs = pow((3.0-sin(phi_rad))/(3.0+sin(phi_rad)),4.0);

    // internal values
    IvNames = "z0", "z1", "evp", "edp";
}

inline void CamClay::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, NIvs);

    // specific void
    v0  = Ini("v0");
    chi = (kap-lam)/v0;

    // internal variables
    Calc_pqt (sta->Sig);
    double M   = Calc_M (t);
    double p0  = p+(q*q)/(p*M*M);
    double OCR = (Ini.HasKey("OCR") ? Ini("OCR") : 1.0);
    if (NewSU)
    {
        sta->Ivs(0) = p0;
        sta->Ivs(1) = OCR*p0;
    }
    else sta->Ivs(0) = OCR*p0;

    // check initial yield function
    double f = YieldFunc (sta->Sig, sta->Ivs);
    if (f>FTol)           throw new Fatal("CamClay:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
    if (NewSU && f<-FTol) throw new Fatal("CamClay:InitIvs: stress point (sig=(%g,%g,%g,%g]) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,f,sta->Ivs(0));
}

inline void CamClay::Gradients (Vec_t const & Sig, Vec_t const & Ivs, bool Potential) const
{
    OctInvs (Sig, p,q,t,th,s, qTol, &dthdsig);
    double M = Calc_M (t);
    Vec_t * VorW = &V;
    if (Potential) VorW = &W;
    else
    {
        Y(0) = -M*M*p;
        Y(1) = 0.0;
        Y(2) = 0.0;
        Y(3) = 0.0;
    }
    double dfdp = M*M*(2.0*p-Ivs(0));
    double m    = -dfdp/Util::SQ3;
    if (q>qTol)
    {
        double dMdt  = 0.25*M*(1.0-wcs)/(1.0+wcs-(1.0-wcs)*t);
        double dtdth = 3.0*cos(3.0*th);
        double dMdth = dMdt*dtdth;
        double dfdth = -2.0*M*dMdth*(Ivs(0)-p);
        (*VorW) = 2.0*s + m*I + dfdth*dthdsig;
    }
    else (*VorW) = m*I;
}

inline void CamClay::Hardening (Vec_t const & Sig, Vec_t const & Ivs) const
{
    H(0) = Ivs(0)*Tra(W)/chi;
    H(1) = 0.0;
    H(2) = 0.0;
    H(3) = 0.0;
    if (NewSU)
    {
        double D = 2.0*Ivs(1)/(Ivs(1)+Ivs(0))-1.0;
        if (D<0.0) D = 0.0;
        H(1) = exp(-BetSU*D)*H(0); //  (D>0.0 ? 0.0 : H(0));
        H(0) = H(0) + AlpSU*(1.0-exp(-BetSU*D));
    }
}

inline double CamClay::YieldFunc (Vec_t const & Sig, Vec_t const & Ivs) const
{
    Calc_pqt (Sig);
    double M = Calc_M (t);
    return q*q + (p-Ivs(0))*p*M*M;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * CamClayMaker(int NDim, SDPair const & Prms, Model const * O) { return new CamClay(NDim,Prms); }

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
