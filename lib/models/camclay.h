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
#include "models/elastoplastic.h"

using std::cout;
using std::endl;

class CamClay : public ElastoPlastic
{
public:
    // Constructor
    CamClay (int NDim, SDPair const & Prms);

    // Derived methods
    void   InitIvs   (SDPair const & Ini, State * Sta)                      const;
    void   Gradients (EquilibState const * Sta, Vec_t & V, Vec_t & Y)       const;
    void   Hardening (EquilibState const * Sta, Vec_t const & W, Vec_t & H) const { H.change_dim(1);  H(0) = Sta->Ivs(0)*Tra(W)/chi; }
    double YieldFunc (EquilibState const * Sta)                             const;
    double FailCrit  (EquilibState const * Sta)                             const;

    // Data
    double         lam;
    double         kap;
    double         nu;
    double         phi;
    double         M;
    Vec_t          I;
    mutable double chi;
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline CamClay::CamClay (int NDim, SDPair const & Prms)
    : ElastoPlastic (NDim,Prms)
{
    lam = Prms("lam");
    kap = Prms("kap");
    nu  = Prms("nu");
    phi = Prms("phi");
    M   = Phi2M(phi);
    I.change_dim(2*NDim);
    if (NDim==2) I = 1.0, 1.0, 1.0, 0.0;
    else         I = 1.0, 1.0, 1.0, 0.0, 0.0, 0.0;
}

inline void CamClay::InitIvs (SDPair const & Ini, State * Sta) const
{
    // initialize state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    sta->Init (Ini, /*NIvs*/1);

    // specific void
    double v0 = Ini("v0");
    chi = (kap-lam)/v0;

    // invariants
    double p,q,t;
    OctInvs (sta->Sig, p,q,t);

    // internal variables
    sta->Ivs(0) = (Ini.HasKey("z0") ? Ini("z0") : p+(q*q)/(p*M*M));

    // check initial yield function
    double f = YieldFunc (sta);
    if (f>1.0e-8) throw new Fatal("CamClay:InitIvs: stress point (sig=(%g,%g,%g,%g], p=%g, q=%g) is outside yield surface (f=%g) with z0=%g",sta->Sig(0),sta->Sig(1),sta->Sig(2),sta->Sig(3)/Util::SQ2,p,q,f,sta->Ivs(0));
}

inline void CamClay::Gradients (EquilibState const * Sta, Vec_t & V, Vec_t & Y) const
{
    // invariants
    double p,q,t;
    Vec_t dev_sig;
    OctInvs (Sta->Sig, p,q,t);
    Dev     (Sta->Sig, dev_sig);

    // gradients
    Y.change_dim (1);
    V    = (M*M*(Sta->Ivs(0)-2.0*p)/(Util::SQ3))*I + 2.0*dev_sig;
    Y(0) = -M*M*p;
}

inline double CamClay::YieldFunc (EquilibState const * Sta) const
{
    double p,q,t;
    OctInvs (Sta->Sig, p,q,t);
    return q*q + (p - Sta->Ivs(0))*p*M*M;
}

inline double CamClay::FailCrit (EquilibState const * Sta) const
{
    double p,q,t;
    OctInvs (Sta->Sig, p,q,t);
    return q - M*p;
}


///////////////////////////////////////////////////////////////////////////////////////// Autoregistration /////


Model * CamClayMaker(int NDim, SDPair const & Prms) { return new CamClay(NDim,Prms); }

int CamClayRegister()
{
    ModelFactory["CamClay"] = CamClayMaker;
    MODEL.Set ("CamClay", (double)MODEL.Keys.Size());
    return 0;
}

int __CamClay_dummy_int = CamClayRegister();

#endif // MECHSYS_CAMCLAY_H
