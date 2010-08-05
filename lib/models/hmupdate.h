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

#ifndef MECHSYS_HMUPDATE_H
#define MECHSYS_HMUPDATE_H

// Std Lib
#include <iostream>

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/hmstate.h>

class HMUpdate
{
public:
    // enum
    enum Scheme_t { ME_t, SingleFE_t }; ///< Integration scheme

    // Constructor
    HMUpdate (Model const * Mdl);

    // Methods
    void Update      (double Dpw, Vec_t const & DEps, State * Sta, Vec_t & DSig) const;
    void TangentIncs (HMState const * sta, double dpw, Vec_t const & deps, Vec_t & dsig, Vec_t & divs) const;

    // Data
    Model const * Mdl;

    // Constants for integration
    Scheme_t Scheme; ///< Scheme: ME_t (Modified-Euler)
    double   STOL;
    double   dTini;
    double   mMin;
    double   mMax;
    size_t   MaxSS;
    bool     CDrift; ///< correct drift ?
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline HMUpdate::HMUpdate (Model const * TheMdl)
    : Mdl    (TheMdl),
      Scheme (ME_t),
      STOL   (1.0e-5),
      dTini  (1.0),
      mMin   (0.1),
      mMax   (10.0),
      MaxSS  (2000),
      CDrift (true)
{
}

inline void HMUpdate::Update (double Dpw, Vec_t const & DEps, State * Sta, Vec_t & DSig) const
{
    // current state
    HMState * sta = static_cast<HMState*>(Sta);
    DSig = sta->Sig; // temporary copy to calculate increment later

    // constants
    size_t ncp = size(sta->Sig); // num of stress components
    size_t niv = size(sta->Ivs); // num of internal variables

    // auxiliar variables
    double dpw;
    Vec_t dsig(ncp);
    Vec_t deps(ncp);
    Vec_t divs(niv);

    if (Scheme==SingleFE_t) // without intersection detection (should be used for linear elasticity only)
    {
        dpw  = Dpw;
        deps = DEps;
        TangentIncs (sta, dpw, deps, dsig, divs);
        sta->pw  += dpw;
        sta->Eps += deps;
        sta->Sig += dsig;
        sta->Ivs += divs;
    }
    else if (Scheme==ME_t)
    {
        // auxiliar variables
        double dpw_1;                                // intermediate increments
        HMState sta_1 (Mdl->NDim);                   // intermediate state
        HMState sta_ME(Mdl->NDim);                   // Modified-Euler state
        Vec_t deps_1(ncp), dsig_1(ncp), divs_1(niv); // intermediate increments
        Vec_t dsig_2(ncp), divs_2(niv);              // ME increments
        Vec_t sig_dif(ncp);                          // ME - FE stress difference

        // loading-unloading ?
        double aint = -1.0; // no intersection
        bool   ldg  = Mdl->LoadCond (sta, DEps, aint); // returns true if there is loading (also when there is intersection)

        // with intersection ?
        if (aint>0.0 && aint<1.0)
        {
            // update to intersection
            dpw  = aint*Dpw;
            deps = aint*DEps;
            TangentIncs (sta, dpw, deps, dsig, divs);
            sta->Eps += deps;
            sta->Sig += dsig;
            sta->Ivs += divs;
            dpw  = fabs(1.0-aint)*Dpw;  // remaining of Dpw to be applied
            deps = fabs(1.0-aint)*DEps; // remaining of DEps to be applied

            // drift correction
            if (CDrift) Mdl->CorrectDrift (sta);
        }
        else
        {
            dpw  = Dpw;  // update with full Dpw
            deps = DEps; // update with full DEps
        }

        // set loading flag (must be after intersection because the TgIncs during intersection must be calc with Ldg=false)
        sta  ->Ldg = ldg;
        sta_1 .Ldg = ldg;
        sta_ME.Ldg = ldg;

        // for each pseudo time T
        double T  = 0.0;
        double dT = dTini;
        size_t k  = 0;
        for (k=0; k<MaxSS; ++k)
        {
            // exit point
            if (T>=1.0) break;

            // FE and ME increments
            dpw_1  = dT*dpw;
            deps_1 = dT*deps;
            TangentIncs (sta, dpw_1, deps_1, dsig_1, divs_1);
            sta_1.pw  = sta->pw  + dpw_1;
            sta_1.Eps = sta->Eps + deps_1;
            sta_1.Sig = sta->Sig + dsig_1;
            sta_1.Ivs = sta->Ivs + divs_1;
            TangentIncs (&sta_1, dpw_1, deps_1, dsig_2, divs_2);
            sta_ME.Sig = sta->Sig + 0.5*(dsig_1+dsig_2);
            sta_ME.Ivs = sta->Ivs + 0.5*(divs_1+divs_2);

            // local error estimate
            sig_dif = sta_ME.Sig - sta_1.Sig;
            double sig_err = Norm(sig_dif)/(1.0+Norm(sta_ME.Sig));
            double ivs_err = 0.0;
            for (size_t i=0; i<niv; ++i) ivs_err += fabs(sta_ME.Ivs(i)-sta_1.Ivs(i))/(1.0+fabs(sta_ME.Ivs(i)));
            double error = sig_err + ivs_err;

            // step multiplier
            double m = (error>0.0 ? 0.9*sqrt(STOL/error) : mMax);

            // update
            if (error<STOL)
            {
                // update state
                T += dT;
                sta->pw  = sta_1 .pw;
                sta->Eps = sta_1 .Eps;
                sta->Sig = sta_ME.Sig;
                sta->Ivs = sta_ME.Ivs;

                // drift correction
                if (CDrift) Mdl->CorrectDrift (sta);

                // limit change on stepsize
                if (m>mMax) m = mMax;
            }
            else if (m<mMin) m = mMin;

            // change next step size
            dT = m * dT;

            // check for last increment
            if (dT>1.0-T) dT = 1.0-T;
        }
        if (k>=MaxSS) throw new Fatal("HMUpdate::Update: Modified-Euler (local) did not converge after %d substeps",k);
    }
    else throw new Fatal("HMUpdate::Update: Scheme is not available yet");

    // return total effective stress increment
    DSig = sta->Sig - DSig;
}

inline void HMUpdate::TangentIncs (HMState const * sta, double dpw, Vec_t const & deps, Vec_t & dsig, Vec_t & divs) const
{
    Mat_t D;
    Vec_t Dw;
    if (size(sta->Ivs)>0)
    {
        Vec_t h, d;
        Mdl->Stiffness (sta, D, Dw, &h, &d);
        divs = h*dot(d,deps);
    }
    else Mdl->Stiffness (sta, D, Dw);
    dsig = D*deps + Dw*dpw;
}

#endif // MECHSYS_HMUPDATE_H