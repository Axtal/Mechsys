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

#ifndef MECHSYS_STRESSUPDATE_H
#define MECHSYS_STRESSUPDATE_H

// Std Lib
#include <iostream>

// MechSys
#include "models/model.h"
#include "models/equilibstate.h"

class StressUpdate
{
public:
    // enum
    enum Scheme_t { ME_t }; ///< Integration scheme

    // Constructor
    StressUpdate (Model const * Mdl);

    // Methods
    void Update      (Vec_t const & DEps, State * Sta, Vec_t & DSig) const;
    void TangentIncs (EquilibState const * sta, Vec_t & deps, Vec_t & dsig, Vec_t & divs) const;

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


inline StressUpdate::StressUpdate (Model const * TheMdl)
    : Mdl    (TheMdl),
      Scheme (ME_t),
      STOL   (1.0e-5),
      dTini  (1.0),
      mMin   (0.1),
      mMax   (10.0),
      MaxSS  (20),
      CDrift (true)
{
}

inline void StressUpdate::Update (Vec_t const & DEps, State * Sta, Vec_t & DSig) const
{
    // current state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    DSig = sta->Sig; // temporary copy to calculate increment later

    // constants
    size_t ncp = sta->Sig.size(); // num of stress components
    size_t niv = sta->Ivs.size(); // num of internal variables

    // auxiliar variables
    Vec_t dsig(ncp);
    Vec_t deps(ncp);
    Vec_t divs(niv);

    if (Scheme==ME_t)
    {
        // auxiliar variables
        EquilibState sta_1 (Mdl->NDim);              // intermediate state
        EquilibState sta_ME(Mdl->NDim);              // Modified-Euler state
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
            deps = aint*DEps;
            TangentIncs (sta, deps, dsig, divs);
            sta->Eps += deps;
            sta->Sig += dsig;
            sta->Ivs += divs;
            deps = fabs(1.0-aint)*DEps; // remaining of DEps to be applied

            // drift correction
            if (CDrift) Mdl->CorrectDrift (sta);
        }
        else deps = DEps; // update with full DEps

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
            deps_1 = dT*deps;
            TangentIncs (sta, deps_1, dsig_1, divs_1);
            sta_1.Eps = sta->Eps + deps_1;
            sta_1.Sig = sta->Sig + dsig_1;
            sta_1.Ivs = sta->Ivs + divs_1;
            TangentIncs (&sta_1, deps_1, dsig_2, divs_2);
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
        if (k>=MaxSS) throw new Fatal("StressUpdate::Update: Modified-Euler (local) did not converge after %d substeps",k);
    }
    else throw new Fatal("StressUpdate::Update: Scheme is not available yet");

    // return total stress increment
    DSig = sta->Sig - DSig;
}

inline void StressUpdate::TangentIncs (EquilibState const * sta, Vec_t & deps, Vec_t & dsig, Vec_t & divs) const
{
    Array<double> h;
    Mat_t         D;
    Vec_t         d;
    Mdl->Stiffness (sta, D, &h, &d);
    dsig = D * deps;
    for (size_t k=0; k<sta->Ivs.size(); ++k) divs(k) = h[k]*dot(d,deps);

    // calculate dez for plane-stress
    if (Mdl->GTy==pse_t) deps(2) = Mdl->CalcDEz(dsig);
}

#endif // MECHSYS_STRESSUPDATE_H
