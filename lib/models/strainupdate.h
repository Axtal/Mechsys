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

#ifndef MECHSYS_STRAINUPDATE_H
#define MECHSYS_STRAINUPDATE_H

// Std Lib
#include <iostream>

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/equilibstate.h>

class StrainUpdate
{
public:
    // enum
    enum Scheme_t { ME_t, SingleFE_t }; ///< Integration scheme

    // Constructor
    StrainUpdate (Model const * Mdl);

    // Methods
    void Update (Vec_t const & DSig, State * Sta, Vec_t & DEps) const;

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


inline StrainUpdate::StrainUpdate (Model const * TheMdl)
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

inline void StrainUpdate::Update (Vec_t const & DSig, State * Sta, Vec_t & DEps) const
{
    // current state
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    DEps = sta->Eps; // temporary copy to calculate increment later

    // constants
    size_t ncp = size(sta->Sig); // num of stress components
    size_t niv = size(sta->Ivs); // num of internal variables

    // auxiliar variables
    Vec_t dsig(ncp);
    Vec_t deps(ncp);
    Vec_t divs(niv);

    if (Scheme==SingleFE_t) // without intersection detection (should be used for linear elasticity only)
    {
        dsig = DSig;
        Mdl->InvTgIncs (sta, dsig, deps, divs);
        sta->Eps += deps;
        sta->Sig += dsig;
        sta->Ivs += divs;
    }
    else if (Scheme==ME_t)
    {
        // auxiliar variables
        EquilibState sta_1 (Mdl->NDim);              // intermediate state
        EquilibState sta_ME(Mdl->NDim);              // Modified-Euler state
        Vec_t deps_1(ncp), dsig_1(ncp), divs_1(niv); // intermediate increments
        Vec_t deps_2(ncp), divs_2(niv);              // ME increments
        Vec_t eps_dif(ncp);                          // ME - FE strain difference

        // update with full DSig
        dsig = DSig;

        // for each pseudo time T
        double T  = 0.0;
        double dT = dTini;
        size_t k  = 0;
        for (k=0; k<MaxSS; ++k)
        {
            // exit point
            if (T>=1.0) break;

            // FE and ME increments
            dsig_1 = dT*dsig;
            Mdl->InvTgIncs (sta, dsig_1, deps_1, divs_1);
            sta_1.Eps = sta->Eps + deps_1;
            sta_1.Sig = sta->Sig + dsig_1;
            sta_1.Ivs = sta->Ivs + divs_1;
            Mdl->InvTgIncs (&sta_1, dsig_1, deps_2, divs_2);
            sta_ME.Eps = sta->Eps + 0.5*(deps_1+deps_2);
            sta_ME.Ivs = sta->Ivs + 0.5*(divs_1+divs_2);

            // local error estimate
            eps_dif = sta_ME.Eps - sta_1.Eps;
            double eps_err = Norm(eps_dif)/(1.0+Norm(sta_ME.Eps));
            double ivs_err = 0.0;
            for (size_t i=0; i<niv; ++i) ivs_err += fabs(sta_ME.Ivs(i)-sta_1.Ivs(i))/(1.0+fabs(sta_ME.Ivs(i)));
            double error = eps_err + ivs_err;

            // step multiplier
            double m = (error>0.0 ? 0.9*sqrt(STOL/error) : mMax);

            // update
            if (error<STOL)
            {
                // update state
                T += dT;
                sta->Eps = sta_ME.Eps;
                sta->Sig = sta_1 .Sig;
                sta->Ivs = sta_ME.Ivs;

                // drift correction
                if (CDrift) Mdl->CorrectDrift (sta);

                // update stress path in model
                Mdl->UpdatePath (sta, Vec_t(0.5*(deps_1+deps_2)), dsig_1);

                // limit change on stepsize
                if (m>mMax) m = mMax;
            }
            else if (m<mMin) m = mMin;

            // change next step size
            dT = m * dT;

            // check for last increment
            if (dT>1.0-T) dT = 1.0-T;
        }
        if (k>=MaxSS) throw new Fatal("StrainUpdate::Update: Modified-Euler (local) did not converge after %d substeps",k);
    }
    else throw new Fatal("StrainUpdate::Update: Scheme is not available yet");

    // return total stress increment
    DEps = sta->Eps - DEps;
}

#endif // MECHSYS_STRAINUPDATE_H
