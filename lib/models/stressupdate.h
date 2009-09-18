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
    enum Scheme_t { FE_t, ME_t }; ///< Integration scheme

    // Constructor
    StressUpdate (Model const * Mdl);

    // Methods
    void Update      (Vec_t const & DEps, State * Sta, Vec_t & DSig) const;
    void TangentIncs (EquilibState const * sta, Vec_t & deps, Vec_t & dsig, Vec_t & divs) const;

    // Data
    Model const * Mdl;

    // Constants for integration
    Scheme_t Scheme; ///< Scheme: FE_t (Forward-Euler), ME_t (Modified-Euler)
    size_t   NDiv;
    double   STOL;
    double   dTini;
    double   mMin;
    double   mMax;
    size_t   maxSS;
    bool     CDrift; ///< correct drift ?
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline StressUpdate::StressUpdate (Model const * TheMdl)
    : Mdl    (TheMdl),
      Scheme (FE_t),
      NDiv   (1),
	  STOL   (1.0e-5),
	  dTini  (1.0),
	  mMin   (0.1),
	  mMax   (10.0),
	  maxSS  (2000),
      CDrift (false)
{
}

inline void StressUpdate::Update (Vec_t const & DEps, State * Sta, Vec_t & DSig) const
{
    EquilibState * sta = static_cast<EquilibState*>(Sta);
    DSig = sta->Sig; // temporary copy to calculate increment later

    Vec_t dsig(sta->Sig.size());
    Vec_t divs(sta->Ivs.size());

    if (Scheme==FE_t)
    {
        Vec_t deps(DEps/NDiv);
        for (size_t i=0; i<NDiv; ++i)
        {
            TangentIncs (sta, deps, dsig, divs);
            sta->Sig += dsig;
            sta->Eps += deps;
            sta->Ivs += divs;
        }
    }
    else throw new Fatal("StressUpdate::Update: Scheme==ME_t is not available yet");

    // return total stress increment
    DSig = sta->Sig - DSig;
}

inline void StressUpdate::TangentIncs (EquilibState const * sta, Vec_t & deps, Vec_t & dsig, Vec_t & divs) const
{
    Mat_t        D;
    Array<Vec_t> d;
    Mdl->Stiffness (sta, D, &d);
    dsig = D * deps;
    for (size_t k=0; k<sta->Ivs.size(); ++k) divs(k) = dot(d[k],deps);

    // calculate dez for plane-stress
    if (Mdl->GTy==pse_t) deps(2) = Mdl->CalcDEz(dsig);
}

#endif // MECHSYS_STRESSUPDATE_H
