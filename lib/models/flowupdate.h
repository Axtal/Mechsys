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

#ifndef MECHSYS_FLOWUPDATE_H
#define MECHSYS_FLOWUPDATE_H

// Std Lib
#include <iostream>

// MechSys
#include <mechsys/models/model.h>
#include <mechsys/models/flowstate.h>

class FlowUpdate
{
public:
    // enum
    enum Scheme_t { FE_t, ME_t }; ///< Integration scheme

    // Constructor
    FlowUpdate (Model const * Mdl);

    // Methods
    void Update (Vec_t const & DGra, State * Sta, Vec_t & DVel) const;

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
};


/////////////////////////////////////////////////////////////////////////////////////////// Implementation /////


inline FlowUpdate::FlowUpdate (Model const * TheMdl)
    : Mdl    (TheMdl),
      Scheme (FE_t),
      NDiv   (1),
	  STOL   (1.0e-5),
	  dTini  (1.0),
	  mMin   (0.1),
	  mMax   (10.0),
	  maxSS  (2000)
{
}

inline void FlowUpdate::Update (Vec_t const & DGra, State * Sta, Vec_t & DVel) const
{
    FlowState * sta = static_cast<FlowState*>(Sta);
    DVel = sta->Vel; // temporary copy to calculate increment later

    Vec_t dvel(size(sta->Vel));
    Vec_t divs(size(sta->Ivs));

    if (Scheme==FE_t)
    {
        Vec_t dgra(DGra/NDiv);
        for (size_t i=0; i<NDiv; ++i)
        {
            Mdl->TgIncs (sta, dgra, dvel, divs);
            sta->Vel += dvel;
            sta->Gra += dgra;
            sta->Ivs += divs;
        }
    }
    else throw new Fatal("FlowUpdate::Update: Scheme==ME_t is not available yet");

    // return total stress increment
    DVel = sta->Vel - DVel;
}

#endif // MECHSYS_FLOWUPDATE_H
